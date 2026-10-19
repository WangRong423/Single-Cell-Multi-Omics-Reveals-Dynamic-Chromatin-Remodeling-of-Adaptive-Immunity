#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Flattened integrated B-cell figure.

Rows: Naive_B, Activated_B
Columns:
  A local triangular Hi-C log2(O/E)
  B DAR-anchored Hi-C O/E profile
  C longitudinal ATAC accessibility
  D ArchR peak-to-gene schematic

Key debug fixes:
  * Correct upper-triangle rotation: V=(X-Y)/2.
  * Use a compact ~240-kb window rather than an ~800-kb window.
  * Candidate_ID is the unique key.
  * Keep public GCB Hi-C labelled as a reference for Activated_B.
  * Plot real genomic separation on the y-axis labels (2*rotated height).
"""

from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import cooler
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
from matplotlib.ticker import FuncFormatter

BASE = Path("/public/users/wangrong/data/public_HiC/analysis/GSE84022")
DIRECT = BASE / "direct_contact"
FIGDIR = BASE / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

REP_FILE = DIRECT / "Bcell_representative_locus_per_subtype.tsv"
ATAC_FILE = DIRECT / "Bcell_top_candidates_ATAC_donor_timepoint.tsv"
NB_COOL = BASE / "cool" / "GSE84022_NB_hg18_10kb.cool"
GCB_COOL = BASE / "cool" / "GSE84022_GCB_hg18_10kb.cool"

OUT_PDF = FIGDIR / "Bcell_integrated_HiC_ATAC_P2G_flat.pdf"
OUT_PNG = FIGDIR / "Bcell_integrated_HiC_ATAC_P2G_flat.png"

RESOLUTION = 10000
BASE_WINDOW_BP = 240000
MARGIN_BP = 30000
HIC_VMAX = 2.5
TIMEPOINTS = [0,1,3,6,14,30,31,33,36,44]
FIG_WIDTH, FIG_HEIGHT = 14.6, 4.65
ASPECT_HIC, ASPECT_PROFILE, ASPECT_ATAC, ASPECT_P2G = 0.43, 0.52, 0.58, 0.43


def first_existing(df, names, required=True):
    for x in names:
        if x in df.columns:
            return x
    if required:
        raise KeyError(f"Missing all aliases: {names}")
    return None


def add_chr(x):
    x = str(x)
    return x if x.startswith("chr") else "chr" + x


for f in [REP_FILE, ATAC_FILE, NB_COOL, GCB_COOL]:
    if not f.exists():
        raise FileNotFoundError(f)

rep = pd.read_csv(REP_FILE, sep="\t")
atac = pd.read_csv(ATAC_FILE, sep="\t")

# Robust aliases for representative-locus fields.
COL = {
    "cid": first_existing(rep, ["Candidate_ID"]),
    "subtype": first_existing(rep, ["B_subtype"]),
    "dar_chr": first_existing(rep, ["DAR_chr"]),
    "dar_start": first_existing(rep, ["DAR_start"]),
    "dar_end": first_existing(rep, ["DAR_end"]),
    "peak_chr18": first_existing(rep, ["chr18_peak", "peak_chr_hg18"]),
    "peak18": first_existing(rep, ["peak_pos_hg18", "Peak_pos_hg18"]),
    "target_chr18": first_existing(rep, ["chr18_target", "target_chr_hg18"]),
    "target18": first_existing(rep, ["target_pos_hg18", "Target_pos_hg18"]),
    "target38": first_existing(rep, ["P2G_target_pos", "Target_pos_hg38", "target_pos_hg38", "P2G_target"], required=False),
    "gene": first_existing(rep, ["Target_gene", "gene", "Gene"], required=False),
    "p2g": first_existing(rep, ["P2G_cor", "value", "P2G_r"], required=False),
    "hic_oe": first_existing(rep, ["HiC_OE"], required=False),
    "hic_pct": first_existing(rep, ["HiC_same_distance_percentile"], required=False),
    "dar_comp": first_existing(rep, ["DAR_comparison", "comparison"], required=False),
    "dar_logfc": first_existing(rep, ["DAR_logFC", "logFC"], required=False),
}

clrs = {
    "Naive_B": cooler.Cooler(str(NB_COOL)),
    "Activated_B": cooler.Cooler(str(GCB_COOL)),
}
ref_label = {
    "Naive_B": "Public naive B Hi-C",
    "Activated_B": "Public GCB Hi-C reference",
}

chr_cache, exp_cache = {}, {}


def get_chr_data(clr, chrom):
    key = (clr.uri, chrom)
    if key not in chr_cache:
        mat = clr.matrix(balance=True, sparse=True).fetch(chrom).tocoo()
        bins = clr.bins().fetch(chrom)
        weights = bins["weight"].to_numpy(dtype=float)
        chr_cache[key] = (mat, weights)
    return chr_cache[key]


def get_expected(clr, chrom, d):
    d = int(abs(d))
    key = (clr.uri, chrom, d)
    if key in exp_cache:
        return exp_cache[key]
    mat, weights = get_chr_data(clr, chrom)
    n = len(weights)
    if d >= n:
        return np.nan
    valid = np.isfinite(weights)
    valid_pairs = valid if d == 0 else (valid[:-d] & valid[d:])
    n_valid = int(valid_pairs.sum())
    if n_valid == 0:
        return np.nan
    mask = (mat.col - mat.row) == d
    vals = mat.data[mask]
    vals = vals[np.isfinite(vals)]
    e = float(vals.sum() / n_valid)
    exp_cache[key] = e
    return e


def determine_window(chrom_length, p1, p2):
    required = abs(p2-p1) + 2*MARGIN_BP
    width = max(BASE_WINDOW_BP, int(np.ceil(required/RESOLUTION))*RESOLUTION)
    mid = int((p1+p2)/2)
    start = max(0, mid-width//2)
    start = (start//RESOLUTION)*RESOLUTION
    end = start+width
    if end > chrom_length:
        end = (chrom_length//RESOLUTION)*RESOLUTION
        start = max(0, end-width)
        start = (start//RESOLUTION)*RESOLUTION
    return int(start), int(end)


def local_oe_matrix(clr, chrom, start, end):
    region = f"{chrom}:{start}-{end}"
    obs = np.asarray(clr.matrix(balance=True, sparse=False).fetch(region), dtype=float)
    n = obs.shape[0]
    oe = np.full_like(obs, np.nan, dtype=float)
    for d in range(n):
        e = get_expected(clr, chrom, d)
        if not np.isfinite(e) or e <= 0:
            continue
        i = np.arange(0, n-d); j = i+d
        vals = obs[i, j]
        oe[i, j] = vals/e
        if d > 0:
            oe[j, i] = vals/e
    return oe


def plot_triangle(ax, oe, start, end, chrom, peak, target, title):
    with np.errstate(divide="ignore", invalid="ignore"):
        z0 = np.log2(oe)
    z0[~np.isfinite(z0)] = np.nan
    upper = np.triu(np.ones(z0.shape, dtype=bool), k=0)
    z = np.ma.masked_where(~upper, z0)
    n = z0.shape[0]
    edges = start + np.arange(n+1)*RESOLUTION
    X, Y = np.meshgrid(edges, edges)
    U = (X+Y)/2.0
    V = (X-Y)/2.0  # correct for upper triangle
    mesh = ax.pcolormesh(U, V, z, shading="auto", cmap="RdBu_r",
                         vmin=-HIC_VMAX, vmax=HIC_VMAX, rasterized=True)
    ax.axvline(peak, ls="--", lw=0.8, color="black")
    ax.axvline(target, ls=":", lw=0.8, color="black")
    cx, cy = (peak+target)/2.0, abs(target-peak)/2.0
    ax.scatter([cx], [cy], s=32, facecolor="#D73027", edgecolor="black", lw=0.6, zorder=10)
    ax.text(cx, cy+4500, "DAR-target", ha="center", va="bottom", fontsize=6.1)
    ax.set_xlim(start, end); ax.set_ylim(0, (end-start)/2.0)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1e6:.2f}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{int(2*y/1000)} kb"))
    ax.set_xlabel(f"{chrom} position (Mb, hg18)", fontsize=7)
    ax.set_ylabel("Genomic separation", fontsize=7)
    ax.set_title(title, fontsize=7.8, loc="left", pad=2)
    ax.tick_params(labelsize=6, length=2.3)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.set_box_aspect(ASPECT_HIC)
    return mesh


def get_anchor_profile(clr, chrom, start, end, anchor):
    binsize = int(clr.binsize)
    b = (int(anchor)//binsize)*binsize
    r_anchor = f"{chrom}:{b}-{b+binsize}"
    r_region = f"{chrom}:{start}-{end}"
    values = np.asarray(clr.matrix(balance=True, sparse=False).fetch(r_anchor, r_region)).reshape(-1)
    centers = start + np.arange(len(values))*binsize + binsize/2
    anchor_i = int((b-start)/binsize)
    oe = np.full(len(values), np.nan)
    for i, v in enumerate(values):
        e = get_expected(clr, chrom, abs(i-anchor_i))
        if np.isfinite(v) and np.isfinite(e) and e > 0:
            oe[i] = v/e
    return centers, oe


def plot_profile(ax, centers, oe, peak, target):
    ax.plot(centers/1e6, oe, lw=1.0)
    ax.axhline(1, ls="--", lw=0.7, color="0.5")
    ax.axvline(peak/1e6, ls="--", lw=0.7, color="black")
    ax.axvline(target/1e6, ls=":", lw=0.7, color="black")
    j = int(np.argmin(np.abs(centers-target)))
    if np.isfinite(oe[j]):
        ax.scatter(centers[j]/1e6, oe[j], s=18, zorder=5)
    ax.set_xlabel("Position (Mb, hg18)", fontsize=7)
    ax.set_ylabel("Hi-C O/E", fontsize=7)
    ax.set_title("DAR-anchored contact", fontsize=7.8, loc="left", pad=2)
    ax.tick_params(labelsize=6, length=2.3)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.set_box_aspect(ASPECT_PROFILE)


def plot_atac(ax, d, row):
    if "Donor" not in d.columns or "Day" not in d.columns or "log2CPM" not in d.columns:
        raise KeyError("ATAC table needs Donor, Day and log2CPM columns.")
    for _, g in d.groupby("Donor"):
        g = g.sort_values("Day")
        ax.plot(g.Day, g.log2CPM, lw=0.55, alpha=0.30, marker="o", ms=2.0)
    s = d.groupby("Day", as_index=False).agg(mean=("log2CPM","mean"), sd=("log2CPM","std"), n=("log2CPM","count"))
    s["sem"] = s["sd"]/np.sqrt(s["n"])
    ax.errorbar(s.Day, s["mean"], yerr=s.sem, marker="o", ms=3.3, lw=1.25, capsize=1.7)
    ax.axvline(30, ls=":", lw=0.7, color="0.5")
    ax.set_xticks(TIMEPOINTS)
    ax.tick_params(axis="x", labelrotation=45, labelsize=5.8, length=2.3)
    ax.tick_params(axis="y", labelsize=6, length=2.3)
    ax.set_xlabel("Day", fontsize=7); ax.set_ylabel(r"$\log_2$(CPM + 1)", fontsize=7)
    title = "Longitudinal ATAC accessibility"
    if COL["dar_comp"] and pd.notna(row[COL["dar_comp"]]):
        title += "\n" + str(row[COL["dar_comp"]])
    if COL["dar_logfc"] and pd.notna(row[COL["dar_logfc"]]):
        title += f" | DAR logFC={float(row[COL['dar_logfc']]):.2f}"
    ax.set_title(title, fontsize=7.8, loc="left", pad=2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.set_box_aspect(ASPECT_ATAC)


def plot_p2g(ax, row):
    p0, p1 = float(row[COL["dar_start"]]), float(row[COL["dar_end"]])
    pmid = (p0+p1)/2.0
    if COL["target38"] is None or pd.isna(row[COL["target38"]]):
        ax.text(0.5,0.5,"hg38 P2G target coordinate\nnot available", ha="center", va="center", transform=ax.transAxes, fontsize=7)
        ax.axis("off"); return
    target = float(row[COL["target38"]])
    left, right = min(p0,target), max(p1,target)
    span = max(right-left, 10000); margin=0.15*span
    ax.plot([left-margin,right+margin],[0,0],lw=0.7,color="black")
    ax.fill_between([p0,p1],[-0.025,-0.025],[0.025,0.025],alpha=0.7)
    ax.vlines(target,-0.045,0.045,lw=1.1)
    middle=(pmid+target)/2.0
    path=MplPath([(pmid,0.025),(middle,0.36),(target,0.025)],[MplPath.MOVETO,MplPath.CURVE3,MplPath.CURVE3])
    ax.add_patch(PathPatch(path,fill=False,lw=1.1))
    gene = str(row[COL["gene"]]) if COL["gene"] and pd.notna(row[COL["gene"]]) else "P2G target"
    ax.text(pmid,-0.075,"DAR",ha="center",va="top",fontsize=6.4)
    ax.text(target,-0.075,gene,ha="center",va="top",fontsize=6.4)
    ann=[]
    if COL["p2g"] and pd.notna(row[COL["p2g"]]): ann.append(f"P2G r={float(row[COL['p2g']]):.2f}")
    if COL["hic_oe"] and pd.notna(row[COL["hic_oe"]]): ann.append(f"Hi-C O/E={float(row[COL['hic_oe']]):.2f}")
    if COL["hic_pct"] and pd.notna(row[COL["hic_pct"]]): ann.append(f"{float(row[COL['hic_pct']]):.1f}th pct.")
    ax.text(0.5,0.91," | ".join(ann),transform=ax.transAxes,ha="center",va="top",fontsize=6.1)
    ax.set_xlim(left-margin,right+margin); ax.set_ylim(-0.13,0.48); ax.set_yticks([])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_: f"{x/1e6:.2f}"))
    ax.set_xlabel(f"{row[COL['dar_chr']]} position (Mb, hg38)",fontsize=7)
    ax.set_title("ArchR peak-to-gene association",fontsize=7.8,loc="left",pad=2)
    ax.tick_params(axis="x",labelsize=6,length=2.3)
    for s in ["top","right","left"]: ax.spines[s].set_visible(False)
    ax.set_box_aspect(ASPECT_P2G)


order = [x for x in ["Naive_B","Activated_B"] if x in set(rep[COL["subtype"]])]
if not order:
    raise ValueError("No Naive_B or Activated_B representative locus found.")

fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
gs = fig.add_gridspec(len(order), 4, width_ratios=[1.60,1.00,1.30,1.18], wspace=0.34, hspace=0.28)
meshes=[]

for ri, subtype in enumerate(order):
    row = rep.loc[rep[COL["subtype"]].eq(subtype)].iloc[0]
    cid = str(row[COL["cid"]])
    clr = clrs[subtype]
    chrom = add_chr(row[COL["peak_chr18"]]); chrom2=add_chr(row[COL["target_chr18"]])
    if chrom != chrom2: raise ValueError(f"{cid}: interchromosomal candidate")
    peak, target = int(float(row[COL["peak18"]])), int(float(row[COL["target18"]]))
    start,end = determine_window(int(clr.chromsizes[chrom]),peak,target)
    oe = local_oe_matrix(clr,chrom,start,end)

    ax1=fig.add_subplot(gs[ri,0]); meshes.append(plot_triangle(ax1,oe,start,end,chrom,peak,target,f"{subtype.replace('_',' ')}\n{ref_label[subtype]}"))
    ax2=fig.add_subplot(gs[ri,1]); centers,profile=get_anchor_profile(clr,chrom,start,end,peak); plot_profile(ax2,centers,profile,peak,target)
    ax3=fig.add_subplot(gs[ri,2]); sub=atac[(atac.Candidate_ID.astype(str)==cid) & (atac.B_subtype==subtype)].copy()
    if len(sub): plot_atac(ax3,sub,row)
    else:
        ax3.text(0.5,0.5,"No ATAC trajectory",ha="center",va="center",transform=ax3.transAxes,fontsize=7); ax3.set_box_aspect(ASPECT_ATAC)
    ax4=fig.add_subplot(gs[ri,3]); plot_p2g(ax4,row)

for x,label in zip([0.022,0.347,0.548,0.795],["A","B","C","D"]):
    fig.text(x,0.974,label,fontsize=10.5,fontweight="bold")
if len(order)>=1: fig.text(0.006,0.72,"Naive B",rotation=90,ha="center",va="center",fontsize=8.5,fontweight="bold")
if len(order)>=2: fig.text(0.006,0.275,"Activated B",rotation=90,ha="center",va="center",fontsize=8.5,fontweight="bold")
if meshes:
    cax=fig.add_axes([0.978,0.58,0.008,0.27]); cb=fig.colorbar(meshes[0],cax=cax); cb.set_label(r"Hi-C log$_2$(O/E)",fontsize=6.8); cb.ax.tick_params(labelsize=5.8,length=2)
fig.subplots_adjust(left=0.034,right=0.962,bottom=0.10,top=0.935)
fig.savefig(OUT_PDF,bbox_inches="tight")
fig.savefig(OUT_PNG,dpi=600,bbox_inches="tight")
plt.close(fig)
print("Saved:", OUT_PDF)
print("Saved:", OUT_PNG)
