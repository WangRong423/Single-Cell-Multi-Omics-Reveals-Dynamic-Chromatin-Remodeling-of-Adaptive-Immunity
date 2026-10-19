#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Score B-cell DAR–P2G candidates using public GSE84022 Hi-C.

Debugged principles retained here:
  * Native Hi-C matrices remain in hg18.
  * Only DAR/target candidate coordinates are lifted to hg18 upstream.
  * Naive_B uses the public naive-B cooler; Activated_B uses the public
    GCB cooler only as an activated/differentiating B-cell reference.
  * Contacts are matrix-balanced first.
  * Direct contact enrichment is O/E relative to valid bin pairs at the
    same chromosome-specific genomic separation.
  * The same-distance percentile is descriptive, not a formal p value.
  * Candidate_ID, not DAR_ID, is used as the unique key because one DAR
    can link to multiple P2G targets.

Expected input:
  direct_contact/Bcell_DAR_P2G_candidates_hg18.tsv

Required semantic fields (column aliases are accepted):
  Candidate_ID, B_subtype, hg18 peak chromosome/position,
  hg18 target chromosome/position, P2G_cor.
"""

from pathlib import Path
import math
import numpy as np
import pandas as pd
import cooler

BASE = Path("/public/users/wangrong/data/public_HiC/analysis/GSE84022")
DIRECT = BASE / "direct_contact"

INPUT = DIRECT / "Bcell_DAR_P2G_candidates_hg18.tsv"
NB_COOL = BASE / "cool" / "GSE84022_NB_hg18_10kb.cool"
GCB_COOL = BASE / "cool" / "GSE84022_GCB_hg18_10kb.cool"

OUT_ALL = DIRECT / "Bcell_DAR_P2G_direct_HiC_scores.tsv"
OUT_BEST = DIRECT / "Bcell_DAR_P2G_direct_HiC_best_per_DAR.tsv"
OUT_TOP = DIRECT / "Bcell_DAR_P2G_direct_HiC_top_candidates.tsv"
OUT_SUMMARY = DIRECT / "Bcell_DAR_P2G_direct_HiC_summary.tsv"


def first_existing(df, names, required=True):
    for x in names:
        if x in df.columns:
            return x
    if required:
        raise KeyError(f"None of these columns were found: {names}")
    return None


def add_chr(x):
    x = str(x)
    return x if x.startswith("chr") else "chr" + x


def robust_log2(x):
    return math.log2(x) if np.isfinite(x) and x > 0 else np.nan


if not INPUT.exists():
    raise FileNotFoundError(INPUT)
for f in (NB_COOL, GCB_COOL):
    if not f.exists():
        raise FileNotFoundError(f)

cand = pd.read_csv(INPUT, sep="\t")

cid_col = first_existing(cand, ["Candidate_ID"])
subtype_col = first_existing(cand, ["B_subtype", "subtype", "celltype"])
peak_chr_col = first_existing(cand, ["chr18_peak", "peak_chr_hg18", "Peak_chr_hg18"])
peak_pos_col = first_existing(cand, ["peak_pos_hg18", "Peak_pos_hg18", "peak_hg18"])
target_chr_col = first_existing(cand, ["chr18_target", "target_chr_hg18", "Target_chr_hg18"])
target_pos_col = first_existing(cand, ["target_pos_hg18", "Target_pos_hg18", "target_hg18"])
p2g_col = first_existing(cand, ["P2G_cor", "value", "P2G_r"], required=False)
dar_col = first_existing(cand, ["DAR_ID"], required=False)

if cand[cid_col].duplicated().any():
    raise ValueError("Candidate_ID must be unique. Do not merge/liftOver by DAR_ID alone.")

clrs = {
    "Naive_B": cooler.Cooler(str(NB_COOL)),
    "Activated_B": cooler.Cooler(str(GCB_COOL)),
}

# Cache chromosome matrices/weights and diagonal statistics.
chr_cache = {}
diag_cache = {}


def get_chr_data(clr, chrom):
    key = (clr.uri, chrom)
    if key in chr_cache:
        return chr_cache[key]
    mat = clr.matrix(balance=True, sparse=True).fetch(chrom).tocoo()
    bins = clr.bins().fetch(chrom)
    weights = bins["weight"].to_numpy(dtype=float)
    start0 = int(bins.iloc[0]["start"])
    binsize = int(clr.binsize)
    chr_cache[key] = (mat, weights, start0, binsize)
    return chr_cache[key]


def diagonal_stats(clr, chrom, d):
    """Return expected mean and empirical nonzero values at diagonal d.

    Missing sparse entries are zeros. Only bin pairs whose balancing weights
    are finite are counted as valid.
    """
    d = int(abs(d))
    key = (clr.uri, chrom, d)
    if key in diag_cache:
        return diag_cache[key]

    mat, weights, _, _ = get_chr_data(clr, chrom)
    n = len(weights)
    if d >= n:
        out = (np.nan, np.array([], dtype=float), 0, 0)
        diag_cache[key] = out
        return out

    valid = np.isfinite(weights)
    valid_pairs = valid if d == 0 else (valid[:-d] & valid[d:])
    n_valid = int(valid_pairs.sum())
    if n_valid == 0:
        out = (np.nan, np.array([], dtype=float), 0, 0)
        diag_cache[key] = out
        return out

    mask = (mat.col - mat.row) == d
    vals = mat.data[mask]
    vals = vals[np.isfinite(vals)]
    # Sparse entries not stored are zeros.
    expected = float(vals.sum() / n_valid)
    n_nonzero = len(vals)
    out = (expected, vals, n_valid, n_nonzero)
    diag_cache[key] = out
    return out


def focal_contact(clr, chrom, p1, p2):
    binsize = int(clr.binsize)
    b1 = (int(p1) // binsize) * binsize
    b2 = (int(p2) // binsize) * binsize
    r1 = f"{chrom}:{b1}-{b1 + binsize}"
    r2 = f"{chrom}:{b2}-{b2 + binsize}"
    x = clr.matrix(balance=True, sparse=False).fetch(r1, r2)
    return float(np.asarray(x)[0, 0]), abs(b2 - b1) // binsize


def empirical_percentile(obs, nonzero_vals, n_valid):
    if not np.isfinite(obs) or n_valid <= 0:
        return np.nan
    # All omitted sparse values are exactly zero. Balanced nonzero entries
    # are normally >=0; include implicit zeros in the percentile.
    n_zero = n_valid - len(nonzero_vals)
    le = n_zero + int(np.sum(nonzero_vals <= obs))
    return 100.0 * le / n_valid

rows = []
for _, r in cand.iterrows():
    subtype = str(r[subtype_col])
    if subtype not in clrs:
        raise ValueError(f"Unsupported B_subtype: {subtype}")
    clr = clrs[subtype]

    chrom1 = add_chr(r[peak_chr_col])
    chrom2 = add_chr(r[target_chr_col])
    p1 = int(float(r[peak_pos_col]))
    p2 = int(float(r[target_pos_col]))

    out = r.to_dict()
    if chrom1 != chrom2:
        out.update({
            "HiC_balanced_observed": np.nan,
            "HiC_expected_same_distance": np.nan,
            "HiC_OE": np.nan,
            "HiC_log2OE": np.nan,
            "HiC_same_distance_percentile": np.nan,
            "Direct_3D_Tier": "interchromosomal_not_scored",
        })
        rows.append(out)
        continue

    obs, d_bins = focal_contact(clr, chrom1, p1, p2)
    exp, vals, n_valid, n_nonzero = diagonal_stats(clr, chrom1, d_bins)
    oe = obs / exp if np.isfinite(obs) and np.isfinite(exp) and exp > 0 else np.nan
    pct = empirical_percentile(obs, vals, n_valid)

    # Descriptive prioritization tiers only; not statistical significance.
    if np.isfinite(oe) and np.isfinite(pct) and oe >= 2 and pct >= 95:
        tier = "Tier1"
    elif np.isfinite(oe) and np.isfinite(pct) and oe >= 1.5 and pct >= 90:
        tier = "Tier2"
    elif np.isfinite(oe) and oe > 1:
        tier = "Tier3"
    else:
        tier = "No_enrichment"

    p2g = abs(float(r[p2g_col])) if p2g_col and pd.notna(r[p2g_col]) else 0.0
    # Transparent ranking score; for locus prioritization only.
    score = (
        (robust_log2(oe) if np.isfinite(robust_log2(oe)) else -5.0)
        + (pct / 100.0 if np.isfinite(pct) else 0.0)
        + p2g
    )

    out.update({
        "HiC_balanced_observed": obs,
        "HiC_expected_same_distance": exp,
        "HiC_distance_bins": int(d_bins),
        "HiC_distance_bp": int(d_bins * clr.binsize),
        "HiC_n_valid_same_distance": n_valid,
        "HiC_n_nonzero_same_distance": n_nonzero,
        "HiC_OE": oe,
        "HiC_log2OE": robust_log2(oe),
        "HiC_same_distance_percentile": pct,
        "Direct_3D_Tier": tier,
        "CandidateScore": score,
    })
    rows.append(out)

res = pd.DataFrame(rows)
res.to_csv(OUT_ALL, sep="\t", index=False)

sort_cols = [c for c in ["HiC_OE", "HiC_same_distance_percentile", "CandidateScore"] if c in res.columns]
res_sorted = res.sort_values(sort_cols, ascending=[False] * len(sort_cols), na_position="last")
res_sorted.head(20).to_csv(OUT_TOP, sep="\t", index=False)

if dar_col:
    best = res_sorted.drop_duplicates(subset=[dar_col], keep="first")
else:
    best = res_sorted.copy()
best.to_csv(OUT_BEST, sep="\t", index=False)

summary = (
    res.groupby(subtype_col, dropna=False)
    .agg(
        n_candidates=(cid_col, "size"),
        n_OE_gt1=("HiC_OE", lambda x: int(np.sum(pd.to_numeric(x, errors="coerce") > 1))),
        median_OE=("HiC_OE", "median"),
        max_OE=("HiC_OE", "max"),
        median_percentile=("HiC_same_distance_percentile", "median"),
    )
    .reset_index()
)
summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

print("Saved:")
for x in [OUT_ALL, OUT_BEST, OUT_TOP, OUT_SUMMARY]:
    print(" ", x)
print("\nSummary:")
print(summary.to_string(index=False))
