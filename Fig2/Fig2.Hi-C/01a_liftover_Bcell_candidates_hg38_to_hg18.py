#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Lift B-cell DAR and P2G target coordinates from hg38 to hg18.

Critical debug fix: use Candidate_ID as the unique key. DAR_ID is not unique
because one DAR can have multiple P2G links.

Set HG38_TO_HG18_CHAIN if the chain file is not on the default path below.
"""

from pathlib import Path
import os, shutil, subprocess, tempfile
import pandas as pd

BASE = Path("/public/users/wangrong/data/public_HiC/analysis/GSE84022/direct_contact")
INPUT = BASE / "Bcell_DAR_P2G_candidates.tsv"
OUTPUT = BASE / "Bcell_DAR_P2G_candidates_hg18.tsv"
CHAIN = Path(os.environ.get("HG38_TO_HG18_CHAIN", "/public/users/wangrong/data/public_HiC/hg38ToHg18.over.chain.gz"))


def first_existing(df, xs):
    for x in xs:
        if x in df.columns: return x
    raise KeyError(f"Missing all aliases: {xs}")

if shutil.which("liftOver") is None:
    raise RuntimeError("UCSC liftOver is not in PATH.")
if not INPUT.exists():
    raise FileNotFoundError(f"Edit INPUT if needed; current path not found: {INPUT}")
if not CHAIN.exists():
    raise FileNotFoundError(f"Set HG38_TO_HG18_CHAIN; chain not found: {CHAIN}")

x = pd.read_csv(INPUT, sep="\t")
if "Candidate_ID" not in x.columns:
    x.insert(0, "Candidate_ID", [f"BcellCand_{i:05d}" for i in range(1, len(x)+1)])
if x["Candidate_ID"].duplicated().any():
    raise ValueError("Candidate_ID is not unique.")

chr_col = first_existing(x, ["DAR_chr","seqnames","Peak_chr"])
start_col = first_existing(x, ["DAR_start","Peak_start"])
end_col = first_existing(x, ["DAR_end","Peak_end"])
tchr_col = first_existing(x, ["Target_chr","P2G_target_chr","DAR_chr","seqnames"])
tpos_col = first_existing(x, ["P2G_target_pos","Target_pos","target_pos_hg38"])

with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    peak_bed = td / "peak.hg38.bed"
    target_bed = td / "target.hg38.bed"
    peak_out = td / "peak.hg18.bed"
    target_out = td / "target.hg18.bed"
    peak_un = td / "peak.unmapped.bed"
    target_un = td / "target.unmapped.bed"

    with peak_bed.open("w") as f:
        for _, r in x.iterrows():
            # BED is 0-based, half-open. Input DAR coordinates are treated as
            # 1-based genomic intervals; subtract 1 from start only.
            f.write(f"{r[chr_col]}\t{int(r[start_col])-1}\t{int(r[end_col])}\t{r['Candidate_ID']}\n")
    with target_bed.open("w") as f:
        for _, r in x.iterrows():
            p = int(r[tpos_col])
            f.write(f"{r[tchr_col]}\t{p-1}\t{p}\t{r['Candidate_ID']}\n")

    subprocess.run(["liftOver", str(peak_bed), str(CHAIN), str(peak_out), str(peak_un)], check=True)
    subprocess.run(["liftOver", str(target_bed), str(CHAIN), str(target_out), str(target_un)], check=True)

    p = pd.read_csv(peak_out, sep="\t", header=None, names=["chr18_peak","start18_peak0","end18_peak","Candidate_ID"])
    t = pd.read_csv(target_out, sep="\t", header=None, names=["chr18_target","start18_target0","end18_target","Candidate_ID"])

    # Require one unique mapping per candidate at both endpoints.
    p = p.groupby("Candidate_ID").filter(lambda g: len(g)==1).copy()
    t = t.groupby("Candidate_ID").filter(lambda g: len(g)==1).copy()
    p["peak_pos_hg18"] = ((p["start18_peak0"] + p["end18_peak"] - 1) // 2).astype(int)
    t["target_pos_hg18"] = t["start18_target0"].astype(int) + 1

    mapped = x.merge(p[["Candidate_ID","chr18_peak","peak_pos_hg18"]], on="Candidate_ID", how="inner")
    mapped = mapped.merge(t[["Candidate_ID","chr18_target","target_pos_hg18"]], on="Candidate_ID", how="inner")

mapped.to_csv(OUTPUT, sep="\t", index=False)
print(f"Input candidates: {len(x)}")
print(f"Uniquely mapped peak+target candidates: {len(mapped)}")
print("Saved:", OUTPUT)
