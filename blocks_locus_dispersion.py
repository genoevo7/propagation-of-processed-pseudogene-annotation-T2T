#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
blocks_locus_dispersion.py

Compute locus dispersion metrics within a specified human T2T (hs1) target window
directly from UCSC chain/net block exports.

Design:
- Runs directly (IDLE-friendly); no argparse or command-line arguments.
- Uses only raw UCSC Table Browser block exports (*.txt).
- Does NOT use chromosome painting, chromosome size files, or merged blocks.
- Quantifies how widely a locus defined on the target (hs1) is dispersed
  across query genomes.

INPUT (files in the same directory):
- chimpanzee.txt
- pygmy chimpanzee.txt
- western lowland gorilla.txt
- Bornean orangutan.txt
- Sumatran orangutan.txt
- siamang.txt
- slow loris.txt
- Ring-tailed lemur.txt
- white-tufted-ear marmoset.txt

Expected UCSC Table Browser block format:
#chrom chromStart chromEnd name score strand tSize qName qSize qStart qEnd chainScore

OUTPUT:
- all_species_locus_by_qname.tsv
  (qName-level locus dispersion metrics)
- all_species_locus_species.tsv
  (species-level summary metrics)

Definitions:
- dupIndex = total_hit_bp / window_len_bp
- mix_fraction = overlap(+,-)_bp / hit_bp
- conflict_fraction = overlap(+,-)_bp / union(+,-)_bp
"""

import os
import sys
import glob
import pandas as pd


# ==========================
# User-defined parameters
# ==========================
TARGET_TNAME  = "chr7"
TARGET_TSTART = 55953761
TARGET_TEND   = 55959097

OUT_PREFIX = "all_species_locus"


def log(msg):
    print(msg, file=sys.stderr)


def find_existing_file(preferred_name: str):
    if os.path.exists(preferred_name):
        return preferred_name
    pat = preferred_name.replace(".tsv", "*.tsv").replace(".txt", "*.txt")
    cands = sorted(glob.glob(pat))
    if cands:
        return cands[0]
    key = os.path.splitext(preferred_name)[0]
    cands2 = sorted(glob.glob(key + "*"))
    return cands2[0] if cands2 else None


# --------------------------
# Interval utilities
# --------------------------
def union_intervals(intervals):
    """Merge overlapping intervals."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = []
    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            merged.append((cs, ce))
            cs, ce = s, e
    merged.append((cs, ce))
    return merged


def intervals_length(intervals):
    return sum(e - s for s, e in intervals)


def intersect_intervals(a, b):
    """Return merged intersections between two merged interval lists."""
    out = []
    i = j = 0
    while i < len(a) and j < len(b):
        s1, e1 = a[i]
        s2, e2 = b[j]
        s = max(s1, s2)
        e = min(e1, e2)
        if e > s:
            out.append((s, e))
        if e1 < e2:
            i += 1
        else:
            j += 1
    return union_intervals(out)


# --------------------------
# Load UCSC blocks
# --------------------------
def load_blocks_raw(path):
    """
    Supports:
    A) Header-present files with columns:
       chrom, chromStart, chromEnd, strand, qName, qSize, qStart, qEnd
    B) Header-less UCSC-style files with fixed column indices
    """
    df0 = pd.read_csv(path, sep="\t", comment="#", header=None, dtype=str)
    first_row = [str(x).strip() for x in df0.iloc[0].tolist()]
    has_header = (
        "chrom" in first_row and
        "chromStart" in first_row and
        "chromEnd" in first_row and
        "qName" in first_row
    )

    if has_header:
        df = pd.read_csv(path, sep="\t", comment="#", header=0, dtype=str)
        need = ["chrom", "chromStart", "chromEnd", "strand", "qName", "qSize", "qStart", "qEnd"]
        for c in need:
            if c not in df.columns:
                raise ValueError(f"Missing column '{c}' in blocks file: {path}")

        blocks = pd.DataFrame({
            "tName":   df["chrom"].str.strip(),
            "tStart":  pd.to_numeric(df["chromStart"], errors="coerce"),
            "tEnd":    pd.to_numeric(df["chromEnd"], errors="coerce"),
            "strand":  df["strand"].str.strip(),
            "qName":   df["qName"].str.strip(),
            "qSize":   pd.to_numeric(df["qSize"], errors="coerce"),
            "qStart":  pd.to_numeric(df["qStart"], errors="coerce"),
            "qEnd":    pd.to_numeric(df["qEnd"], errors="coerce"),
        })
    else:
        if df0.shape[1] < 11:
            raise ValueError(f"Blocks file has too few columns (need >=11): {path}")
        blocks = pd.DataFrame({
            "tName":   df0.iloc[:, 0].str.strip(),
            "tStart":  pd.to_numeric(df0.iloc[:, 1], errors="coerce"),
            "tEnd":    pd.to_numeric(df0.iloc[:, 2], errors="coerce"),
            "strand":  df0.iloc[:, 5].str.strip(),
            "qName":   df0.iloc[:, 7].str.strip(),
            "qSize":   pd.to_numeric(df0.iloc[:, 8], errors="coerce"),
            "qStart":  pd.to_numeric(df0.iloc[:, 9], errors="coerce"),
            "qEnd":    pd.to_numeric(df0.iloc[:,10], errors="coerce"),
        })

    blocks = blocks.dropna(subset=["tName", "tStart", "tEnd", "qName", "qSize", "qStart", "qEnd"])
    blocks = blocks.astype({
        "tStart": int, "tEnd": int,
        "qSize": int, "qStart": int, "qEnd": int
    })

    blocks["tS"] = blocks[["tStart", "tEnd"]].min(axis=1)
    blocks["tE"] = blocks[["tStart", "tEnd"]].max(axis=1)
    blocks["qS"] = blocks[["qStart", "qEnd"]].min(axis=1)
    blocks["qE"] = blocks[["qStart", "qEnd"]].max(axis=1)

    blocks.loc[~blocks["strand"].isin(["+", "-"]), "strand"] = "mix"
    return blocks


# --------------------------
# Crop blocks to target window
# --------------------------
def crop_to_window(blocks, tname, wstart, wend):
    wS, wE = min(wstart, wend), max(wstart, wend)
    sub = blocks[blocks["tName"] == tname].copy()
    if sub.empty:
        return sub.assign(qS_crop=pd.Series(dtype=int), qE_crop=pd.Series(dtype=int),
                          tS_crop=pd.Series(dtype=int), tE_crop=pd.Series(dtype=int))

    ovS = sub["tS"].clip(lower=wS, upper=wE)
    ovE = sub["tE"].clip(lower=wS, upper=wE)
    sub = sub.loc[ovE > ovS].copy()
    if sub.empty:
        return sub.assign(qS_crop=pd.Series(dtype=int), qE_crop=pd.Series(dtype=int),
                          tS_crop=pd.Series(dtype=int), tE_crop=pd.Series(dtype=int))

    sub["tS_crop"] = ovS.loc[sub.index].astype(int)
    sub["tE_crop"] = ovE.loc[sub.index].astype(int)

    sub["t_len"] = (sub["tE"] - sub["tS"]).astype(int)
    sub["q_len"] = (sub["qE"] - sub["qS"]).astype(int)

    def map_q_crop(row):
        tS, tE = row["tS"], row["tE"]
        qS, qE = row["qS"], row["qE"]
        tSc, tEc = row["tS_crop"], row["tE_crop"]
        t_len = max(1, row["t_len"])
        q_len = max(0, row["q_len"])

        left_off = max(0, tSc - tS)
        right_off = max(0, tE - tEc)
        if q_len <= 0 or tEc <= tSc:
            return None

        scale = q_len / t_len
        qs = qS + int(round(left_off * scale))
        qe = qE - int(round(right_off * scale))

        qs = max(qS, min(qE, qs))
        qe = max(qS, min(qE, qe))
        return (qs, qe) if qe > qs else None

    mapped = sub.apply(map_q_crop, axis=1)
    sub = sub.loc[mapped.notnull()].copy()
    sub["qS_crop"] = mapped.map(lambda x: x[0]).astype(int)
    sub["qE_crop"] = mapped.map(lambda x: x[1]).astype(int)
    return sub


# --------------------------
# Metric computation
# --------------------------
def compute_window_metrics(cropped, species_label, window_len):
    rows = []

    if cropped.empty:
        return (
            pd.DataFrame(columns=[
                "species","qName","chr_len_bp","window_len_bp",
                "hit_bp","union_pm_bp","overlap_pm_bp",
                "dupIndex","mix_fraction","conflict_fraction"
            ]),
            pd.DataFrame([{
                "species": species_label,
                "window_len_bp": window_len,
                "dispersion_qnames": 0,
                "total_hit_bp": 0,
                "total_union_pm_bp": 0,
                "total_overlap_pm_bp": 0,
                "dupIndex": 0.0,
                "mix_fraction": 0.0,
                "conflict_fraction": 0.0,
            }])
        )

    for qname, sub in cropped.groupby("qName", sort=False):
        chr_len = int(pd.Series(sub["qSize"]).value_counts().index[0])

        all_u = union_intervals(list(zip(sub["qS_crop"], sub["qE_crop"])))
        hit_bp = intervals_length(all_u)

        plus_u = union_intervals(list(zip(
            sub.loc[sub["strand"] == "+", "qS_crop"],
            sub.loc[sub["strand"] == "+", "qE_crop"]
        )))
        minus_u = union_intervals(list(zip(
            sub.loc[sub["strand"] == "-", "qS_crop"],
            sub.loc[sub["strand"] == "-", "qE_crop"]
        )))

        union_pm = union_intervals(plus_u + minus_u)
        overlap_bp = intervals_length(intersect_intervals(plus_u, minus_u))

        rows.append({
            "species": species_label,
            "qName": qname,
            "chr_len_bp": chr_len,
            "window_len_bp": window_len,
            "hit_bp": hit_bp,
            "union_pm_bp": intervals_length(union_pm),
            "overlap_pm_bp": overlap_bp,
            "dupIndex": hit_bp / window_len if window_len else 0.0,
            "mix_fraction": overlap_bp / hit_bp if hit_bp else 0.0,
            "conflict_fraction": overlap_bp / intervals_length(union_pm) if union_pm else 0.0,
        })

    by_q = pd.DataFrame(rows)

    sp = pd.DataFrame([{
        "species": species_label,
        "window_len_bp": window_len,
        "dispersion_qnames": by_q["qName"].nunique(),
        "total_hit_bp": by_q["hit_bp"].sum(),
        "total_union_pm_bp": by_q["union_pm_bp"].sum(),
        "total_overlap_pm_bp": by_q["overlap_pm_bp"].sum(),
        "dupIndex": by_q["hit_bp"].sum() / window_len if window_len else 0.0,
        "mix_fraction": by_q["overlap_pm_bp"].sum() / by_q["hit_bp"].sum() if by_q["hit_bp"].sum() else 0.0,
        "conflict_fraction": by_q["overlap_pm_bp"].sum() / by_q["union_pm_bp"].sum() if by_q["union_pm_bp"].sum() else 0.0,
    }])

    return by_q, sp


def main():
    wS, wE = min(TARGET_TSTART, TARGET_TEND), max(TARGET_TSTART, TARGET_TEND)
    window_len = wE - wS

    jobs = [
        ("Chimpanzee",                "chimpanzee.txt"),
        ("Pygmy chimpanzee",          "pygmy chimpanzee.txt"),
        ("Western lowland gorilla",   "western lowland gorilla.txt"),
        ("Bornean orangutan",         "Bornean orangutan.txt"),
        ("Sumatran orangutan",        "Sumatran orangutan.txt"),
        ("Siamang",                   "siamang.txt"),
        ("Slow loris",                "slow loris.txt"),
        ("Ring-tailed lemur",         "Ring-tailed lemur.txt"),
        ("White-tufted-ear marmoset", "white-tufted-ear marmoset.txt"),
    ]

    all_by = []
    all_sp = []

    for label, fname in jobs:
        path = find_existing_file(fname)
        if not path:
            log(f"[SKIP] {fname} not found")
            continue

        blocks = load_blocks_raw(path)
        cropped = crop_to_window(blocks, TARGET_TNAME, wS, wE)
        by_q, sp = compute_window_metrics(cropped, label, window_len)

        all_by.append(by_q)
        all_sp.append(sp)

    if all_by:
        pd.concat(all_by, ignore_index=True).to_csv(
            f"{OUT_PREFIX}_by_qname.tsv", sep="\t", index=False
        )
    if all_sp:
        pd.concat(all_sp, ignore_index=True).to_csv(
            f"{OUT_PREFIX}_species.tsv", sep="\t", index=False
        )


if __name__ == "__main__":
    main()
