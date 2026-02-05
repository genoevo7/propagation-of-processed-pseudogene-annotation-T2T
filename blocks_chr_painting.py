#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
blocks_chr_painting.py

Chromosome painting-style visualization of query-side block dispersion, using:
- UCSC Table Browser chain/net block exports (*.txt) for intervals (qName/qStart/qEnd/strand),
- NCBI Datasets "Chromosomes" TSV exports for per-chromosome lengths and labels.

Key features:
- Robust to RefSeq vs GenBank accession mismatches by auto-selecting the best keyspace.
- Handles cases where chromosome labels (e.g., chr7) are duplicated across haplotypes/patches
  (e.g., hap1/hap2) by automatically switching to a chrLabel keyspace so the plot is produced.
- Safe against empty blocks after filtering/merging (no crashes; skips gracefully).

Inputs (same folder; file names can be adjusted in `jobs`):
- <species>.txt                      (UCSC blocks export; contains qName/qStart/qEnd/strand...)
- <species>_chr_size.tsv             (NCBI Datasets Chromosomes TSV; has "Seq length",
                                     "Chromosome name", and RefSeq/GenBank accessions)

Outputs (per species):
- <stem>_chr_painting.pdf
- <stem>_chr_painting.png
- <stem>_accession_to_chr.tsv        (mapping table derived from the sizes TSV)
"""

import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, MultipleLocator


# ==========================
# STYLE
# ==========================
BG_ALPHA   = 0.02
HIT_ALPHA  = 1.0

COLOR_PLUS  = "#0047FF"   # vivid blue
COLOR_MINUS = "#E60000"   # vivid red
COLOR_MIX   = "#444444"

CHR_LW = 1.2
HIT_LW = 1.8

CHR_FONTSIZE   = 10.5
AXIS_FONTSIZE  = 11
TITLE_FONTSIZE = 13

ONLY_HIT_CHROMS = False
MERGE_INTERVALS = True
DROP_MT = True

MAJOR_MB = 50
MINOR_MB = 10

LABEL_X_FRACTION = 0.004

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"


# ==========================
# OPTIONAL GENE MARKERS
# ==========================
GENE_LOCS = {
    "Pygmy chimpanzee chromosome painting": {
        "SEPTIN14": ("chr6", 60635115, 60707848),
        "CIC":      ("chr20", 48373383, 48400734),
    },
    "Chimpanzee chromosome painting": {
        "SEPTIN14": ("chr6", 61565981, 61648745),
        "CIC":      ("chr20", 45184628, 45211932),
    },
    "Western lowland gorilla chromosome painting": {
        "SEPTIN14": ("chr6", 63974164, 64042674),
        "CIC":      ("chr20", 53876205, 53903501),
    },
    "Sumatran orangutan chromosome painting": {
        "SEPTIN14": ("chr6", 27604232, 27675065),
        "CIC":      ("chr20", 48189850, 48217143),
    },
    "Bornean orangutan chromosome painting": {
        "SEPTIN14": ("chr6", 25771729, 25850920),
        "CIC":      ("chr20", 48107503, 48134773),
    },
    "Siamang chromosome painting": {
        "SEPTIN14": ("chr9", 74547767, 74639845),
        "CIC":      ("chr17", 20182190, 20209520),
    },
    "White-tufted-ear marmoset chromosome painting": {
        "SEPTIN14": ("chr2", 512771, 558032),
        "CIC":      ("chr22", 35885382, 35913396),
    },
    "Ring-tailed lemur chromosome painting": {
        "SEPTIN14": ("chr2", 34400839, 34499849),
        "CIC":      ("chr19", 27902469, 27928691),
    },
    "Slow loris chromosome painting": {
        "SEPTIN14": ("chr12", 81025679, 81085424),
        "CIC":      ("chr10", 121666506, 121697202),
    },
}


def log(msg):
    print(msg, file=sys.stderr)


def normalize_chr_label(x):
    x = str(x).strip()
    if x.lower().startswith("chr"):
        return x
    if x.isdigit() or x in {"X", "Y"}:
        return "chr" + x
    return x


def is_mt_label(lab):
    return str(lab).upper() in {"CHRMT", "MT", "M", "CHRM"}


def is_primary_chr_label(lab: str) -> bool:
    s = str(lab).strip()
    if not s.lower().startswith("chr"):
        return False
    tail = s[3:]
    return tail.isdigit() or tail in {"X", "Y"}


def strip_version(acc: str) -> str:
    s = str(acc).strip()
    if "." in s:
        left, right = s.rsplit(".", 1)
        if right.isdigit():
            return left
    return s


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
# LOAD SIZES (robust to hap duplicates)
# --------------------------
def load_sizes(path):
    df = pd.read_csv(path, sep="\t", dtype=str)

    need = ["Seq length", "Chromosome name"]
    for c in need:
        if c not in df.columns:
            raise ValueError(f"Sizes TSV missing column '{c}' in {path}")

    has_gb = "GenBank seq accession" in df.columns
    has_rs = "RefSeq seq accession" in df.columns
    if not has_gb and not has_rs:
        raise ValueError(
            f"Sizes TSV missing BOTH 'GenBank seq accession' and 'RefSeq seq accession' in {path}"
        )

    df["Seq length"] = pd.to_numeric(df["Seq length"], errors="coerce")
    df = df.dropna(subset=["Seq length"])
    df["Seq length"] = df["Seq length"].astype(int)

    df["chrLabel"] = df["Chromosome name"].apply(normalize_chr_label)

    if DROP_MT:
        df = df[~df["chrLabel"].apply(is_mt_label)]

    df = df[df["chrLabel"].apply(is_primary_chr_label)]

    if df.empty:
        map_df = pd.DataFrame(columns=[
            "Chromosome_name", "chrLabel", "seqLen",
            "GenBank", "GenBank_stripped", "RefSeq", "RefSeq_stripped"
        ])
        return ({}, {}, {}, {}, {}, {}, {}, {}, map_df)

    # If chrLabel is duplicated (hap1/hap2 etc.), switch to chrLabel keyspace.
    use_chrlabel_keyspace = df["chrLabel"].duplicated().any()

    if has_gb:
        df["gb"] = df["GenBank seq accession"].astype(str).str.strip()
        df["gb_str"] = df["gb"].apply(strip_version)
    else:
        df["gb"] = ""
        df["gb_str"] = ""

    if has_rs:
        df["rs"] = df["RefSeq seq accession"].astype(str).str.strip()
        df["rs_str"] = df["rs"].apply(strip_version)
    else:
        df["rs"] = ""
        df["rs_str"] = ""

    map_df = df[["Chromosome name", "chrLabel", "Seq length", "gb", "gb_str", "rs", "rs_str"]].copy()
    map_df.columns = ["Chromosome_name", "chrLabel", "seqLen", "GenBank", "GenBank_stripped", "RefSeq", "RefSeq_stripped"]

    if not use_chrlabel_keyspace:
        # Accession keyspace
        size_map_gb, label_map_gb = {}, {}
        for _, r in df.iterrows():
            if r["gb"] and str(r["gb"]).lower() != "nan":
                size_map_gb[r["gb"]] = int(r["Seq length"])
                label_map_gb[r["gb"]] = r["chrLabel"]

        size_map_rs, label_map_rs = {}, {}
        for _, r in df.iterrows():
            if r["rs"] and str(r["rs"]).lower() != "nan":
                size_map_rs[r["rs"]] = int(r["Seq length"])
                label_map_rs[r["rs"]] = r["chrLabel"]

        rs_str_to_rs, gb_str_to_gb = {}, {}
        for _, r in df.iterrows():
            if r["rs_str"] and str(r["rs_str"]).lower() != "nan" and r["rs_str"] not in rs_str_to_rs:
                rs_str_to_rs[r["rs_str"]] = r["rs"]
            if r["gb_str"] and str(r["gb_str"]).lower() != "nan" and r["gb_str"] not in gb_str_to_gb:
                gb_str_to_gb[r["gb_str"]] = r["gb"]

        len_to_gb = {}
        if size_map_gb:
            tmp = {}
            for acc, L in size_map_gb.items():
                tmp.setdefault(L, []).append(acc)
            for L, accs in tmp.items():
                if len(accs) == 1:
                    len_to_gb[L] = accs[0]

        len_to_rs = {}
        if size_map_rs:
            tmp = {}
            for acc, L in size_map_rs.items():
                tmp.setdefault(L, []).append(acc)
            for L, accs in tmp.items():
                if len(accs) == 1:
                    len_to_rs[L] = accs[0]

        return (size_map_gb, label_map_gb, gb_str_to_gb, len_to_gb,
                size_map_rs, label_map_rs, rs_str_to_rs, len_to_rs,
                map_df)

    # chrLabel keyspace mode
    chr_len = df.groupby("chrLabel")["Seq length"].max().to_dict()

    size_map_gb = dict(chr_len)
    size_map_rs = dict(chr_len)

    label_map_gb = {k: k for k in chr_len.keys()}
    label_map_rs = {k: k for k in chr_len.keys()}

    gb_str_to_gb = {}
    rs_str_to_rs = {}
    for _, r in df.iterrows():
        lab = r["chrLabel"]
        if r.get("gb_str") and str(r["gb_str"]).lower() != "nan":
            gb_str_to_gb.setdefault(r["gb_str"], lab)
        if r.get("rs_str") and str(r["rs_str"]).lower() != "nan":
            rs_str_to_rs.setdefault(r["rs_str"], lab)

    len_to_gb = {}
    tmp = {}
    for lab, L in chr_len.items():
        tmp.setdefault(int(L), []).append(lab)
    for L, labs in tmp.items():
        if len(labs) == 1:
            len_to_gb[int(L)] = labs[0]
    len_to_rs = dict(len_to_gb)

    return (size_map_gb, label_map_gb, gb_str_to_gb, len_to_gb,
            size_map_rs, label_map_rs, rs_str_to_rs, len_to_rs,
            map_df)


# --------------------------
# LOAD BLOCKS (header/no-header)
# --------------------------
def load_blocks(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        first = f.readline().strip()

    has_header = (
        first.startswith("#chrom")
        or ("\tqName\t" in first)
        or first.lstrip().startswith("chrom\t")
    )

    if has_header:
        df = pd.read_csv(path, sep="\t", header=0, dtype=str)
        df.columns = [c.lstrip("#").strip() for c in df.columns]

        need_cols = {"qName", "qSize", "qStart", "qEnd", "strand"}
        if not need_cols.issubset(set(df.columns)):
            raise ValueError(
                f"Blocks header TSV missing columns {sorted(list(need_cols - set(df.columns)))} in {path}"
            )

        blocks = pd.DataFrame({
            "qName":  df["qName"].astype(str).str.strip(),
            "qSize":  pd.to_numeric(df["qSize"], errors="coerce"),
            "qStart": pd.to_numeric(df["qStart"], errors="coerce"),
            "qEnd":   pd.to_numeric(df["qEnd"], errors="coerce"),
            "strand": df["strand"].astype(str).str.strip()
        })
    else:
        df = pd.read_csv(path, sep="\t", comment="#", header=None)
        if df.shape[1] < 11:
            raise ValueError(f"Blocks file has too few columns (need >=11): {path}")

        blocks = pd.DataFrame({
            "qName":  df.iloc[:, 7].astype(str).str.strip(),
            "qSize":  pd.to_numeric(df.iloc[:, 8], errors="coerce"),
            "qStart": pd.to_numeric(df.iloc[:, 9], errors="coerce"),
            "qEnd":   pd.to_numeric(df.iloc[:,10], errors="coerce"),
            "strand": df.iloc[:, 5].astype(str).str.strip()
        })

    blocks = blocks.dropna(subset=["qName", "qStart", "qEnd"])
    if blocks.empty:
        # Return an empty DF while preserving expected columns
        blocks["qName_str"] = []
        blocks["qS"] = []
        blocks["qE"] = []
        return blocks.iloc[0:0]

    blocks["qStart"] = pd.to_numeric(blocks["qStart"], errors="coerce")
    blocks["qEnd"]   = pd.to_numeric(blocks["qEnd"], errors="coerce")
    blocks = blocks.dropna(subset=["qStart", "qEnd"])
    if blocks.empty:
        blocks["qName_str"] = []
        blocks["qS"] = []
        blocks["qE"] = []
        return blocks.iloc[0:0]

    blocks["qStart"] = blocks["qStart"].astype(int)
    blocks["qEnd"]   = blocks["qEnd"].astype(int)
    blocks["qS"] = blocks[["qStart", "qEnd"]].min(axis=1)
    blocks["qE"] = blocks[["qStart", "qEnd"]].max(axis=1)

    blocks.loc[~blocks["strand"].isin(["+", "-"]), "strand"] = "mix"
    blocks["qName_str"] = blocks["qName"].apply(strip_version)

    blocks["qSize"] = pd.to_numeric(blocks["qSize"], errors="coerce")
    blocks = blocks.dropna(subset=["qSize"])
    if blocks.empty:
        return blocks.iloc[0:0]

    blocks["qSize"] = blocks["qSize"].astype(int)
    return blocks


# --------------------------
# MERGE INTERVALS (safe on empty DF)
# --------------------------
def merge_intervals(df):
    if df is None or df.empty:
        return df.copy()

    out = []
    for q, sub in df.groupby("qName", sort=False):
        sub = sub.sort_values(["qS", "qE"])
        cs, ce = None, None
        strands = set()
        qstr = None
        qsize = None

        for _, r in sub.iterrows():
            s, e, st = int(r["qS"]), int(r["qE"]), r["strand"]
            if qstr is None:
                qstr = r.get("qName_str", strip_version(q))
            if qsize is None:
                qsize = int(r.get("qSize", 0)) if pd.notna(r.get("qSize", None)) else None

            if cs is None:
                cs, ce = s, e
                strands = {st}
            elif s <= ce:
                ce = max(ce, e)
                strands.add(st)
            else:
                out.append((q, qstr, qsize, cs, ce, strands))
                cs, ce = s, e
                strands = {st}

        if cs is not None:
            out.append((q, qstr, qsize, cs, ce, strands))

    rows = []
    for q, qstr, qsize, s, e, stset in out:
        strand = "+" if stset == {"+"} else "-" if stset == {"-"} else "mix"
        rows.append({"qName": q, "qName_str": qstr, "qSize": qsize, "qS": s, "qE": e, "strand": strand})

    cols = ["qName", "qName_str", "qSize", "qS", "qE", "strand"]
    if not rows:
        return df.iloc[0:0][cols].copy()

    return pd.DataFrame(rows, columns=cols)


def _add_gene_markers(ax, title, size_map, label_map, order, y_pos, bar_h):
    if title not in GENE_LOCS:
        return

    chrlabel_to_key = {}
    for k in order:
        lab = label_map.get(k, k)
        if lab not in chrlabel_to_key:
            chrlabel_to_key[lab] = k

    info = GENE_LOCS[title]

    def _mid(chrlab, s, e):
        if chrlab not in chrlabel_to_key:
            return None
        key = chrlabel_to_key[chrlab]
        L = size_map.get(key, None)
        if L is None:
            return None
        m = int((int(s) + int(e)) / 2)
        m = max(0, min(m, L))
        return key, m

    if "SEPTIN14" in info:
        ch, s, e = info["SEPTIN14"]
        got = _mid(ch, s, e)
        if got is not None:
            key, m = got
            y = y_pos[key] + bar_h/2
            ax.scatter([m], [y], marker="o", s=120, facecolors="none",
                       edgecolors="black", linewidths=2.2, zorder=10)

    if "CIC" in info:
        ch, s, e = info["CIC"]
        got = _mid(ch, s, e)
        if got is not None:
            key, m = got
            y = y_pos[key] + bar_h/2
            ax.scatter([m], [y], marker="x", s=180, c="black",
                       linewidths=2.4, zorder=10)


def plot_species(size_map, label_map, order, blocks, title, out_pdf, out_png):
    max_size = max(size_map.values())
    bar_h = 0.65
    y_gap = 0.25
    n = len(order)

    fig, ax = plt.subplots(figsize=(14, max(6, 0.18*n + 1.5)), dpi=300)
    y_pos = {c: (n-1-i)*(bar_h+y_gap) for i, c in enumerate(order)}
    label_x = -LABEL_X_FRACTION * max_size

    for c in order:
        y = y_pos[c]
        L = size_map[c]
        ax.add_patch(Rectangle((0, y), L, bar_h, fill=False, lw=CHR_LW, edgecolor="black"))
        ax.add_patch(Rectangle((0, y), L, bar_h, alpha=BG_ALPHA, lw=0))
        ax.text(label_x, y + bar_h/2, label_map.get(c, c),
                ha="right", va="center", fontsize=CHR_FONTSIZE, fontweight="bold")

    for _, r in blocks.iterrows():
        c = r["qName"]
        if c not in size_map:
            continue
        s = max(0, min(int(r["qS"]), size_map[c]))
        e = max(0, min(int(r["qE"]), size_map[c]))
        if e <= s:
            continue

        color = COLOR_PLUS if r["strand"] == "+" else COLOR_MINUS if r["strand"] == "-" else COLOR_MIX

        ax.add_patch(Rectangle((s, y_pos[c]), e-s, bar_h,
                               facecolor=color, edgecolor=color,
                               linewidth=HIT_LW, alpha=HIT_ALPHA))

    _add_gene_markers(ax, title, size_map, label_map, order, y_pos, bar_h)

    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{int(x/1e6)}"))
    ax.xaxis.set_major_locator(MultipleLocator(MAJOR_MB * 1e6))
    ax.xaxis.set_minor_locator(MultipleLocator(MINOR_MB * 1e6))

    ax.set_xlim(label_x*1.1, max_size)
    ax.set_ylim(-0.5, max(y_pos.values()) + bar_h + 0.5)
    ax.set_yticks([])
    ax.set_xlabel("Position (Mb)", fontsize=AXIS_FONTSIZE, fontweight="bold")
    ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight="bold")

    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_pdf)
    log(f"[DONE] {out_pdf}")
    plt.savefig(out_png, dpi=600)
    log(f"[DONE] {out_png}")
    plt.close()


def safe_stem(title: str) -> str:
    s = title.lower().replace("chromosome painting", "").strip()
    s = s.replace(" ", "_").replace("-", "_")
    while "__" in s:
        s = s.replace("__", "_")
    return s.strip("_")


def choose_keyspace(blocks, size_map_gb, gb_str_to_gb, size_map_rs, rs_str_to_rs):
    b_full = set(blocks["qName"])
    b_str  = set(blocks["qName_str"])

    gb_full = set(size_map_gb.keys())
    rs_full = set(size_map_rs.keys())
    gb_str  = set(gb_str_to_gb.keys())
    rs_str  = set(rs_str_to_rs.keys())

    gb_score = len(b_full & gb_full) + len(b_str & gb_str)
    rs_score = len(b_full & rs_full) + len(b_str & rs_str)

    return "refseq" if rs_score > gb_score else "genbank"


def remap_blocks_to_sizes(blocks, mode,
                          size_map_gb, gb_str_to_gb, len_to_gb,
                          size_map_rs, rs_str_to_rs, len_to_rs):
    if mode == "refseq":
        size_map = size_map_rs
        str_to_full = rs_str_to_rs
        len_to_full = len_to_rs
    else:
        size_map = size_map_gb
        str_to_full = gb_str_to_gb
        len_to_full = len_to_gb

    def fix(q, qs, qsize):
        if q in size_map:
            return q
        if qs in str_to_full and str_to_full[qs] in size_map:
            return str_to_full[qs]
        if qsize in len_to_full and len_to_full[qsize] in size_map:
            return len_to_full[qsize]
        return None

    blocks = blocks.copy()
    blocks["qName_fixed"] = [
        fix(q, qs, int(qsz))
        for q, qs, qsz in zip(blocks["qName"], blocks["qName_str"], blocks["qSize"])
    ]
    blocks = blocks.dropna(subset=["qName_fixed"])
    if blocks.empty:
        return blocks.iloc[0:0]

    blocks["qName_fixed"] = blocks["qName_fixed"].astype(str)
    blocks["qName"] = blocks["qName_fixed"]
    blocks = blocks.drop(columns=["qName_fixed"])
    return blocks


def run_one(title, sizes_name, blocks_name):
    sizes_path = find_existing_file(sizes_name)
    blocks_path = find_existing_file(blocks_name)

    if not sizes_path:
        log(f"[SKIP] sizes not found: {sizes_name}")
        return
    if not blocks_path:
        log(f"[SKIP] blocks not found: {blocks_name}")
        return

    try:
        (size_map_gb, label_map_gb, gb_str_to_gb, len_to_gb,
         size_map_rs, label_map_rs, rs_str_to_rs, len_to_rs,
         map_df) = load_sizes(sizes_path)

        if (not size_map_gb) and (not size_map_rs):
            raise ValueError("After filtering to primary chromosomes, sizes are empty (no chrN/X/Y rows).")

        blocks = load_blocks(blocks_path)

        if MERGE_INTERVALS:
            blocks = merge_intervals(blocks)

        if blocks is None or blocks.empty:
            log(f"[SKIP] {title}: blocks empty after load/merge (no intervals). file={blocks_path}")
            return

        if "qName" not in blocks.columns:
            log(f"[SKIP] {title}: blocks missing qName columns={list(blocks.columns)} file={blocks_path}")
            return

    except Exception as e:
        log(f"[SKIP] load failed for {title}: {e}")
        return

    mode = choose_keyspace(blocks, size_map_gb, gb_str_to_gb, size_map_rs, rs_str_to_rs)

    if mode == "refseq":
        log(f"[INFO] {title}: using RefSeq accessions (NC_...)")
        blocks = remap_blocks_to_sizes(
            blocks, "refseq",
            size_map_gb, gb_str_to_gb, len_to_gb,
            size_map_rs, rs_str_to_rs, len_to_rs
        )
        size_map = size_map_rs
        label_map = label_map_rs
        if not size_map:
            log(f"[WARN] {title}: RefSeq size_map empty; fallback to GenBank")
            mode = "genbank"

    if mode == "genbank":
        log(f"[INFO] {title}: using GenBank accessions (CM/... or chrLabel)")
        blocks = remap_blocks_to_sizes(
            blocks, "genbank",
            size_map_gb, gb_str_to_gb, len_to_gb,
            size_map_rs, rs_str_to_rs, len_to_rs
        )
        size_map = size_map_gb
        label_map = label_map_gb

    if blocks is None or blocks.empty:
        log(f"[SKIP] {title}: blocks empty after remap.")
        return

    matched = len(set(blocks["qName"]) & set(size_map.keys()))
    log(f"[INFO] {title}: matched={matched}")
    if matched == 0:
        log(f"[ERROR] still mismatch for {title}")
        log("  - example block qName (top10): " + ", ".join(list(blocks["qName"].unique())[:10]))
        log("  - sizes keys (top10): " + ", ".join(list(list(size_map.keys())[:10])))
        return

    order = list(size_map.keys())
    if ONLY_HIT_CHROMS:
        hit = set(blocks["qName"])
        order = [c for c in order if c in hit]

    stem = safe_stem(title)
    out_pdf = f"{stem}_chr_painting.pdf"
    out_png = f"{stem}_chr_painting.png"
    out_map = f"{stem}_accession_to_chr.tsv"

    map_df.to_csv(out_map, sep="\t", index=False)
    log(f"[DONE] {out_map}")

    plot_species(size_map, label_map, order, blocks, title, out_pdf, out_png)


def main():
    jobs = [
        ("Chimpanzee chromosome painting",                "chimpanzee_chr_size.tsv",                 "chimpanzee.txt"),
        ("Pygmy chimpanzee chromosome painting",          "pygmy chimpanzee_chr_size.tsv",           "pygmy chimpanzee.txt"),
        ("Western lowland gorilla chromosome painting",   "western lowland gorilla_chr_size.tsv",    "western lowland gorilla.txt"),
        ("Bornean orangutan chromosome painting",         "Bornean orangutan_chr_size.tsv",          "Bornean orangutan.txt"),
        ("Sumatran orangutan chromosome painting",        "Sumatran orangutan_chr_size.tsv",         "Sumatran orangutan.txt"),
        ("Siamang chromosome painting",                   "siamang_chr_size.tsv",                    "siamang.txt"),
        ("Slow loris chromosome painting",                "slow loris_chr_size.tsv",                 "slow loris.txt"),
        ("Ring-tailed lemur chromosome painting",         "Ring-tailed lemur_chr_size.tsv",          "Ring-tailed lemur.txt"),
        ("White-tufted-ear marmoset chromosome painting", "white-tufted-ear marmoset_chr_size.tsv",  "white-tufted-ear marmoset.txt"),
    ]

    for title, sizes_name, blocks_name in jobs:
        run_one(title, sizes_name, blocks_name)

    log("[ALL DONE]")


if __name__ == "__main__":
    main()
