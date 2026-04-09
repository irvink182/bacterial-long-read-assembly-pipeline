#!/usr/bin/env python3

import pandas as pd
import argparse
from pathlib import Path


def load_map(map_file):
    df = pd.read_csv(map_file, sep="\t", header=None, names=["old", "new"])
    df["old"] = df["old"].astype(str).str.strip()
    df["new"] = df["new"].astype(str).str.strip()
    return dict(zip(df["old"], df["new"]))


def load_meta(meta_file):
    if not Path(meta_file).exists() or Path(meta_file).stat().st_size == 0:
        return {}

    df = pd.read_csv(meta_file, sep="\t", header=None, names=["contig", "rotated", "gene"], engine="python")
    df["gene"] = df["gene"].fillna("NA")
    df["contig"] = df["contig"].astype(str).str.strip()
    return dict(zip(df["contig"], df["gene"]))


def load_mosdepth(mos_file):
    df = pd.read_csv(mos_file, sep="\t")
    
    df = df[~df["chrom"].str.contains("_region")]
    df = df[df["chrom"] != "total"]
    df["chrom"] = df["chrom"].astype(str).str.strip()

    return df[["chrom", "length", "mean"]]


def classify(row, min_cov, high_cov):
    cov = row["coverage"]
    gene = row["rotated_gene"]

    if cov < min_cov:
        return "low_cov"
    elif "dnaA" in gene:
        return "chromosome"
    elif "repA" in gene:
        return "plasmid_candidate"
    elif cov > high_cov:
        return "high_cov"
    else:
        return "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--mosdepth", required=True)
    ap.add_argument("--min-cov", type=float, required=True)
    ap.add_argument("--out", required=True)

    args = ap.parse_args()

    mapping = load_map(args.map)
    mapping = {str(k).strip(): str(v).strip() for k, v in mapping.items()}
    metadata = load_meta(args.meta)
    mos = load_mosdepth(args.mosdepth)

    rows = []

    for _, r in mos.iterrows():
        old = str(r["chrom"]).strip()

        if old not in mapping:
            continue

        new = mapping[old]

        gene = metadata.get(new, "NA")

        rows.append({
            "sample": args.sample,
            "contig": new,
            "length": r["length"],
            "coverage": r["mean"],
            "rotated_gene": gene
        })

    df = pd.DataFrame(rows)

    if df.empty:
        print("[WARNING] Empty contig stats")
        df.to_csv(args.out, sep="\t", index=False)
        return

    high_cov = args.min_cov * 10

    df["type"] = df.apply(lambda x: classify(x, args.min_cov, high_cov), axis=1)

    df = df.sort_values("length", ascending=False)

    df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()