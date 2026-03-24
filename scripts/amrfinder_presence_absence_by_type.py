#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd

def main():
    ap = argparse.ArgumentParser(
        description="Build presence/absence matrices per Type using AMRFinder Element symbol."
    )
    ap.add_argument("-i", "--input", required=True, help="AMRFinder combined TSV.")
    ap.add_argument("-o", "--outprefix", required=True,
                    help="Output prefix. Creates <prefix>.AMR.tsv, <prefix>.STRESS.tsv, etc.")
    ap.add_argument("--types", default="AMR,STRESS,VIRULENCE",
                    help="Comma-separated list of Types to export (default: AMR,STRESS,VIRULENCE).")
    ap.add_argument("--min-prevalence", type=int, default=1,
                    help="Keep symbols present in at least N samples (default: 1).")
    args = ap.parse_args()

    types = [t.strip() for t in args.types.split(",") if t.strip()]

    df = pd.read_csv(args.input, sep="\t", dtype=str)

    # Remove duplicated header row(s)
    if "Name" in df.columns:
        df = df[df["Name"] != "Name"].copy()

    # Strip whitespace
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()

    # Normalize missing values
    for c in ["Type", "Element symbol", "Name"]:
        if c in df.columns:
            df[c] = df[c].replace({"nan": np.nan, "": np.nan})

    required = {"Name", "Type", "Element symbol"}
    if not required.issubset(df.columns):
        raise SystemExit(f"Missing required columns: {sorted(required - set(df.columns))}")

    df = df[df["Name"].notna() & df["Type"].notna() & df["Element symbol"].notna()].copy()

    # Ensure uniqueness at (sample, symbol, type)
    df = df.drop_duplicates(subset=["Name", "Type", "Element symbol"])

    all_samples = sorted(df["Name"].unique().tolist())

    for t in types:
        sub = df[df["Type"] == t].copy()
        if sub.empty:
            print(f"[WARN] No rows for Type={t}. Skipping.")
            continue

        # Presence/absence: 1 if symbol present in sample
        mat = (
            sub.assign(present=1)
               .pivot_table(index="Name", columns="Element symbol", values="present",
                            aggfunc="max", fill_value=0)
        )

        # Reindex to include all samples (fill missing with 0)
        mat = mat.reindex(index=all_samples, fill_value=0)

        # Filter by prevalence (#samples with 1)
        prevalence = (mat > 0).sum(axis=0)
        mat = mat.loc[:, prevalence >= args.min_prevalence]

        out = f"{args.outprefix}.{t}.presence_absence.tsv"
        mat.reset_index().rename(columns={"Name": "sample"}).to_csv(out, sep="\t", index=False)
        print(f"Wrote {out}  (samples={mat.shape[0]}, symbols={mat.shape[1]})")

if __name__ == "__main__":
    main()
