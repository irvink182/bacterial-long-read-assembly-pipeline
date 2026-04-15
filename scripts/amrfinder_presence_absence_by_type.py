#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

def resolve_amrfinder_columns(df):
    # SAMPLE/NAME
    if "Name" in df.columns:
        sample_col = "Name"
    elif "Sample" in df.columns:
        sample_col = "Sample"
    else:
        raise ValueError("Missing sample column (Name/Sample)")

    # GENE/ELEMENT
    if "Element symbol" in df.columns:
        gene_col = "Element symbol"
    elif "Gene symbol" in df.columns:
        gene_col = "Gene symbol"
    else:
        raise ValueError("Missing gene column (Element symbol / Gene symbol)")

    # TYPE/ELEMENT TYPE
    if "Type" in df.columns:
        type_col = "Type"
    elif "Element type" in df.columns:
        type_col = "Element type"
    else:
        raise ValueError("Missing type column (Type / Element type)")

    # CLASS
    class_col = "Class" if "Class" in df.columns else None

    return sample_col, gene_col, type_col, class_col

def build_matrix(df, sample_col, gene_col):
    if df.empty:
        return pd.DataFrame(columns=["sample"])

    matrix = (
        df.assign(value=1)
        .pivot_table(index=sample_col, columns=gene_col, values="value", fill_value=0)
        .astype(int)
        .reset_index()
        .rename(columns={sample_col: "sample"})
    )

    return matrix

def main():
    ap = argparse.ArgumentParser(
        description="Build presence/absence matrices per AMRFinder type."
    )
    ap.add_argument("-i", "--input", required=True, help="AMRFinder combined TSV.")
    ap.add_argument("-o", "--outprefix", required=True,
                    help="Output prefix. Creates <prefix>.AMR.tsv, <prefix>.STRESS.tsv, etc.")
    ap.add_argument("--types", default="AMR,STRESS,VIRULENCE",
                    help="Comma-separated list of Types to export (default: AMR,STRESS,VIRULENCE).")
    ap.add_argument("--min-prevalence", type=int, default=1,
                    help="Keep symbols present in at least N samples (default: 1).")
    args = ap.parse_args()

    # Load dataframe
    df = pd.read_csv(args.input, sep="\t", dtype=str)

    # Remove duplicated header row(s)
    if "Name" in df.columns:
        df = df[df["Name"] != "Name"].copy()

    # Strip whitespace
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()

    # Resolve columns names
    sample_col, gene_col, type_col, class_col = resolve_amrfinder_columns(df)

    print(f"[INFO] Using columns: sample={sample_col}, gene={gene_col}, type={type_col}")

    # Normalize missing values
    df[gene_col] = df[gene_col].replace({"nan": np.nan, "": np.nan})
    df = df[df[gene_col].notna()].copy()

    # Normalize type
    df[type_col] = df[type_col].str.upper()

    type_list = [t.strip().upper() for t in args.types.split(",")]

    for t in type_list:
        print(f"[INFO] Processing type: {t}")

        df_t = df[df[type_col] == t].copy()

        if df_t.empty:
            print(f"[WARN] No genes found for type {t}")
            out = pd.DataFrame(columns=["sample"])
        else:
            # Filter by prevalence
            counts = df_t.groupby(gene_col)[sample_col].nunique()
            keep = counts[counts >= args.min_prevalence].index
            df_t = df_t[df_t[gene_col].isin(keep)]

            out = build_matrix(df_t, sample_col, gene_col)

        out_file = f"{args.outprefix}.{t}.tsv"
        out.to_csv(out_file, sep="\t", index=False)

        print(f"[INFO] Written: {out_file}")

if __name__ == "__main__":
    main()
