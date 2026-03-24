#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd

def uniq_join(series: pd.Series) -> str:
    vals = [x for x in series.dropna().astype(str) if x and x.lower() != "nan"]
    return ";".join(sorted(set(vals)))

def main():
    ap = argparse.ArgumentParser(
        description="Summarize AMRFinder results using UNIQUE Element symbols, split by Type and/or Class."
    )
    ap.add_argument("-i", "--input", required=True, help="AMRFinder combined TSV.")
    ap.add_argument("-o", "--output", required=True, help="Output TSV.")

    ap.add_argument("--include-type-counts", action="store_true", default=True,
                    help="Include unique Element symbol counts per Type (default: on).")
    ap.add_argument("--include-class-counts", action="store_true",
                    help="Include unique Element symbol counts per Class.")
    ap.add_argument("--include-type-lists", action="store_true",
                    help="Include lists of unique Element symbols per Type.")
    ap.add_argument("--include-class-lists", action="store_true",
                    help="Include lists of unique Element symbols per Class.")

    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str)

    # Remove duplicated header row(s)
    if "Name" in df.columns:
        df = df[df["Name"] != "Name"].copy()

    # Strip whitespace
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()

    # Normalize NA
    for c in ["Type", "Class", "Element symbol"]:
        if c in df.columns:
            df[c] = df[c].replace({"nan": np.nan, "": np.nan})

    # Keep rows with Element symbol
    if "Element symbol" not in df.columns or "Name" not in df.columns:
        raise SystemExit("Missing required columns: Name and/or Element symbol")

    df = df[df["Element symbol"].notna()].copy()

    # ---- Total unique symbols per sample ----
    total_symbols = (
        df.groupby("Name")["Element symbol"]
        .nunique()
        .rename("total_unique_element_symbols")
    )

    out_parts = [total_symbols]

    # ---- Unique symbol counts per Type ----
    type_counts = None
    if args.include_type_counts and "Type" in df.columns:
        type_counts = (
            df.groupby(["Name", "Type"])["Element symbol"]
            .nunique()
            .unstack(fill_value=0)
        )
        type_counts.columns = [f"n_unique_symbols_type__{c}" for c in type_counts.columns.astype(str)]
        out_parts.append(type_counts)

    # ---- Unique symbol counts per Class ----
    class_counts = None
    if args.include_class_counts and "Class" in df.columns:
        class_counts = (
            df.groupby(["Name", "Class"])["Element symbol"]
            .nunique()
            .unstack(fill_value=0)
        )
        class_counts.columns = [f"n_unique_symbols_class__{c}" for c in class_counts.columns.astype(str)]
        out_parts.append(class_counts)

    # ---- Lists per Type ----
    if args.include_type_lists and "Type" in df.columns:
        type_lists = (
            df.groupby(["Name", "Type"])["Element symbol"]
            .apply(uniq_join)
            .unstack(fill_value="")
        )
        type_lists.columns = [f"element_symbols__{c}" for c in type_lists.columns.astype(str)]

        # Fill empty with "none" where corresponding count == 0
        if type_counts is not None:
            for col in type_lists.columns:
                type_name = col.replace("element_symbols__", "")
                count_col = f"n_unique_symbols_type__{type_name}"
                if count_col in type_counts.columns:
                    zero_mask = type_counts[count_col] == 0
                    type_lists.loc[zero_mask, col] = "none"
        else:
            # If counts not computed, replace empty with none for consistency
            type_lists = type_lists.replace("", "none")

        out_parts.append(type_lists)

    # ---- Lists per Class ----
    if args.include_class_lists and "Class" in df.columns:
        class_lists = (
            df.groupby(["Name", "Class"])["Element symbol"]
            .apply(uniq_join)
            .unstack(fill_value="")
        )
        class_lists.columns = [f"element_symbols_class__{c}" for c in class_lists.columns.astype(str)]

        # Fill empty with "none" where corresponding count == 0
        if class_counts is not None:
            for col in class_lists.columns:
                class_name = col.replace("element_symbols_class__", "")
                count_col = f"n_unique_symbols_class__{class_name}"
                if count_col in class_counts.columns:
                    zero_mask = class_counts[count_col] == 0
                    class_lists.loc[zero_mask, col] = "none"
        else:
            class_lists = class_lists.replace("", "none")

        out_parts.append(class_lists)

    summary = (
        pd.concat(out_parts, axis=1)
        .reset_index()
        .rename(columns={"Name": "sample"})
    )

    summary.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
