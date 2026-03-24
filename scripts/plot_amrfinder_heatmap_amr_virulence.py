#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

# Make SVG text editable in vector editors (Inkscape/Illustrator)
matplotlib.rcParams["svg.fonttype"] = "none"


def read_matrix(path: str) -> pd.DataFrame:
    """
    Reads a presence/absence matrix TSV with a 'sample' column and gene columns (0/1).
    Returns a DataFrame indexed by sample.
    """
    df = pd.read_csv(path, sep="\t")
    if "sample" not in df.columns:
        raise ValueError(f"Expected a 'sample' column in {path}")
    df = df.set_index("sample")

    # Force numeric (0/1)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    return df


def select_and_order(mat: pd.DataFrame, top_n: int, order: str) -> pd.DataFrame:
    """
    Selects top_n genes and orders them either by prevalence (desc) or alphabetically.
    If top_n == 0 -> keep all genes.
    """
    if mat.shape[1] == 0:
        return mat

    if top_n == 0 or top_n > mat.shape[1]:
        top_n = mat.shape[1]

    if order == "prevalence":
        prev = mat.sum(axis=0).sort_values(ascending=False)
        cols = prev.head(top_n).index.tolist()
    elif order == "alphabetical":
        cols = sorted(mat.columns)[:top_n]
    else:
        raise ValueError("order must be 'alphabetical' or 'prevalence'")

    return mat.loc[:, cols]


def main():
    ap = argparse.ArgumentParser(
        description="Plot exploratory combined heatmap (presence/absence) for AMR + VIRULENCE matrices."
    )
    ap.add_argument("--amr", required=True, help="AMR presence/absence TSV (sample + gene columns 0/1).")
    ap.add_argument("--virulence", required=True, help="VIRULENCE presence/absence TSV (sample + gene columns 0/1).")
    ap.add_argument("-o", "--output", required=True,
                    help="Output basename or filename. Will write <base>.png and <base>.svg.")

    ap.add_argument("--top-amr", type=int, default=20,
                    help="Number of AMR genes to display (default: 20; use 0 to include all).")
    ap.add_argument("--top-vir", type=int, default=20,
                    help="Number of VIRULENCE genes to display (default: 20; use 0 to include all).")

    ap.add_argument("--order", default="alphabetical",
                    choices=["alphabetical", "prevalence"],
                    help="Gene ordering method within each block (default: alphabetical).")

    ap.add_argument("--figwidth", type=float, default=14, help="Figure width (inches).")
    ap.add_argument("--figheight", type=float, default=8, help="Figure height (inches).")
    ap.add_argument("--title", default=None,
                    help="Optional plot title (default: no title).")
    ap.add_argument("--no-sample-labels", action="store_true",
                    help="Hide sample labels on Y axis (useful for many samples).")

    ap.add_argument("--png-dpi", type=int, default=300, help="DPI for PNG output (default: 300).")

    args = ap.parse_args()

    # Determine base name (strip extension if provided)
    base, ext = os.path.splitext(args.output)
    out_base = base if ext else args.output

    # Load matrices
    amr = read_matrix(args.amr)
    vir = read_matrix(args.virulence)

    # Keep shared samples only
    shared_samples = sorted(set(amr.index) & set(vir.index))
    if len(shared_samples) == 0:
        raise ValueError("No shared samples between AMR and VIRULENCE matrices.")
    amr = amr.loc[shared_samples]
    vir = vir.loc[shared_samples]

    # Select + order genes
    amr_top = select_and_order(amr, args.top_amr, args.order)
    vir_top = select_and_order(vir, args.top_vir, args.order)

    # Combine side-by-side
    combined = pd.concat([amr_top, vir_top], axis=1)

    # ---- Plot ----
    plt.figure(figsize=(args.figwidth, args.figheight))
    plt.imshow(combined, aspect="auto")
    plt.colorbar(label="Presence (1) / Absence (0)")

    # Y labels (samples)
    if not args.no_sample_labels:
        plt.yticks(range(len(combined.index)), combined.index)
    else:
        plt.yticks([])

    # X labels (genes)
    plt.xticks(range(len(combined.columns)), combined.columns, rotation=90)

    # Separator + section labels
    if len(amr_top.columns) > 0 and len(vir_top.columns) > 0:
        sep_x = len(amr_top.columns) - 0.5
        plt.axvline(x=sep_x)
        plt.text(len(amr_top.columns) / 2, -2, "AMR", ha="center")
        plt.text(len(amr_top.columns) + (len(vir_top.columns) / 2), -2, "VIRULENCE", ha="center")

    if args.title:
        plt.title(args.title)
    plt.tight_layout()

    # ---- Save both PNG and SVG ----
    png_path = f"{out_base}.png"
    svg_path = f"{out_base}.svg"

    plt.savefig(png_path, dpi=args.png_dpi)
    plt.savefig(svg_path, format="svg")

    print("Wrote:", png_path)
    print("Wrote:", svg_path)


if __name__ == "__main__":
    main()
