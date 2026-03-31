#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd


def read_tsv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def normalize_sample_col(df: pd.DataFrame) -> pd.DataFrame:
    for col in ["Sample", "sample", "Sample_ID", "sample_id"]:
        if col in df.columns:
            return df.rename(columns={col: "Sample"})
    raise ValueError("No sample column found. Expected one of: Sample, sample, Sample_ID, sample_id")


def normalize_sample_name(sample: str) -> str:
    if not isinstance(sample, str):
        return ""
    sample = sample.strip()

    # sufijos frecuentes de lecturas
    for suffix in [".fastplong", ".porechop_filtlong", ".clean", ".fastq.gz", ".fastq", ".fq.gz", ".fq"]:
        if sample.endswith(suffix):
            sample = sample[: -len(suffix)]

    # casos encadenados
    for token in [".fastplong.clean", ".porechop_filtlong.clean"]:
        if token in sample:
            sample = sample.split(token)[0]

    # sufijos tipo contig
    sample = sample.replace(".fa", "").replace(".fasta", "")
    sample = sample.replace(".assembly", "")
    sample = sample.replace(".contigs", "")
    sample = sample.strip()

    return sample


def standardize_sample_column(df: pd.DataFrame) -> pd.DataFrame:
    df = normalize_sample_col(df)
    df["Sample"] = df["Sample"].astype(str).map(normalize_sample_name)
    df["Sample"] = df["Sample"].str.replace(r"_contig_\d+$", "", regex=True)
    df["Sample"] = df["Sample"].str.strip()
    return df


def first_non_empty(series: pd.Series) -> str:
    for x in series:
        s = str(x).strip()
        if s and s.lower() != "nan":
            return s
    return ""


def split_unique_join(values: Iterable[str], sep: str = ";") -> str:
    items: list[str] = []
    seen = set()
    for v in values:
        if not isinstance(v, str):
            continue
        v = v.strip()
        if not v or v == "-" or v.lower() == "nan":
            continue
        for part in v.split(sep):
            p = part.strip()
            if not p or p == "-" or p.lower() == "nan":
                continue
            if p not in seen:
                seen.add(p)
                items.append(p)
    return ";".join(items)


def count_semicolon_items(value: str) -> int:
    if not isinstance(value, str):
        return 0
    value = value.strip()
    if not value or value == "-" or value.lower() == "nan":
        return 0
    return len([x for x in value.split(";") if x.strip()])


def parse_samples(samples_path: str | Path) -> pd.DataFrame:
    df = read_tsv(samples_path)
    df = standardize_sample_column(df)

    rename_map = {}
    if "asm_type" not in df.columns and "ASM_TYPE" in df.columns:
        rename_map["ASM_TYPE"] = "asm_type"
    if "expected_genome_size" not in df.columns and "genome_size" in df.columns:
        rename_map["genome_size"] = "expected_genome_size"
    if rename_map:
        df = df.rename(columns=rename_map)

    wanted = ["Sample"]
    for c in ["asm_type", "expected_genome_size"]:
        if c in df.columns:
            wanted.append(c)

    out = df[wanted].copy()
    out = out.groupby("Sample", as_index=False).agg(first_non_empty)
    return out


def parse_reads_qc(reads_qc_dir: str | Path, samples: list[str]) -> pd.DataFrame:
    rows = []

    for sample in samples:
        f = Path(reads_qc_dir) / f"{sample}.clean.txt"

        data = {
            "Sample": sample,
            "clean_read_count": "",
            "clean_total_bases": "",
            "clean_mean_length": "",
            "clean_median_length": "",
            "clean_read_N50": "",
            "clean_mean_qscore": "",
        }

        if f.exists():
            with open(f) as fh:
                for line in fh:
                    line = line.strip()

                    if line.startswith("Number of reads:"):
                        data["clean_read_count"] = line.split(":")[1].strip()

                    elif line.startswith("Total bases:"):
                        data["clean_total_bases"] = line.split(":")[1].strip()

                    elif line.startswith("Mean read length:"):
                        data["clean_mean_length"] = line.split(":")[1].strip()

                    elif line.startswith("Median read length:"):
                        data["clean_median_length"] = line.split(":")[1].strip()

                    elif line.startswith("Read length N50:"):
                        data["clean_read_N50"] = line.split(":")[1].strip()

                    elif line.startswith("Mean read quality:"):
                        data["clean_mean_qscore"] = line.split(":")[1].strip()

        rows.append(data)

    return pd.DataFrame(rows)


def parse_quast(quast_dir: str | Path, samples: list[str]) -> pd.DataFrame:
    rows = []

    for sample in samples:
        report = Path(quast_dir) / sample / "transposed_report.tsv"

        if report.exists():
            df = pd.read_csv(report, sep="\t")

            # QUAST transposed: una fila con métricas como columnas
            row = df.iloc[0].to_dict()

            row["Sample"] = normalize_sample_name(sample)

            rows.append(row)

        else:
            # muestra sin QUAST
            rows.append({
                "Sample": sample,
                "assembly_size": "",
                "contigs": "",
                "N50": "",
                "L50": "",
                "GC_percent": "",
            })

    df = pd.DataFrame(rows)

    # Renombrar columnas si existen
    rename_map = {
        "# contigs": "contigs",
        "Total length": "assembly_size",
        "Total length (>= 0 bp)": "assembly_size",
        "GC (%)": "GC_percent",
        "GC%": "GC_percent",
    }

    df = df.rename(columns=rename_map)

    # Asegurar columnas
    wanted = ["Sample", "assembly_size", "contigs", "N50", "L50", "GC_percent"]
    for c in wanted:
        if c not in df.columns:
            df[c] = ""

    return df[wanted]

def parse_checkm2(checkm2_path: str | Path) -> pd.DataFrame:
    df = read_tsv(checkm2_path)
    df = standardize_sample_column(df)

    rename_map = {}
    for src, dst in {
        "Completeness": "checkm2_completeness",
        "Completeness_General": "checkm2_completeness",
        "completeness": "checkm2_completeness",
        "Contamination": "checkm2_contamination",
        "contamination": "checkm2_contamination",
    }.items():
        if src in df.columns:
            rename_map[src] = dst

    df = df.rename(columns=rename_map)

    for col in ["checkm2_completeness", "checkm2_contamination"]:
        if col not in df.columns:
            df[col] = ""

    return df[["Sample", "checkm2_completeness", "checkm2_contamination"]].groupby("Sample", as_index=False).agg(first_non_empty)


def parse_taxonomy(taxonomy_path: str | Path) -> pd.DataFrame:
    df = read_tsv(taxonomy_path)
    df = standardize_sample_column(df)

    # elimina posibles filas repetidas de header incrustadas
    df = df[df["Sample"].astype(str).str.strip().ne("Sample")].copy()

    preferred = {
        "final_species": [
            "final_species",
            "primary_species_call",
            "selected_species",
            "taxonomy_species",
            "consensus_species",
        ],
        "final_genus": [
            "final_genus",
            "skani_genus",
            "selected_genus",
            "taxonomy_genus",
            "consensus_genus",
            "sourmash_genus",
        ],
        "ANI": [
            "ANI",
            "ani",
            "skani_ani",
            "sourmash_ani_est",
            "sylph_adjusted_ani",
        ],
        "confidence_overall": [
            "confidence_overall",
            "Confidence_overall",
            "overall_confidence",
            "overall_label",
            "overall_score",
        ],
        "skani_species": ["skani_species"],
        "sourmash_species": ["sourmash_species"],
    }

    out = pd.DataFrame({"Sample": df["Sample"]})
    for dst, candidates in preferred.items():
        src = next((c for c in candidates if c in df.columns), None)
        out[dst] = df[src] if src else ""

    return out.groupby("Sample", as_index=False).agg(first_non_empty)


def parse_mlst(mlst_path: str | Path) -> pd.DataFrame:
    rows = []

    with open(mlst_path) as f:
        header = f.readline().rstrip("\n").split("\t")

        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue

            parts = line.split("\t")

            # Need at least Sample / Scheme / ST
            if len(parts) < 3:
                continue

            sample = parts[0].strip()
            scheme = parts[1].strip()
            st = parts[2].strip()

            if not sample or sample in ("Sample", "sample", "Sample_ID"):
                continue

            rows.append({
                "Sample": normalize_sample_name(sample),
                "mlst_scheme": scheme,
                "mlst_ST": st
            })

    df = pd.DataFrame(rows)

    if df.empty:
        return pd.DataFrame(columns=["Sample", "mlst_scheme", "mlst_ST"])

    return df.groupby("Sample", as_index=False).agg(first_non_empty)


def parse_bakta_stats(bakta_dir: str | Path, samples: list[str]) -> pd.DataFrame:
    rows = []

    for sample in samples:
        f = Path(bakta_dir) / f"{sample}.txt"

        data = {
            "Sample": sample,
            "bakta_cds": "",
            "bakta_tRNA": "",
            "bakta_rRNA": "",
            "bakta_pseudogenes": "",
        }

        if f.exists():
            with open(f) as fh:
                for line in fh:
                    line = line.strip()

                    if line.startswith("CDSs:"):
                        data["bakta_cds"] = line.split(":")[1].strip()

                    elif line.startswith("tRNAs:"):
                        data["bakta_tRNA"] = line.split(":")[1].strip()

                    elif line.startswith("rRNAs:"):
                        data["bakta_rRNA"] = line.split(":")[1].strip()

                    elif line.startswith("pseudogenes:"):
                        data["bakta_pseudogenes"] = line.split(":")[1].strip()

        rows.append(data)

    return pd.DataFrame(rows)


def parse_amrfinder(amr_combined_path: str | Path) -> pd.DataFrame:
    df = read_tsv(amr_combined_path)

    sample_col = next((c for c in ["Sample", "sample", "name", "Name"] if c in df.columns), None)
    if sample_col is None:
        raise ValueError("Could not find sample column in AMRFinder combined table")

    df = df.rename(columns={sample_col: "Sample"})
    df["Sample"] = df["Sample"].astype(str).map(normalize_sample_name)

    type_col = next((c for c in ["Element type", "Type", "type"] if c in df.columns), None)
    if type_col is None:
        out = (
            df[["Sample"]]
            .drop_duplicates()
            .assign(amr_gene_count="", virulence_gene_count="", stress_gene_count="")
        )
        return out

    df = df.rename(columns={type_col: "ElementType"})
    t = df["ElementType"].astype(str).str.upper()

    out = pd.DataFrame({"Sample": sorted(df["Sample"].dropna().unique())})
    out = out.merge(df.loc[t == "AMR"].groupby("Sample").size().rename("amr_gene_count"), on="Sample", how="left")
    out = out.merge(df.loc[t == "VIRULENCE"].groupby("Sample").size().rename("virulence_gene_count"), on="Sample", how="left")
    out = out.merge(df.loc[t == "STRESS"].groupby("Sample").size().rename("stress_gene_count"), on="Sample", how="left")
    out = out.fillna(0)

    for c in ["amr_gene_count", "virulence_gene_count", "stress_gene_count"]:
        out[c] = out[c].astype(int)

    return out


def parse_plasmidfinder(plasmidfinder_summary_path: str | Path) -> pd.DataFrame:
    df = read_tsv(plasmidfinder_summary_path)
    df = standardize_sample_column(df)

    inc_col = "inc_types" if "inc_types" in df.columns else None
    if inc_col is None:
        possible = [c for c in df.columns if "inc" in c.lower()]
        inc_col = possible[0] if possible else None

    out = pd.DataFrame({"Sample": df["Sample"]})
    out["plasmidfinder_replicons"] = df[inc_col] if inc_col else ""
    out["plasmidfinder_hits"] = out["plasmidfinder_replicons"].apply(count_semicolon_items)
    out["plasmidfinder_has_plasmid"] = out["plasmidfinder_hits"].gt(0)

    return out.groupby("Sample", as_index=False).agg({
        "plasmidfinder_replicons": split_unique_join,
        "plasmidfinder_hits": "max",
        "plasmidfinder_has_plasmid": "max",
    })


def parse_mobtyper(mobtyper_combined_path: str | Path) -> pd.DataFrame:
    df = read_tsv(mobtyper_combined_path)

    sample_col = next((c for c in ["sample_id", "Sample", "sample"] if c in df.columns), None)
    if sample_col is None:
        raise ValueError("Could not find sample column in MOB-typer combined table")

    df = df.rename(columns={sample_col: "sample_id_raw"})
    df["Sample"] = df["sample_id_raw"].astype(str).map(normalize_sample_name)
    df["Sample"] = df["Sample"].str.replace(r"_contig_\d+$", "", regex=True)

    def get_col(name: str) -> pd.Series:
        return df[name] if name in df.columns else pd.Series([""] * len(df), index=df.index)

    df["rep_types"] = get_col("rep_type(s)")
    df["predicted_mobility"] = get_col("predicted_mobility")
    df["relaxase_types"] = get_col("relaxase_type(s)")
    df["mpf_types"] = get_col("mpf_type")
    if "mpf_type(s)" in df.columns:
        df["mpf_types"] = df["mpf_type(s)"]

    positive = df["rep_types"].astype(str).str.strip().replace({"": "-", "nan": "-"}).ne("-")
    df["mob_positive"] = positive

    rows = []
    for sample, sub in df.groupby("Sample", dropna=False):
        rows.append(
            {
                "Sample": sample,
                "mobtyper_positive_contigs": int(sub["mob_positive"].sum()),
                "mobtyper_rep_types": split_unique_join(sub["rep_types"]),
                "mobtyper_predicted_mobility": split_unique_join(sub["predicted_mobility"]),
                "mobtyper_relaxase_types": split_unique_join(sub["relaxase_types"]),
                "mobtyper_mpf_types": split_unique_join(sub["mpf_types"]),
                "mobtyper_has_plasmid": bool(sub["mob_positive"].any()),
            }
        )

    return pd.DataFrame(rows)


def add_estimated_coverage(df: pd.DataFrame) -> pd.DataFrame:
    def safe_float(x):
        try:
            return float(str(x).replace(",", ""))
        except Exception:
            return None

    coverage = []
    for _, row in df.iterrows():
        bases = safe_float(row.get("clean_total_bases", ""))
        gsize = safe_float(row.get("expected_genome_size", ""))
        if bases is None or gsize is None or gsize == 0:
            coverage.append("")
        else:
            coverage.append(round(bases / gsize, 2))
    df["estimated_coverage"] = coverage
    return df

def add_input_presence_flag(df: pd.DataFrame, col: str, flag_name: str) -> pd.DataFrame:
    if col not in df.columns:
        df[flag_name] = False
        return df

    series = df[col]

    # 🔥 FIX: si es DataFrame (columnas duplicadas), coger la primera
    if isinstance(series, pd.DataFrame):
        series = series.iloc[:, 0]

    df[flag_name] = (
        series.astype(str)
        .str.strip()
        .replace({"": pd.NA, "NA": pd.NA, "None": pd.NA, "nan": pd.NA})
        .notna()
    )

    return df

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--reads-qc", required=True)
    ap.add_argument("--quast", required=True)
    ap.add_argument("--checkm2", required=True)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--mlst", required=True)
    ap.add_argument("--bakta", required=True)
    ap.add_argument("--amr", required=True)
    ap.add_argument("--plasmidfinder", required=True)
    ap.add_argument("--mobtyper", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-xlsx", default="")
    args = ap.parse_args()

    samples = parse_samples(args.samples)
    sample_list = samples["Sample"].tolist()
    reads_qc = parse_reads_qc(args.reads_qc, sample_list)
    quast = parse_quast(args.quast, sample_list)
    checkm2 = parse_checkm2(args.checkm2)
    taxonomy = parse_taxonomy(args.taxonomy)
    mlst = parse_mlst(args.mlst)
    bakta = parse_bakta_stats(args.bakta, sample_list)
    amr = parse_amrfinder(args.amr)
    plasmidfinder = parse_plasmidfinder(args.plasmidfinder)
    mobtyper = parse_mobtyper(args.mobtyper)

    df = samples.copy()
    for part in [reads_qc, quast, checkm2, taxonomy, mlst, bakta, amr, plasmidfinder, mobtyper]:
        df = df.merge(part, on="Sample", how="left")
        df = df.loc[:, ~df.columns.duplicated()]
        df = add_input_presence_flag(df, "assembly_size", "has_quast")
        df = add_input_presence_flag(df, "checkm2_completeness", "has_checkm2")
        df = add_input_presence_flag(df, "final_species", "has_taxonomy")
        df = add_input_presence_flag(df, "mlst_scheme", "has_mlst")
        df = add_input_presence_flag(df, "bakta_cds", "has_bakta")
        df = add_input_presence_flag(df, "amr_gene_count", "has_amr")
        df = add_input_presence_flag(df, "plasmidfinder_hits", "has_plasmidfinder")
        df = add_input_presence_flag(df, "mobtyper_positive_contigs", "has_mobtyper")

    df = add_estimated_coverage(df)

    final_cols = [
        "Sample",
        "asm_type",
        "expected_genome_size",
        "clean_read_count",
        "clean_total_bases",
        "clean_mean_length",
        "clean_median_length",
        "clean_read_N50",
        "clean_mean_qscore",
        "estimated_coverage",
        "assembly_size",
        "contigs",
        "N50",
        "L50",
        "GC_percent",
        "checkm2_completeness",
        "checkm2_contamination",
        "final_species",
        "final_genus",
        "ANI",
        "confidence_overall",
        "skani_species",
        "sourmash_species",
        "mlst_scheme",
        "mlst_ST",
        "bakta_cds",
        "bakta_tRNA",
        "bakta_rRNA",
        "bakta_pseudogenes",
        "amr_gene_count",
        "virulence_gene_count",
        "stress_gene_count",
        "plasmidfinder_hits",
        "plasmidfinder_replicons",
        "plasmidfinder_has_plasmid",
        "mobtyper_positive_contigs",
        "mobtyper_rep_types",
        "mobtyper_predicted_mobility",
        "mobtyper_relaxase_types",
        "mobtyper_mpf_types",
        "mobtyper_has_plasmid",
        "has_quast",
        "has_checkm2",
        "has_taxonomy",
        "has_mlst",
        "has_bakta",
        "has_amr",
        "has_plasmidfinder",
        "has_mobtyper",
    ]

    for c in final_cols:
        if c not in df.columns:
            df[c] = ""

    df = df[final_cols].copy()

    # reemplazar celdad vacias por NA
    df = df.replace("", pd.NA)
    df = df.fillna("NA")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)

    if args.out_xlsx:
        xlsx_path = Path(args.out_xlsx)
        xlsx_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_excel(xlsx_path, index=False)


if __name__ == "__main__":
    main()