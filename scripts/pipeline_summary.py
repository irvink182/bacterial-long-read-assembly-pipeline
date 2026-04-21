#!/usr/bin/env python3

import pandas as pd
import argparse

#HINTS
def get_hint(step: str, reason: str) -> str:
    step = str(step).lower()
    reason = str(reason).lower()

    # --- ASSEMBLY ---
    if step == "assembly":
        if "low_reads" in reason:
            return "Check read depth or filtering thresholds (step_01)"
        if "flye_failed" in reason:
            return "Inspect assembly logs (Flye), possible low coverage or corrupted reads"

    # --- TAXONOMY ---
    if step == "taxonomy":
        if "missing_assembly" in reason:
            return "Assembly failed, check previous step"
        if "sourmash_failed" in reason:
            return "Check sourmash signatures/searches and/or genome quality"
        if "skani_failed" in reason:
            return "Check skani searches and/or genome quality"

    # --- ANNOTATION ---
    if step == "annotation":
        if "bakta_failed" in reason:
            return "Check Bakta logs and database availability"

    # --- CHARACTERIZATION ---
    if step == "characterization":
        if "amrfinder_failed" in reason:
            return "Check AMRFinder DB and input format"
        if "plasmidfinder_failed" in reason:
            return "Check PlasmidFinder DB and FASTA format"
        if "mobtyper_failed" in reason:
            return "Check MOB-suite input files"

    # --- GENERIC ---
    if "missing" in reason:
        return "Upstream step likely failed"

    return "Check logs for details"

#HTML REPORT
def build_html_report(df_status, summary, failed_samples, skipped_samples, output_path,
    pipeline_name, pipeline_version, author, institution):

    # -----------------------------
    # Run summary automatic
    # -----------------------------
    # unique samples
    all_sample_ids = df_status[df_status["sample"] != "GLOBAL"]["sample"].unique()
    total_samples = len(all_sample_ids)

    # Identify FAILED samples
    samples_with_failure = df_status[
        (df_status["sample"] != "GLOBAL") & 
        (df_status["status"] == "FAILED") & 
        (df_status["step"] != "cleanup")
        ]["sample"].unique()

    failed = len(samples_with_failure)

    # Identify COMPLETED samples
    completed = df_status[
    (df_status["step"] == "characterization") & 
    (df_status["status"] == "COMPLETED")
    ]["sample"].nunique()

    # Identify SKIPPED samples
    skipped = total_samples - completed - failed

    from datetime import datetime
    today = datetime.now().strftime("%Y-%m-%d")

    # -----------------------------
    # footer (author metadata)
    # -----------------------------

    footer = f"""
    <hr>
    <p style="font-size: 12px; color: #666;">
    <b>{pipeline_name}</b><br>
    Version: {pipeline_version}<br>
    Author: {author}<br>
    Institution: {institution}<br>
    Date: {today}
    </p>
    """

    # -----------------------------
    # Pipeline overview (static)
    # -----------------------------
    pipeline_overview = """
    <h2>Pipeline Overview</h2>
    <ul>
      <li><b>Step 01 – Trimming:</b> fastplong / porechop + filtlong</li>
      <li><b>Step 02 – Assembly:</b> Flye + Medaka</li>
      <li><b>Step 03 – QC:</b> QUAST + CheckM2</li>
      <li><b>Step 04 – Taxonomy:</b> Sourmash + Skani + MLST</li>
      <li><b>Step 05 – Annotation:</b> Bakta</li>
      <li><b>Step 06 – Characterization:</b> AMRFinder, PlasmidFinder, MOB-suite</li>
      <li><b>Step 07 – Final summary:</b> Merge reports</li>
      <li><b>Step 08 – Cleanup:</b> Remove intermediate and temporary files</li>
    </ul>
    """

    # -----------------------------
    # HTML
    # -----------------------------
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Pipeline Report</title>
        <style>
            body {{ font-family: Arial; margin: 40px; }}
            h1, h2 {{ color: #2c3e50; }}
            table {{ border-collapse: collapse; width: 100%; margin-bottom: 30px; }}
            th, td {{ border: 1px solid #ccc; padding: 6px; }}
            th {{ background-color: #f4f4f4; }}
            .ok {{ color: green; font-weight: bold; }}
            .fail {{ color: red; font-weight: bold; }}
            .skip {{ color: orange; font-weight: bold; }}
            .low-rate {{ background-color: #ffcccc; }}
            .high-rate {{ background-color: #ccffcc; }}
            .progress-bar {{
            background: #eee; 
            border-radius: 4px; 
            width: 100px; 
            height: 12px; 
            display: inline-block; 
            }}
        </style>
    </head>

    <body>

    <h1>Pipeline Report</h1>

    <!-- ========================= -->
    <!-- Run Summary -->
    <!-- ========================= -->
    <h2>Run Summary</h2>

    <ul>
      <li><b>Date:</b> {today}</li>
      <li><b>Total samples:</b> {total_samples}</li>
      <li class="ok">Completed: {completed}</li>
      <li class="fail">Failed: {failed}</li>
      <li class="skip">Skipped: {skipped}</li>
    </ul>

    <!-- ========================= -->
    <!-- Pipeline Overview -->
    <!-- ========================= -->
    {pipeline_overview}

    <!-- ========================= -->
    <!-- Tables -->
    <!-- ========================= -->

    <h2>Step Summary</h2>
    {summary.to_html(index=False)}

    <h2>Failed Samples</h2>
    {failed_samples.to_html(index=False) if not failed_samples.empty else "<p>None</p>"}

    <h2>Skipped Samples</h2>
    {skipped_samples.to_html(index=False) if not skipped_samples.empty else "<p>None</p>"}

    {footer}

    </body>
    </html>
    """

    with open(output_path, "w") as f:
        f.write(html)

#MAIN
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", required=True, help="pipeline_status.tsv")
    ap.add_argument("--out_tsv", default="", help="output TSV summary")
    ap.add_argument("--out_excel", default="", help="output Excel summary")
    ap.add_argument("--html", default="")
    ap.add_argument("--pipeline-name", default="Pipeline")
    ap.add_argument("--pipeline-version", default="unknown")
    ap.add_argument("--author", default="")
    ap.add_argument("--institution", default="")
    args = ap.parse_args()

    df = pd.read_csv(args.status, sep="\t")

    # ----------------------------------
    # Clean
    # ----------------------------------
    df["sample"] = df["sample"].astype(str)
    df["step"] = df["step"].astype(str)
    df["status"] = df["status"].astype(str)

    # ----------------------------------
    # Separate GLOBAL vs samples
    # ----------------------------------
    df_samples = df[df["sample"] != "GLOBAL"].copy()
    df_global = df[df["sample"] == "GLOBAL"].copy()

    # ----------------------------------
    # Keep LAST status per sample/step
    # ----------------------------------

    # 1. Primero preparamos los datos base
    df_samples = df_samples.sort_values("timestamp")
    df_samples = df_samples[df_samples["step"].notna() & (df_samples["step"] != "nan")]

    # 2. ASIGNAMOS el valor a last_status (Esto evita el UnboundLocalError)
    last_status = df_samples.groupby(["sample", "step"], as_index=False).tail(1)

    # step order
    step_order = [
    "trimming",
    "assembly",
    "assembly_qc",
    "taxonomy",
    "annotation",
    "characterization",
    "final_summary",
    "cleanup",
    "run_summary"
    ]

    all_samples = df_samples["sample"].unique()
    full_index = pd.MultiIndex.from_product(
        [all_samples, step_order],
        names=["sample", "step"]
    )

    # Reindex
    last_status = last_status.set_index(["sample", "step"]).reindex(full_index).reset_index()

    # reorder steps
    last_status["step"] = pd.Categorical(last_status["step"], categories=step_order, ordered=True)
    last_status = last_status.sort_values(["sample", "step"])

    # 
    last_status["status"] = last_status["status"].fillna("SKIPPED")
    last_status["reason"] = last_status.groupby("sample")["reason"].ffill()
    
    # Fulfill empty/NaN values
    last_status["reason"] = last_status["reason"].fillna("missing_assembly")

    # Add hints for FAILED and SKIPPED samples
    failed_mask = last_status["status"] == "FAILED"
    skipped_mask = last_status["status"] == "SKIPPED"

    # Update table
    failed_samples = last_status[failed_mask][["sample", "step", "reason"]].copy()
    skipped_samples = last_status[skipped_mask][["sample", "step", "reason"]].copy()

    # apply hints
    failed_samples["hint"] = failed_samples.apply(lambda x: get_hint(x["step"], x["reason"]), axis=1)
    skipped_samples["hint"] = skipped_samples.apply(lambda x: get_hint(x["step"], x["reason"]), axis=1)

    # ----------------------------------
    # Summary per step
    # ----------------------------------
    # Define specific functional/global steps
    steps_to_exclude = ["final_summary", "run_summary", "cleanup"]

    # Filter last_status to statistic summary
    df_for_summary = last_status[~last_status["step"].isin(steps_to_exclude)].copy()

    summary = (
        df_for_summary
        .groupby(["step", "status"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )

    for col in ["COMPLETED", "FAILED", "SKIPPED"]:
        if col not in summary.columns:
            summary[col] = 0

    summary = summary.reset_index()
    summary.columns.name = None

    summary["step"] = pd.Categorical(summary["step"], categories=step_order, ordered=True)
    summary = summary.sort_values("step")

    valid_status = ["COMPLETED", "FAILED", "SKIPPED"]
    summary = summary[["step"] + [c for c in summary.columns if c in valid_status]]


    # Add totals
    summary["total"] = summary.drop(columns=["step"]).sum(axis=1)

    # Add success rate
    if "COMPLETED" in summary.columns:
        summary["success_rate_%"] = (
            summary["COMPLETED"] / summary["total"] * 100
        ).round(2)
    else:
        summary["success_rate_%"] = 0

    # Remove GLOBAL steps from individual reports
    steps_to_exclude = ["final_summary", "run_summary", "cleanup"]
    last_status = last_status[~((last_status["sample"] != "GLOBAL") & (last_status["step"].isin(steps_to_exclude)))]

    # ----------------------------------
    # Failed samples
    # ----------------------------------
    failed = last_status[last_status["status"] == "FAILED"]
    failed_samples = failed[["sample", "step", "reason"]].copy()
    failed_samples["hint"] = failed_samples.apply(
        lambda x: get_hint(x["step"], x["reason"]),
        axis=1
    )

    # ----------------------------------
    # Skipped samples
    # ----------------------------------
    skipped = last_status[last_status["status"] == "SKIPPED"]
    skipped_samples = skipped[["sample", "step", "reason"]].copy()
    skipped_samples["hint"] = skipped_samples.apply(
        lambda x: get_hint(x["step"], x["reason"]),
        axis=1
    )

    # ----------------------------------
    # Print summary
    # ----------------------------------
    print("\n=== PIPELINE SUMMARY ===\n")
    print(summary.to_string(index=False))

    print("\n=== FAILED SAMPLES ===\n")
    if failed_samples.empty:
        print("None")
    else:
        print(failed_samples.to_string(index=False))

    print("\n=== SKIPPED SAMPLES ===\n")
    if skipped_samples.empty:
        print("None")
    else:
        print(skipped_samples.to_string(index=False))

    # ----------------------------------
    # Save TSV
    # ----------------------------------
    if args.out_tsv:
        summary.to_csv(args.out_tsv, sep="\t", index=False)

    # ----------------------------------
    # Save Excel
    # ----------------------------------
    if args.out_excel:
        # IMPORTANTE: Usamos args.out_excel para respetar la ruta de Bash
        with pd.ExcelWriter(args.out_excel, engine="xlsxwriter") as writer:
            summary.to_excel(writer, sheet_name="Summary", index=False)
            failed_samples.to_excel(writer, sheet_name="Failed", index=False)
            skipped_samples.to_excel(writer, sheet_name="Skipped", index=False)

            workbook  = writer.book
            worksheet = writer.sheets["Summary"]

            # Formato de porcentaje
            percent_fmt = workbook.add_format({'num_format': '0.00'})

            if "success_rate_%" in summary.columns:
                col_idx = summary.columns.get_loc("success_rate_%")
                # Ajustamos un poco el ancho (18) para que quepa bien el título
                worksheet.set_column(col_idx, col_idx, 18, percent_fmt)

    # ----------------------------------
    # Save HTML report
    # ----------------------------------
    if args.html:
        build_html_report(
            df,
            summary,
            failed_samples,
            skipped_samples,
            args.html,
            args.pipeline_name,
            args.pipeline_version,
            args.author,
            args.institution
        )


if __name__ == "__main__":
    main()