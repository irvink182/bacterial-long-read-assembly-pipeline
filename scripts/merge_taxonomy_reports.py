#!/usr/bin/env python3

import argparse
import csv
import glob
import math
import os
import re
from typing import Dict, List, Optional, Tuple

ACC_RE = re.compile(r'GC[AF]_\d+\.\d+')
GENUS_RE = re.compile(r'^[A-Z][a-zA-Z_-]+$')

# ----------------------------
# helpers
# ----------------------------

def clip(x, lo, hi):
    if x is None:
        return None
    return max(lo, min(hi, x))

def to_float(x):
    try:
        if x is None or x == "":
            return None
        return float(x)
    except Exception:
        return None

def to_int(x):
    try:
        if x is None or x == "":
            return None
        return int(float(x))
    except Exception:
        return None

def label_from_score(s: Optional[float], high: float, medium: float) -> str:
    if s is None:
        return "NA"
    if s >= high:
        return "HIGH"
    if s >= medium:
        return "MEDIUM"
    return "LOW"

def extract_accession(text: str) -> Optional[str]:
    if text is None:
        return None
    m = ACC_RE.search(str(text))
    return m.group(0) if m else None

def species_from_prefixed_name(text: str) -> str:
    """
    Handles:
      GCF_000742135.1 Klebsiella pneumoniae
      NZ_KN046818.1 Klebsiella pneumoniae strain ...
      QZVU01000097.1 Arthrobacter frigidicola strain ...
      Klebsiella pneumoniae
    Returns:
      Genus species
    """
    if text is None:
        return ""
    parts = str(text).strip().split()
    if len(parts) < 2:
        return str(text).strip()

    # contig/WGS-like accession followed by genus + species
    if re.match(r'^[A-Z0-9_]+\.\d+$', parts[0]) and len(parts) >= 3 and GENUS_RE.match(parts[1]):
        return " ".join(parts[1:3])

    # standard GCF/GCA/NZ prefix
    if parts[0].startswith(("GCF_", "GCA_", "NZ_")) and len(parts) >= 3:
        return " ".join(parts[1:3])

    return " ".join(parts[:2])

def genus_from_species(species: str) -> str:
    if not species:
        return ""
    parts = str(species).split()
    return parts[0] if parts else ""

def read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)

def read_csv(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        return list(reader)

def basename_sample(path: str) -> str:
    return os.path.basename(path).split(".")[0]

def glob_map(pattern: str) -> Dict[str, str]:
    d = {}
    for path in glob.glob(pattern):
        d[basename_sample(path)] = path
    return d

def read_sample_list(path: str) -> Dict[str, Dict[str, str]]:
    """
    Reads a tab-separated sample sheet with header.
    Expected at minimum:
      Sample_ID or Sample
    Optional:
      expected_genome_size
    """
    result = {}

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Empty or invalid sample list: {path}")

        fieldnames = [x.strip() for x in reader.fieldnames]

        if "Sample_ID" in fieldnames:
            sample_col = "Sample_ID"
        elif "Sample" in fieldnames:
            sample_col = "Sample"
        else:
            raise ValueError(
                f"Sample list must contain 'Sample_ID' or 'Sample' column. Found: {fieldnames}"
            )

        genome_col = "expected_genome_size" if "expected_genome_size" in fieldnames else None

        for row in reader:
            sample = (row.get(sample_col) or "").strip()
            if not sample:
                continue
            if sample in ("Sample_ID", "Sample"):
                continue

            expected = (row.get(genome_col) or "").strip() if genome_col else ""

            result[sample] = {
                "Sample": sample,
                "expected_genome_size": expected
            }

    return result

# ----------------------------
# score functions
# ----------------------------

def skani_score_comment(ani, af_ref, af_q, ani_high, af_high, af_medium, high_thr, med_thr):
    if ani is None or af_ref is None or af_q is None:
        return None, "NA", "missing skani metrics"

    af_avg = (af_ref + af_q) / 2.0

    # hard fail conditions
    if ani < ani_high:
        return 20.0, "LOW", f"ANI={ani:.2f}; ANI<{ani_high}"
    if af_ref < 10 or af_q < 10:
        return 10.0, "LOW", f"ANI={ani:.2f}; AF_ref={af_ref:.2f}; AF_query={af_q:.2f}; very_low_AF"

    if af_ref >= af_high and af_q >= af_high:
        score = min(100.0, 85.0 + max(0, ani - ani_high) * 3.0)
        return score, "HIGH", f"ANI={ani:.2f}; AF_ref={af_ref:.2f}; AF_query={af_q:.2f}; high_ani_high_af"

    if af_ref >= af_medium and af_q >= af_medium:
        score = min(84.0, 65.0 + max(0, ani - ani_high) * 2.0)
        return score, "MEDIUM", f"ANI={ani:.2f}; AF_ref={af_ref:.2f}; AF_query={af_q:.2f}; high_ani_moderate_af"

    return 35.0, "LOW", f"ANI={ani:.2f}; AF_ref={af_ref:.2f}; AF_query={af_q:.2f}; low_AF"

def sylph_score_comment(tax_ab, adj_ani, cov, tax_high, tax_med, ani_high, cov_high, cov_low, high_thr, med_thr):
    if tax_ab is None:
        return None, "NA", "missing sylph metrics"

    tax_ab = clip(tax_ab, 0, 100)
    adj_ani = clip(adj_ani, 0, 100) if adj_ani is not None else None
    cov = cov if cov is not None else None
    if cov is not None and cov < cov_low:
        comment = f"Top1={tax_ab:.1f}%"
        if adj_ani is not None:
            comment += f"; AdjANI={adj_ani:.2f}"
        comment += f"; cov={cov:.1f}; low_read_support"
        return 35.0, "LOW", comment

    if tax_ab >= tax_high and (adj_ani is not None and adj_ani >= ani_high) and (cov is not None and cov >= cov_high):
        return 95.0, "HIGH", f"Top1={tax_ab:.1f}%; AdjANI={adj_ani:.2f}; cov={cov:.1f}; single_dominant_species"

    if tax_ab >= tax_med and (adj_ani is not None and adj_ani >= ani_high):
        comment = f"Top1={tax_ab:.1f}%"
        if cov is not None:
            comment += f"; cov={cov:.1f}"
        comment += "; possible_minor_contamination"
        return 70.0, "MEDIUM", comment

    comment = f"Top1={tax_ab:.1f}%"
    if adj_ani is not None:
        comment += f"; AdjANI={adj_ani:.2f}"
    if cov is not None:
        comment += f"; cov={cov:.1f}"
    comment += "; possible_mixed_sample"
    return 35.0, "LOW", comment

def checkm2_score_comment(comp, cont, asm_type, high_thr, med_thr,
                          comp_high_long, comp_med_long, cont_high_long, cont_med_long,
                          comp_high_short, comp_med_short, cont_high_short, cont_med_short):
    if comp is None or cont is None:
        return None, "NA", "missing checkm2 metrics"

    if asm_type in ("long", "hybrid"):
        if cont >= 50:
            return 5.0, "LOW", f"comp={comp:.1f}; cont={cont:.1f}; extremely_high_contamination"
        if cont >= 20:
            return 35.0, "LOW", f"comp={comp:.1f}; cont={cont:.1f}; high_contamination"
        if cont > cont_high_long or comp < comp_high_long:
            return 70.0, "MEDIUM", f"comp={comp:.1f}; cont={cont:.1f}; moderate_quality"
        return 95.0, "HIGH", f"comp={comp:.1f}; cont={cont:.1f}; high_quality_assembly"

    # short reads
    if cont >= 50:
        return 5.0, "LOW", f"comp={comp:.1f}; cont={cont:.1f}; extremely_high_contamination"
    if cont >= 20:
        return 35.0, "LOW", f"comp={comp:.1f}; cont={cont:.1f}; high_contamination"
    if cont > cont_high_short:
        return 70.0, "MEDIUM", f"comp={comp:.1f}; cont={cont:.1f}; moderate_contamination; low_completeness_expected_for_short_reads"
    return 90.0, "HIGH", f"comp={comp:.1f}; cont={cont:.1f}; contamination_low"

def sourmash_consistency_score_comment(consistent, rank, high_thr, med_thr):
    if consistent is None:
        return None, "NA", "missing sourmash/skani consistency"
    if consistent:
        score = 100 if rank == 1 else max(70, 100 - (rank - 1) * 3)
        label = label_from_score(score, high_thr, med_thr)
        return score, label, f"skani_best_in_sourmash_top20 rank={rank}"
    score = 20
    label = label_from_score(score, high_thr, med_thr)
    return score, label, "skani_best_not_in_sourmash_top20"

def overall_score_comment(skani_score, sylph_score, checkm2_score, sourmash_score,
                          w_skani, w_sylph, w_checkm2, w_sourmash,
                          overall_high, overall_medium, require_skani_pass, skani_label,
                          genus_consensus=True, species_consensus=True,
                          sylph_mixed_flag=False, checkm2_contamination=None):
    items = []
    if skani_score is not None:
        items.append(("skani", skani_score, w_skani))
    if sylph_score is not None:
        items.append(("sylph", sylph_score, w_sylph))
    if checkm2_score is not None and w_checkm2 > 0:
        items.append(("checkm2", checkm2_score, w_checkm2))
    if sourmash_score is not None:
        items.append(("sourmash", sourmash_score, w_sourmash))

    if not items:
        return None, "NA", "no_scores_available"

    wsum = sum(w for _, _, w in items)
    score = sum(s * w for _, s, w in items) / wsum
    label = label_from_score(score, overall_high, overall_medium)

    comments = []

    if require_skani_pass and skani_label == "LOW":
        label = "LOW"
        score = min(score, overall_medium - 1)
        comments.append("skani_low")

    if checkm2_contamination is not None:
        if checkm2_contamination >= 50:
            label = "LOW"
            score = min(score, overall_medium - 1)
            comments.append("high_checkm2_contamination")
        elif checkm2_contamination >= 20 and label == "HIGH":
            label = "MEDIUM"
            score = min(score, overall_high - 1)
            comments.append("moderate_checkm2_contamination")

    if sylph_mixed_flag:
        if label == "HIGH":
            label = "MEDIUM"
            score = min(score, overall_high - 1)
        comments.append("possible_mixed_sample")

    if genus_consensus is False:
        if label == "HIGH":
            label = "MEDIUM"
            score = min(score, overall_high - 1)
        comments.append("genus_discordance")

    if species_consensus is False and sylph_mixed_flag:
        label = "LOW"
        score = min(score, overall_medium - 1)
        comments.append("species_discordance")

    if not comments:
        comments.append("all_tools_consistent")

    return score, label, "; ".join(comments)

def primary_call_status_comment(row):
    checkm2_cont = to_float(row.get("checkm2_contamination"))
    sylph_mixed = row.get("sylph_mixed_flag", False)
    overall = row.get("overall_label", "NA")
    skani_label = row.get("skani_label", "NA")
    species_consensus = row.get("species_consensus", False)
    genus_consensus = row.get("genus_consensus", False)

    # mixed first
    if sylph_mixed or (checkm2_cont is not None and checkm2_cont >= 50):
        return "mixed", "mixed_sample_detected"

    # uncertain if skani itself is weak
    if skani_label == "LOW":
        return "uncertain", "weak_skani_support"

    # confirmed
    if overall == "HIGH":
        return "confirmed", "high_confidence_species_call"

    # provisional
    if overall == "MEDIUM" and species_consensus and genus_consensus:
        return "provisional", "consistent_species_call_with_moderate_support"

    # uncertain fallback
    return "uncertain", "discordant_or_limited_support"

# ----------------------------
# parsing per tool
# ----------------------------

def parse_skani(path: str, min_af_candidate: float = 10.0) -> Dict[str, object]:
    rows = read_tsv(path)
    if not rows:
        return {}

    # parse numeric values first
    parsed = []
    for r in rows:
        ani = to_float(r.get("ANI"))
        af_ref = to_float(r.get("Align_fraction_ref"))
        af_query = to_float(r.get("Align_fraction_query"))
        parsed.append((r, ani, af_ref, af_query))

    # sort by ANI descending
    parsed.sort(key=lambda x: x[1] if x[1] is not None else -1, reverse=True)

    # keep original top1 for debug/reference
    raw_best = parsed[0][0]
    raw_best_name = raw_best.get("Ref_name", "")
    raw_best_species = species_from_prefixed_name(raw_best_name)
    raw_best_acc = extract_accession(raw_best.get("Ref_file", "")) or extract_accession(raw_best_name)

    # filter out spurious low-AF hits
    filtered = [
        (r, ani, af_ref, af_query)
        for (r, ani, af_ref, af_query) in parsed
        if af_ref is not None and af_query is not None
        and af_ref >= min_af_candidate
        and af_query >= min_af_candidate
    ]

    if filtered:
        best, ani, af_ref, af_query = filtered[0]
        selection_reason = f"selected_best_hit_after_AF_filter>={min_af_candidate}"
    else:
        # fallback to raw top1 if everything fails filter
        best, ani, af_ref, af_query = parsed[0]
        selection_reason = f"fallback_raw_top1_no_hit_passed_AF_filter>={min_af_candidate}"

    ref_file = best.get("Ref_file", "")
    ref_name = best.get("Ref_name", "")
    species = species_from_prefixed_name(ref_name)
    genus = genus_from_species(species)
    ref_acc = extract_accession(ref_file) or extract_accession(ref_name)

    return {
        "skani_species": species,
        "skani_genus": genus,
        "skani_ref_accession": ref_acc or "",
        "skani_ref_name": ref_name,
        "ANI": ani,
        "AF_ref": af_ref,
        "AF_query": af_query,
        "skani_selection_reason": selection_reason,
        "skani_raw_top1_species": raw_best_species,
        "skani_raw_top1_accession": raw_best_acc or "",
    }

def parse_sourmash(path: str, skani_ref_accession: str, topn: int) -> Dict[str, object]:
    rows = read_csv(path)
    if not rows:
        return {}
    rows.sort(key=lambda r: to_float(r.get("similarity")) or -1, reverse=True)
    top1 = rows[0]

    top1_name = top1.get("name", "")
    top1_species = species_from_prefixed_name(top1_name)
    top1_genus = genus_from_species(top1_species)

    topn_rows = rows[:topn]
    rank = None
    if skani_ref_accession:
        for i, r in enumerate(topn_rows, start=1):
            acc = extract_accession(r.get("name", "")) or r.get("name", "").split(" ")[0]
            if acc == skani_ref_accession:
                rank = i
                break

    consistent = (rank is not None)

    return {
        "sourmash_species": top1_species,
        "sourmash_genus": top1_genus,
        "sourmash_similarity": to_float(top1.get("similarity")),
        "sourmash_ani_est": to_float(top1.get("ani")),
        "skani_best_in_sourmash_top20": consistent,
        "sourmash_rank_of_skani_best": rank if rank is not None else "",
    }

def parse_sylph(path: str) -> Dict[str, object]:
    rows = read_tsv(path)
    if not rows:
        return {}

    rows.sort(key=lambda r: to_float(r.get("Taxonomic_abundance")) or -1, reverse=True)

    top1 = rows[0]
    top2 = rows[1] if len(rows) > 1 else {}

    top1_name = top1.get("Contig_name", "") or top1.get("Name", "")
    top2_name = top2.get("Contig_name", "") or top2.get("Name", "")

    top1_species = species_from_prefixed_name(top1_name)
    top2_species = species_from_prefixed_name(top2_name) if top2 else ""
    top1_genus = genus_from_species(top1_species)
    top2_genus = genus_from_species(top2_species)

    top1_tax = to_float(top1.get("Taxonomic_abundance"))
    top2_tax = to_float(top2.get("Taxonomic_abundance")) if top2 else None

    mixed_flag = False
    if top1_tax is not None and top1_tax < 85:
        mixed_flag = True
    if top2_tax is not None and top2_tax >= 10:
        mixed_flag = True

    return {
        "sylph_species": top1_species,
        "sylph_genus": top1_genus,
        "sylph_taxonomic_abundance": top1_tax,
        "sylph_sequence_abundance": to_float(top1.get("Sequence_abundance")),
        "sylph_adjusted_ani": to_float(top1.get("Adjusted_ANI")),
        "sylph_true_cov": to_float(top1.get("True_cov")),
        "sylph_top2_species": top2_species,
        "sylph_top2_genus": top2_genus,
        "sylph_top2_taxonomic_abundance": top2_tax if top2_tax is not None else "",
        "sylph_mixed_flag": mixed_flag,
    }

def parse_checkm2(path: str) -> Dict[str, object]:
    rows = read_tsv(path)
    if not rows:
        return {}
    row = rows[0]
    return {
        "checkm2_completeness": to_float(row.get("Completeness_General")),
        "checkm2_contamination": to_float(row.get("Contamination")),
        "size": to_int(row.get("Genome_Size")),
        "n_contigs": to_int(row.get("Total_Contigs")),
    }

# ----------------------------
# main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()

    # mandatory
    ap.add_argument("--asm-type", required=True, choices=["short", "long", "hybrid"])
    ap.add_argument("--sample-list", required=True)
    ap.add_argument("--skani-files", required=True)
    ap.add_argument("--sylph-files", required=True)
    ap.add_argument("--sourmash-files", required=True)
    ap.add_argument("--out", required=True)

    # optional
    ap.add_argument("--checkm2-files", default=None)
    ap.add_argument("--no-checkm2", action="store_true")

    # sourmash
    ap.add_argument("--sourmash-topn", type=int, default=20)

    # skani thresholds
    ap.add_argument("--skani-ani-high", type=float, default=95.0)
    ap.add_argument("--skani-af-high", type=float, default=50.0)
    ap.add_argument("--skani-af-medium", type=float, default=30.0)
    ap.add_argument("--skani-min-af-candidate", type=float, default=10.0)

    # sylph thresholds
    ap.add_argument("--sylph-tax-high", type=float, default=90.0)
    ap.add_argument("--sylph-tax-medium", type=float, default=70.0)
    ap.add_argument("--sylph-ani-high", type=float, default=95.0)
    ap.add_argument("--sylph-cov-high", type=float, default=10.0)
    ap.add_argument("--sylph-cov-low", type=float, default=5.0)

    # checkm2 thresholds
    ap.add_argument("--checkm2-comp-high-long", type=float, default=90.0)
    ap.add_argument("--checkm2-comp-med-long", type=float, default=80.0)
    ap.add_argument("--checkm2-cont-high-long", type=float, default=5.0)
    ap.add_argument("--checkm2-cont-med-long", type=float, default=10.0)

    ap.add_argument("--checkm2-comp-high-short", type=float, default=70.0)
    ap.add_argument("--checkm2-comp-med-short", type=float, default=50.0)
    ap.add_argument("--checkm2-cont-high-short", type=float, default=5.0)
    ap.add_argument("--checkm2-cont-med-short", type=float, default=10.0)

    # weights
    ap.add_argument("--weight-skani", type=float, default=0.60)
    ap.add_argument("--weight-sylph", type=float, default=0.25)
    ap.add_argument("--weight-checkm2", type=float, default=0.15)
    ap.add_argument("--weight-sourmash", type=float, default=0.05)

    # labels
    ap.add_argument("--score-high", type=float, default=85.0)
    ap.add_argument("--score-medium", type=float, default=60.0)
    ap.add_argument("--overall-high", type=float, default=85.0)
    ap.add_argument("--overall-medium", type=float, default=60.0)

    # behavior
    ap.add_argument("--require-skani-pass", action="store_true")
    ap.add_argument("--write-comments", action="store_true")

    args = ap.parse_args()

    if args.asm_type == "short":
        args.no_checkm2 = True

    # reweight if no checkm2
    if args.no_checkm2:
        args.weight_skani = 0.65
        args.weight_sylph = 0.30
        args.weight_sourmash = 0.05
        args.weight_checkm2 = 0.0

    sample_meta = read_sample_list(args.sample_list)
    skani_map = glob_map(args.skani_files)
    sourmash_map = glob_map(args.sourmash_files)
    sylph_map = glob_map(args.sylph_files)
    checkm2_map = {} if args.no_checkm2 or not args.checkm2_files else glob_map(args.checkm2_files)

    rows_out = []

    for sample in sorted(sample_meta.keys()):
        row = {
            "Sample": sample,
            "asm_type": args.asm_type,
            "expected_genome_size": sample_meta[sample].get("expected_genome_size", ""),
            "size": "",
            "n_contigs": "",
            "skani_species": "",
            "skani_genus": "",
            "skani_ref_accession": "",
            "skani_ref_name": "",
            "skani_selection_reason": "",
            "skani_raw_top1_species": "",
            "skani_raw_top1_accession": "",
            "sourmash_species": "",
            "sourmash_genus": "",
            "sylph_species": "",
            "sylph_genus": "",
            "sylph_top2_species": "",
            "sylph_top2_genus": "",
            "sylph_top2_taxonomic_abundance": "",
            "sylph_mixed_flag": "",
            "ANI": "",
            "AF_ref": "",
            "AF_query": "",
            "sourmash_similarity": "",
            "sourmash_ani_est": "",
            "sylph_taxonomic_abundance": "",
            "sylph_sequence_abundance": "",
            "sylph_adjusted_ani": "",
            "sylph_true_cov": "",
            "checkm2_completeness": "",
            "checkm2_contamination": "",
            "primary_species_call": "",
            "primary_call_status": "",
            "primary_call_comment": "",
        }

        # skani
        skani_data = parse_skani(
        skani_map[sample],
        min_af_candidate=args.skani_min_af_candidate
        ) if sample in skani_map else {}
        row.update(skani_data)

        # sourmash
        sourmash_data = parse_sourmash(
            sourmash_map[sample],
            row.get("skani_ref_accession", ""),
            args.sourmash_topn
        ) if sample in sourmash_map else {}
        row.update(sourmash_data)

        # sylph
        sylph_data = parse_sylph(sylph_map[sample]) if sample in sylph_map else {}
        row.update(sylph_data)

        # checkm2
        if not args.no_checkm2 and sample in checkm2_map:
            checkm2_data = parse_checkm2(checkm2_map[sample])
            row.update(checkm2_data)

        # consensus
        genus_set = {g for g in [row.get("skani_genus"), row.get("sourmash_genus"), row.get("sylph_genus")] if g}
        species_set = {s for s in [row.get("skani_species"), row.get("sourmash_species"), row.get("sylph_species")] if s}
        row["genus_consensus"] = (len(genus_set) == 1 and len(genus_set) > 0)
        row["species_consensus"] = (len(species_set) == 1 and len(species_set) > 0)

        # skani score
        sk_score, sk_label, sk_comment = skani_score_comment(
            to_float(row.get("ANI")),
            to_float(row.get("AF_ref")),
            to_float(row.get("AF_query")),
            args.skani_ani_high,
            args.skani_af_high,
            args.skani_af_medium,
            args.score_high,
            args.score_medium,
        )
        row["skani_score"] = round(sk_score, 2) if sk_score is not None else ""
        row["skani_label"] = sk_label
        row["skani_comment"] = sk_comment if args.write_comments else ""

        # sylph score
        sy_score, sy_label, sy_comment = sylph_score_comment(
            to_float(row.get("sylph_taxonomic_abundance")),
            to_float(row.get("sylph_adjusted_ani")),
            to_float(row.get("sylph_true_cov")),
            args.sylph_tax_high,
            args.sylph_tax_medium,
            args.sylph_ani_high,
            args.sylph_cov_high,
            args.sylph_cov_low,
            args.score_high,
            args.score_medium,
        )
        row["sylph_score"] = round(sy_score, 2) if sy_score is not None else ""
        row["sylph_label"] = sy_label
        row["sylph_comment"] = sy_comment if args.write_comments else ""

        # checkm2 score
        cm_score, cm_label, cm_comment = checkm2_score_comment(
            to_float(row.get("checkm2_completeness")),
            to_float(row.get("checkm2_contamination")),
            args.asm_type,
            args.score_high,
            args.score_medium,
            args.checkm2_comp_high_long,
            args.checkm2_comp_med_long,
            args.checkm2_cont_high_long,
            args.checkm2_cont_med_long,
            args.checkm2_comp_high_short,
            args.checkm2_comp_med_short,
            args.checkm2_cont_high_short,
            args.checkm2_cont_med_short,
        ) if not args.no_checkm2 else (None, "NA", "disabled_for_short_reads")
        row["checkm2_score"] = round(cm_score, 2) if cm_score is not None else ""
        row["checkm2_label"] = cm_label
        row["checkm2_comment"] = cm_comment if args.write_comments else ""

        # sourmash consistency score
        sm_score, sm_label, sm_comment = sourmash_consistency_score_comment(
            row.get("skani_best_in_sourmash_top20"),
            to_int(row.get("sourmash_rank_of_skani_best")),
            args.score_high,
            args.score_medium,
        )
        row["sourmash_score"] = round(sm_score, 2) if sm_score is not None else ""
        row["sourmash_label"] = sm_label
        row["sourmash_comment"] = sm_comment if args.write_comments else ""

        # overall
        ov_score, ov_label, ov_comment = overall_score_comment(
            sk_score, sy_score, cm_score, sm_score,
            args.weight_skani, args.weight_sylph, args.weight_checkm2, args.weight_sourmash,
            args.overall_high, args.overall_medium,
            args.require_skani_pass, sk_label,
            row.get("genus_consensus"),
            row.get("species_consensus"),
            row.get("sylph_mixed_flag", False),
            to_float(row.get("checkm2_contamination"))
        )
        row["overall_score"] = round(ov_score, 2) if ov_score is not None else ""
        row["overall_label"] = ov_label
        row["overall_comment"] = ov_comment if args.write_comments else ""

        # Primary call_status
        status, status_comment= primary_call_status_comment(row)
        row["primary_call_status"] = status
        row["primary_call_comment"] = status_comment

        if status == "mixed":
            row["primary_species_call"] = "mixed"
        elif status == "uncertain":
            row["primary_species_call"] = "uncertain"
        else:
            row["primary_species_call"] = row.get("skani_species", "")

        rows_out.append(row)

    fieldnames = [
        "Sample", "asm_type", "expected_genome_size", "size", "n_contigs",
        "primary_species_call", "primary_call_status", "primary_call_comment",
        "skani_species", "skani_genus", "skani_ref_accession","skani_ref_name",
        "skani_selection_reason", "skani_raw_top1_species", "skani_raw_top1_accession",
        "sourmash_species", "sourmash_genus",
        "sylph_species", "sylph_genus", "sylph_top2_species", "sylph_top2_genus",
        "sylph_top2_taxonomic_abundance", "sylph_mixed_flag",
        "ANI", "AF_ref", "AF_query",
        "sourmash_similarity", "sourmash_ani_est",
        "skani_best_in_sourmash_top20", "sourmash_rank_of_skani_best",
        "sylph_taxonomic_abundance", "sylph_sequence_abundance", "sylph_adjusted_ani", "sylph_true_cov",
        "checkm2_completeness", "checkm2_contamination",
        "skani_score", "skani_label", "skani_comment",
        "sylph_score", "sylph_label", "sylph_comment",
        "checkm2_score", "checkm2_label", "checkm2_comment",
        "sourmash_score", "sourmash_label", "sourmash_comment",
        "genus_consensus", "species_consensus",
        "overall_score", "overall_label", "overall_comment",
    ]

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Wrote {args.out} with {len(rows_out)} samples")

if __name__ == "__main__":
    main()