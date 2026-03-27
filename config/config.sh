#!/usr/bin/env bash

####################################
# GENERAL
####################################

THREADS="${THREADS:-16}"
LOG_DIR="${LOG_DIR:-logs}"

####################################
# INPUT
####################################

RAW_READS_DIR="${RAW_READS_DIR:-data/raw_reads}"
SAMPLE_LIST_DEFAULT="${SAMPLE_LIST_DEFAULT:-samplesheet.tsv}"

####################################
# MAIN RESULTS DIRECTORIES
####################################

RESULTS_DIR="${RESULTS_DIR:-results}"

READS_RESULTS_DIR="${READS_RESULTS_DIR:-${RESULTS_DIR}/reads}"
ASSEMBLY_RESULTS_DIR="${ASSEMBLY_RESULTS_DIR:-${RESULTS_DIR}/assemblies}"
QC_RESULTS_DIR="${QC_RESULTS_DIR:-${RESULTS_DIR}/quality_check}"
TAXONOMY_RESULTS_DIR="${TAXONOMY_RESULTS_DIR:-${RESULTS_DIR}/taxonomy}"
ANNOTATION_RESULTS_DIR="${ANNOTATION_RESULTS_DIR:-${RESULTS_DIR}/annotation}"
CHARACTERIZATION_RESULTS_DIR="${CHARACTERIZATION_RESULTS_DIR:-${RESULTS_DIR}/characterization}"

####################################
# STEP 01 - READS OUTPUT DIRECTORIES
####################################

PORECHOP_READS_DIR="${PORECHOP_READS_DIR:-${READS_RESULTS_DIR}/porechop_reads}"
CLEAN_READS_DIR="${CLEAN_READS_DIR:-${READS_RESULTS_DIR}/clean_reads}"

####################################
# STEP 01 - READS QC DIRECTORIES
####################################

READS_QC_DIR="${READS_QC_DIR:-${READS_RESULTS_DIR}/qc}"

FASTPLONG_QC_DIR="${FASTPLONG_QC_DIR:-${READS_QC_DIR}/fastplong}"
PORECHOP_QC_DIR="${PORECHOP_QC_DIR:-${READS_QC_DIR}/porechop}"
FILTLONG_QC_DIR="${FILTLONG_QC_DIR:-${READS_QC_DIR}/filtlong}"

NANOSTAT_RAW_DIR="${NANOSTAT_RAW_DIR:-${READS_QC_DIR}/nanostat_raw}"
NANOSTAT_CLEAN_DIR="${NANOSTAT_CLEAN_DIR:-${READS_QC_DIR}/nanostat_clean}"
QC_RAW_DIR="${NANOSTAT_RAW_DIR}"
QC_CLEAN_DIR="${NANOSTAT_CLEAN_DIR}"


MULTIQC_DIR="${MULTIQC_DIR:-${READS_QC_DIR}/multiqc}"
MULTIQC_RAW_DIR="${MULTIQC_RAW_DIR:-${MULTIQC_DIR}/raw}"
MULTIQC_CLEAN_DIR="${MULTIQC_CLEAN_DIR:-${MULTIQC_DIR}/clean}"

####################################
# STEP 01 - SYLPH OUTPUT DIRECTORIES
####################################

SYLPH_DIR="${SYLPH_DIR:-${READS_RESULTS_DIR}/sylph}"
SYLPH_SKETCH_DIR="${SYLPH_SKETCH_DIR:-${SYLPH_DIR}/sylph_sketches}"
SYLPH_TABLE_DIR="${SYLPH_TABLE_DIR:-${SYLPH_DIR}/sylph_tables}"

####################################
# STEP 02 - ASSEMBLY OUTPUT DIRECTORIES
####################################

FLYE_DIR="${FLYE_DIR:-${ASSEMBLY_RESULTS_DIR}/flye_asm}"
MEDAKA_DIR="${MEDAKA_DIR:-${ASSEMBLY_RESULTS_DIR}/medaka_asm}"
FINAL_ASM_DIR="${FINAL_ASM_DIR:-${ASSEMBLY_RESULTS_DIR}/final_asm}"

####################################
# STEP 03 - ASSEMBLY QC OUTPUT DIRECTORIES
####################################

QUAST_DIR="${QUAST_DIR:-${QC_RESULTS_DIR}/quast}"
QUAST_SUMMARY_DIR="${QUAST_SUMMARY_DIR:-${QUAST_DIR}/summary}"

CHECKM2_DIR="${CHECKM2_DIR:-${QC_RESULTS_DIR}/checkm2}"
CHECKM2_INDIVIDUAL_DIR="${CHECKM2_INDIVIDUAL_DIR:-${CHECKM2_DIR}/individual}"
CHECKM2_SUMMARY_DIR="${CHECKM2_SUMMARY_DIR:-${CHECKM2_DIR}/summary}"

####################################
# STEP 04 - TAXONOMY OUTPUT DIRECTORIES
####################################

TAXONOMY_RESULTS_DIR="${TAXONOMY_RESULTS_DIR:-results/taxonomy}"

SOURMASH_DIR="${SOURMASH_DIR:-${TAXONOMY_RESULTS_DIR}/sourmash}"
SOURMASH_SIG_DIR="${SOURMASH_SIG_DIR:-${SOURMASH_DIR}/signatures}"
SOURMASH_SEARCH_DIR="${SOURMASH_SEARCH_DIR:-${SOURMASH_DIR}/searches}"

SKANI_DIR="${SKANI_DIR:-${TAXONOMY_RESULTS_DIR}/skani}"
SKANI_SEARCH_DIR="${SKANI_SEARCH_DIR:-${SKANI_DIR}/searches}"

MLST_DIR="${MLST_DIR:-${TAXONOMY_RESULTS_DIR}/mlst}"
MLST_INDIVIDUAL_DIR="${MLST_INDIVIDUAL_DIR:-${MLST_DIR}/individual}"
MLST_SUMMARY_DIR="${MLST_SUMMARY_DIR:-${MLST_DIR}/summary}"

TAXONOMY_MERGED_DIR="${TAXONOMY_MERGED_DIR:-${TAXONOMY_RESULTS_DIR}/merged}"

####################################
# STEP 05 - ANNOTATION OUTPUT DIRECTORIES
####################################

ANNOTATION_RESULTS_DIR="${ANNOTATION_RESULTS_DIR:-results/annotation}"

BAKTA_DIR="${BAKTA_DIR:-${ANNOTATION_RESULTS_DIR}/bakta}"
ANNOTATION_SUMMARY_DIR="${ANNOTATION_SUMMARY_DIR:-${ANNOTATION_RESULTS_DIR}/summary}"

####################################
# STEP 06 - CHARACTERIZATION OUTPUT DIRECTORIES
####################################

CHARACTERIZATION_RESULTS_DIR="${CHARACTERIZATION_RESULTS_DIR:-results/characterization}"

AMRFINDER_DIR="${AMRFINDER_DIR:-${CHARACTERIZATION_RESULTS_DIR}/amrfinderplus}"
AMRFINDER_INDIVIDUAL_DIR="${AMRFINDER_INDIVIDUAL_DIR:-${AMRFINDER_DIR}/individual}"
AMRFINDER_SUMMARY_DIR="${AMRFINDER_SUMMARY_DIR:-${AMRFINDER_DIR}/summary}"

PLASMIDFINDER_DIR="${PLASMIDFINDER_DIR:-${CHARACTERIZATION_RESULTS_DIR}/plasmidfinder}"
PLASMIDFINDER_INDIVIDUAL_DIR="${PLASMIDFINDER_INDIVIDUAL_DIR:-${PLASMIDFINDER_DIR}/individual}"
PLASMIDFINDER_SUMMARY_DIR="${PLASMIDFINDER_SUMMARY_DIR:-${PLASMIDFINDER_DIR}/summary}"

MOBSUITE_DIR="${MOBSUITE_DIR:-${CHARACTERIZATION_RESULTS_DIR}/mobsuite}"
MOBSUITE_INDIVIDUAL_DIR="${MOBSUITE_INDIVIDUAL_DIR:-${MOBSUITE_DIR}/individual}"
MOBSUITE_SUMMARY_DIR="${MOBSUITE_SUMMARY_DIR:-${MOBSUITE_DIR}/summary}"


#################################################
# STEP 07 — FINAL SUMMARY
#################################################

FINAL_SUMMARY_DIR="${FINAL_SUMMARY_DIR:-results/final_summary}"

####################################
# FILE PREFIXES
####################################
#READS
PORECHOP_PREFIX="${PORECHOP_PREFIX:-porechop}"
FILTLONG_PREFIX="${FILTLONG_PREFIX:-clean}"

#ASSEMBLY
FLYE_PREFIX="${FLYE_PREFIX:-onlyflye}"
MEDAKA_PREFIX="${MEDAKA_PREFIX:-medaka}"
FINAL_ASM_PREFIX="${FINAL_ASM_PREFIX:-final}"

#QC QUAST & CHECKM2
QUAST_PREFIX="${QUAST_PREFIX:-quast}"
CHECKM2_PREFIX="${CHECKM2_PREFIX:-checkm2}"

#TAXONOMY
SOURMASH_PREFIX="${SOURMASH_PREFIX:-sourmash}"
SKANI_PREFIX="${SKANI_PREFIX:-skani}"
MLST_PREFIX="${MLST_PREFIX:-mlst}"

#ANNOTATION
BAKTA_PREFIX="${BAKTA_PREFIX:-bakta}"

#CHARACTERIZATION
AMRFINDER_PREFIX="${AMRFINDER_PREFIX:-amrfinder}"
PLASMIDFINDER_PREFIX="${PLASMIDFINDER_PREFIX:-plasmidfinder}"
MOBSUITE_PREFIX="${MOBSUITE_PREFIX:-mobtyper}"

####################################
# PARAMETERS/FLAGS
####################################

####################################
# TRIMMING STRATEGY
####################################

TRIMMER="${TRIMMER:-fastplong}"   # fastplong | porechop_filtlong

####################################
# READ FILTERING
####################################

FASTPLONG_MIN_LENGTH="${FASTPLONG_MIN_LENGTH:-1000}"

PORECHOP_DISCARD_MIDDLE="${PORECHOP_DISCARD_MIDDLE:-true}"

FILTLONG_MIN_LENGTH="${FILTLONG_MIN_LENGTH:-500}"
FILTLONG_MIN_MEAN_Q="${FILTLONG_MIN_MEAN_Q:-75}"
FILTLONG_KEEP_PERCENT="${FILTLONG_KEEP_PERCENT:-95}"

####################################
# ASSEMBLY AND QC
####################################

FLYE_ASM_COVERAGE="${FLYE_ASM_COVERAGE:-50}"
FLYE_ITERATIONS="${FLYE_ITERATIONS:-3}"

MEDAKA_MODEL="${MEDAKA_MODEL:-r1041_e82_400bps_sup_v5.2.0}"
MEDAKA_BATCH_SIZE="${MEDAKA_BATCH_SIZE:-50}"

QUAST_MIN_CONTIG="${QUAST_MIN_CONTIG:-500}"

####################################
# TAXONOMY
####################################

SOURMASH_KSIZE="${SOURMASH_KSIZE:-31}"
SOURMASH_SCALED="${SOURMASH_SCALED:-1000}"

SKANI_TOPN="${SKANI_TOPN:-10}"
SKANI_BOTH_MIN_AF="${SKANI_BOTH_MIN_AF:-0.4}"

####################################
# ANNOTATION
####################################

BAKTA_GENUS="${BAKTA_GENUS:-}"
BAKTA_SPECIES="${BAKTA_SPECIES:-}"
BAKTA_TRANSLATION_TABLE="${BAKTA_TRANSLATION_TABLE:-11}"

####################################
# ANNOTATION
###################################

AMRFINDER_COVERAGE_MIN="${AMRFINDER_COVERAGE_MIN:-0.8}"

####################################
# THREADS
####################################

####################################
# NANOSTAT / MULTIQC
####################################

NANOSTAT_THREADS="${NANOSTAT_THREADS:-12}"

####################################
# SYLPH
####################################

SYLPH_THREADS="${SYLPH_THREADS:-12}"

####################################
# ASSEMBLY
####################################

FLYE_THREADS="${FLYE_THREADS:-12}"
MEDAKA_THREADS="${MEDAKA_THREADS:-12}"

PYPOLCA_THREADS="${PYPOLCA_THREADS:-10}"

####################################
# ASSEMBLY QC
####################################

QUAST_THREADS="${QUAST_THREADS:-8}"
CHECKM2_THREADS="${CHECKM2_THREADS:-12}"

####################################
# TAXONOMY
####################################

TOTAL_CORES=${TOTAL_CORES:-32}
THREADS_PER_JOB=${THREADS_PER_JOB:-8}

MLST_JOBS="${MLST_JOBS:-4}"

ASM_TYPE_FOR_MERGE="${ASM_TYPE_FOR_MERGE:-long}"
FINAL_SPECIES_TABLE="${FINAL_SPECIES_TABLE:-${TAXONOMY_MERGED_DIR}/final_species_table.tsv}"

####################################
# ANNOTATION
####################################

BAKTA_THREADS="${BAKTA_THREADS:-12}"

####################################
# CHARACTERIZATION
####################################

AMRFINDER_THREADS="${AMRFINDER_THREADS:-12}"
PLASMIDFINDER_THREADS="${PLASMIDFINDER_THREADS:-4}"
MOBSUITE_THREADS="${MOBSUITE_THREADS:-12}"

####################################
# CONDA ENVIRONMENTS
####################################

#### Main environments
ENV_TRIMMING="${ENV_TRIMMING:-env_trim_longreads}"
ENV_SYLPH="${ENV_SYLPH:-sylph}"
ENV_ASSEMBLY="${ENV_ASSEMBLY:-env_assembly}"
ENV_QUAST="${ENV_QUAST:-env_quast}"
ENV_CHECKM2="${ENV_CHECKM2:-checkm2}"
ENV_TAXONOMY="${ENV_TAXONOMY:-env_taxonomy}"
ENV_ANNOTATION="${ENV_ANNOTATION:-env_annotation}"
ENV_CHARACTERIZATION="${ENV_CHARACTERIZATION:-env_characterization}"

#ENV_AMRFINDER="${ENV_AMRFINDER:-epitools_amrfinder}"
#ENV_MOBSUITE="${ENV_MOBSUITE:-mob_suite}"
#ENV_PLASMIDFINDER="${ENV_PLASMIDFINDER:-epitools}"

ENV_ANALYSIS="${ENV_ANALYSIS:-python_r_utils_env}"

#### Helper scripts
RSCRIPT_PLOT="${RSCRIPT_PLOT:-scripts/abundance_species_plot.R}"
MERGE_TAXONOMY_BIN="${MERGE_TAXONOMY_BIN:-scripts/merge_taxonomy_reports.py}"
#SCRIPT_SUMMARY="${SCRIPT_SUMMARY:-amrfinder_unique_symbols_by_type_and_class.py}"
#SCRIPT_MATRIX="${SCRIPT_MATRIX:-amrfinder_presence_absence_by_type.py}"
#SCRIPT_HEATMAP="${SCRIPT_HEATMAP:-plot_amrfinder_heatmap_amr_virulence.py}"
#SCRIPT_PLASMIDFINDER_SUMMARY="${SCRIPT_PLASMIDFINDER_SUMMARY:-summarize_plasmidfinder.py}"
SCRIPT_MERGE_FINAL_REPORT="${SCRIPT_MERGE_FINAL_REPORT:-scripts/merge_final_report.py}"