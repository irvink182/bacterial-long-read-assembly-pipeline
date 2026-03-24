# Bacterial WGS long-read / hybrid pipeline (v1.0)
A modular pipeline for bacterial whole-genome sequencing (WGS) analysis from long-read data, with optional short-read polishing for hybrid workflows.

## Overview
This pipeline processes bacterial sequencing data from raw reads to a final integrated summary table including:

- Sequencing reads QC
- Assembly
- Assembly QC
- Taxonomy
- Annotation
- Genomic characterization
- Final Summary Table

The workflow is currently implemented as a series of modular Bash scripts (`step_01` to `step_07`) invoked through a main launcher script (`run_all_pipeline.sh`).

Two analysis modes are supported:

- `long`: long-read assembly only
- `hybrid`: long-read assembly followed by additional short-read polishing

Please note that `hybrid` mode refers to:

- long-read assembly
- ONT correction with `Medaka`
- short-read polishing with `pypolca`

It does **not** perform true hybrid *de novo* assembly with a dedicated hybrid assembler.

## Workflow design

### `long` mode
```text
long reads
  -> Trimming (porechop + filtlong or fastplong) / read QC
  -> Flye
  -> Medaka
  -> Assembly QC (quast & CheckM2)
  -> Taxonomy (Sylph, Sourmash & skani)
  -> Annotation (Bakta)
  -> Genomic characterization (amrfinder, mob_typer, plasmidfinder)
  -> Final Summary Table
```
### `hybrid` mode
```text
long reads + short reads
  -> Trimming (porechop + filtlong or fastplong) / read QC
  -> Flye
  -> Medaka
  -> pypolca (for short-read polishing)
  -> Assembly QC (quast & CheckM2)
  -> Taxonomy (Sylph, Sourmash & skani)
  -> Annotation (Bakta)
  -> Genomic characterization (amrfinder, mob_typer, plasmidfinder)
  -> Final Summary Table
```
## Proposed layout
```text
Project/
├── config/
│   ├── config.sh
│   └── databases.config
├── data/
│   └── raw_reads/
├── envs/
│   ├── env_analysis.yml
│   ├── env_annotation.yml
│   ├── env_assembly.yml
│   ├── env_characterization.yml
│   ├── env_checkm2.yml
│   ├── env_quast.yml
│   ├── env_sylph.yml
│   ├── env_taxonomy.yml
│   └── env_trim_longreads.yml
├── logs/
├── results/
│   ├── annotation/
│   ├── assemblies/
│   ├── characterization/
│   ├── final_summary/
│   ├── quality_check/
│   ├── reads/
│   └── taxonomy/
├── scripts/
│   ├── abundance_species_plot.R
│   ├── amrfinder_presence_absence_by_type.py
│   ├── amrfinder_unique_symbols_by_type_and_class.py
│   ├── merge_final_report.py
│   ├── merge_taxonomy_reports.py
│   ├── plot_amrfinder_heatmap_amr_virulence.py
│   ├── summarize_plasmidfinder.py
│   └── utils.sh
├── steps/
│   ├── step_01_trimming.sh
│   ├── step_02_assembly.sh
│   ├── step_03_assembly_qc.sh
│   ├── step_04_taxonomy.sh
│   ├── step_05_annotation.sh
│   ├── step_06_characterization.sh
│   └── step_07_final_summary.sh
├── run_all_pipeline.sh
├── run_pipeline.sh
├── samplesheet.tsv
└── README.md
```
## Installation
### conda distribution
Before installing the piepñine, make sure that you have installed `Conda` distribution in your system. 
If `conda` is not available on your system, you need to install the distribution first.

On Linux:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh
```
`conda` Official documentation:
- [Conda installation guide](https://docs.conda.io/docs/user-guide/install/download.html)
- [Installing on Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

### pipeline installation
Clone the repository and create the required conda environments from the `envs/` directory.

Example:

```bash
conda env create -f envs/env_trim_longreads.yml
conda env create -f envs/env_assembly.yml
conda env create -f envs/env_quast.yml
conda env create -f envs/env_checkm2.yml
conda env create -f envs/env_taxonomy.yml
conda env create -f envs/env_annotation.yml
conda env create -f envs/env_characterization.yml
conda env create -f envs/env_analysis.yml
```

## Configuration
Pipeline internal settings (paths, prefixes, etc) are stored in:

- `config/config.sh`: general pipeline parameters and directory settings
- `config/databases.config`: paths to external databases

Please note that before running the pipeline, make sure all database paths are correctly defined for your system.

## Running the pipeline
To run the complete pipeline (`all seven steps`), use:
```bash
bash run_all_pipeline.sh --samples samplesheet.tsv --trimmer porechop_filtlong
```
Please note that if you want to use `fastplong` instead of porechop & filtlong, you have to run `--trimmer fastplong` option:
```bash
bash run_all_pipeline.sh --samples samplesheet.tsv --trimmer fastplong
```

You can also run a specific step of the pipeline.
For example, if you want to perform only the `trimming` step:
```bash
bash run_pipeline.sh --step trimming --samples samplesheet.tsv --trimmer porechop_filtlong
```

For example, if you want to perform only the `taxonomy` step:
```bash
bash run_pipeline.sh --step taxonomy --samples samplesheet.tsv
```

## Pipeline steps

### `step_01_trimming.sh`
Long-read preprocessing:
- read trimming with `fastplong` or `porechop + filtlong`
- raw & cleaned read QC with `NanoStat`
- cleaned read taxonomic profiling with `sylph`

### `step_02_assembly.sh`
Assembly and polishing:
- long-read assembly with `Flye`
- ONT correction with `Medaka`
- optional short-read polishing with `pypolca` when `asm_type=hybrid`
- Final assembly FASTA file

### `step_03_assembly_qc.sh`
Assembly quality control:
- Assembly stats with `QUAST`
- Assembly completeness and contamination with `CheckM2`

### `step_04_taxonomy.sh`
Combination of different tools for taxonomic assignment:
- `skani` (based on GTDB r226 database)
- `sourmash` (based on GTDB r226 database)
- `sylph` (based on GTDB r226 database)
- merged taxonomy summary table

### `step_05_annotation.sh`
Genome annotation:
- `Bakta`

### `step_06_characterization.sh`
Genomic characterization:
- Resistance determinants: `AMRFinderPlus`
- in silico identification and typing of plasmids: `PlasmidFinder`
- Plasmid typing: `MOB-typer`
- PubMLST typing schemes: `MLST`

### `step_07_final_summary.sh`
Final reporting:
- merge of all summary tables
- Final report per-sample

### Input samplesheet

The pipeline uses a tab-delimited samplesheet (for example, a `.tsv` extension) with header

### Required columns

| Column | Description |
|---|---|
| `Sample_ID` | Sample ID |
| `asm_type` | Analysis mode: `long` or `hybrid` |
| `expected_genome_size` | Expected genome size in bp |
| `long_reads` | Long-read FASTQ file |
| `short_r1` | Short-read R1 FASTQ file |
| `short_r2` | Short-read R2 FASTQ file |

### Example samplesheet
| Sample_ID | asm_type | expected_genome_size | long_reads | short_r1 | short_r2 |
|---|---|---:|---|---|---|
| isolate1 | hybrid | 6000000 | isolate1.fastq.gz | isolate1.R1.fastq.gz | isolate1.R2.fastq.gz |
| isolate_2 | long | 5300000 | isolate_2.fastq.gz | NA | NA |

Example as TSV:

```tsv
Sample_ID	asm_type	expected_genome_size	long_reads	short_r1	short_r2
isolate1	hybrid	5000000	isolate1.fastq.gz	isolate1.R1.fastq.gz	isolate1.R2.fastq.gz
isolate_2	long	5300000	isolate_2.fastq.gz	NA	NA
```

## Expected outputs

Main output directories are written to `results/`:

- `results/reads/`: cleaned reads (fastq files) and read QC tables
- `results/assemblies/`: Final assembly files (`Medaka` or polished `pypolca`)
- `results/quality_check/`: QUAST and CheckM2 results
- `results/taxonomy/`: taxonomic assignment outputs
- `results/annotation/`: Bakta annotations
- `results/characterization/`: AMR, plasmid, & MLST typing results
- `results/final_summary/`: Final merged summary report


## Software requirements

### Read QC and trimming
- [NanoStat](https://github.com/wdecoster/nanostat)
- [fastplong](https://github.com/OpenGene/fastplong)
- [Porechop](https://github.com/rrwick/Porechop)
- [Filtlong](https://github.com/rrwick/Filtlong)
- [sylph](https://github.com/bluenote-1577/sylph)

### Assembly and polishing
- [Flye](https://github.com/mikolmogorov/Flye)
- [Medaka](https://github.com/nanoporetech/medaka)
- [pypolca](https://github.com/gbouras13/pypolca)

### Assembly QC
- [QUAST](https://github.com/ablab/quast)
- [CheckM2](https://github.com/chklovski/CheckM2)

### Taxonomy
- [skani](https://github.com/bluenote-1577/skani)
- [sourmash](https://github.com/sourmash-bio/sourmash)
- [sylph](https://github.com/bluenote-1577/sylph)

### Annotation
- [Bakta](https://github.com/oschwengers/bakta)

### Characterization
- [AMRFinderPlus](https://github.com/ncbi/amr)
- [PlasmidFinder](https://github.com/genomicepidemiology/plasmidfinder)
- [MOB-suite](https://github.com/phac-nml/mob-suite)
- [mlst](https://github.com/tseemann/mlst)
