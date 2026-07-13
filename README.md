# Spatiotemporal dynamics of multi-kingdom microbial communities in hospital sinks

This repository contains the analysis scripts used in the study
**"Spatiotemporal dynamics of multi-kingdom microbial communities in hospital sinks."**

It includes the bash-based amplicon and shotgun metagenomics preprocessing
pipeline (`scripts/`) and, once added, the R scripts used for downstream
statistical analysis and figure generation (`R/`).

## Repository structure

```
.
├── scripts/
│   ├── 01_primer_trimming.sh          # Adapter/primer removal (16S, 18S, ITS) with Cutadapt
│   ├── 02_quality_trimming.sh         # Quality/length trimming with Trimmomatic
│   ├── 03_host_read_removal.sh        # Removal of human reads with Bowtie2
│   └── 04_taxonomic_classification.sh # Taxonomic classification with Kraken2
└── R/
    ├── 01_import_kraken.R             # Parses Kraken2 reports into a combined table
    ├── 02_build_taxonomy_tables.R     # Genus-level abundance tables per marker
    ├── 03_rarefaction.R               # Rarefaction curves and rarefied OTU/ASV tables
    ├── 04_build_taxa_tables.R         # Taxonomy reference tables (lineage + OTU_ID)
    └── 05_microbial_analysis.R        # Composition, alpha/beta diversity, PERMANOVA
```

The numbering above reflects a suggested processing order; adjust it if your
pipeline runs the steps in a different sequence (e.g., quality trimming
before primer removal, or host read removal before/after taxonomic
classification, depending on the dataset — amplicon vs. shotgun).

## Requirements

| Tool | Used by |
|---|---|
| [Cutadapt](https://cutadapt.readthedocs.io/) | `01_primer_trimming.sh` |
| [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) + Java | `02_quality_trimming.sh` |
| [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) | `03_host_read_removal.sh` |
| [Kraken2](https://ccb.jhu.edu/software/kraken2/) (+ [GNU parallel](https://www.gnu.org/software/parallel/), optional) | `04_taxonomic_classification.sh` |

All bash scripts expect paired-end FASTQ files named with the Illumina-style
convention `<sample>_R1_001.fastq.gz` / `<sample>_R2_001.fastq.gz`.

The R scripts require R (≥ 4.2 recommended) with the following packages:
`tidyverse`, `stringr`, `doParallel`, `dbplyr`, `vegan`, `patchwork`,
`ggtext`, `microbiome`, `RColorBrewer`, `pairwiseAdonis`, `ggpubr`,
`ggvenn`, `phyloseq`, `fantaxtic`, `scales`.

Each R script has a `work_dir` variable near the top — set it to the root
of your project directory before running. The R scripts expect the
following folder layout under `work_dir` (created by the earlier steps):

```
work_dir/
├── Control_sink/{16S,18S,ITS}/Analysis/Kraken/reports/   # Kraken2 reports
└── {16S,18S,ITS}_SINKS/
    ├── tables/     # Intermediate and final tables
    ├── figures/    # Marker-level figures
    ├── tables_area/  # Area-stratified tables (16S only, see note below)
    └── figures_area/ # Area-stratified figures (16S only, see note below)
```

## Usage

### 1. Primer trimming (amplicon data: 16S / 18S / ITS)

```bash
./scripts/01_primer_trimming.sh -i raw_fastq/ -o primer_trimmed/ -m 16S
./scripts/01_primer_trimming.sh -i raw_fastq/ -o primer_trimmed/ -m 18S
./scripts/01_primer_trimming.sh -i raw_fastq/ -o primer_trimmed/ -m ITS
```

This single script replaces three previously separate, near-identical
scripts (one per marker gene); the `-m` flag selects the adapter and primer
sequences for the target marker. The underlying Cutadapt command and the
sequences themselves are unchanged from the original per-marker scripts.

### 2. Quality trimming

```bash
./scripts/02_quality_trimming.sh -i primer_trimmed/ -o quality_trimmed/
```

Produces `paired/`, `unpaired/`, and `logs/` subfolders in the output
directory. Trimming parameters (head crop, sliding window, minimum length,
minimum average quality) are fixed in the script to match the settings used
in this study.

### 3. Host (human) read removal

```bash
./scripts/03_host_read_removal.sh -x /path/to/human_bowtie2_index -i quality_trimmed/paired/ -o host_filtered/ -t 8
```

Aligns reads to a human reference genome index and retains only read pairs
that do not align (non-host reads).

### 4. Taxonomic classification

```bash
./scripts/04_taxonomic_classification.sh -i host_filtered/ -o kraken2_out/ -d /path/to/kraken2_db -t 4 -p 1
```

Runs Kraken2 per sample, writing classification output and reports;
`-p` enables parallel processing of samples (requires GNU parallel).

### 5. Downstream analysis (R)

Run the R scripts in order from within the `R/` directory (or adjust the
`source()`/paths accordingly):

```r
source("01_import_kraken.R")        # defines import_kraken()
source("02_build_taxonomy_tables.R") # -> tax_tab_{16S,18S,ITS}.tsv
source("03_rarefaction.R")           # -> rarefacted_table_{16S,18S,ITS}.tsv
source("04_build_taxa_tables.R")     # -> TAXA_{16S,18S,ITS}.tsv
source("05_microbial_analysis.R")    # composition, diversity, PERMANOVA
```

`05_microbial_analysis.R` additionally expects a metadata table
(`rare_meta_table_<marker>.tsv`) and a manually curated taxonomy table
(`TAXA_fullcurated_<marker>.tsv`) and OTU table
(`otutable_fullcurated_<marker>.tsv`) per marker; these correspond to the
outputs of steps 1–4 after any manual curation described in the
manuscript's methods.

## Notes on this repository

- All scripts were adapted from the original lab pipeline for public
  release: comments and messages were translated to English and
  documentation headers were added for clarity, but the underlying
  commands and parameters (Cutadapt adapters/primers, Trimmomatic
  parameters, Bowtie2 options, Kraken2 options) are unchanged from the
  scripts used to generate the results in the manuscript.
- `01_primer_trimming.sh` consolidates three original per-marker scripts
  into one; this is an interface simplification only (see above).
- `04_build_taxa_tables.R` consolidates three original per-marker blocks
  (identical logic, different file paths) into a loop; the function and
  its outputs are unchanged.
- All hard-coded personal file paths (e.g. `/home/ivan/...`) were replaced
  with a `work_dir` variable that must be set by whoever runs the scripts.
- `05_microbial_analysis.R` and `03_rarefaction.R` were **not** restructured
  into loops/functions, despite containing substantial repetition across
  markers (16S/18S/ITS) and, in `05_microbial_analysis.R`, across areas
  (A/B/C). The repeated blocks differ in enough small ways (thresholds,
  regex patterns, column selections) that collapsing them safely needs a
  careful, dedicated pass — happy to do this as a follow-up if wanted.

### Open questions (pending confirmation, not yet changed)

1. In `05_microbial_analysis.R`, within each area sub-analysis (A/B/C,
   16S only), the Bray-Curtis pairwise PERMANOVA results are saved to
   `tables_area/`, but the Jaccard results are saved to `tables/` instead
   — this pattern repeats in all three sub-blocks and looks like a
   copy-paste inconsistency rather than an intentional choice.
2. The area-stratified (A/B/C) sub-analysis in `05_microbial_analysis.R`
   is only present for the 16S dataset — there is no 18S/ITS equivalent
   in this script. Confirm whether this is intentional.
3. In `03_rarefaction.R`, the estimated-richness bar plot is built
   differently for 16S/ITS (values taken by position, assuming row order
   matches) versus 18S (values taken via an explicit join by `Sample`).
   The join-based approach is safer; confirm whether to standardize all
   three to it.

## Citation

If you use these scripts, please cite:

> [Iván Linares-Ambohades, Natalia Guerra-Pinto, Sandra Mingo-Ramirez, Silvia Serrano-Calleja, Francisco Amaro, Ana Alastruey-Izquierdo, María Cruz Soriano, Raúl de Pablo6, Val F. Lanza, Rafael Cantón, Fernando Baquero, Teresa M. Coque1, Ana Elena Pérez-Cobas]. Spatiotemporal dynamics of multi-kingdom microbial
> communities in hospital sinks. [Journal, year, DOI — to be added upon
> publication]
