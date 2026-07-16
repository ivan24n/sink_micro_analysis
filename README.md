# Spatiotemporal dynamics of multi-kingdom microbial communities in hospital sinks

This repository contains the analysis scripts used in the study
**"Spatiotemporal dynamics of multi-kingdom microbial communities in hospital sinks."**

It includes the bash-based preprocessing pipeline for 16S, ITS, and 18S amplicon metagenomics (`scripts/`)
and the R scripts used for downstream statistical analysis and figure generation (`R/`).

## Repository structure

```
.
├── scripts/
│   ├── 01_primer_trimming.sh          # Adapter/primer removal (16S, 18S, ITS) with Cutadapt
│   ├── 02_quality_trimming.sh         # Quality/length trimming with Trimmomatic
│   ├── 03_host_read_removal.sh        # Removal of human reads with Bowtie2
│   └── 04_taxonomic_classification.sh # Taxonomic classification with Kraken2
├── Databases/
│   ├── build_16S_database.sh          # Kraken2 DB for 16S (pre-built RDP index)
│   ├── build_18S_database.sh          # Kraken2 DB for 18S (EukRibo + SILVA + PR2 + protozoa)
│   └── build_ITS_database.sh          # Kraken2 DB for ITS (UNITE + NCBI fungi)
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
classification, depending on the dataset — amplicon vs. shotgun). The
scripts in `db_build/` are a one-time setup step: build each Kraken2
database once, then reuse it across runs of
`scripts/04_taxonomic_classification.sh`.

## Requirements

| Tool | Used by |
|---|---|
| [Cutadapt](https://cutadapt.readthedocs.io/) | `01_primer_trimming.sh` |
| [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) + Java | `02_quality_trimming.sh` |
| [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) | `03_host_read_removal.sh` |
| [Kraken2](https://ccb.jhu.edu/software/kraken2/) (+ [GNU parallel](https://www.gnu.org/software/parallel/), optional) | `04_taxonomic_classification.sh` |
| `kraken2-build`, `wget`, `tar` | `db_build/*.sh` |

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
    ├── tables/                # General (non-stratified) tables
    │   ├── area/               # Area-grouped statistics (alpha/beta diversity)
    │   └── period/             # Period-grouped statistics (alpha/beta diversity)
    └── figures/                # General figures
        ├── area/
        └── period/
```

This structure is applied consistently across all three markers.

## Usage

### 0. Build the Kraken2 databases (one-time setup)

```bash
cd db_build/
./build_16S_database.sh   # -> kraken2_db_16S/
./build_18S_database.sh   # -> kraken2_db_18S/
./build_ITS_database.sh   # -> kraken2_db_ITS/  (requires manually downloading
                          #    the UNITE archive first, see script header)
```

Each script downloads/builds its database in the current working
directory, under a folder named `kraken2_db_{16S,18S,ITS}` — these
folder names/paths are what you then pass to the `-d` flag of
`scripts/04_taxonomic_classification.sh`. This only needs to be run once;
the resulting database folders can be reused across all sequencing runs.

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
./scripts/04_taxonomic_classification.sh -i host_filtered/16S/ -o kraken2_out/16S/ -d /path/to/kraken2_db_16S -t 4 -p 1
./scripts/04_taxonomic_classification.sh -i host_filtered/18S/ -o kraken2_out/18S/ -d /path/to/kraken2_db_18S -t 4 -p 1 -c 0.1
```

Runs Kraken2 per sample, writing classification output and reports;
`-p` enables parallel processing of samples (requires GNU parallel).
`-c` sets Kraken2's `--confidence` threshold (0-1); for 18S rRNA
classification, a higher value than the default is recommended, since
eukaryotic reads are more prone to spurious low-confidence hits across
kingdoms. The example above uses `0.1` as a placeholder — replace it with
the value used in the manuscript's methods.

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
- `db_build/*.sh` were adapted the same way (English, generic/documented,
  no logic changes) with two safe additions: a check that the UNITE
  archive exists before `build_ITS_database.sh` proceeds (with download
  instructions, since UNITE has no stable direct-download URL), and an
  unused leftover variable removed from `build_18S_database.sh`. All
  three now build into consistently named folders
  (`kraken2_db_16S/18S/ITS`) matching the `-d` flag used in
  `04_taxonomic_classification.sh`.
- All hard-coded personal file paths (e.g. `/home/ivan/...`) were replaced
  with a `work_dir` variable that must be set by whoever runs the scripts.
- `05_microbial_analysis.R` and `03_rarefaction.R` were **not** restructured
  into loops/functions, despite containing repetition across markers
  (16S/18S/ITS). The repeated blocks differ in enough small ways
  (thresholds, regex patterns, column selections) that collapsing them
  safely would need a careful, dedicated pass — happy to do this as a
  follow-up if wanted.
- The per-area (A/B/C) sub-analysis section that previously existed at the
  end of `05_microbial_analysis.R` (16S only) has been removed at the
  authors' request, as it was an extra analysis not needed for the
  manuscript.
- The "Period Analysis" section (16S only) has also been removed at the
  authors' request: this was a further breakdown of each sampling period
  ("Empty" and "Multifunctional") into area-level sample subgroups, and
  was likewise deemed an extra analysis not needed for the manuscript.
- Output tables and figures in `05_microbial_analysis.R` now follow one
  consistent folder scheme per marker (`tables/area/`, `tables/period/`,
  `figures/area/`, `figures/period/`), instead of the original's three
  inconsistent naming conventions (`tables/..._area_16S.tsv`,
  `tables_area/...`, `tables_unit/...`).
- Fixed a copy-paste bug found in the original script while reorganizing
  these paths: in the per-area Bray/Jaccard pairwise PERMANOVA saves,
  Jaccard results were written to the general `tables/` folder instead of
  the area-specific one.
- Per the authors' request, the 18S dataset now only computes Jaccard
  beta-diversity (no Bray-Curtis): the corresponding Bray-Curtis code,
  plots, and PERMANOVA tables for 18S were removed, and the 18S composite
  summary figures (`figures/area/plot_18S.svg`,
  `figures/period/plot_18S.svg`) use a 6-panel layout instead of the
  7-panel layout used for 16S/ITS (which still include both metrics).
- `03_rarefaction.R`: the estimated-richness bar plot is now built the
  same way for all three markers — joining explicitly by `Sample` rather
  than assuming matching row order (as the original 18S block already
  did; the 16S/ITS blocks were updated to match).

## Citation

If you use these scripts, please cite:

> Linares-Ambohades I, Guerra-Pinto N, Mingo-Ramirez S, Serrano-Calleja S,
> Amaro F, Alastruey-Izquierdo A, Soriano MC, de Pablo R, Lanza VF,
> Cantón R, Baquero F, Coque TM, Pérez-Cobas AE. Spatiotemporal dynamics
> of multi-kingdom microbial communities in hospital sinks.
> [Journal, year, DOI — to be added upon publication]

**Authors and affiliations**

Iván Linares-Ambohades¹,², Natalia Guerra-Pinto¹,³, Sandra Mingo-Ramirez¹,
Silvia Serrano-Calleja¹, Francisco Amaro⁴, Ana Alastruey-Izquierdo⁵,
María Cruz Soriano⁶,⁷, Raúl de Pablo⁶,⁷, Val F. Lanza³,ˣ, Rafael Cantón¹,³,
Fernando Baquero¹,ˣ, Teresa M. Coque¹,³*, Ana Elena Pérez-Cobas¹,³*

1. Department of Microbiology, Ramón y Cajal Institute for Health Research (IRYCIS), Ramón y Cajal University Hospital, Madrid, Spain
2. Escuela de Doctorado, Universidad Autónoma de Madrid, Madrid, Spain
3. CIBER in Infectious Diseases (CIBERINFEC), Madrid, Spain
4. University Complutense of Madrid (UCM), Madrid, Spain
6. Intensive Medicine, Ramón y Cajal University Hospital and Ramón y Cajal Health Research Institute (IRYCIS), Madrid, Spain
7. University of Alcalá (UAH), Madrid, Spain
8. Translational Genomics (NGS) and Bioinformatics Unit, Ramón y Cajal Health Research Institute (IRYCIS), Madrid, Spain

*Correspondence:*
Teresa M. Coque — mariateresa.coque@salud.madrid.org, teresacoque@gmail.com
Ana Elena Pérez-Cobas — anaelena.perez@salud.madrid.org

*Department of Microbiology, Ramón y Cajal Institute for Health Research
(IRYCIS), Ramón y Cajal University Hospital, Carretera de Colmenar, km.
9.1, Madrid 28034, Spain.*
