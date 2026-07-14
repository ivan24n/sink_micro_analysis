##############################################################################
# 02_build_taxonomy_tables.R
#
# Builds per-marker genus-level taxonomic abundance tables (16S, 18S, ITS)
# from Kraken2 reports, using import_kraken() to combine per-sample reports
# and then cleaning/collapsing lineages down to genus level.
#
# Part of the analysis pipeline for:
#   "Spatiotemporal dynamics of multi-kingdom microbial communities in
#    hospital sinks"
#
# Requires 01_import_kraken.R (in this repository) to be sourced first.
##############################################################################

library(tidyverse)
library(stringr)

# Assumes 01_import_kraken.R is in the same directory as this script.
source("01_import_kraken.R")

# Set this to the root of the project directory, i.e. the folder that
# contains the Control_sink/ subfolder with the per-marker Kraken2 reports.
work_dir <- "path/to/project"

#' Clean a semicolon-delimited lineage string
#'
#' Drops the root/domain-level entry, drops any sub-rank levels not marked
#' with a recognized single-letter-plus-digits prefix, and optionally
#' strips the remaining rank prefixes (e.g. "G__").
#'
#' @param taxa_string Semicolon-delimited lineage string.
#' @param remove_main_prefixes If TRUE, strip rank prefixes (e.g. "G__")
#'   from the retained levels.
#' @return Cleaned, semicolon-delimited lineage string.
clean_taxa <- function(taxa_string, remove_main_prefixes = FALSE) {
  taxa_levels <- strsplit(taxa_string, ";")[[1]]
  taxa_levels <- trimws(taxa_levels)
  cleaned_levels <- character()
  if (length(taxa_levels) > 0) {
    taxa_levels <- taxa_levels[-1]
  }
  for (level in taxa_levels) {
    if (grepl("^[A-Z]+\\d+__", level)) {
      next
    }
    if (remove_main_prefixes) {
      level <- gsub("^[A-Z]+__", "", level)
    }
    cleaned_levels <- c(cleaned_levels, level)
  }
  cleaned_taxa <- paste(cleaned_levels, collapse = ";")
  return(cleaned_taxa)
}

#' Extract the main taxonomic ranks from a lineage string
#'
#' Keeps only the domain/phylum/class/order/family/genus levels (identified
#' by their "D__"/"P__"/.../"G__" prefixes), in that order, with the
#' prefixes stripped.
#'
#' @param tax_string Semicolon-delimited lineage string.
#' @return Semicolon-delimited string of the six main taxonomic ranks.
extract_taxonomy <- function(tax_string) {
  levels <- strsplit(tax_string, ";")[[1]]
  levels <- trimws(levels)

  main_prefixes <- c("D__", "P__", "C__", "O__", "F__", "G__")
  cleaned_levels <- character()

  for (prefix in main_prefixes) {
    match <- levels[grepl(paste0("^", prefix), levels)]
    if (length(match) > 0) {
      cleaned <- sub("^.*?__", "", match[1])
      cleaned_levels <- c(cleaned_levels, cleaned)
    }
  }

  return(paste(cleaned_levels, collapse = ";"))
}


# ---- 16S ----
files_16S <- dir(file.path(work_dir, "Control_sink/16S/Analysis/Kraken/reports/"),
                  pattern = "_kraken2_report.txt", full.names = TRUE)
table_16S <- import_kraken(files_16S, threads = 20)
write.table(table_16S,
            file.path(work_dir, "Control_sink/16S/Analysis/Kraken/tax_abundance.tsv"),
            col.names = TRUE)

filtab_16S <- table_16S %>%
  mutate(Sample = gsub("_L001_kraken2_report.txt", "", Sample)) %>%
  filter(TaxRank == "G" & depth > 6) %>%
  select(Fullname2, Reads, Sample)

final_table_16S <- filtab_16S %>%
  pivot_wider(names_from = Sample, values_from = Reads, values_fill = 0) %>%
  as.data.frame()

final_table_16S <- final_table_16S %>%
  mutate(total = rowSums(across(-1))) %>%
  select(Fullname2, total, everything()) %>%
  filter(total > 1)

final_table_16S$Fullname2 <- sapply(final_table_16S$Fullname2, clean_taxa,
                                    remove_main_prefixes = TRUE)

final_table_16S <- final_table_16S %>%
  group_by(Fullname2) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  as.data.frame()

rownames(final_table_16S) <- final_table_16S$Fullname2
final_table_16S <- final_table_16S %>%
  arrange(desc(total)) %>%
  select(-total, -Fullname2)

write.table(final_table_16S,
            file.path(work_dir, "Control_sink/16S/Analysis/tables/tax_tab_16S.tsv"),
            sep = "\t", col.names = NA, row.names = TRUE)

# ---- 18S ----
files_18S <- dir(file.path(work_dir, "Control_sink/18S/Analysis/Kraken/reports/"),
                  pattern = "_kraken2_report.txt", full.names = TRUE)
table_18S <- import_kraken(files_18S, threads = 20)
write.table(table_18S,
            file.path(work_dir, "Control_sink/18S/Analysis/Kraken/reports/tax_abundance.tsv"),
            col.names = TRUE, row.names = FALSE, sep = "\t")

filtab_18S <- table_18S %>%
  mutate(Sample = gsub("_L001_kraken2_report.txt", "", Sample)) %>%
  filter(TaxRank == "G" & depth > 6 & !str_detect(Fullname2, "D__Bacteria|D__Archaea")) %>%
  select(Fullname2, Reads, Sample)

final_table_18S <- filtab_18S %>%
  pivot_wider(names_from = Sample, values_from = Reads, values_fill = 0) %>%
  as.data.frame()

final_table_18S <- final_table_18S %>%
  mutate(total = rowSums(across(-1))) %>%
  select(total, everything())

final_table_18S$Fullname2 <- sapply(final_table_18S$Fullname2, extract_taxonomy)
final_table_18S <- final_table_18S %>%
  group_by(Fullname2) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  as.data.frame()

rownames(final_table_18S) <- final_table_18S$Fullname2
final_table_18S <- final_table_18S %>%
  arrange(desc(total)) %>%
  select(-total, -Fullname2)

write.table(final_table_18S,
            file.path(work_dir, "Control_sink/18S/Analysis/tables/tax_tab_18S.tsv"),
            sep = "\t", col.names = NA, row.names = TRUE)

# ---- ITS ----
files_ITS <- dir(file.path(work_dir, "Control_sink/ITS/Analysis/Kraken/reports/"),
                  pattern = "_kraken2_report.txt", full.names = TRUE)
table_ITS <- import_kraken(files_ITS, threads = 20)
write.table(table_ITS,
            file.path(work_dir, "Control_sink/ITS/Analysis/Kraken/tax_abundance.tsv"),
            col.names = TRUE)

filtab_ITS <- table_ITS %>%
  mutate(Sample = gsub("_L001_kraken2_report.txt", "", Sample)) %>%
  filter(TaxRank == "G" & depth > 6) %>%
  select(Fullname2, Reads, Sample)

final_table_ITS <- filtab_ITS %>%
  pivot_wider(names_from = Sample, values_from = Reads, values_fill = 0) %>%
  as.data.frame()

final_table_ITS$Fullname2 <- sapply(final_table_ITS$Fullname2, clean_taxa,
                                    remove_main_prefixes = TRUE)

final_table_ITS$Fullname2 <- gsub("Fungi;", "", final_table_ITS$Fullname2)

final_table_ITS <- final_table_ITS %>%
  group_by(Fullname2) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  as.data.frame()

final_table_ITS <- final_table_ITS %>%
  mutate(total = rowSums(across(-1))) %>%
  select(total, everything()) %>%
  filter(total > 1)

rownames(final_table_ITS) <- final_table_ITS$Fullname2
final_table_ITS <- final_table_ITS %>%
  arrange(desc(total)) %>%
  select(-total, -Fullname2)

write.table(final_table_ITS,
            file.path(work_dir, "Control_sink/ITS/Analysis/tables/tax_tab_ITS.tsv"),
            sep = "\t", col.names = NA, row.names = TRUE)
