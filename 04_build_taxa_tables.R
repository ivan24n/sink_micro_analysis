##############################################################################
# 04_build_taxa_tables.R
#
# Builds the taxonomy reference table for each marker (16S, 18S, ITS) from
# the rarefied OTU/ASV table: splits the lineage string into individual
# ranks and assigns a stable OTU_ID to each row.
#
# Part of the analysis pipeline for:
#   "Spatiotemporal dynamics of multi-kingdom microbial communities in
#    hospital sinks"
##############################################################################

library(dbplyr)

# Set this to the root of the project directory, i.e. the folder that
# contains the 16S_SINKS/, 18S_SINKS/, and ITS_SINKS/ subfolders.
work_dir <- "path/to/project"

#' Build a taxonomy reference table from a rarefied OTU/ASV table
#'
#' @param rute_OTU_table Path to a rarefied OTU/ASV table (tab-separated,
#'   first column a semicolon-delimited lineage string, remaining columns
#'   per-sample counts).
#' @return A data frame with one row per taxon, an assigned OTU_ID, the
#'   original lineage string, and the lineage split into Domain, Phylum,
#'   Class, Order, Family, and Genus columns.
create_taxa_table <- function(rute_OTU_table) {
  OTU_table <- read.table(rute_OTU_table,
                              header = TRUE, row.names = 1, check.names = FALSE,
                              sep = "\t")

  TAXA <- OTU_table %>%
    rownames_to_column(var = "Taxa") %>%
    mutate(seqs = rowSums(across(-1))) %>%
    arrange(desc(seqs)) %>%
    separate(Taxa, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
             sep = ";", extra = "merge", remove = FALSE) %>%
    mutate(OTU_ID = paste0("Otu_", sprintf("%04d", row_number()))) %>%
    select(OTU_ID, Taxa, Domain, Phylum, Class, Order, Family, Genus)

  return(TAXA)
}

markers <- c("16S", "18S", "ITS")

for (marker in markers) {
  input_path <- file.path(work_dir, paste0(marker, "_SINKS/tables/rarefacted_table_", marker, ".tsv"))
  output_path <- file.path(work_dir, paste0(marker, "_SINKS/tables/TAXA_", marker, ".tsv"))

  taxa_table <- create_taxa_table(input_path)

  write.table(taxa_table,
              output_path,
              sep = "\t",
              row.names = FALSE)
}
