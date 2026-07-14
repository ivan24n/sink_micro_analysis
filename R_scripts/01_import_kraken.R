##############################################################################
# 01_import_kraken.R
#
# Defines import_kraken(), which reads a list of Kraken2 report files
# (standard 6-column Kraken2 report format) and combines them into a single
# long-format table with full lineage strings reconstructed from the
# report's indentation-based taxonomy tree.
#
# Part of the analysis pipeline for:
#   "Spatiotemporal dynamics of multi-kingdom microbial communities in
#    hospital sinks"
##############################################################################

#' Import and combine Kraken2 reports
#'
#' Reads a set of Kraken2 report files and returns a single combined table
#' with, for every taxon in every sample, its relative abundance, read
#' counts, taxonomic rank, and full lineage (by name, by rank-prefixed name,
#' and by NCBI taxon ID) reconstructed from the report's indentation.
#'
#' @param file_list Character vector of paths to Kraken2 report files
#'   (`*_kraken2_report.txt`, tab-separated, no header: proportion, reads,
#'   assigned reads, rank code, NCBI taxon ID, indented name).
#' @param threads Number of parallel workers used to read the files
#'   (via doParallel/foreach).
#'
#' @return A tibble combining all input reports, with one row per taxon per
#'   sample and columns: Prop, Reads, Assigned, TaxRank, NCBI, Name, Sample,
#'   depth, Name2, Name3, Fullname, Fullname2, FullNCBI.
#'
#' @export
import_kraken <- function(file_list, threads)
{
  require(doParallel)

  res <- data.frame() %>% as_tibble

  cl <- makeCluster(threads)
  registerDoParallel(cl)
  res <- foreach(f = file_list, .combine = rbind) %dopar% {
    {
      table <- read.delim(f, sep = "\t", header = F, stringsAsFactors = F)
      colnames(table) <- c("Prop", "Reads", "Assigned", "TaxRank", "NCBI", "Name")
      table$Sample <- basename(f)

      # Taxonomic depth is inferred from the number of leading two-space
      # indents in the report's Name column.
      table$depth <- nchar(gsub("\\S.*", "", table$Name)) / 2
      table$Name2 <- gsub("^ *", "", table$Name)
      table$depth <- table$depth + 1
      table$Name3 <- paste0(table$TaxRank, "__", table$Name2)

      # Reconstruct the full lineage for each row by tracking the last seen
      # taxon at each depth level (name, rank-prefixed name, and NCBI ID).
      last <- c()
      for (i in 1:nrow(table)) {
        last[table$depth[i]] <- table$Name2[i]
        table$Fullname[i] <- paste(last[1:table$depth[i]], sep = ";", collapse = ";")
      }
      last <- c()
      for (i in 1:nrow(table)) {
        last[table$depth[i]] <- table$Name3[i]
        table$Fullname2[i] <- paste(last[1:table$depth[i]], sep = ";", collapse = ";")
      }
      last <- c()
      for (i in 1:nrow(table)) {
        last[table$depth[i]] <- table$NCBI[i]
        table$FullNCBI[i] <- paste(last[1:table$depth[i]], sep = ";", collapse = ";")
      }

      res <- rbind(res, table)
    }
  }
  res <- res %>% distinct()
  stopCluster(cl)
  return(res)
}
