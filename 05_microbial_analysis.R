##############################################################################
# 05_microbial_analysis.R
#
# Downstream statistical analysis and figure generation for 16S rRNA
# (bacteria/archaea), 18S rRNA (eukaryotes), and ITS (fungi) amplicon data.
#
# Covers: relative-abundance/composition plots, core microbiome, alpha- and
# beta-diversity (by area and by sampling period), and pairwise PERMANOVA
# (Bray-Curtis and Jaccard for 16S/ITS; Jaccard only for 18S, per study
# design).
#
# Part of the analysis pipeline for:
#   "Spatiotemporal dynamics of multi-kingdom microbial communities in
#    hospital sinks"
#
# Input: rarefied OTU/ASV tables, taxonomy tables, and metadata produced by
# 01_import_kraken.R -> 02_build_taxonomy_tables.R -> 03_rarefaction.R ->
# 04_build_taxa_tables.R.
#
# NOTE: this script is intentionally left with its original repeated
# per-marker (16S/18S/ITS) blocks rather than being refactored into
# loops/functions, to avoid changing analysis results during the
# publication cleanup pass. See the repository README for a proposed
# follow-up refactor.
##############################################################################

####LOAD LIBRARY####
library(tidyverse)
library(ggtext)
library(vegan)
library(microbiome)
library(RColorBrewer)
library(patchwork)
library(pairwiseAdonis)
library(ggpubr)
library(ggvenn)
library(phyloseq)
library(fantaxtic)
library(scales)

####LOAD DATA####
set.seed(24111997)

# Set this to the root of the project directory, i.e. the folder that
# contains the 16S_SINKS/, 18S_SINKS/, and ITS_SINKS/ subfolders.
work_dir <- "path/to/project"
setwd(work_dir)

#####Palette colors#####
area_colors <- c("#3f88c5", "#E57A44", "#70A520")

period_colors <- c("#5d2e8c","#C6AC18")

comp_colors <- c("#009292FF", "#004949FF", "#006DDBFF", "#B6DBFFFF",
                 "#6DB6FFFF", "#F1985D", "#FF7000FF", "#920000FF",
                 "#DC5D5D", "#B87EF2","#490092FF","#755B8E") 

#Metadata
metadata_16S <- read.table("16S_SINKS/tables/rare_meta_table_16S.tsv",
                           header = TRUE, check.names = FALSE) 

metadata_18S <- read.table("18S_SINKS/tables/rare_meta_table_18S.tsv",
                           header = TRUE, check.names = FALSE) 

metadata_ITS <- read.table("ITS_SINKS/tables/rare_meta_table_ITS.tsv",
                           header = TRUE, check.names = FALSE) 

#Taxonomy
taxa_16S <- read_tsv("16S_SINKS/tables/TAXA_fullcurated_16S.tsv")
taxa_18S <- read_tsv("18S_SINKS/tables/TAXA_fullcurated_18S_test.tsv")
taxa_ITS <- read_tsv("ITS_SINKS/tables/TAXA_fullcurated_ITS.tsv") 

#Otu tables
OTU_table_16S <- read.table("16S_SINKS/tables/otutable_fullcurated_16S.tsv", 
                            header = TRUE, check.names = FALSE,
                            sep = "\t")

OTU_table_18S <- read.table("18S_SINKS/tables/otutable_fullcurated_18S_test.tsv", 
                            header = TRUE, check.names = FALSE,
                            sep = "\t")

OTU_table_ITS <- read.table("ITS_SINKS/tables/otutable_fullcurated_ITS.tsv", 
                            header = TRUE, check.names = FALSE,
                            sep = "\t")

#DF FULL DATA
df_16S <- OTU_table_16S %>% 
  pivot_longer(-Taxa, names_to = "Sample.id", values_to = "Counts") %>% 
  inner_join(metadata_16S, ., by = "Sample.id") %>% 
  inner_join(., taxa_16S, by = "Taxa") %>% 
  group_by(Sample.id) %>% 
  mutate(rel_abund = 100*Counts/sum(Counts)) %>% 
  ungroup() %>%
  select(-Taxa) %>% 
  pivot_longer(cols = c("Domain", "Phylum", "Class", "Order", "Family", "Genus",
                        "OTU_ID"),
               names_to = "Tax_Rank",
               values_to = "Taxon") %>% 
  as_tibble()

sample_depth_16S <- sum(OTU_table_16S$`SMPL16S_02-22_10_C_M`)
write.table(df_16S, "16S_SINKS/tables/full_df_16S.tsv",
            row.names = FALSE)

df_18S <- OTU_table_18S %>% 
  pivot_longer(-Taxa, names_to = "Sample.id", values_to = "Counts") %>% 
  inner_join(metadata_18S, ., by = "Sample.id") %>% 
  inner_join(., taxa_18S, by = "Taxa") %>% 
  group_by(Sample.id) %>% 
  mutate(rel_abund = 100*Counts/sum(Counts)) %>% 
  ungroup() %>%
  dplyr::select(-Taxa) %>% 
  pivot_longer(cols = c("Domain", "Phylum", "Class", "Order", "Family","Genus",
                        "OTU_ID"),
               names_to = "Tax_Rank",
               values_to = "Taxon") %>% 
  as_tibble()

sample_depth_18S <- sum(OTU_table_18S$`SMPL18S_02-22_11_C_M`)
write.table(df_18S, "18S_SINKS/tables/full_df_18S.tsv",
            row.names = FALSE)

df_ITS <- OTU_table_ITS %>% 
  pivot_longer(-Taxa, names_to = "Sample.id", values_to = "Counts") %>% 
  inner_join(metadata_ITS, ., by = "Sample.id") %>% 
  inner_join(., taxa_ITS, by = "Taxa") %>% 
  group_by(Sample.id) %>% 
  mutate(rel_abund = 100*Counts/sum(Counts)) %>% 
  ungroup() %>%
  dplyr::select(-Taxa) %>% 
  pivot_longer(cols = c("Domain", "Phylum", "Class", "Order", "Family", "Genus",
                        "OTU_ID"),
               names_to = "Tax_Rank",
               values_to = "Taxon") %>% 
  as_tibble()

sample_depth_ITS <- sum(OTU_table_ITS$`SMPL-ITS_03-23_1_A`)
write.table(df_ITS, "ITS_SINKS/tables/full_df_ITS.tsv",
            row.names = FALSE)

#### COMPOSITION ####

#Table_filter
otu_counts_filt <- function(otu_table, tax_table, otu_col = "Taxa", 
                            genus_col = "Genus", min_prevalence = 0, 
                            min_abundance_pct = 0) {
  
  combined_data <- otu_table %>%
    inner_join(tax_table %>% select(all_of(c(otu_col, genus_col))), 
               by = setNames(otu_col, otu_col)) %>%
    filter(!is.na(.data[[genus_col]]))
  
  genus_counts <- combined_data %>%
    group_by(.data[[genus_col]]) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  
  count_matrix <- as.matrix(genus_counts %>% select(where(is.numeric)))
  total_study_reads <- sum(count_matrix, na.rm = TRUE)
  prevalence_vec <- rowSums(count_matrix > 0) / ncol(count_matrix)
  global_abundance_vec <- (rowSums(count_matrix) / total_study_reads) * 100
  
  keep_genus <- genus_counts[[genus_col]][prevalence_vec >= min_prevalence & 
                                            global_abundance_vec >= min_abundance_pct]
  
  final_counts <- genus_counts %>%
    filter(.data[[genus_col]] %in% keep_genus) %>%
    rename(OTU_ID = !!sym(genus_col))
  
  return(final_counts)
}

#Relative abundance table
otu_reltab_creator <- function(otu_table, tax_table, otu_col = "Taxa", 
                               genus_col = "Genus", min_prevalence = 0, 
                               min_abundance_pct = 0) {
  
  combined_data <- otu_table %>%
    inner_join(tax_table %>% select(all_of(c(otu_col, genus_col))), 
               by = setNames(otu_col, otu_col)) %>%
    filter(!is.na(.data[[genus_col]]))
  
  genus_counts <- combined_data %>%
    group_by(.data[[genus_col]]) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  
  count_matrix <- as.matrix(genus_counts %>% select(where(is.numeric)))
  total_study_reads <- sum(count_matrix, na.rm = TRUE)
  prevalence_vec <- rowSums(count_matrix > 0) / ncol(count_matrix)
  global_abundance_vec <- (rowSums(count_matrix) / total_study_reads) * 100
  
  keep_genus <- genus_counts[[genus_col]][prevalence_vec >= min_prevalence & 
                                            global_abundance_vec >= min_abundance_pct]
  
  genus_filtered <- genus_counts %>%
    filter(.data[[genus_col]] %in% keep_genus)
  
  final_table <- genus_filtered %>%
    mutate(across(where(is.numeric), ~ {
      s <- sum(.x, na.rm = TRUE)
      if(s > 0) .x / s else 0
    })) %>%
    rename(OTU_ID = !!sym(genus_col))
  
  final_table <- final_table %>%
    mutate(mean_rel = rowMeans(select(., where(is.numeric)))) %>%
    arrange(desc(mean_rel)) %>%
    select(-mean_rel)
  
  return(final_table)
}

#Core function
core_microbiome <- function(abundance_table, prevalence_threshold = 0.7) {
  abundance_norm <- abundance_table
  abundance_norm[-1] <- sweep(
    abundance_norm[-1],
    2,
    colSums(abundance_norm[-1]),
    FUN = "/"
  ) * 100
  
  abundance_norm %>%
    rowwise() %>%
    mutate(
      Prevalence = sum(c_across(-Taxa) > 0) / (ncol(.) - 1),
      Mean_Abundance = mean(c_across(-Taxa))
    ) %>%
    ungroup() %>%
    filter(Prevalence >= prevalence_threshold) %>%
    select(Taxa, Prevalence, Mean_Abundance) %>%
    arrange(desc(Mean_Abundance)) %>%
    separate(
      Taxa,
      c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
      sep = ";"
    )
}

#Nested composition
plot_nested_abundance <- function(metadata, OTU_table, taxa, plot_title) {
  
  metadata_data <- metadata %>% column_to_rownames("Sample.id")
  OTU_table_data <- OTU_table %>% column_to_rownames("Taxa")
  taxa_data <- taxa %>% column_to_rownames("Taxa") %>% select(-OTU_ID) %>% 
    as.matrix()
  OTU_table_data <- OTU_table_data[,rownames(metadata_data)]
  OTU_table_data <- OTU_table_data[rownames(taxa_data),]
  
  
  OTU = otu_table(OTU_table_data, taxa_are_rows = TRUE)
  TAX = tax_table(taxa_data)
  sampledata = sample_data(metadata_data)
  
  pseq <- phyloseq(OTU, TAX, sampledata)
  
  
  top_nested <- nested_top_taxa(pseq,
                                top_tax_level = "Phylum",
                                nested_tax_level = "Family",
                                n_top_taxa = 3, 
                                n_nested_taxa = 3)
  
  plot <- plot_nested_bar(ps_obj = top_nested$ps_obj,
                  top_level = "Phylum",
                  nested_level = "Family") +
    scale_y_continuous(labels = label_percent(accuracy = 1),
                       expand = c(0, 0)) + 
    theme(
      axis.text.x = element_blank(),    
      axis.ticks.x = element_blank(),   
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
    ) +
    labs(y = "Relative Abundance (%)", x = "Samples",
         title = plot_title) +
    facet_wrap(~ Box, scales = "free_x", nrow = 1)
  
  return(plot)
}
  
#Composition plot function by group
composition_analysis <- function(df, tax_rank_selected,
                                 grouping_variable,
                                 minimum_percentage) {
  
  rel_abund <- df %>% 
    filter(Tax_Rank == tax_rank_selected) %>% 
    group_by(!!sym(grouping_variable), Taxon, Sample.id) %>% 
    summarise(rel_abund = sum(rel_abund), .groups = "drop") %>% 
    group_by(!!sym(grouping_variable), Taxon) %>% 
    summarise(mean_rel_abund = mean(rel_abund), .groups = "drop")
  
  pool <- rel_abund %>% 
    group_by(Taxon) %>% 
    summarise(pool = max(mean_rel_abund) < minimum_percentage,
              mean = mean(mean_rel_abund),
              .groups = "drop") 
  
  plot_df <- inner_join(rel_abund, pool, by = "Taxon") %>% 
    mutate(Taxon = if_else(pool, "Other", Taxon)) %>%
    group_by(!!sym(grouping_variable), Taxon) %>%
    summarise(mean_rel_abund = sum(mean_rel_abund),
              mean = min(mean),
              .groups = "drop") %>%
    group_by(Taxon) %>%
    mutate(global_abund = mean(mean_rel_abund)) %>%
    ungroup() %>%
    mutate(Taxon = if_else(!(Taxon %in% 
                               names(sort(tapply(global_abund, Taxon, mean), decreasing = TRUE)[1:15])) &
                             Taxon != "Other",
                           "Other", Taxon)) %>%
    group_by(!!sym(grouping_variable), Taxon) %>%
    summarise(mean_rel_abund = sum(mean_rel_abund),
              .groups = "drop") %>%
    mutate(Taxon = factor(Taxon)) %>%
    mutate(Taxon = fct_relevel(Taxon, "Other",
                               after = length(unique(Taxon)) %/% 1))
  
  comp_palette <- c("#69EBD0","#009292FF", "#004949FF","#006DDBFF", "#B6DBFFFF",
                    "#6DB6FFFF", "#F1985D", "#FF7000FF","#DD1C1A", "#920000FF",
                    "#DC5D5D", "#D7B8F3","#490092FF","#9546E3","#635380") 
  
  taxon_levels <- levels(plot_df$Taxon)
  legend_breaks <- c(setdiff(taxon_levels, "Other"), "Other") 
  legend_values <- setNames(comp_palette[1:(length(legend_breaks) - 1)],
                            legend_breaks[-length(legend_breaks)])
  legend_values["Other"] <- "#e5e5e5"  
  
  plot <- ggplot(plot_df, aes(x = !!sym(grouping_variable),
                              y = mean_rel_abund, fill = Taxon)) +
    geom_col(width = 0.75) +
    scale_x_discrete(expand = c(0.25, 0)) +
    scale_fill_manual(name = str_to_title(tax_rank_selected),
                      breaks = legend_breaks,  
                      values = legend_values) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL,
         y = "Relative Abundance (%)",
         title = "Composition") +
    theme_classic2() +
    theme(legend.text = element_text(face = "italic", size = 11),
          legend.title = element_text(size = 12),
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 10),
          legend.key.spacing.y = unit(2, "pt"))
  
  return(plot)
}

#Venn diagram function
create_venn_diagram <- function(df, group_col, color_palette) {
  otu_list <- df %>%
    filter(Tax_Rank == "Genus" & Counts > 0) %>%
    select(!!sym(group_col), Taxon) %>%
    group_by(!!sym(group_col)) %>%
    summarise(otus = list(unique(Taxon))) %>%
    pull(otus, name = !!sym(group_col))
  
  diagram <- ggvenn(
    otu_list,
    show_counts = TRUE,
    show_percentage = TRUE,
    fill_color = color_palette,
    fill_alpha = 0.5,
    stroke_size = 0,
    set_name_size = 4,
    text_size = 3
  ) +
    labs(title = "Genus Overlap") +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      text = element_text(color = "black")
    )
  
  return(diagram)
}

##### GENERAL COMPOSITION #####

#16S
otu_reltab_16S <- otu_reltab_creator(OTU_table_16S, taxa_16S)

write.table(otu_reltab_16S, "Correlations/SPARCC/SparCC-master/RNA16S/otu_table_relative_16S_allsamps.txt", sep ="\t",
            row.names = FALSE)

comp_box_16S <- plot_nested_abundance(metadata_16S,
                                      OTU_table_16S,
                                      taxa_16S,
                                      "Microbial Community Structure per Room: 16S rRNA Amplicon Data")

core_16 <- core_microbiome(OTU_table_16S, 0.8)


write.table(core_16, "16S_SINKS/tables/core_tax_16S.tsv",
            row.names = FALSE)

ggsave("16S_SINKS/figures/comp_box_16S.svg", plot = comp_box_16S,
       width = 10,
       height = 10,
       units = "in",
       dpi = 100)


#18S
comp_box_18S <- plot_nested_abundance(metadata_18S,
                                      OTU_table_18S,
                                      taxa_18S,
                                      "Microbial Community Structure per Room: 18S rRNA Amplicon Data")

core_18 <- core_microbiome(OTU_table_18S, 0.5)

write.table(core_18, "18S_SINKS/tables/core_tax_18S.tsv",
            row.name = FALSE)

ggsave("18S_SINKS/figures/comp_box_18S.svg", plot = comp_box_18S,
       width = 10,
       height = 10,
       units = "in",
       dpi = 100)

#ITS
otu_reltab_ITS <- otu_reltab_creator(OTU_table_ITS, taxa_ITS)

write.table(otu_reltab_ITS, "Correlations/SPARCC/SparCC-master/ITS/otu_table_relative_ITS_allsamps.txt", sep ="\t",
            row.names = FALSE)

comp_box_ITS <- plot_nested_abundance(metadata_ITS,
                                      OTU_table_ITS,
                                      taxa_ITS,
                                      "Fungal Community Structure per Room: ITS Amplicon Sequencing Data")

core_ITS <- core_microbiome(OTU_table_ITS, 0.7)

write.table(core_ITS, "ITS_SINKS/tables/core_tax_ITS.tsv",
            row.names = FALSE)

ggsave("ITS_SINKS/figures/comp_box_ITS.svg", plot = comp_box_ITS,
       width = 10,
       height = 10,
       units = "in",
       dpi = 100)

##### Area #####

#16S
venn_area_16S <- create_venn_diagram(df_16S, "Area", area_colors)


comp_area_16S <- composition_analysis(df = df_16S, tax_rank_selected = "Genus",
                                      grouping_variable = "Area",
                                      minimum_percentage = 2)
#18S
venn_area_18S <- create_venn_diagram(df_18S, "Area", area_colors)


comp_area_18S <- composition_analysis(df = df_18S, tax_rank_selected = "Genus",
                                      grouping_variable = "Area",
                                      minimum_percentage = 5)
#ITS
venn_area_ITS <- create_venn_diagram(df_ITS, "Area", area_colors)


comp_area_ITS <- composition_analysis(df = df_ITS, tax_rank_selected = "Genus",
                                      grouping_variable = "Area",
                                      minimum_percentage = 2)

##### Unit #####

#16S
venn_unit_16S <- create_venn_diagram(df_16S, "Period", period_colors)



comp_unit_16S <- composition_analysis(df = df_16S, tax_rank_selected = "Genus",
                                      grouping_variable = "Period",
                                      minimum_percentage = 2)

#18S
venn_unit_18S <- create_venn_diagram(df_18S, "Period", period_colors)



comp_unit_18S <- composition_analysis(df = df_18S, tax_rank_selected = "Genus",
                                      grouping_variable = "Period",
                                      minimum_percentage = 2)

#### ALPHA-DIVERSITY ####

#16S
short_df_16S <- OTU_table_16S %>%
  inner_join(taxa_16S, by = "Taxa") %>% 
  select(OTU_ID, where(is.numeric)) %>% 
  column_to_rownames("OTU_ID") %>% 
  t()

order <- match(metadata_16S$Sample.id, rownames(short_df_16S))
short_df_16S <- short_df_16S[order,] %>% as.matrix()

div_16S <- microbiome::alpha(t(short_df_16S)) %>%
  rownames_to_column(var = "Sample.id") %>%
  inner_join(., metadata_16S, by = "Sample.id") %>% 
  rename(
    Sp_obs = observed,
    Shannon = diversity_shannon,
    Simpson = diversity_gini_simpson,
    InvSimpson = diversity_inverse_simpson,
    Chao1 = chao1,
    Berger_Parker = dominance_dbp
  )
write.table(short_df_16S, "16S_SINKS/tables/otu_matrix.tsv",
            col.names = NA, row.names = TRUE)

#18S
short_df_18S <- OTU_table_18S %>%
  inner_join(taxa_18S, by = "Taxa") %>% 
  select(OTU_ID, where(is.numeric)) %>% 
  column_to_rownames("OTU_ID") %>% 
  t()

order <- match(metadata_18S$Sample.id, rownames(short_df_18S))
short_df_18S <- short_df_18S[order,] %>% as.matrix()

div_18S <- microbiome::alpha(t(short_df_18S)) %>%
  rownames_to_column(var = "Sample.id") %>%
  inner_join(., metadata_18S, by = "Sample.id") %>% 
  rename(
    Sp_obs = observed,
    Shannon = diversity_shannon,
    Simpson = diversity_gini_simpson,
    InvSimpson = diversity_inverse_simpson,
    Chao1 = chao1,
    Berger_Parker = dominance_dbp
  )

write.table(short_df_18S, "18S_SINKS/tables/otu_matrix.tsv", row.names = TRUE,
            col.names = NA)

#ITS
short_df_ITS <- OTU_table_ITS %>%
  inner_join(taxa_ITS, by = "Taxa") %>% 
  select(OTU_ID, where(is.numeric)) %>% 
  column_to_rownames("OTU_ID") %>% 
  t()

order <- match(metadata_ITS$Sample.id, rownames(short_df_ITS))
short_df_ITS <- short_df_ITS[order,] %>% as.matrix()

div_ITS <- microbiome::alpha(t(short_df_ITS)) %>%
  rownames_to_column(var = "Sample.id") %>%
  inner_join(., metadata_ITS, by = "Sample.id") %>% 
  rename(
    Sp_obs = observed,
    Shannon = diversity_shannon,
    Simpson = diversity_gini_simpson,
    InvSimpson = diversity_inverse_simpson,
    Chao1 = chao1,
    Berger_Parker = dominance_dbp
  )

write.table(short_df_ITS, "ITS_SINKS/tables/otu_matrix.tsv", row.names = TRUE,
            col.names = NA)
              
alpha_index <- c("Sp_obs", "Chao1", "Shannon", "Simpson", "InvSimpson",
                 "Berger_Parker") 

#Alpha-Plot function
plot_alpha_diversity <- function(data, x_col, y_col,
                               title = "Diversity Plot",
                               color_palette = NULL) {

  data[[x_col]] <- as.factor(data[[x_col]])
  
  unique_categories <- levels(data[[x_col]])
  

  if (is.null(color_palette)) {
    color_palette <- rainbow(length(unique_categories))
  }

  p <- ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]],
                        fill = .data[[x_col]], color = .data[[x_col]])) +
    geom_violin(trim = FALSE, alpha = 0.4, scale = "width") +
    geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.4,
                 show.legend = FALSE) +
    geom_jitter(position = position_jitter(0.2), size = 0.65,
                show.legend = FALSE) +
    labs(title = title,
         y = "Index value",
         x = NULL) +
    scale_x_discrete(expand = expansion(0.2, 0.2)) +
    scale_fill_manual(values = color_palette) +
    scale_color_manual(values = color_palette) +
    theme_light() +
    theme(legend.position = "none",
          legend.text = element_text(size=12),
          legend.title = element_text(size=13),
          legend.spacing.y = unit(0.5, "cm"),
          plot.title = element_text(size = 13, hjust = 0.5),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(10,10,10,10))
  
  return(p)
}

#Wilcoxon significance table function
pairwise_wilcox_table <- function(data, value_cols, group_col, 
                                  p_adjust = "bonferroni") {
  
  data[[group_col]] <- as.factor(data[[group_col]])
  
  result_list <- list()
  
  for (value_col in value_cols) {

    wilcox_results <- pairwise.wilcox.test(
      data[[value_col]], data[[group_col]], 
      p.adjust.method = p_adjust, exact = FALSE
    )

    p_values_df <- as.data.frame(as.table(wilcox_results$p.value))

    p_values_df$Var1 <- as.character(p_values_df$Var1)
    p_values_df$Var2 <- as.character(p_values_df$Var2)

    p_values_df$Metric <- value_col

    p_values_df$Significance <- ifelse(
      p_values_df$Freq < 0.001, "***",
      ifelse(p_values_df$Freq < 0.01, "**",
             ifelse(p_values_df$Freq < 0.05, "*", "ns")
      )
    )
    
    colnames(p_values_df) <- c("Group1", "Group2", 
                               "P_Value", "Metric", "Significance")
    
    p_values_df <- p_values_df[p_values_df$Group1 != p_values_df$Group2, ]
    
    result_list[[value_col]] <- p_values_df
  }
  
  final_table <- do.call(rbind, result_list)
  
  final_table <- final_table[, c("Metric", "Group1", 
                                 "Group2", "P_Value", "Significance")]
  
  return(final_table)
}


#####Area#####

#16S
p_shannon_area_16S <- plot_alpha_diversity(div_16S, "Area", "Shannon",
                   color_palette = area_colors,
                   title = "Shannon")

p_chao1_area_16S <- plot_alpha_diversity(div_16S, "Area", "Chao1",
                   color_palette = area_colors,
                   title = "Chao-1")


p_simpson_area_16S <- plot_alpha_diversity(div_16S, "Area", "Simpson",
                   color_palette = area_colors,
                   title = "Simpson")

p_invsimpson_area_16S <- plot_alpha_diversity(div_16S, "Area", "InvSimpson",
                   color_palette = area_colors,
                   title = "Inv-Simpson")

p_sobs_area_16S <- plot_alpha_diversity(div_16S, "Area", "Sp_obs",
                    color_palette = area_colors,
                    title = "Observed Species")

p_bdp_area_16S <- plot_alpha_diversity(div_16S, "Area", "Berger_Parker",
                    color_palette = area_colors,
                    title = "Berger-Parker")


alpha_area_sign_16S <- pairwise_wilcox_table(div_16S, alpha_index, "Area")

write.table(alpha_area_sign_16S, "16S_SINKS/tables/area/a-sign_16S.tsv",
            row.names = FALSE)

write.table(div_16S, "16S_SINKS/tables/a-div_16S.tsv",
            row.names = FALSE)

#18S
p_shannon_area_18S <- plot_alpha_diversity(div_18S, "Area", "Shannon",
                                           color_palette = area_colors,
                                           title = "Shannon")

p_chao1_area_18S <- plot_alpha_diversity(div_18S, "Area", "Chao1",
                                         color_palette = area_colors,
                                         title = "Chao-1")

p_simpson_area_18S <- plot_alpha_diversity(div_18S, "Area", "Simpson",
                                           color_palette = area_colors,
                                           title = "Simpson")

p_invsimpson_area_18S <- plot_alpha_diversity(div_18S, "Area", "InvSimpson",
                                              color_palette = area_colors,
                                              title = "Inv-Simpson")

p_sobs_area_18S <- plot_alpha_diversity(div_18S, "Area", "Sp_obs",
                                        color_palette = area_colors,
                                        title = "Observed Species")

p_bdp_area_18S <- plot_alpha_diversity(div_18S, "Area", "Berger_Parker",
                                       color_palette = area_colors,
                                       title = "Berger-Parker")


alpha_area_sign_18S <- pairwise_wilcox_table(div_18S, alpha_index, "Area")

write.table(alpha_area_sign_18S, "18S_SINKS/tables/area/a-sign_18S.tsv",
            row.names = FALSE)

write.table(div_18S, "18S_SINKS/tables/a-div_18S.tsv",
            row.names = FALSE)

#ITS
p_shannon_area_ITS <- plot_alpha_diversity(div_ITS, "Area", "Shannon",
                                           color_palette = area_colors,
                                           title = "Shannon")

p_chao1_area_ITS <- plot_alpha_diversity(div_ITS, "Area", "Chao1",
                                         color_palette = area_colors,
                                         title = "Chao-1")

p_simpson_area_ITS <- plot_alpha_diversity(div_ITS, "Area", "Simpson",
                                           color_palette = area_colors,
                                           title = "Simpson")

p_invsimpson_area_ITS <- plot_alpha_diversity(div_ITS, "Area", "InvSimpson",
                                              color_palette = area_colors,
                                              title = "Inv-Simpson")

p_sobs_area_ITS <- plot_alpha_diversity(div_ITS, "Area", "Sp_obs",
                                        color_palette = area_colors,
                                        title = "Observed Species")

p_bdp_area_ITS <- plot_alpha_diversity(div_ITS, "Area", "Berger_Parker",
                                       color_palette = area_colors,
                                       title = "Berger-Parker")


alpha_area_sign_ITS <- pairwise_wilcox_table(div_ITS, alpha_index, "Area")

write.table(alpha_area_sign_ITS, "ITS_SINKS/tables/area/a-sign_ITS.tsv",
            row.names = FALSE)

write.table(div_ITS, "ITS_SINKS/tables/a-div_ITS.tsv",
            row.names = FALSE)

#####Unit#####

#16S
p_shannon_unit_16S <- plot_alpha_diversity(div_16S, "Period", "Shannon",
                                    title = "Shannon",
                                    color_palette =period_colors)

p_chao1_unit_16S <- plot_alpha_diversity(div_16S, "Period", "Chao1",
                                    title = "Chao-1",
                                   color_palette =period_colors)

p_simpson_unit_16S <- plot_alpha_diversity(div_16S, "Period", "Simpson",
                                    title = "Simpson",
                                    color_palette =period_colors)

p_invsimpson_unit_16S <- plot_alpha_diversity(div_16S, "Period", "InvSimpson",
                                    title = "Inverted Simpson",
                                    color_palette =period_colors)

p_sobs_unit_16S <- plot_alpha_diversity(div_16S, "Period", "Sp_obs",
                                    title = "Observed Species",
                                    color_palette =period_colors)

p_bdp_unit_16S <- plot_alpha_diversity(div_16S, "Period", "Berger_Parker",
                                    title = "Berger-Parker",
                                    color_palette =period_colors)


alpha_unit_sign_16S <- pairwise_wilcox_table(div_16S, alpha_index, "Period")

write.table(alpha_unit_sign_16S, "16S_SINKS/tables/period/a-sign_16S.tsv",
            row.names = FALSE)



#18S
p_shannon_unit_18S <- plot_alpha_diversity(div_18S, "Period", "Shannon",
                                           title = "Shannon",
                                           color_palette =period_colors)

p_chao1_unit_18S <- plot_alpha_diversity(div_18S, "Period", "Chao1",
                                         title = "Chao-1",
                                         color_palette =period_colors)

p_simpson_unit_18S <- plot_alpha_diversity(div_18S, "Period", "Simpson",
                                           title = "Simpson",
                                           color_palette =period_colors)

p_invsimpson_unit_18S <- plot_alpha_diversity(div_18S, "Period", "InvSimpson",
                                              title = "Inverted Simpson",
                                              color_palette =period_colors)

p_sobs_unit_18S <- plot_alpha_diversity(div_18S, "Period", "Sp_obs",
                                        title = "Observed Species",
                                        color_palette =period_colors)

p_bdp_unit_18S <- plot_alpha_diversity(div_18S, "Period", "Berger_Parker",
                                       title = "Berger-Parker",
                                       color_palette =period_colors)


alpha_unit_sign_18S <- pairwise_wilcox_table(div_18S, alpha_index, "Period")

write.table(alpha_unit_sign_18S, "18S_SINKS/tables/period/a-sign_18S.tsv",
            row.names = FALSE)

#### BETA-DIVERSITY ####

#Beta NMDS plot
plot_beta_diversity <- function(abundance_matrix, metadata, 
                                grouping_var, method = "bray", 
                                sample_depth = NULL, color_palette = NULL) {
  
  if (tolower(method) %in% c("jacc", "jaccard")) {
    dist_matrix <- vegdist(abundance_matrix, method = "jaccard", 
                           binary = TRUE, sample = sample_depth)
    method_label <- "Jaccard"
  } else {
    dist_matrix <- vegdist(abundance_matrix, method = "bray", 
                           sample = sample_depth)
    method_label <- "Bray-Curtis"
  }
  
  NMDS_result <- metaMDS(dist_matrix) %>% 
    scores() %>% 
    as_tibble(rownames = "Sample.id")
  
  NMDS_meta <- NMDS_result %>% 
    inner_join(metadata, by = "Sample.id")
  
  centroid <- NMDS_meta %>% 
    group_by(.data[[grouping_var]]) %>% 
    summarise(axis1 = mean(NMDS1), axis2 = mean(NMDS2))
  
  if (is.null(color_palette)) {
    color_palette <- scales::hue_pal()(length(unique(NMDS_meta[[grouping_var]])))
    names(color_palette) <- unique(NMDS_meta[[grouping_var]])
  }
  
  ggplot(NMDS_meta, aes(x = NMDS1, y = NMDS2, 
                        color = .data[[grouping_var]], 
                        fill = .data[[grouping_var]])) +
    stat_ellipse(geom = "polygon", level = 0.8, alpha = 0.2, show.legend = FALSE) +
    geom_point(size = 1.5) +
    geom_point(data = centroid, aes(x = axis1, y = axis2), 
               shape = 22, size = 4, color = "black", show.legend = FALSE) + 
    coord_fixed() +
    labs(
      title = method_label,
      x = "NMDS1",
      y = "NMDS2"
    ) +
    theme_light() +
    theme(
      plot.title = element_text(size = 13, hjust = 0.5),
      legend.position = "none",
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    ) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette)
}

#Saver adonis pairwise
save_pairwise_adonis_tsv <- function(pairwise_result,
                                     base_filename = "pairwise_adonis_results") {
  if (!is.list(pairwise_result)) {
    stop("The input must be a list (result of pairwise.adonis2).")
  }
  list_df_with_comparison <- list()
  for (element_name in names(pairwise_result)) {
    if (element_name == "parent_call") next
    table <- pairwise_result[[element_name]]
    if (is.data.frame(table)) {
      df_with_comparison <- as.data.frame(table) %>%
        tibble::rownames_to_column(var = "Statistic") %>%
        dplyr::mutate(Comparison = element_name)
      list_df_with_comparison[[element_name]] <- df_with_comparison
    }
  }
  combined_results <- dplyr::bind_rows(list_df_with_comparison)
  tsv_filename <- paste0(base_filename, ".tsv")
  write.table(combined_results, file = tsv_filename, sep = "\t",
              quote = FALSE, row.names = FALSE)
}

##### Area #####

#16S
p_bray_area_16S <- plot_beta_diversity(short_df_16S, metadata_16S, "Area",
                                   method = "bray",
                                   color_palette = area_colors)

p_jacc_area_16S <-plot_beta_diversity(short_df_16S, metadata_16S, "Area",
                                  method = "jacc",
                                  color_palette = area_colors)

dist_bray_16S <- vegdist(short_df_16S, sample = sample_depth_16S, method = "bray")

dist_jacc_16S <- vegdist(short_df_16S, sample = sample_depth_16S, method = "jacc", binary = TRUE)

area_bray_adonispw_16S <- pairwise.adonis2(dist_bray_16S ~ Area, metadata_16S,
                                 p.adjust.m = "bonferroni")

area_jacc_adonispw_16S <- pairwise.adonis2(dist_jacc_16S ~ Area, metadata_16S,
                                           p.adjust.m = "bonferroni")

save_pairwise_adonis_tsv(area_bray_adonispw_16S, base_filename = "16S_SINKS/tables/area/adonispw_bray_16S")

save_pairwise_adonis_tsv(area_jacc_adonispw_16S, base_filename = "16S_SINKS/tables/area/adonispw_jacc_16S")


area_bray_adonis_16S <- adonis2(dist_bray_16S ~ Area, metadata_16S, permutations = 999)

area_jacc_adonis_16S <- adonis2(dist_jacc_16S ~ Area, metadata_16S, permutations = 999)

write.table(area_bray_adonis_16S, "16S_SINKS/tables/area/adonis_bray_16S.tsv",
            col.names = NA, row.names = TRUE)

write.table(area_jacc_adonis_16S, "16S_SINKS/tables/area/adonis_jacc_16S.tsv",
            col.names = NA, row.names = TRUE)

#18S
# Note: for the 18S dataset, only Jaccard beta-diversity is computed
# (Bray-Curtis is intentionally omitted per study design).
dist_jacc_18S <- vegdist(short_df_18S, sample = sample_depth_18S, method = "jacc", binary = TRUE)

p_jacc_area_18S <-plot_beta_diversity(short_df_18S, metadata_18S, "Area",
                                      method = "jacc",
                                      color_palette = area_colors)

area_jacc_adonispw_18S <- pairwise.adonis2(dist_jacc_18S ~ Area, metadata_18S,
                                      p.adjust.m = "bonferroni")

save_pairwise_adonis_tsv(area_jacc_adonispw_18S, base_filename = "18S_SINKS/tables/area/adonispw_jacc_18S")

area_jacc_adonis_18S <- adonis2(dist_jacc_18S ~ Area, metadata_18S, permutations = 999)

write.table(area_jacc_adonis_18S, "18S_SINKS/tables/area/adonis_jacc_18S.tsv",
            col.names = NA, row.names = TRUE)


#ITS
p_bray_area_ITS <- plot_beta_diversity(short_df_ITS, metadata_ITS, "Area",
                                       method = "bray",
                                       color_palette = area_colors)

p_jacc_area_ITS <-plot_beta_diversity(short_df_ITS, metadata_ITS, "Area",
                                      method = "jacc",
                                      color_palette = area_colors)

dist_bray_ITS <- vegdist(short_df_ITS, sample = sample_depth_ITS, method = "bray")

dist_jacc_ITS <- vegdist(short_df_ITS, sample = sample_depth_ITS, method = "jacc", binary = TRUE)

area_bray_adonispw_ITS <- pairwise.adonis2(dist_bray_ITS ~ Area, metadata_ITS,
                                      p.adjust.m = "bonferroni")

area_jacc_adonispw_ITS <- pairwise.adonis2(dist_jacc_ITS ~ Area, metadata_ITS,
                                           p.adjust.m = "bonferroni")

save_pairwise_adonis_tsv(area_bray_adonispw_ITS, base_filename = "ITS_SINKS/tables/area/adonispw_bray_ITS")

save_pairwise_adonis_tsv(area_jacc_adonispw_ITS, base_filename = "ITS_SINKS/tables/area/adonispw_jacc_ITS")


area_bray_adonis_ITS <- adonis2(dist_bray_ITS ~ Area, metadata_ITS, permutations = 999)

area_jacc_adonis_ITS <- adonis2(dist_jacc_ITS ~ Area, metadata_ITS, permutations = 999)

write.table(area_bray_adonis_ITS, "ITS_SINKS/tables/area/adonis_bray_ITS.tsv",
            col.names = NA, row.names = TRUE)

write.table(area_jacc_adonis_ITS, "ITS_SINKS/tables/area/adonis_jacc_ITS.tsv",
            col.names = NA, row.names = TRUE)


##### Unit #####

#16S
p_bray_unit_16S <- plot_beta_diversity(short_df_16S, metadata_16S, "Period",
                                   method = "bray",
                                   color_palette = period_colors)

p_jacc_unit_16S <- plot_beta_diversity(short_df_16S, metadata_16S, "Period",
                                   method = "jacc",
                                   color_palette = period_colors)

unit_bray_adonis_16S <- adonis2(dist_bray_16S ~ Period, metadata_16S, permutations = 999)

unit_jacc_adonis_16S <- adonis2(dist_jacc_16S ~ Period, metadata_16S, permutations = 999)

write.table(unit_bray_adonis_16S, "16S_SINKS/tables/period/adonis_bray_16S.tsv",
            col.names = NA, row.names = TRUE)

write.table(unit_jacc_adonis_16S, "16S_SINKS/tables/period/adonis_jacc_16S.tsv",
            col.names = NA, row.names = TRUE)

#18S
# Note: for the 18S dataset, only Jaccard beta-diversity is computed
# (Bray-Curtis is intentionally omitted per study design).
p_jacc_unit_18S <- plot_beta_diversity(short_df_18S, metadata_18S, "Period",
                                       method = "jacc",
                                       color_palette = period_colors)

unit_jacc_adonis_18S <- adonis2(dist_jacc_18S ~ Period, metadata_18S, permutations = 999)

write.table(unit_jacc_adonis_18S, "18S_SINKS/tables/period/adonis_jacc_18S.tsv",
            col.names = NA, row.names = TRUE)

  #### Final figures ####

lay <- "
112233
445566
445566
777766
777766
"

#16S
p_alpha_area_16S <- (p_shannon_area_16S + p_chao1_area_16S + p_bdp_area_16S +
                  plot_layout(guides = "collect", axis_titles = "collect")) +
                  plot_annotation(
                     title = "α-Diversity Indexes",
                     theme = theme(plot.title = element_text(face = "bold",
                                                             hjust = 0.5))) &
                   theme(legend.position = "right",
                         legend.title.position = "top",
                         axis.text.x = element_blank(),
                         legend.key.spacing.y = unit(5, "pt"),
                         legend.key.size = unit(0.25, "cm"),
                         legend.box.spacing = unit(0, "pt")
                         )


p_beta_area_16S <- (p_bray_area_16S + p_jacc_area_16S) +
  plot_layout(guides = "collect", axis_titles = "collect") +
  plot_annotation(
    title = "β-Diversity",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ) &
  theme(
    legend.position = "right",
    legend.title.position = "top",
    axis.text.x = element_blank(),
    legend.key.spacing.y = unit(5, "pt"),
    legend.key.size = unit(0.25, "cm"),
    legend.box.spacing = unit(0, "pt")
  )


p_alpha_unit_16S <- (p_shannon_unit_16S + p_chao1_unit_16S + p_sobs_unit_16S +
                   plot_layout(guides = "collect", axis_titles = "collect")) +
  plot_annotation(
    title = "α-Diversity Indexes",
    theme = theme(plot.title = element_text(face = "bold",
                                            hjust = 0.5))) &
  theme(legend.position = "right",
        legend.title.position = "top",
        axis.text.x = element_blank(),
        legend.key.spacing.y = unit(5, "pt"),
        legend.key.size = unit(0.25, "cm"),
        legend.box.spacing = unit(0, "pt")
  )

p_beta_unit_16S <- (p_bray_unit_16S + p_jacc_unit_16S +
                      plot_layout(guides = "collect", axis_titles = "collect")) +
  plot_annotation(
    title = "β-Diversity",
    theme = theme(plot.title = element_text(face = "bold",
                                            hjust = 0.5))) &
  theme(legend.position = "right",
        legend.title.position = "top",
        axis.text.x = element_blank(),
        legend.key.spacing.y = unit(5, "pt"),
        legend.key.size = unit(0.25, "cm"),
        legend.box.spacing = unit(0, "pt")
  )

p_final_area_16S <- p_alpha_area_16S + p_bray_area_16S + p_jacc_area_16S + 
  comp_area_16S + venn_area_16S +
  plot_layout(design = lay) +
  plot_annotation(title = "Spatial Variation in Diversity and Composition in 16S rRNA",
                  theme = theme(plot.title = element_text(face = "bold",
                                                          hjust = 0.5,
                                                          size = 14)))

p_final_unit_16S <- p_alpha_unit_16S + p_bray_unit_16S + p_jacc_unit_16S + 
  comp_unit_16S + venn_unit_16S +
  plot_layout(design = lay) +
  plot_annotation(title = "Temporal Variation in Diversity and Composition in 16S rRNA",
                  theme = theme(plot.title = element_text(face = "bold",
                                                          hjust = 0.5,
                                                          size = 14)))

ggsave("16S_SINKS/figures/area/plot_16S.svg", plot = p_final_area_16S,
       width = 10,
       height = 10,
       units = "in",
       dpi = 100)

ggsave("16S_SINKS/figures/period/plot_16S.svg", plot = p_final_unit_16S,
       width = 10,
       height = 10,
       units = "in",
       dpi = 100)

#18S
p_alpha_area_18S <- (p_shannon_area_18S + p_chao1_area_18S + p_bdp_area_18S +
                       plot_layout(guides = "collect", axis_titles = "collect")) +
  plot_annotation(
    title = "α-Diversity Indexes",
    theme = theme(plot.title = element_text(face = "bold",
                                            hjust = 0.5))) &
  theme(legend.position = "right",
        legend.title.position = "top",
        axis.text.x = element_blank(),
        legend.key.spacing.y = unit(5, "pt"),
        legend.key.size = unit(0.25, "cm"),
        legend.box.spacing = unit(0, "pt")
  )

p_alpha_unit_18S <- (p_shannon_unit_18S + p_chao1_unit_18S + p_bdp_unit_18S +
                       plot_layout(guides = "collect", axis_titles = "collect")) +
  plot_annotation(
    title = "α-Diversity Indexes",
    theme = theme(plot.title = element_text(face = "bold",
                                            hjust = 0.5))) &
  theme(legend.position = "right",
        legend.title.position = "top",
        axis.text.x = element_blank(),
        legend.key.spacing.y = unit(5, "pt"),
        legend.key.size = unit(0.25, "cm"),
        legend.box.spacing = unit(0, "pt")
  )

# 18S only has Jaccard beta-diversity (see note above), so its composite
# figures use a 6-panel layout (one fewer panel than 16S/ITS, which also
# include a Bray-Curtis panel).
lay_18S <- "
112233
444455
444455
666655
666655
"

p_final_area_18S <- p_alpha_area_18S + p_jacc_area_18S +
  comp_area_18S + venn_area_18S +
  plot_layout(design = lay_18S) +
  plot_annotation(title = "Spatial Variation in Diversity and Composition in 18S rRNA",
                  theme = theme(plot.title = element_text(face = "bold",
                                                          hjust = 0.5,
                                                          size = 14)))

p_final_unit_18S <- p_alpha_unit_18S + p_jacc_unit_18S +
  comp_unit_18S + venn_unit_18S +
  plot_layout(design = lay_18S) +
  plot_annotation(title = "Temporal Variation in Diversity and Composition in 18S rRNA",
                  theme = theme(plot.title = element_text(face = "bold",
                                                          hjust = 0.5,
                                                          size = 14)))

ggsave("18S_SINKS/figures/area/plot_18S.svg", plot = p_final_area_18S,
       width = 10,
       height = 10,
       units = "in",
       dpi = 100)

ggsave("18S_SINKS/figures/period/plot_18S.svg", plot = p_final_unit_18S,
       width = 10,
       height = 10,
       units = "in",
       dpi = 100)

#ITS
p_alpha_area_ITS <- (p_shannon_area_ITS + p_chao1_area_ITS + p_bdp_area_ITS +
                       plot_layout(guides = "collect", axis_titles = "collect")) +
  plot_annotation(
    title = "α-Diversity Indexes",
    theme = theme(plot.title = element_text(face = "bold",
                                            hjust = 0.5))) &
  theme(legend.position = "right",
        legend.title.position = "top",
        axis.text.x = element_blank(),
        legend.key.spacing.y = unit(5, "pt"),
        legend.key.size = unit(0.25, "cm"),
        legend.box.spacing = unit(0, "pt")
  )

p_beta_area_ITS <- (p_bray_area_ITS + p_jacc_area_ITS +
                      plot_layout(guides = "collect", axis_titles = "collect")) +
  plot_annotation(
    title = "β-Diversity",
    theme = theme(plot.title = element_text(face = "bold",
                                            hjust = 0.5))) &
  theme(legend.position = "right",
        legend.title.position = "top",
        axis.text.x = element_blank(),
        legend.key.spacing.y = unit(5, "pt"),
        legend.key.size = unit(0.25, "cm"),
        legend.box.spacing = unit(0, "pt")
  )

p_final_area_ITS <- p_alpha_area_ITS + p_bray_area_ITS + p_jacc_area_ITS + 
  comp_area_ITS + venn_area_ITS +
  plot_layout(design = lay) +
  plot_annotation(title = "Spatial Variation in Diversity and Composition in ITS",
                  theme = theme(plot.title = element_text(face = "bold",
                                                          hjust = 0.5,
                                                          size = 14)))


ggsave("ITS_SINKS/figures/area/plot_ITS.svg", plot = p_final_area_ITS,
       width = 10,
       height = 10,
       units = "in",
       dpi = 100)

