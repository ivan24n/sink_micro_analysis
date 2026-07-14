##############################################################################
# 03_rarefaction.R
#
# Rarefaction analysis for each marker (16S, 18S, ITS): sequencing depth and
# richness summaries per sample, rarefaction curves, and generation of a
# rarefied OTU/ASV table at a fixed subsampling depth per marker.
#
# Part of the analysis pipeline for:
#   "Spatiotemporal dynamics of multi-kingdom microbial communities in
#    hospital sinks"
#
# NOTE: the 16S, 18S, and ITS blocks below are structurally very similar
# but are kept as separate blocks (as in the original script) rather than
# merged into a single function, since some parameters differ by marker
# (rarefaction depth, sample-ID naming pattern). The estimated-richness bar
# plot is now built the same (safer) way for all three markers: values are
# joined explicitly by Sample rather than assumed to be in matching row
# order.
##############################################################################

library(tidyverse)
library(vegan)
library(patchwork)
set.seed(12)

# Set this to the root of the project directory, i.e. the folder that
# contains the 16S_SINKS/, 18S_SINKS/, and ITS_SINKS/ subfolders.
work_dir <- "path/to/project"
setwd(work_dir)

                      ######## RNA 16S #######

#### Load OTU table: 16S ####
otu_table_16S <- read.table("16S_SINKS/tables/tax_tab_16S.tsv",
                            header = TRUE,
                            row.names = 1,
                            check.names = FALSE,
                            sep = "\t")


#### Summarize sequencing depth and richness: 16S ####
total_sum_16S <- data.frame(Sample = names(otu_table_16S),
                            Total = colSums(otu_table_16S),
                            Richness = colSums(otu_table_16S != 0),
                            row.names = NULL)

top_seqs_16S <- total_sum_16S %>%
  arrange(desc(Total)) %>%
  rename(n_seqs = Total)

ggplot(total_sum_16S, aes(x = Total)) +
  geom_histogram()

seqs_samp_16S <- ggplot(total_sum_16S, aes(x = Sample, y = Total)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title =
  "Comparison total sequences, diversity and estimated diversity after rarefy",
    x = NULL,
    y = "n_seqs") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, max(total_sum_16S$Total), by = 50000))

rich_samp_16S <- ggplot(total_sum_16S, aes(x = Sample, y = Richness)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "Richness") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, max(total_sum_16S$Richness), by = 100))

#### Rarefy: 16S ####
otu_wider_16S <- otu_table_16S %>%
  t() %>%
  as.data.frame()

min_seqs_16S <- min(total_sum_16S$Total)

otu_wider_16S <- otu_wider_16S %>%
  rownames_to_column(var = "Sample") %>%
  rowwise() %>%
  filter(sum(c_across(-Sample)) >= min_seqs_16S) %>%
  ungroup() %>%
  column_to_rownames(var = "Sample")

total_sum_16S <- total_sum_16S %>%
  filter(Total >= min_seqs_16S)


rich_rare_16S <- rarefy(otu_wider_16S, min_seqs_16S) %>%
  as_tibble(rownames = "Sample") %>%
  select(Sample, Vegan_est = value)

# Join explicitly by Sample rather than relying on row order, to avoid any
# risk of misalignment between total_sum_16S and rich_rare_16S.
est_rich_data_16S <- total_sum_16S %>%
  inner_join(rich_rare_16S, by = "Sample")

est_rich_samp_16S <- ggplot(est_rich_data_16S, aes(x = Sample,
                                               y = Vegan_est)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "Estimated Richness") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(breaks = seq(0, max(total_sum_16S$Richness), by = 100))

combined_plot_16S <- seqs_samp_16S / rich_samp_16S / est_rich_samp_16S

rare_curves_16S <- rarecurve(otu_wider_16S, step = 100)

df_rare_curves_16S <- map_dfr(rare_curves_16S, bind_rows) %>%
  bind_cols(Sample = rownames(otu_wider_16S), .) %>%
  pivot_longer(-Sample) %>%
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_remove(name, "N"))) %>%
  rename(estimated_spe = value) %>%
  extract(
    col = Sample,
    into = c("period", "box", "area", "fase"),
    regex = "SMPL16S_([^_]+)_([^_]+)_(.)_(.)",
    remove = FALSE) %>%
  select(-name)

df_rare_curves_16S$box <- factor(as.numeric(df_rare_curves_16S$box))

plot_rc_16S <- ggplot(df_rare_curves_16S,
                     aes(x = n_seqs, y = estimated_spe, group = Sample)) +
  geom_line(linewidth = 0.4, alpha = 0.7, colour = "steelblue") +
  geom_vline(xintercept = min_seqs_16S, color = "darkorange",
             linetype = "dashed", linewidth = 0.8) +
  coord_cartesian(xlim = c(0, 100000)) +
  labs(title = "Rarefaction curves of RNA 16S sequences",
       x = "Number of Sequences",
       y = "Estimated Richness (Species)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.spacing = unit(1, "lines")) +
  scale_x_continuous(
    breaks = seq(0, 100000, by = 25000),
    labels = scales::comma) +
  facet_wrap(~ box)

rare_tab16S <- rrarefy(otu_wider_16S, sample = min_seqs_16S)
rare_tab16S <- t(rare_tab16S)
rare_tab16S <- as.data.frame(rare_tab16S)
rare_tab16S <- rare_tab16S %>%
  rownames_to_column(var = "Taxa") %>%
  mutate(seqs = rowSums(across(-1))) %>%
  filter(seqs > 1) %>%
  arrange(desc(seqs)) %>%
  select(-seqs)

#### Save figures and tables: 16S ####

write_tsv(total_sum_16S, "16S_SINKS/tables/total_seqs_16S.tsv")
write.table(rare_tab16S, "16S_SINKS/tables/rarefacted_table_16S.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE)

ggsave("16S_SINKS/tables/rare_curves_16S.png", plot_rc_16S)
ggsave("16S_SINKS/tables/reads_rich_16S.png", combined_plot_16S)

                      ######## RNA 18S #######

#### Load OTU table: 18S ####
otu_table_18S <- read.table("18S_SINKS/tables/tax_tab_18S.tsv",
                            header = TRUE,
                            row.names = 1,
                            check.names = FALSE,
                            sep = "\t")


#### Summarize sequencing depth and richness: 18S ####
total_sum_18S <- data.frame(Sample = names(otu_table_18S),
                            Total = colSums(otu_table_18S),
                            Richness = colSums(otu_table_18S != 0),
                            row.names = NULL)

top_seqs_18S <- total_sum_18S %>%
  arrange(desc(Total)) %>%
  rename(n_seqs = Total)

ggplot(total_sum_18S, aes(x = Total)) +
  geom_histogram(fill = "steelblue", binwidth = 8915)


seqs_samp_18S <- ggplot(total_sum_18S, aes(x = Sample, y = Total)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title =
  "Comparison total sequences, diversity and estimated diversity after rarefy",
       x = NULL,
       y = "n_seqs") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, max(total_sum_18S$Total), by = 100000))

rich_samp_18S <- ggplot(total_sum_18S, aes(x = Sample, y = Richness)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "Richness") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, max(total_sum_18S$Richness), by = 100))

#### Rarefy: 18S ####
otu_wider_18S <- otu_table_18S %>%
  t() %>%
  as.data.frame()

# Fixed rarefaction depth chosen for the 18S dataset (see manuscript
# methods for the rationale behind this value).
min_seqs_18S <- 10760

otu_wider_18S <- otu_wider_18S %>%
  rownames_to_column(var = "Sample") %>%
  rowwise() %>%
  filter(sum(c_across(-Sample)) >= min_seqs_18S) %>%
  ungroup() %>%
  column_to_rownames(var = "Sample")

total_sum_18S <- total_sum_18S %>%
  filter(Total >= min_seqs_18S)

rich_rare_18S <- rarefy(otu_wider_18S, min_seqs_18S) %>%
  as_tibble(rownames = "Sample") %>%
  select(Sample, Vegan_est = value)

est_rich_samp_18S <- ggplot(rich_rare_18S, aes(x = Sample,
                                               y = Vegan_est)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "Estimated Richness") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

combined_plot_18S <- seqs_samp_18S / rich_samp_18S / est_rich_samp_18S

rare_curves_18S <- rarecurve(otu_wider_18S, step = 100)

df_rare_curves_18S <- map_dfr(rare_curves_18S, bind_rows) %>%
  bind_cols(Sample = rownames(otu_wider_18S), .) %>%
  pivot_longer(-Sample) %>%
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_remove(name, "N"))) %>%
  rename(estimated_spe = value) %>%
  extract(
    col = Sample,
    into = c("period", "box", "area", "fase"),
    regex = "SMPL18S_([^_]+)_([^_]+)_(.)_(.)",
    remove = FALSE) %>%
  select(-name)

df_rare_curves_18S$box <- factor(as.numeric(df_rare_curves_18S$box))

plot_rc_18S <- ggplot(df_rare_curves_18S,
                      aes(x = n_seqs, y = estimated_spe, group = Sample)) +
  geom_line(linewidth = 0.4, alpha = 0.7, colour = "steelblue") +
  geom_vline(xintercept = min_seqs_18S, color = "darkorange",
             linetype = "dashed", linewidth = 0.8) +
  coord_cartesian(xlim = c(0, 100000)) +
  labs(title = "Rarefaction curves of RNA 18S sequences",
       x = "Number of Sequences",
       y = "Estimated Richness (Species)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.spacing = unit(1, "lines")) +
  scale_x_continuous(
    breaks = seq(0, 100000, by = 25000),
    labels = scales::comma
  ) +
  facet_wrap(~ box)

rare_tab18S <- rrarefy(otu_wider_18S, sample = min_seqs_18S)
rare_tab18S <- as.data.frame(t(rare_tab18S))
rare_tab18S <- as.data.frame(rare_tab18S)
rare_tab18S <- rare_tab18S %>%
  rownames_to_column(var = "Taxa") %>%
  mutate(seqs = rowSums(across(-1))) %>%
  arrange(desc(seqs)) %>%
  select(-seqs)

#### Save figures and tables: 18S ####

write_tsv(total_sum_18S, "18S_SINKS/tables/total_seqs_18S.tsv")
write.table(rare_tab18S, "18S_SINKS/tables/rarefacted_table_18S.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE)

ggsave("18S_SINKS/figures/rare_curves_18S.png", plot_rc_18S)
ggsave("18S_SINKS/figures/reads_rich_18S.png", combined_plot_18S)

                      ######## RNA ITS #######

#### Load OTU table: ITS ####

otu_table_ITS <- read.table("ITS_SINKS/tables/tax_tab_ITS.tsv",
                            header = TRUE,
                            row.names = 1,
                            check.names = FALSE,
                            sep = "\t")


#### Summarize sequencing depth and richness: ITS ####
total_sum_ITS <- data.frame(Sample = names(otu_table_ITS),
                            Total = colSums(otu_table_ITS),
                            Richness = colSums(otu_table_ITS != 0),
                            row.names = NULL)

top_seqs_ITS <- total_sum_ITS %>%
  arrange(desc(Total)) %>%
  rename(n_seqs = Total)

ggplot(total_sum_ITS, aes(x = Total)) +
  geom_histogram()

seqs_samp_ITS <- ggplot(total_sum_ITS, aes(x = Sample, y = Total)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title =
         "Comparison total sequences, diversity and estimated diversity after rarefy",
       x = NULL,
       y = "n_seqs") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, max(total_sum_ITS$Total), by = 50000))

rich_samp_ITS <- ggplot(total_sum_ITS, aes(x = Sample, y = Richness)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "Richness") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, max(total_sum_ITS$Richness), by = 100))

#### Rarefy: ITS ####
otu_wider_ITS <- otu_table_ITS %>%
  t() %>%
  as.data.frame()

# Fixed rarefaction depth chosen for the ITS dataset (see manuscript
# methods for the rationale behind this value).
min_seqs_ITS <- 9027

otu_wider_ITS <- otu_wider_ITS %>%
  rownames_to_column(var = "Sample") %>%
  rowwise() %>%
  filter(sum(c_across(-Sample)) >= min_seqs_ITS) %>%
  ungroup() %>%
  column_to_rownames(var = "Sample")

total_sum_ITS <- total_sum_ITS %>%
  filter(Total >= min_seqs_ITS)

rich_rare_ITS <- rarefy(otu_wider_ITS, min_seqs_ITS) %>%
  as_tibble(rownames = "Sample") %>%
  select(Sample, Vegan_est = value)

# Join explicitly by Sample rather than relying on row order, to avoid any
# risk of misalignment between total_sum_ITS and rich_rare_ITS.
est_rich_data_ITS <- total_sum_ITS %>%
  inner_join(rich_rare_ITS, by = "Sample")

est_rich_samp_ITS <- ggplot(est_rich_data_ITS, aes(x = Sample,
                                               y = Vegan_est)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "Estimated Richness") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(breaks = seq(0, max(total_sum_ITS$Richness), by = 100))

combined_plot_ITS <- seqs_samp_ITS / rich_samp_ITS / est_rich_samp_ITS

rare_curves_ITS <- rarecurve(otu_wider_ITS, step = 100)

df_rare_curves_ITS <- map_dfr(rare_curves_ITS, bind_rows) %>%
  bind_cols(Sample = rownames(otu_wider_ITS), .) %>%
  pivot_longer(-Sample) %>%
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_remove(name, "N"))) %>%
  rename(estimated_spe = value) %>%
  extract(
    col = Sample,
    into = c("period", "box", "area"),
    regex = "SMPL-ITS_([^_]+)_([^_]+)_(.)",
    remove = FALSE) %>%
  select(-name)

df_rare_curves_ITS$box <- factor(as.numeric(df_rare_curves_ITS$box))

plot_rc_ITS <- ggplot(df_rare_curves_ITS,
                      aes(x = n_seqs, y = estimated_spe, group = Sample)) +
  geom_line(linewidth = 0.4, alpha = 0.7, colour = "steelblue") +
    geom_vline(xintercept = min_seqs_ITS, color = "darkorange",
             linetype = "dashed", linewidth = 0.8) +
  coord_cartesian(xlim = c(0, 100000)) +
  labs(title = "Rarefaction curves of RNA ITS sequences",
       x = "Number of Sequences",
       y = "Estimated Richness (Species)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.spacing = unit(1, "lines")) +
  scale_x_continuous(
    breaks = seq(0, 100000, by = 25000),
    labels = scales::comma
  ) +
  facet_wrap(~ box)

rare_tabITS <- rrarefy(otu_wider_ITS, sample = min_seqs_ITS)
rare_tabITS <- as.data.frame(t(rare_tabITS))
rare_tabITS <- as.data.frame(rare_tabITS)
rare_tabITS <- rare_tabITS %>%
  rownames_to_column(var = "Taxa") %>%
  mutate(seqs = rowSums(across(-1))) %>%
  arrange(desc(seqs)) %>%
  select(-seqs)

#### Save figures and tables: ITS ####

write_tsv(total_sum_ITS, "ITS_SINKS/tables/total_seqs_ITS.tsv")
write.table(rare_tabITS, "ITS_SINKS/tables/rarefacted_table_ITS.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE)

ggsave("ITS_SINKS/figures/rare_curves_ITS.png", plot_rc_ITS)
ggsave("ITS_SINKS/figures/reads_rich_ITS.png", combined_plot_ITS)
