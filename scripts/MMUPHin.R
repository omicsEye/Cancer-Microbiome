library(MMUPHin)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(omicsArt)
library(dplyr)
library(stringr)
library(omicsArt)
library(ComplexHeatmap)
library(colorRamp2)
library(readxl)

setwd("/Users/xinyang/Library/CloudStorage/Box-Box/Cancer_Mircrobiome/MetaAnalysis")

Daver <- read.delim(
  "data/Davar_merged_table.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
species_Daver <- Daver
species_Daver <- species_Daver[grepl("\\|s__[^\\|]+$", rownames(species_Daver)),]
rownames(species_Daver) <- gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "s__", rownames(species_Daver))

cols_temp_Daver <- str_split(colnames(species_Daver), "_")
colnames(species_Daver) <- sapply(cols_temp_Daver, "[[", 2)
write.table(species_Daver, file = "Metaphlan4_Daver2021.tsv", col.names = NA, sep="\t")

Frankel <- read.delim(
   #"data/Frankel2_merged_abundance_table.txt",
  "data/Frankel_merged_table.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
species_Frankel <- Frankel
species_Frankel <- species_Frankel[grepl("\\|s__[^\\|]+$", rownames(species_Frankel)),]
rownames(species_Frankel) <- gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "s__", rownames(species_Frankel))
cols_temp_Frankel <- str_split(colnames(species_Frankel), "_")
colnames(species_Frankel) <- sapply(cols_temp_Frankel, "[[", 2)
# species_Frankel
write.table(species_Frankel, file = "Metaphlan4_Frankel2017.tsv", col.names = NA, sep="\t")

Baruch <- read.delim(
  "data/Baruch_merged_table.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
species_Baruch <- Baruch
species_Baruch <- species_Baruch[grepl("\\|s__[^\\|]+$", rownames(species_Baruch)),]
rownames(species_Baruch) <- gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "s__", rownames(species_Baruch))
cols_temp_Baruch <- str_split(colnames(species_Baruch), "_")
colnames(species_Baruch) <- sapply(cols_temp_Baruch, "[[", 2)
write.table(species_Baruch, file = "Metaphlan4_Baruch2021.tsv", col.names = NA, sep="\t")

Lee <- read.delim(
  "data/Lee_merged_table.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
species_Lee <- Lee
species_Lee <- species_Lee[grepl("\\|s__[^\\|]+$", rownames(species_Lee)),]
rownames(species_Lee) <- gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "s__", rownames(species_Lee))
# cols_temp_Lee <- str_split(colnames(species_Lee), "_")
# colnames(species_Lee) <- sapply(cols_temp_Lee, "[[", 2)
write.table(species_Lee, file = "Metaphlan4_Lee2022.tsv", col.names = NA, sep="\t")

Mat <- read.delim(
  "data/Mat_merged_table.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
species_Mat <- Mat
species_Mat <- species_Mat[grepl("\\|s__[^\\|]+$", rownames(species_Mat)),]
rownames(species_Mat) <- gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "s__", rownames(species_Mat))
cols_temp_Mat <- str_split(colnames(species_Mat), "_")
colnames(species_Mat) <- sapply(cols_temp_Mat, "[[", 2)
write.table(species_Mat, file = "Metaphlan4_Mat2018.tsv", col.names = NA, sep="\t")

Sep <- read.delim(
  "data/Sep_merged_table.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
species_Sep <- Sep
species_Sep <- species_Sep[grepl("\\|s__[^\\|]+$", rownames(species_Sep)),]
rownames(species_Sep) <- gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "s__", rownames(species_Sep))
cols_temp_Sep <- str_split(colnames(species_Sep), "_")
colnames(species_Sep) <- sapply(cols_temp_Sep, "[[", 2)
write.table(species_Sep, file = "Metaphlan4_Sep2021.tsv", col.names = NA, sep="\t")

species_Gop <- read.delim(
  "data/Gop_merged_abundance_table_species.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
cols_temp_Gop <- str_split(colnames(species_Gop), "_")
colnames(species_Gop) <- sapply(cols_temp_Gop, "[[", 2)
write.table(species_Gop, file = "Metaphlan4_Gop2018.tsv", col.names = NA, sep="\t")

Peter <- read.delim(
  "data/Peter_merged_table.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

species_Peter <- Peter
species_Peter <- species_Peter[grepl("\\|s__[^\\|]+$", rownames(species_Peter)),]
rownames(species_Peter) <- gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "s__", rownames(species_Peter))
cols_temp_Peter <- str_split(colnames(species_Peter), "_")
colnames(species_Peter) <- sapply(cols_temp_Peter, "[[", 2)
write.table(species_Peter, file = "Metaphlan4_Peter2019.tsv", col.names = NA, sep="\t")


library(purrr)
library(tidyverse)
df_species_Lee <- species_Lee 
species <- rownames(df_species_Lee)
rownames(df_species_Lee) <- NULL
df_species_Lee <- cbind(species,df_species_Lee)

df_species_Frankel <- species_Frankel 
species <- rownames(df_species_Frankel)
rownames(df_species_Frankel) <- NULL
df_species_Frankel <- cbind(species,df_species_Frankel)

df_species_Daver <- species_Daver 
species <- rownames(df_species_Daver)
rownames(df_species_Daver) <- NULL
df_species_Daver <- cbind(species,df_species_Daver)

df_species_Sep <- species_Sep 
species <- rownames(df_species_Sep)
rownames(df_species_Sep) <- NULL
df_species_Sep <- cbind(species,df_species_Sep)

df_species_Gop <- species_Gop 
species <- rownames(df_species_Gop)
rownames(df_species_Gop) <- NULL
df_species_Gop <- cbind(species,df_species_Gop)

df_species_Mat <- species_Mat 
species <- rownames(df_species_Mat)
rownames(df_species_Mat) <- NULL
df_species_Mat <- cbind(species,df_species_Mat)

df_species_Baruch <- species_Baruch 
species <- rownames(df_species_Baruch)
rownames(df_species_Baruch) <- NULL
df_species_Baruch <- cbind(species,df_species_Baruch)

####_________________________HOLD______________________
# df_species_Peter <- species_Peter 
# species <- rownames(df_species_Peter)
# rownames(df_species_Peter) <- NULL
# df_species_Peter <- cbind(species,df_species_Peter)
####_________________________HOLD______________________


list_of_dfs <- list(df_species_Lee, df_species_Frankel, df_species_Daver, 
                    df_species_Sep, df_species_Gop, df_species_Mat, 
                    df_species_Baruch)

result <- purrr::reduce(list_of_dfs, ~full_join(.x, .y, by="species"))
result2 <- result[,-1]
rownames(result2) <- result[,1]


# data("CRC_abd", "CRC_meta")
# CRC_abd[1:5, 1, drop = FALSE]
# CRC_meta[1, 1:5]
# table(CRC_meta$studyID)

meta <- read.delim(
  "data/meta.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
)

table(meta$dataset)

meta2 <- meta[,-2]
rownames(meta2) <- meta[,2]
result2[is.na(result2)] <- 0
result2 <- result2[,1:ncol(result2)]/100

result2 <- result2[,order(colnames(result2))]
meta2 <- meta2[order(row.names(meta2)),]

fit_adjust_batch2 <- adjust_batch(feature_abd = result2,
                                 batch = "dataset",
                                 covariates = "response",
                                 data = meta2,
                                 control = list(verbose = FALSE))
###After batch effect adjust matrix
fit_adjust_batch2_adj <- fit_adjust_batch2$feature_abd_adj
write.table(fit_adjust_batch2_adj, file = "After_batch_effect_adjust_mergedTable.tsv", col.names = NA, sep="\t")

library(vegan, quietly = TRUE)
D_before <- vegdist(t(result2))
D_after <- vegdist(t(fit_adjust_batch2_adj))

set.seed(1)
fit_adonis_before <- adonis(D_before ~ dataset, data = meta2)
fit_adonis_after <- adonis(D_after ~ dataset, data = meta2)

print(fit_adonis_before$aov.tab)
print(fit_adonis_after$aov.tab)


fit_lm_meta <- lm_meta(feature_abd = fit_adjust_batch2_adj,
                       batch = "dataset",
                       exposure = "response",
                       covariates = NULL,
                       data = meta2,
                       control = list(verbose = FALSE,
                                      analysis_method="CPLM"))
fit_lm_meta[["control"]][["forest_plot"]]
meta_fits <- fit_lm_meta$meta_fits
meta_fits %>% 
  filter(qval.fdr < 0.25) %>% 
  arrange(coef) %>% 
  mutate(feature = factor(feature, levels = feature)) %>% 
  ggplot(aes(y = coef, x = feature)) +
  geom_bar(stat = "identity") +
  coord_flip()


# control_meta <- subset(meta2, response == "R")
# control_abd_adj <- fit_adjust_batch2_adj[, rownames(control_meta)]
# 
# D_control <- vegdist(t(control_abd_adj))
# fit_discrete <- discrete_discover(D = D_control,
#                                   batch = "dataset",
#                                   data = control_meta,
#                                   control = list(k_max = 8,
#                                                  verbose = FALSE))
# 
# study_id = "Lee_2022"
# 
# internal <- data.frame(K = 2:7,
#                        statistic = fit_discrete$internal_mean[, study_id],
#                        se = fit_discrete$internal_se[, study_id],
#                        type = "internal")
# external <- data.frame(K = 2:7,
#                        statistic = fit_discrete$external_mean[, study_id],
#                        se = fit_discrete$external_se[, study_id],
#                        type = "external")
# rbind(internal, external) %>% 
#   ggplot(aes(x = K, y = statistic, color = type)) +
#   geom_point(position = position_dodge(width = 0.5)) + 
#   geom_line(position = position_dodge(width = 0.5)) +
#   geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
#                 position = position_dodge(width = 0.5), width = 0.5) +
#   ggtitle("Evaluation of discrete structure in control stool microbiome (Lee_2022)")
# 
# fit_continuous <- continuous_discover(feature_abd = control_abd_adj,
#                                       batch = "dataset",
#                                       data = control_meta,
#                                       control = list(var_perc_cutoff = 0.5,
#                                                      verbose = FALSE))
