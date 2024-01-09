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
  "data/BGC_Davar_Merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
bgc_Daver <- Daver

Frankel <- read.delim(
  "data/BGC_Frankel_Merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
bgc_Frankel <- Frankel

Baruch <- read.delim(
  "data/BGC_Baruch_merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
bgc_Baruch <- Baruch

Lee <- read.delim(
  "data/BGC_Lee_merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
bgc_Lee <- Lee

Mat <- read.delim(
  "data/BGC_Matson_Merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
bgc_Mat <- Mat

Sep <- read.delim(
  "data/BGC_Spencer_Merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
bgc_Sep <- Sep

bgc_Gop <- read.delim(
  "data/BGC_Gop_Merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)


library(purrr)
library(tidyverse)
df_bgc_Lee <- bgc_Lee 
BGC <- rownames(df_bgc_Lee)
rownames(df_bgc_Lee) <- NULL
df_bgc_Lee <- cbind(BGC,df_bgc_Lee)

df_bgc_Frankel <- bgc_Frankel 
BGC <- rownames(df_bgc_Frankel)
rownames(df_bgc_Frankel) <- NULL
df_bgc_Frankel <- cbind(BGC,df_bgc_Frankel)

df_bgc_Daver <- bgc_Daver 
BGC <- rownames(df_bgc_Daver)
rownames(df_bgc_Daver) <- NULL
df_bgc_Daver <- cbind(BGC,df_bgc_Daver)

df_bgc_Sep <- bgc_Sep 
BGC <- rownames(df_bgc_Sep)
rownames(df_bgc_Sep) <- NULL
df_bgc_Sep <- cbind(BGC,df_bgc_Sep)

df_bgc_Gop <- bgc_Gop 
BGC <- rownames(df_bgc_Gop)
rownames(df_bgc_Gop) <- NULL
df_bgc_Gop <- cbind(BGC,df_bgc_Gop)

df_bgc_Mat <- bgc_Mat 
BGC <- rownames(df_bgc_Mat)
rownames(df_bgc_Mat) <- NULL
df_bgc_Mat <- cbind(BGC,df_bgc_Mat)

df_bgc_Baruch <- bgc_Baruch 
BGC <- rownames(df_bgc_Baruch)
rownames(df_bgc_Baruch) <- NULL
df_bgc_Baruch <- cbind(BGC,df_bgc_Baruch)


list_of_dfs <- list(df_bgc_Lee, df_bgc_Frankel, df_bgc_Daver, 
                    df_bgc_Sep, df_bgc_Gop, df_bgc_Mat, df_bgc_Baruch)

result <- purrr::reduce(list_of_dfs, ~full_join(.x, .y, by="BGC"))
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

fit_adjust_batch_bgc <- adjust_batch(feature_abd = result2,
                                  batch = "dataset",
                                  covariates = "response",
                                  data = meta2,
                                  control = list(verbose = FALSE))
fit_adjust_batch_bgc_adj <- fit_adjust_batch_bgc$feature_abd_adj
write.table(fit_adjust_batch_bgc_adj, file = "After_batch_effect_adjust_BGCmergedTable.tsv", col.names = NA, sep="\t")

library(vegan, quietly = TRUE)
D_before <- vegdist(t(result2))
D_after <- vegdist(t(fit_adjust_batch_bgc_adj))

set.seed(1)
fit_adonis_before <- adonis(D_before ~ dataset, data = meta2)
fit_adonis_after <- adonis(D_after ~ dataset, data = meta2)

print(fit_adonis_before$aov.tab)
print(fit_adonis_after$aov.tab)

fit_adjust_batch_bgc_adj_log <- fit_adjust_batch_bgc_adj ** (1/2)


fit_lm_meta_bgc <- lm_meta(feature_abd = fit_adjust_batch_bgc_adj_log,
                       batch = "dataset",
                       exposure = "response",
                       covariates = NULL,
                       data = meta2,
                       control = list(verbose = FALSE,
                                      analysis_method="CPLM"))

fit_lm_meta_bgc[["control"]][["forest_plot"]]
meta_fits_bgc <- fit_lm_meta_bgc$meta_fits
meta_fits_bgc %>% 
  filter(qval.fdr < 0.05) %>% 
  arrange(coef) %>% 
  mutate(feature = factor(feature, levels = feature)) %>% 
  ggplot(aes(y = coef, x = feature)) +
  geom_bar(stat = "identity") +
  coord_flip()

meta_fits_bgc_selected = meta_fits_bgc %>% 
  filter(qval.fdr < 0.05) 
meta_fits_bgc_selected$feature

