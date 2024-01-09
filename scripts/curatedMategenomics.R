library(curatedMetagenomicData)
library(curatedMetagenomicAnalyses)
library(dplyr)
library(tidyverse)
library(DT)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)

# BiocManager::install("waldronlab/curatedMetagenomicDataAnalyses", dependencies = TRUE)
tse_LeeKA_2022 <- curatedMetagenomicData("LeeKA_2022.relative_abundance", dryrun = FALSE)
LeeKA_2022_relativeAbundance <- assays(tse_LeeKA_2022[["2022-04-13.LeeKA_2022.relative_abundance"]])[[1]]
LeeKA_2022_relativeAbundance
species_Lee <- LeeKA_2022_relativeAbundance
species_Lee <- species_Lee[grepl("\\|s__[^\\|]+$", rownames(species_Lee)),]
rownames(species_Lee) <- gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "s__", rownames(species_Lee))
write.table(species_Lee, file = "CMD_Lee2022.tsv", col.names = NA, sep="\t")


tse_GopalakrishnanV_2018 <- curatedMetagenomicData("GopalakrishnanV_2018.relative_abundance", dryrun = FALSE)
GopalakrishnanV_2018_relativeAbundance <- assays(tse_GopalakrishnanV_2018[["2021-10-14.GopalakrishnanV_2018.relative_abundance"]])[[1]]
species_Gop <- GopalakrishnanV_2018_relativeAbundance
species_Gop <- species_Gop[grepl("\\|s__[^\\|]+$", rownames(species_Gop)),]
rownames(species_Gop) <- gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "s__", rownames(species_Gop))
write.table(species_Gop, file = "CMD_Gop2018.tsv", col.names = NA, sep="\t")

tse_MatsonV_2018 <- curatedMetagenomicData("MatsonV_2018.relative_abundance", dryrun = FALSE)
MatsonV_2018_relativeAbundance <- assays(tse_MatsonV_2018[["2021-10-14.MatsonV_2018.relative_abundance"]])[[1]]
species_Mat <- MatsonV_2018_relativeAbundance
species_Mat <- species_Mat[grepl("\\|s__[^\\|]+$", rownames(species_Mat)),]
rownames(species_Mat) <- gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "s__", rownames(species_Mat))
write.table(species_Mat, file = "CMD_Mat2018.tsv", col.names = NA, sep="\t")


tse_FrankelAE_2017 <- curatedMetagenomicData("FrankelAE_2017.relative_abundance", dryrun = FALSE)
FrankelAE_2017_relativeAbundance <- assays(tse_FrankelAE_2017[["2022-04-13.FrankelAE_2017.relative_abundance"]])[[1]]
species_Frankel <- FrankelAE_2017_relativeAbundance
species_Frankel <- species_Frankel[grepl("\\|s__[^\\|]+$", rownames(species_Frankel)),]
rownames(species_Frankel) <- gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "s__", rownames(species_Frankel))
write.table(species_Frankel, file = "CMD_Frankel2017.tsv", col.names = NA, sep="\t")

metadata_Frankel2017 <- sampleMetadata |>
  filter(study_name == "FrankelAE_2017")

tse_PetersBA_2019 <- curatedMetagenomicData("PetersBA_2019.relative_abundance", dryrun = FALSE)
metadata_PetersBA_2019 <- sampleMetadata |>
  filter(study_name == "PetersBA_2019")


test <- sampleMetadata |>
  filter(study_name == "PetersBA_2019") |>
  select(where(~ !any(is.na(.x)))) |>
  datatable(options = list(dom = "t"), extensions = "Responsive")
