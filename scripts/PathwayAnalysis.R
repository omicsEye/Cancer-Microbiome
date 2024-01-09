library(iGSEA)
library(Maaslin2)

setwd("/Users/xinyang/Library/CloudStorage/Box-Box/Cancer_Mircrobiome/MetaAnalysis")
input_data <- read.delim(
  "PathwayAbundanceData/Spencer2021_PathcoverageMergedTable.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
)
pathway_Spencer <- input_data
cols_temp_Spencer <- str_split(colnames(pathway_Spencer), "_")
colnames(pathway_Spencer) <- sapply(cols_temp_Spencer, "[[", 1)

input_metadata <- read.delim(
  "PathwayAbundanceData/Spencer2021_metadata.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
)
fit_data <- Maaslin2(
  pathway_Spencer, input_metadata,'Spencer_output', transform = "AST",
  normalization = 'NONE',
  reference = 'response,R',
  standardize = FALSE)


