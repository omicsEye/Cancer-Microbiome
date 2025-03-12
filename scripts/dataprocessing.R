sp <- read.delim(
  "data/maaslin2_output/significant_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
sp_abd <- read.delim(
  "N_NR_sorted_mergedTable.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

significant_sp <- as.matrix(rownames(sp))
sp_abd_significant <- subset(sp_abd,rownames(sp_abd) %in% significant_sp)
write.table(sp_abd_significant, file = "N_NR_sp_mergedTable.tsv", col.names = NA, sep="\t")

bgc <- read.delim(
  "data/maaslin2_bgc_output/significant_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
bgc_abd <- read.delim(
  "N_NR_sorted_bgc_mergedTable.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)


significant_bgc <- as.matrix(rownames(bgc))
bgc_abd_significant <- subset(bgc_abd,rownames(bgc_abd) %in% significant_bgc)
write.table(bgc_abd_significant, file = "N_NR_bgc_mergedTable.tsv", col.names = NA, sep="\t")
