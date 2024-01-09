library(Maaslin2)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
# load libraries
library(ggplot2)
library(ggplotify)
library(pheatmap)
library(patchwork)
setwd("/Users/xinyang/Library/CloudStorage/Box-Box/Cancer_Mircrobiome/MetaAnalysis")


### METAPHLAN - species ###
All <- read.delim(
  "After_batch_effect_adjust_mergedTable.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

All_bgc <- read.delim(
  "After_batch_effect_adjust_BGCmergedTable.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

All_meta <- read.delim(
  "data/meta.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
)

All_meta <- subset(All_meta,select = c("sampleid","response"))
All_meta <- All_meta[order(All_meta$response,decreasing = TRUE),]
sort_by_response_group <- c(All_meta[,1])
All <- All[sort_by_response_group]
All_bgc <- All_bgc[sort_by_response_group]
rownames(All_meta) = All_meta[,1]
All_meta <- subset(All_meta,select = c("response"))

write.table(All, file = "N_NR_sorted_mergedTable.tsv", col.names = NA, sep="\t")
write.table(All_bgc, file = "N_NR_sorted_bgc_mergedTable.tsv", col.names = NA, sep="\t")

Maaslin2::Maaslin2(
  input_data  = t(All),
  input_metadata =  All_meta,
  output =  'data/maaslin2_output',
  max_significance= 0.1,
  analysis_method = 'LM',
  standardize = FALSE,
  transform = 'LOG',
  normalization = 'NONE',
  heatmap_first_n = 20,
  reference = c("response,R"))

### FDR concept
masslin2_stats  <- read.delim(
  "data/maaslin2_output/significant_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = 1
)


top <- 30
sp_top_lst <- as.matrix(sort(rowSums(All),decreasing = TRUE)[1:top])
All_sp_top <- subset(All,rownames(All) %in% rownames(sp_top_lst))

All_sp_top <- All_sp_top[order(rowSums(All_sp_top),decreasing = TRUE),]
write.table(All_sp_top, file = "Top30_N_NR_sorted_mergedTable.tsv", col.names = NA, sep="\t")
# sort_by_response_group <- c(All_meta[,1])
# sort_by_response_group <- rowSums(All_sp_top100)
fakepcl <- list(meta=All_meta, x=as.matrix(t(All_sp_top)),
                    ns=dim(All_sp_top)[2], nf=dim(All_sp_top)[1])
heat_plot <- omicsArt:::pcl.heatmap(fakepcl, sqrtspace = T, gamma = 2, meta= T, 
                                            show_colnames = F, show_rownames = T,
                                            treeheight_row = 0, treeheight_col= 5)
ggsave(filename='data/All_heat_plot_top30.png', plot=heat_plot, width = 11, height = 18, units = "in", dpi = 300)
library("viridis") 
new_log <- All_sp_top ** (1/2)
pic <- pheatmap(as.matrix(All_sp_top),
                #annotation_col = All_meta,
                cluster_cols = FALSE,
                cluster_rows = FALSE,
                show_colnames = FALSE,
                legend_breaks = seq(0.1,0.8,by=0.1),
                #color = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "gray95")))(100),
                colorRampPalette(c("gray97","goldenrod1","firebrick"))(100),
                treeheight_row = 0,
                treeheight_col = 5,
                main = "Top30 Specises")
ggsave(filename='data/p_heat_plot_top30.png', plot=pic, width = 11, height = 18, units = "in", dpi = 300)

### seqSight - BGC ###

Maaslin2::Maaslin2(
  input_data  = t(All_bgc),
  input_metadata =  All_meta,
  output =  'data/maaslin2_bgc_output',
  max_significance= 0.1,
  analysis_method = 'LM',
  standardize = FALSE,
  transform = 'LOG',
  normalization = 'NONE',
  heatmap_first_n = 20,
  reference = c("response,R"))

### FDR concept
masslin2_stats  <- read.delim(
  "data/maaslin2_bgc_output/significant_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = 1
)

filtered_All_bgc_data <- All_bgc[rowSums(log(All_bgc+1)[])>0,]
bgc_top_lst <- as.matrix(sort(rowSums(filtered_All_bgc_data),decreasing = TRUE)[1:50])
All_BGC_top50 <- filtered_All_bgc_data[rownames(filtered_All_bgc_data) %in% rownames(bgc_top_lst),]
write.table(All_BGC_top50, file = "Top50_N_NR_sorted_BGC_mergedTable.tsv", col.names = NA, sep="\t")

new_fakepcl <- list(meta=All_meta, x=as.matrix(t(All_BGC_top50)),
                    ns=dim(All_BGC_top50)[2], nf=dim(All_BGC_top50)[1])

new_BGC_heat_plot <- omicsArt:::pcl.heatmap(new_fakepcl, sqrtspace = T, gamma = 2, meta= T, 
                                            show_colnames = F, show_rownames = T,
                                            treeheight_row = 0, treeheight_col= 5)
ggsave(filename='data/All_heat_plot_bgc_filtered.png', plot=new_BGC_heat_plot, width = 11, height = 18, units = "in", dpi = 300)

# fakepcl_bgc <- list(meta=All_meta, x=as.matrix(t(All_bgc)),
#                 ns=dim(All_bgc)[2], nf=dim(All_bgc)[1])
# 
# heat_plot_bgc <- omicsArt:::pcl.heatmap(fakepcl_bgc, sqrtspace = T, gamma = 2, meta= T, 
#                                     show_colnames = F, show_rownames = T,
#                                     treeheight_row = 0, treeheight_col= 5)

ggsave(filename='data/All_heat_plot_bgc.png', plot=heat_plot_bgc, width = 11, height = 18, units = "in", dpi = 300)


pic_bgc <- pheatmap(as.matrix(All_BGC_top50),
                annotation_col = All_meta,
                cluster_cols = FALSE,
                cluster_rows = FALSE,
                show_colnames = FALSE,
                legend_breaks = seq(-0.1,0.8,by=0.1),
                color = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "white")))(100),
                treeheight_row = 2,
                treeheight_col = 3,
                main = "Top30 BGC")
ggsave(filename='data/bgc_heat_plot_top30.png', plot=pic_bgc, width = 18, height = 18, units = "in", dpi = 300)



#####NEW TEST for species########
Top30_R <- read.delim(
  "data/Top3_R_abd_sorted_numeric.tsv",
  sep = ',',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Top30_R <- Top30_R ** (1/2)
new_Rpic <- pheatmap(as.matrix(Top30_R),
                annotation_col = All_meta,
                border_color = NA,
                cellwidth = 3,
                cellheight = 20,
                cluster_cols = FALSE,
                cluster_rows = FALSE,
                show_colnames = FALSE,
                show_rownames = FALSE,
                legend_breaks = seq(0.1,1,by=0.1),
                #color = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "gray95")))(100),
                # colorRampPalette(c("gray97","goldenrod1","firebrick"))(100),
                colorRampPalette(c("gray97","darkblue","red"))(100),
                treeheight_row = 0,
                treeheight_col = 5,
                main = "R-Top30Specises")
ggsave(filename='data/R_new_pic_sp.png', plot=new_Rpic, width = 20, height = 18, units = "in", dpi = 300)


Top30_NR <- read.delim(
  "data/Top3_NR_abd_sorted_numeric_ascending.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Top30_NR <- Top30_NR ** (1/2)
newNR_pic <- pheatmap(as.matrix(Top30_NR),
                    annotation_col = All_meta,
                    annotation_legend = FALSE,
                    border_color = NA,
                    cellwidth = 3,
                    cellheight = 20,
                    cluster_cols = FALSE,
                    cluster_rows = FALSE,
                    show_colnames = FALSE,
                    legend_breaks = seq(0.1,1,by=0.1),
                    legend = FALSE,
                    colorRampPalette(c("gray97","darkblue","red"))(100),
                    main = "NR-Top30Specises")
ggsave(filename='data/NR_new_pic.png', plot=newNR_pic, width = 20, height = 18, units = "in", dpi = 300)

p5 <- as.ggplot(new_Rpic)
p6 <- as.ggplot(newNR_pic)

# use patchwork to arrange them together
combine3 <- p6 + p5
ggsave(filename='data/combine3.png', plot=combine3, width = 32, height = 20, units = "in", dpi = 300)

#####NEW TEST for bgc########
Top50_R_bgc <- read.delim(
  "data/Top50_R_bgc_sorted_numeric_descending.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
new_pic_BGC <- pheatmap(as.matrix(Top50_R_bgc),
                    annotation_col = All_meta,
                    border_color = NA,
                    cellwidth = 3,
                    cellheight = 20,
                    cluster_cols = FALSE,
                    cluster_rows = FALSE,
                    show_colnames = FALSE,
                    show_rownames = FALSE,
                    legend_breaks = seq(0,0.01,by=0.001),
                    colorRampPalette(c("gray97","darkblue","red"))(100),
                    treeheight_row = 0,
                    treeheight_col = 5,
                    main = "R-Top50BGC")
ggsave(filename='data/R_new_pic_bgc.png', plot=new_pic_BGC, width = 20, height = 18, units = "in", dpi = 300)

Top50_NR_bgc <- read.delim(
  "data/Top50_NR_bgc_sorted_numeric_ascending.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
newNR_pic_bgc <- pheatmap(as.matrix(Top50_NR_bgc),
                      annotation_col = All_meta,
                      annotation_legend = FALSE,
                      border_color = NA,
                      cellwidth = 3,
                      cellheight = 20,
                      cluster_cols = FALSE,
                      cluster_rows = FALSE,
                      show_colnames = FALSE,
                      legend_breaks = seq(0.1,0.8,by=0.1),
                      legend = TRUE,
                      #color = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "gray95")))(100),
                      colorRampPalette(c("gray97","darkblue","red"))(100),
                      main = "NR-Top50BGC")
ggsave(filename='data/NR_new_pic_bgc.png', plot=newNR_pic_bgc, width = 20, height = 18, units = "in", dpi = 300)

p7 <- as.ggplot(new_pic_BGC)
p8 <- as.ggplot(newNR_pic_bgc)

# use patchwork to arrange them together
combine4 <- p8 + p7
ggsave(filename='data/combine4.png', plot=combine4, width = 31, height = 20, units = "in", dpi = 300)

### Maaslin2 significant sp results###

significant_R_sp <- read.delim(
  "data/R_significant_sp_abd_descending.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
significant_R_sp_pic <- pheatmap(as.matrix(significant_R_sp),
                        annotation_col = All_meta,
                        border_color = NA,
                        cellwidth = 3,
                        cellheight = 20,
                        cluster_cols = FALSE,
                        cluster_rows = FALSE,
                        show_colnames = FALSE,
                        show_rownames = FALSE,
                        legend_breaks = seq(0.1,0.8,by=0.1),
                        colorRampPalette(c("gray98","darkblue","red"))(100),
                        # colorRampPalette(c("gray97","darkblue","red"))(100),
                        treeheight_row = 0,
                        treeheight_col = 5,
                        main = "R-significant")
ggsave(filename='data/significant_R_sp_pic.png', plot=significant_R_sp_pic, width = 20, height = 25, units = "in", dpi = 300)

significant_NR_sp <- read.delim(
  "data/NR_significant_sp_abd_ascending.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
significant_NR_sp_pic <- pheatmap(as.matrix(significant_NR_sp),
                                  annotation_col = All_meta,
                                  annotation_legend = FALSE,
                                  border_color = NA,
                                  cellwidth = 3,
                                  cellheight = 20,
                                  cluster_cols = FALSE,
                                  cluster_rows = FALSE,
                                  show_colnames = FALSE,
                                  legend_breaks = seq(0.1,0.8,by=0.1),
                                  legend = TRUE,
                                  colorRampPalette(c("gray98","darkblue","red"))(100),
                                  main = "NR-significant")
ggsave(filename='data/significant_NR_sp_pic.png', plot=significant_NR_sp_pic, width = 20, height = 25, units = "in", dpi = 300)
# significant_R_bgc_pic + significant_NR_bgc_pic
# save plots
p3 <- as.ggplot(significant_R_sp_pic)
p4 <- as.ggplot(significant_NR_sp_pic)

# use patchwork to arrange them together
combine2 <- p4 + p3
ggsave(filename='data/combine2.png', plot=combine2, width = 30, height = 40, units = "in", dpi = 300)

### Maaslin2 significant bgc results###
significant_R_bgc <- read.delim(
  "data/R_significant_bgc_no_zero_descending.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
significant_R_bgc_pic <- pheatmap(as.matrix(significant_R_bgc),
                                 annotation_col = All_meta,
                                 border_color = NA,
                                 cellwidth = 3,
                                 cellheight = 20,
                                 cluster_cols = FALSE,
                                 cluster_rows = FALSE,
                                 show_colnames = FALSE,
                                 show_rownames = FALSE,
                                 legend_breaks = seq(0,0.01,by=0.001),
                                 colorRampPalette(c("yellow","darkblue","red"))(100),
                                 main = "R-bgc-significant")
ggsave(filename='data/significant_R_bgc_pic_nonzero.png', plot=significant_R_bgc_pic, width = 13, height = 18, units = "in", dpi = 300)

significant_NR_bgc <- read.delim(
  "data/NR_significant_bgc_no_zero_ascending.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
significant_NR_bgc_pic <- pheatmap(as.matrix(significant_NR_bgc),
                                  annotation_col = All_meta,
                                  annotation_legend = FALSE,
                                  border_color = NA,
                                  cellwidth = 3,
                                  cellheight = 20,
                                  cluster_cols = FALSE,
                                  cluster_rows = FALSE,
                                  show_colnames = FALSE,
                                  legend_breaks = seq(0,0.01,by=0.001),
                                  legend = FALSE,
                                  colorRampPalette(c("yellow","darkblue","red"))(100),
                                  main = "NR-bgc-significant")
ggsave(filename='data/significant_NR_bgc_pic_nonzero.png', plot=significant_NR_bgc_pic, width = 13, height = 18, units = "in", dpi = 300)


# significant_R_bgc_pic + significant_NR_bgc_pic
# save plots
p1 <- as.ggplot(significant_R_bgc_pic)
p2 <- as.ggplot(significant_NR_bgc_pic)

# use patchwork to arrange them together
combine <- p2 + p1
ggsave(filename='data/combine.png', plot=combine, width = 13, height = 18, units = "in", dpi = 300)
 