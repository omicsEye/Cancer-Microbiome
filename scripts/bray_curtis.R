library(vegan)
library(phyloseq)
library(ggplot2)
library(plyr)
library(ade4)

# Frankel Study

Frankel_M4 <- read.delim(
  "Metaphlan4/Metaphlan4_Frankel2017.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Frankel_M3 <- read.delim(
  "CMD/CMD_Frankel2017.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

# calculate Bray-Curtis distance using the vegan package
Frankel_M4_dist_bc <- as.matrix(vegdist(t(Frankel_M4), method = "bray")) 
Frankel_M4_dist_bc

Frankel_M3_dist_bc <- as.matrix(vegdist(t(Frankel_M3), method = "bray")) 
Frankel_M3_dist_bc

Frankel_M4_dist_bc <- as.dist(Frankel_M4_dist_bc)
Frankel_M3_dist_bc <- as.dist(Frankel_M3_dist_bc)

# mantel.rtest(station.dists, ozone.dists, nrepet = 9999)
# plot(Frankel_r1 <- mantel.rtest(Frankel_M4_dist_bc,Frankel_M3_dist_bc), main = "Frankel Mantel's test")

plot(Frankel_r1 <- mantel.rtest(Frankel_M4_dist_bc,Frankel_M3_dist_bc, nrepet = 9999), main = "Frankel Mantel's test")
Frankel_r1

# Gop Study

Gop_M4 <- read.delim(
  "Metaphlan4/Metaphlan4_Gop2018.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Gop_M3 <- read.delim(
  "CMD/CMD_Gop2018.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

# calculate Bray-Curtis distance using the vegan package
Gop_M4_dist_bc <- as.matrix(vegdist(t(Gop_M4), method = "bray")) 
Gop_M4_dist_bc

Gop_M3_dist_bc <- as.matrix(vegdist(t(Gop_M3), method = "bray")) 
Gop_M3_dist_bc

Gop_M4_dist_bc <- as.dist(Gop_M4_dist_bc)
Gop_M3_dist_bc <- as.dist(Gop_M3_dist_bc)
plot(Gop_r1 <- mantel.rtest(Gop_M4_dist_bc,Gop_M3_dist_bc, nrepet = 9999), main = "Gop Mantel's test")
Gop_r1


# Lee Study

Lee_M4 <- read.delim(
  "Metaphlan4/Metaphlan4_Lee2022.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Lee_M3 <- read.delim(
  "CMD/CMD_Lee2022.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

# calculate Bray-Curtis distance using the vegan package
Lee_M4_dist_bc <- as.matrix(vegdist(t(Lee_M4), method = "bray")) 
Lee_M4_dist_bc

Lee_M3_dist_bc <- as.matrix(vegdist(t(Lee_M3), method = "bray")) 
Lee_M3_dist_bc

Lee_M4_dist_bc <- as.dist(Lee_M4_dist_bc)
Lee_M3_dist_bc <- as.dist(Lee_M3_dist_bc)
plot(Lee_r1 <- mantel.rtest(Lee_M4_dist_bc,Lee_M3_dist_bc, nrepet = 9999), main = "Lee Mantel's test")
Lee_r1

# Matson Study

Matson_M4 <- read.delim(
  "Metaphlan4/Metaphlan4_Mat2018.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Matson_M3 <- read.delim(
  "CMD/CMD_Mat2018.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

# calculate Bray-Curtis distance using the vegan package
Matson_M4_dist_bc <- as.matrix(vegdist(t(Matson_M4), method = "bray")) 
Matson_M4_dist_bc

Matson_M3_dist_bc <- as.matrix(vegdist(t(Matson_M3), method = "bray")) 
Matson_M3_dist_bc

Matson_M4_dist_bc <- as.dist(Matson_M4_dist_bc)
Matson_M3_dist_bc <- as.dist(Matson_M3_dist_bc)
plot(Mat_r1 <- mantel.rtest(Matson_M4_dist_bc,Matson_M3_dist_bc, nrepet = 9999), main = "Matson Mantel's test")
Mat_r1

# Daver
Daver_M4 <- read.delim(
  "Metaphlan4/Metaphlan4_Daver2021.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

# Spenncer
Spe_M4 <- read.delim(
  "Metaphlan4/Metaphlan4_Spe2021.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

# Baruch
Baruch_M4 <- read.delim(
  "Metaphlan4/Metaphlan4_Baruch2021.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)



Frankel_M4$avg_of_ra <- rep(c(rowSums(Frankel_M4)/(ncol(Frankel_M4)-1)))
Gop_M4$avg_of_ra <- rep(c(rowSums(Gop_M4)/(ncol(Gop_M4)-1)))
Lee_M4$avg_of_ra <- rep(c(rowSums(Lee_M4)/(ncol(Lee_M4)-1)))
Matson_M4$avg_of_ra <- rep(c(rowSums(Matson_M4)/(ncol(Matson_M4)-1)))
Baruch_M4$avg_of_ra <- rep(c(rowSums(Baruch_M4)/(ncol(Baruch_M4)-1)))
Spe_M4$avg_of_ra <- rep(c(rowSums(Spe_M4)/(ncol(Spe_M4)-1)))
Daver_M4$avg_of_ra <- rep(c(rowSums(Daver_M4)/(ncol(Daver_M4)-1)))


Frankel_rowname <- rownames(subset(Frankel_M4, avg_of_ra >0.4))
Gop_rowname <- rownames(subset(Gop_M4, avg_of_ra >0.4))
Lee_rowname <- rownames(subset(Lee_M4, avg_of_ra >0.4))
Matson_rowname <- rownames(subset(Matson_M4, avg_of_ra >0.4))

Baruch_rowname <- rownames(subset(Baruch_M4, avg_of_ra >0.4))
Spe_rowname <- rownames(subset(Spe_M4, avg_of_ra >0.4))
Daver_rowname <- rownames(subset(Daver_M4, avg_of_ra >0.4))

common_species <- Reduce(intersect, list(Frankel_rowname,Gop_rowname,Lee_rowname,
                                         Matson_rowname, Baruch_rowname,
                                         Spe_rowname,Daver_rowname))
common_species


# Add species to individual BGC table which share the same sample id
# Frankel
Frankel_BGC <- read.delim(
  "Data/BGC_Frankel_Merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Frankel_BGC[is.na(Frankel_BGC)] <- 0
# Frankel_species_common <- Frankel_M4[rownames(Frankel_M4) %in% common_species,
#                                      !(names(Frankel_M4) %in% "avg_of_ra")]/100  
# Frankel_BGC_species <- rbind(Frankel_BGC,Frankel_species_common)
# Frankel_BGC_species_dist_bc <- as.matrix(vegdist(Frankel_BGC_species, method = "bray")) 
# library(ComplexHeatmap)
# Heatmap(Frankel_BGC_species_dist_bc,heatmap_legend_param = list(at = seq(0,1, 0.1)))
# 
# 
# selected_dist_bc <- Frankel_BGC_species_dist_bc[,c("s__Faecalibacterium_prausnitzii",
#                                "s__Akkermansia_muciniphila",
#                                "s__Phocaeicola_vulgatus",
#                                "s__Bacteroides_uniformis",
#                                "s__Alistipes_putredinis",
#                                "s__Eubacterium_rectale",
#                                "s__Prevotella_copri_clade_A",
#                                "s__Ruminococcus_bicirculans",
#                                "s__Phocaeicola_dorei")]
# subset(selected_dist_bc, avg_of_ra >0.4)
# Heatmap(selected_dist_bc)
# plot(selected_dist_bc)


# Frankel + bgc_Annotation
library(taxonomizr)

Frankel_BGC_Taxa <- read.delim(
  "Data/BGC_Frankel_TaxaInfo.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Frankel_BGC_RA <- Frankel_BGC
Frankel_col_taxa <- as.list(Frankel_BGC_Taxa$TaxaId)
li2 <- list()
li3 <- list()
for (id in Frankel_col_taxa){
  linage <- getTaxonomy(id,'/Users/xinyang/accessionTaxa.sql')
  li2 <- append(li2,linage[,"phylum"])
  li3 <- append(li3,linage[,"family"])
}

Frankel_BGC_RA$phylum <- unlist(li2)
Frankel_BGC_RA$family <- unlist(li3)
# Frankel_BGC_RA <- as.data.frame(lapply(Frankel_BGC_RA, unlist))
# rownames(Frankel_BGC_RA) <- Frankel_BGC_RA$Family
print(unique(Frankel_BGC_RA$family))

Frankel_BGC_RA <- Frankel_BGC_RA[order(Frankel_BGC_RA$phylum),]
Frankel_BGC_RA_family <- subset(Frankel_BGC_RA,select = -c(40))
Frankel_BGC_RA_phylum <- subset(Frankel_BGC_RA,select = -c(41))
Frankel_BGC_RA_data <- subset(Frankel_BGC_RA,select = -c(41,40))

Frankel_BGC_RA_dist_bc <- as.matrix(vegdist(Frankel_BGC_RA_data),method = "bray")
library(ComplexHeatmap)
Heatmap(Frankel_BGC_RA_dist_bc,
        heatmap_legend_param = list(at = seq(0,1, 0.1)),
        row_title_gp = gpar(fontsize = 1),
        column_title_gp = gpar(fontsize = 1))

# library(devtools)
# install_github('omicsEye/omicsArt', force = TRUE)
library(omicsArt)
Frankel_BGC_RA_metadata <- Frankel_BGC_RA[,c("phylum","family")]

pcoa_plots_Frankel_BGC <- omicsArt::ordplots(Frankel_BGC_RA_dist_bc, Frankel_BGC_RA_metadata, 
                                 output = 'Data/', outputname = 'pcoa_BGC', method = 'pcoa')

tsne_plots_Frankel_BGC <- omicsArt::ordplots(Frankel_BGC_RA_dist_bc, Frankel_BGC_RA_metadata, 
                                 output = 'Data/', outputname = 'tsne_BGC', method = 'tsne')
# tsne_plots_Frankel_BGC_Heatmap <- omicsArt::heatmap()
#   ordplots(Frankel_BGC_RA_dist_bc, Frankel_BGC_RA_metadata, 
#          output = 'Data/', outputname = 'tsne_BGC', method = 'tsne')

pcoa_plots_Frankel_taxa <- omicsArt::ordplots(Frankel_M4_dist_bc, Frankel_BGC_RA_metadata, 
                                              output = 'Data/', outputname = 'pcoa_taxa', method = 'pcoa')
tsne_plots_Frankel_taxa <- omicsArt::ordplots(Frankel_M4_dist_bc, Frankel_BGC_RA_metadata, 
                                             output = 'Data/', outputname = 'tsne_taxa', method = 'tsne')

# tsne_plots_Frankel_taxa_Heatmap <- omicsArt::heatmap()
# ordplots(Frankel_M4_dist_bc, Frankel_BGC_RA_metadata, 
#            output = 'Data/', outputname = 'tsne_BGC', method = 'tsne')
Frankel_BGC_RA_data

filtered_Frankel_BGC_RA_data <- Frankel_BGC_RA_data[rowSums(log(Frankel_BGC_RA_data+1)[])>0,]

### draft to keep only thos BGCs that has max value atleast 10% 
threshold = 0.1
max_values = apply(filtered_Frankel_BGC_RA_data, 1, max)

new_data = filtered_Frankel_BGC_RA_data[max_values >= threshold, ]


new_metadata <- Frankel_BGC_RA_metadata[rownames(new_data) %in% rownames(new_data),]

new_fakepcl <- list(meta=new_metadata, x=as.matrix(new_data),
                    ns=dim(new_data)[2], nf=dim(new_data)[1])

new_BGC_heat_plot <- omicsArt:::pcl.heatmap(new_fakepcl, sqrtspace = T, gamma = 2, meta= T, 
                                            show_colnames = F, show_rownames = T,
                                            treeheight_row = 0, treeheight_col= 5 )


ggsave(filename='data/Frankel_BGC_heatmap.png', plot=new_BGC_heat_plot, width = 11, height = 18, units = "in", dpi = 300)
####
filtered_Frankel_BGC_RA_metadata <- Frankel_BGC_RA_metadata[rownames(Frankel_BGC_RA_metadata) 
                                                            %in% rownames(filtered_Frankel_BGC_RA_data),]

fakepcl <- list(meta=filtered_Frankel_BGC_RA_metadata, x=as.matrix(filtered_Frankel_BGC_RA_data),
                  ns=dim(filtered_Frankel_BGC_RA_data)[2], nf=dim(filtered_Frankel_BGC_RA_data)[1])
BGC_heat_plot <- omicsArt:::pcl.heatmap(fakepcl, sqrtspace = T, gamma = 2, meta= T, 
                                              show_colnames = F, show_rownames = T,
                                              treeheight_row = 0, treeheight_col= 5 )
ggsave(filename='data/BGC_heatmap.png', plot=BGC_heat_plot, width = 11, height = 18, units = "in", dpi = 300)

Full_meta <- read.delim(
  "data/Full_Metadata",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
)
Frankel_Taxa_metadata <- subset(Full_meta,dataset == "Frankel_2017",select = c("sampleid","response","Sex"))
rownames(Frankel_Taxa_metadata) = Frankel_Taxa_metadata[,1]
Frankel_Taxa_metadata <- Frankel_Taxa_metadata[,-1]
Frankel_Taxa_metadata <- Frankel_Taxa_metadata[order(Frankel_Taxa_metadata$response,Frankel_Taxa_metadata$Sex),]
Frankel_M4_sorted <- Frankel_M4[,rownames(Frankel_Taxa_metadata)]
# calculate Bray-Curtis distance using the vegan package
Frankel_M4_dist_bc_sorted <- as.matrix(vegdist(t(Frankel_M4_sorted), method = "bray")) 
Frankel_M4_dist_bc_sorted

speices_heat_plot2 <- pheatmap(Frankel_M4_dist_bc_sorted,
         annotation_row = Frankel_Taxa_metadata,
         cutree_rows = 2,
         cluster_rows = FALSE
         )

fakepcl_taxa <- list(meta=Frankel_Taxa_metadata, x=as.matrix(Frankel_M4_dist_bc_sorted),
                ns=dim(data)[2], nf=dim(data)[1])
speices_heat_plot <- omicsArt:::pcl.heatmap(fakepcl_taxa, sqrtspace = T, gamma = 2, meta= T, 
                                            show_colnames = F, show_rownames = T,
                                            treeheight_row = 0, treeheight_col= 5,
                                            cluster_cols = FALSE)

ggsave(filename='data/species_heatmap.png', plot=speices_heat_plot, width = 11, height = 18, units = "in", dpi = 300)


library(Maaslin2)
Maaslin2::Maaslin2(
  input_data  = t(Frankel_BGC_RA_without_phylum),
  input_metadata =  Frankel_Taxa_metadata,
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
  "data/maaslin2_output/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = 1
)


h_pval <- hist(masslin2_stats$pval,breaks=seq(0,1,0.05))

polygon(c(0,0.05,0.05,0),c(0,0,h_pval$counts[1],h_pval$counts[1]),col="grey")

h_qval <- hist(masslin2_stats$qval,breaks=seq(0,1,0.1))
polygon(c(0,0.1,0.1,0),c(0,0,h_qval$counts[1],h_qval$counts[1]),col="grey")

# taxa<-getTaxonomy(c("1955450"),'/Users/xinyang/accessionTaxa.sql')
# print(taxa)

# Gop
Gop_BGC <- read.delim(
  "Data/BGC_Gop_Merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Gop_BGC[is.na(Gop_BGC)] <- 0
Gop_species_common <- Gop_M4[rownames(Gop_M4) %in% common_species, ]  
Gop_BGC_species <- rbind(Gop_BGC,Gop_species_common)
Gop_BGC_species_dist_bc <- as.matrix(vegdist(Gop_BGC_species, method = "bray")) 
Gop_BGC_species_dist_bc

# Matson
Matson_BGC <- read.delim(
  "Data/BGC_Matson_Merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

Matson_BGC[is.na(Matson_BGC)] <- 0
Matson_species_common <- Matson_M4[rownames(Matson_M4) %in% common_species, ]  
Matson_BGC_species <- rbind(Matson_BGC,Matson_species_common)
Matson_BGC_species_dist_bc <- as.matrix(vegdist(Matson_BGC_species, method = "bray")) 
Matson_BGC_species_dist_bc


# Lee
Lee_BGC <- read.delim(
  "Data/BGC_Lee_Merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Lee_BGC[is.na(Lee_BGC)] <- 0
Lee_species_common <- Lee_M4[rownames(Lee_M4) %in% common_species, ]  
Lee_BGC_species <- rbind(Lee_BGC,Lee_species_common)
Lee_BGC_species_dist_bc <- as.matrix(vegdist(Lee_BGC_species, method = "bray")) 
Lee_BGC_species_dist_bc

# Davar
Davar_BGC <- read.delim(
  "Data/BGC_Davar_Merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Davar_BGC[is.na(Davar_BGC)] <- 0
Davar_species_common <- Daver_M4[rownames(Daver_M4) %in% common_species, ]  
Davar_BGC_species <- rbind(Davar_BGC,Davar_species_common)
Davar_BGC_species_dist_bc <- as.matrix(vegdist(Lee_BGC_species, method = "bray")) 
Davar_BGC_species_dist_bc

# Spencer
Spe_BGC <- read.delim(
  "Data/BGC_Spencer_Merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Spe_BGC[is.na(Spe_BGC)] <- 0
Spe_species_common <- Spe_M4[rownames(Spe_M4) %in% common_species, ]  
Spe_BGC_species <- rbind(Spe_BGC,Spe_species_common)
Spe_BGC_species_dist_bc <- as.matrix(vegdist(Spe_BGC_species, method = "bray")) 
Spe_BGC_species_dist_bc

# Baruch
Baruch_BGC <- read.delim(
  "Data/BGC_Baruch_Merged.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
Baruch_BGC[is.na(Baruch_BGC)] <- 0
Baruch_species_common <- Baruch_M4[rownames(Baruch_M4) %in% common_species, ]  
Baruch_BGC_species <- rbind(Baruch_M4,Baruch_species_common)
Baruch_BGC_species_dist_bc <- as.matrix(vegdist(Baruch_BGC_species, method = "bray")) 
Baruch_BGC_species_dist_bc







