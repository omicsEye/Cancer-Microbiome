library(tidyverse)
library(vegan)
library(ggpubr)
setwd("/Users/xinyang/Library/CloudStorage/Box-Box/Cancer_Mircrobiome/MetaAnalysis")
Frankel1 <- read.delim(
  "data/species_Frankel1.tsv",
  sep = ',',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

Frankel2 <- read.delim(
  "data/species_Frankel2.tsv",
  sep = ',',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

df <- Frankel1 %>% 
  rownames_to_column(var='species') %>% 
  gather(sample, 'db1', -species) %>% 
  mutate(sample = unlist(strsplit(sample, '_'))[c(F,T)]) %>%
  inner_join(
    Frankel2 %>% 
      rownames_to_column(var='species') %>% 
      gather(sample, 'db2', -species) %>% 
      mutate(sample = unlist(strsplit(sample, '_'))[c(F,T)])
)

linear_m = lm(log(db1+1)~log(db2+1)-1, data = df)
summary(linear_m)

rsq = summary(linear_m)$adj.r.squared

scatter = ggplot(df,aes(log(db1+1), log(db2+1))) +
  geom_point() +
  geom_smooth(method='lm') +
  xlab('Log(Metaphlan4)') +
  ylab('Log(Metaphlan2)') +
  annotate('text',x = 0.1, y = 4,
           label = paste('R2 = ', round(rsq, 2)))+
  theme_classic()
scatter
ggsave(filename = 'scatterplot.png', plot = scatter,
       width = 7.2, height = 5, units = 'in', dpi = 300)


ord2 <- decorana(Frankel2)

common_name <- intersect(rownames(Frankel1),rownames(Frankel2))

Frankel1 <- Frankel1[which(rownames(Frankel1)%in%common_name),]
Frankel2 <- Frankel2[which(rownames(Frankel2)%in%common_name),]
combine_Frankel <- cbind(Frankel1, Frankel2)

dist()


