library(tidyverse)
library(vegan)
library(ggpubr)

setwd("~/Downloads/")
Frankel1 <- read.delim(
  "species_Frankel1.tsv",
  sep = ',',
  header = TRUE,
  fill = FALSE,
  comment.char = "",
  check.names = FALSE,
  row.names = 1
)

Frankel2 <- read.delim(
  "species_Frankel2.tsv",
  sep = ',',
  header = TRUE,
  fill = FALSE,
  comment.char = "",
  check.names = FALSE,
  row.names = 1
)

df <- Frankel1 %>%
  rownames_to_column(var = 'species') %>%
  gather(sample, 'db1', -species) %>%
  mutate(sample = unlist(strsplit(sample, '_'))[c(FALSE, TRUE)]) %>%
  inner_join(
    Frankel2 %>%
      rownames_to_column(var = 'species') %>%
      gather(sample, 'db2', -species) %>%
      mutate(sample = unlist(strsplit(sample, '_'))[c(FALSE, TRUE)])
  )

df <- df %>%
  mutate(x_val = log(db1 + 1),
         y_val = log(db2 + 1),
         group = case_when(
           x_val == 0 & y_val != 0 ~ "green",  # only x is zero -> green
           y_val == 0 & x_val != 0 ~ "red",    # only y is zero -> red
           TRUE ~ "gray"                      # all other cases (including both zero) -> gray
         ))

df$group <- factor(df$group, levels = c("green", "red", "gray"))

green_count <- sum(df$group == "green")
red_count   <- sum(df$group == "red")
gray_count  <- sum(df$group == "gray")

gray_df <- df %>% filter(group == "gray")
gray_linear_m <- lm(x_val ~ y_val - 1, data = gray_df)
gray_rsq <- summary(gray_linear_m)$adj.r.squared


scatter <- ggplot(df, aes(x = x_val, y = y_val, color = group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', data = gray_df, color = "black") +
  xlab('Log(MetaPhlAn4)') +
  ylab('Log(MetaPhlAn2)') +
  annotate('text', x = 1.3, y = 4,
           label = paste('R2 (MetaPhlAn2 vs. MetaPhlAn4) =', round(gray_rsq, 2))) +
  theme_classic() +
  omicsArt::theme_omicsEye() +
  scale_color_manual(values = c("green" = "green", "red" = "red", "gray" = "gray"),
                     labels = c(
                       paste0("Detected by MetaPhlAn2 (n = ", green_count, ")"),
                       paste0("Detected by MetaPhlAn4 (n = ", red_count, ")"),
                       paste0("Detected by MetaPhlAn2&4 (n = ", gray_count, ")")
                     )) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_classic() +
  theme(legend.position = c(1.15, 0.01),        
        legend.justification = c(1.15, 0.01),         
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA))         

scatter

ggsave(filename = 'scatterplot.png', plot = scatter,
       width = 4, height = 3, units = 'in', dpi = 300)