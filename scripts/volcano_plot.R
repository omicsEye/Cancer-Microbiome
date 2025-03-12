library(omicsArt)
library(ggplot2)

## Load files
pathway_all_results <- read.delim(
  "~/pathway_batch_REML.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = 1
)
pathway_sig_results <- read.delim(
  "~/pathway_sig_results.csv",
  sep = ',',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = 1
)
# Load necessary libraries
# install.packages("ggplot2")
# install.packages("ggrepel")
library(ggplot2)
library(ggrepel)

# Read the data from the TSV file
# data <- read.delim("Cancer-Microbiome/maaslin2_pathway_output/all_results.tsv", header = TRUE)
data <- pathway_all_results
data$std_coef <- data$coef / data$stderr

# Ensure the data contains a column named 'pathway' for pathway names
if (!"feature" %in% colnames(data)) {
  stop("The data must contain a column named 'feature' for pathway names.")
}

# Define the volcano_plot function
volcano_plot <- function(stats_table,
                         threshold = 0.05,
                         method = 'nominal',
                         pvalue_col = "pval",
                         fdr_col = "qval",
                         orderby = 'std_coef',
                         y_label = "-log(p-val)",
                         x_label = 'Coefficient') {
  colnames(stats_table) <- gsub(orderby, "std_coef", colnames(stats_table))
  colnames(stats_table) <- gsub(pvalue_col, "P.Value", colnames(stats_table))
  colnames(stats_table) <- gsub(fdr_col, "fdr", colnames(stats_table))
  
  if (method == 'nominal')
    stats_table$fdr <- stats_table$P.Value
  
  stats_table$feature <- stats_table$feature
  stats_table$feature[stats_table$fdr >= threshold] <- NA
  
  p <- ggplot(stats_table, aes(
    x = std_coef,
    y = -log10(P.Value),
    label = feature
  )) +
    scale_fill_gradient(low = "lightgray", high = "navy") +
    scale_color_gradient(low = "lightgray", high = "navy") +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", show.legend = FALSE) +
    geom_point(
      data = subset(stats_table, fdr <= threshold & std_coef > 0.0),
      fill = "green",
      color = 'black',
      alpha = .5,
      shape = 21,
      size = 3,
      stroke = 0.05
    ) +
    geom_point(
      data = subset(stats_table, fdr <= threshold & std_coef < 0.0),
      fill = 'red',
      color = 'black',
      alpha = .5,
      shape = 21,
      size = 3,
      stroke = 0.05
    ) +
    geom_point(
      data = subset(stats_table, fdr > threshold),
      fill =  'gray',
      color = "black",
      alpha = 0.2,
      shape = 21,
      size = 1,
      stroke = 0.05
    ) +
    geom_vline(xintercept = 0, size = 0.1) +
    geom_hline(yintercept = 0, size = 0.1) +
    geom_hline(
      yintercept = -log10(threshold),
      linetype = "dashed",
      size = 0.1,
      color = "red"
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      size = 0.1,
      color = "red"
    ) +
    theme_linedraw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 6), 
          axis.text.y = element_text(size = 6)) +
    xlab(x_label) +
    ylab(y_label) +
    annotate(
      "text",
      x = min(stats_table$std_coef, na.rm = TRUE) + .1,
      y = min(-log10(stats_table$P.Value), na.rm = TRUE) + .75,
      label = "Significant up",
      size = 3,
      color = "black",
      hjust = 0
    ) +
    annotate(
      "point",
      x = min(stats_table$std_coef, na.rm = TRUE),
      y = min(-log10(stats_table$P.Value), na.rm = TRUE) + .75,
      color = "green"
    ) +
    annotate(
      "text",
      x = min(stats_table$std_coef, na.rm = TRUE) + .1,
      y = min(-log10(stats_table$P.Value), na.rm = TRUE) + .5,
      label = "Significant down",
      size = 3,
      color = "black",
      hjust = 0
    ) +
    annotate(
      "point",
      x = min(stats_table$std_coef, na.rm = TRUE),
      y = min(-log10(stats_table$P.Value), na.rm = TRUE) + .5,
      color = "red"
    ) +
    annotate(
      "text",
      x = min(stats_table$std_coef, na.rm = TRUE),
      y = min(-log10(stats_table$P.Value), na.rm = TRUE) + 1.0,
      label = paste(method, " threshold: ", threshold, sep = ""),
      size = 3,
      color = "black",
      hjust = 0
    ) +
    geom_text_repel(
      size = 2,
      force = 1,
      box.padding = 0.5,  # Increase padding around text labels
      point.padding = 0.3,  # Increase padding around points
      max.overlaps = 6,  # Control the maximum number of overlapping labels to avoid
      nudge_y = 0.2,  # Slightly nudge labels to avoid overlap
      nudge_x = 0.1,  # Slightly nudge labels to avoid overlap
      fontface = "italic"
    )
  return (p)
}

# Call the function to create the volcano plot
volcano <- volcano_plot(data, threshold = 0.05, method = 'nominal', pvalue_col = "pval", fdr_col = "qval", orderby = 'std_coef')
ggsave(filename='~volcano_plot.pdf',plot=volcano, width = 7.2, height = 4.0, units = "in", dpi = 350)
