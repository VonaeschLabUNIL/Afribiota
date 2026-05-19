###########################################################
# Title: Analysis Of Metabolomics Data for the Study
# "Lithocholic acid boosts the growth of butyrate-producing bacteria and is decreased in the feces of stunted children" 
# Author: Mariia Beliaeva (mariia.beliaeva@embl.de)
# Description: R code for data analysis and plots
###########################################################

### Import libraries ----
library(tidyverse)
library(scales)
library(viridis)

## for plotting ##
scientific_10 <- function(x) {   
  parse(text = gsub("e\\+*", " %*% 10^", scientific_format()(x)))}

### READ RAW DATA ----

# Set paths to working directory and data folders
wd_path <- getwd()
setwd(wd_path)

df <- read_csv(paste0(wd_path, "/0_raw_clean_data.csv"))

df$sample <- factor(df$sample, levels = c("control","bacterial culture","supernatant"))
df$condition <- factor(df$condition, levels = c("media","S. flexneri","F. nucleatum",
                                                "S. salivarius", "R. intestinalis",
                                                "C. comes", "G. elegans"))

### PLOT FIGURES ----

#### Taurochenodeoxycholic Acid ----

for (m in unique(df$media)) {
  pdf(paste0(wd_path, "/plots/3_TaurochenodeoxycholicAcid_", m, ".pdf"), width = 10, height = 4)
  
  p <- df %>%
    filter(media == m) %>%
    filter(compound == "Taurochenodeoxycholic Acid") %>%
    filter(compound_analysed %in% c("Taurochenodeoxycholic Acid", "Chenodeoxycholic Acid")) %>%
    group_by(condition, sample, compound) %>%
    ggplot(aes(x = condition, y = area, fill = sample)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    scale_y_continuous(label = scientific_10, limits = c(0, NA)) +
    labs(x = "", y = "Peak Area") +
    facet_wrap(~compound_analysed, scales = "free") +
    theme_bw() +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 20, colour = "black"),
      axis.text.y = element_text(size = 14, vjust = 0.5, hjust = 1,
                                 face = "italic", colour = "black"),
      axis.text.x = element_text(size = 14, angle = 20, vjust = 0.5, hjust = 0.5,
                                 face = "italic", colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 20, colour = "black"),
      legend.position = "right",
      strip.text = element_text(size = 12, colour = "black"),
      plot.title = element_text(colour = "black")) +
    ggtitle(paste("samples treated with Taurochenodeoxycholic Acid in", m))
  
  print(p)
  dev.off()
}

#### Taurolithocholic Acid ----
for (m in unique(df$media)) {
  pdf(paste0(wd_path, "/plots/3_TaurolithocholicAcid_", m, ".pdf"), width = 10, height = 4)
  
  p <- df %>%
    filter(media == m) %>%
    filter(compound == "Taurolithocholic Acid") %>%
    filter(compound_analysed %in% c("Taurolithocholic Acid", "Lithocholic Acid")) %>%
    group_by(condition, sample, compound) %>%
    ggplot(aes(x = condition, y = area, fill = sample)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    scale_y_continuous(label = scientific_10, limits = c(0, NA)) +
    labs(x = "", y = "Peak Area") +
    facet_wrap(~compound_analysed, scales = "free") +
    theme_bw() +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 20, colour = "black"),
      axis.text.y = element_text(size = 14, vjust = 0.5, hjust = 1,
                                 face = "italic", colour = "black"),
      axis.text.x = element_text(size = 14, angle = 20, vjust = 0.5, hjust = 0.5,
                                 face = "italic", colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 20, colour = "black"),
      legend.position = "right",
      strip.text = element_text(size = 12, colour = "black"),
      plot.title = element_text(colour = "black")) +
    ggtitle(paste("samples treated with Taurolithocholic Acid in", m))
  
  print(p)
  dev.off()
}

#### Chenodeoxycholic Acid ----
for (m in unique(df$media)) {
  pdf(paste0(wd_path, "/plots/3_ChenodeoxycholicAcid_", m, ".pdf"), width = 10, height = 4)
  
  p <- df %>%
    filter(media == m) %>%
    filter(compound == "Chenodeoxycholic Acid") %>%
    filter(compound_analysed %in% c("Taurochenodeoxycholic Acid", "Chenodeoxycholic Acid")) %>%
    group_by(condition, sample, compound) %>%
    ggplot(aes(x = condition, y = area, fill = sample)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    scale_y_continuous(label = scientific_10, limits = c(0, NA)) +
    labs(x = "", y = "Peak Area") +
    facet_wrap(~compound_analysed, scales = "free") +
    theme_bw() +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 20, colour = "black"),
      axis.text.y = element_text(size = 14, vjust = 0.5, hjust = 1,
                                 face = "italic", colour = "black"),
      axis.text.x = element_text(size = 14, angle = 20, vjust = 0.5, hjust = 0.5,
                                 face = "italic", colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 20, colour = "black"),
      legend.position = "right",
      strip.text = element_text(size = 12, colour = "black"),
      plot.title = element_text(colour = "black")) +
    ggtitle(paste("samples treated with Chenodeoxycholic Acid in", m))
  
  print(p)
  dev.off()
}

#### Lithocholic Acid ----
for (m in unique(df$media)) {
  pdf(paste0(wd_path, "/plots/3_LithocholicAcid_", m, ".pdf"), width = 10, height = 4)
  
  p <- df %>%
    filter(media == m) %>%
    filter(compound == "Lithocholic Acid") %>%
    filter(compound_analysed %in% c("Taurolithocholic Acid", "Lithocholic Acid")) %>%
    group_by(condition, sample, compound) %>%
    ggplot(aes(x = condition, y = area, fill = sample)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    scale_y_continuous(label = scientific_10, limits = c(0, NA)) +
    labs(x = "", y = "Peak Area") +
    facet_wrap(~compound_analysed, scales = "free") +
    theme_bw() +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 20, colour = "black"),
      axis.text.y = element_text(size = 14, vjust = 0.5, hjust = 1,
                                 face = "italic", colour = "black"),
      axis.text.x = element_text(size = 14, angle = 20, vjust = 0.5, hjust = 0.5,
                                 face = "italic", colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 20, colour = "black"),
      legend.position = "right",
      strip.text = element_text(size = 12, colour = "black"),
      plot.title = element_text(colour = "black")) +
    ggtitle(paste("samples treated with Lithocholic Acid in", m))
  
  print(p)
  dev.off()
}

