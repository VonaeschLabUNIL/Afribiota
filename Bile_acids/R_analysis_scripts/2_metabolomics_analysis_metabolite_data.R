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
library(rstatix)
library(ggpubr)

### DATA ANALYSIS ----
#### READ DATA ----

# Set paths to working directory and data folders
wd_path <- getwd()
setwd(wd_path)

# This is the cleaned data without outliers, 
# which were identified in script 1 based on parent molecules and excluded  
df_full <- read_csv(paste0(wd_path, "/output/1_all_data_wo_outliers.csv"))

df_full$sample <- factor(df_full$sample, levels = c("control","bacterial culture","supernatant"))
df_full$condition <- factor(df_full$condition, levels = c("media","S. flexneri","F. nucleatum",
                                                          "S. salivarius", "R. intestinalis",
                                                          "C. comes", "G. elegans"))

#### EXTRACT METABOLITES ONLY ----
# Include hydrolysis products (conversion of Tauro-BA into BA)
# Exclude data where mean of area per combination == 5000 (no metabolite detected) 
# This excludes S. flexneri, F. nucleatum and S. salivarius - no hydrolysis metabolites detected there
df_metabolites <- df_full %>% 
  filter((compound == "Taurochenodeoxycholic Acid" & compound_analysed == "Chenodeoxycholic Acid") |
           (compound == "Taurolithocholic Acid" & compound_analysed == "Lithocholic Acid")) %>% 
  group_by(media, compound, compound_analysed) %>% 
  mutate(mean_area = mean(area)) %>% 
  filter(mean_area != 5000) %>% 
  ungroup() %>% 
  select(-mean_area) %>% 
  write_csv(paste0(wd_path, "/output/2_metabolite_data_wo_outliers.csv"))

#### CALCULATE P-VALUES ----
df <- df_metabolites %>% 
  mutate(condition = as.factor(condition),
         sample = as.factor(sample),
         position = as.factor(position),
         area = as.numeric(area)) 

df_controls <- df %>% 
  filter(sample == "control")

df_species <- df %>% 
  filter(sample != "control")

species <- as.character(unique(df$condition)[-1])
compounds <- unique(df$compound)

# Test whether each compound is different to media control and between supernatant and whole culture
df_anova <- data.frame()
df_anova_list <- list()
df_anova_list[[1]] <- df_anova

for (s in species) {
  # Find correct media for all species
  m <- df_species %>% 
    filter(condition%in%c(s)) %>% 
    pull(media) %>% 
    unique()
  
  for (c in compounds) {
    df_pval <- df_species %>% 
      filter(compound == c, media == m, condition == s)
    
    tmp <- df_controls %>% 
      filter(compound == c, media == m) %>% 
      bind_rows(df_pval)
    
    # If there is no variation (all areas == 5000, no compound detected), insert 1 as PVAL
    # Control vs bacterial culture
    control_bacteria_ttest <- tryCatch({tmp %>%
        filter(sample %in% c("bacterial culture", "control")) %>%
        mutate(sample = droplevels(sample)) %>%
        t_test(area ~ sample) %>%
        pull(p)}, error = function(e) {1})
    
    # Control vs supernatant
    control_supernatant_ttest <- tryCatch({tmp %>%
        filter(sample %in% c("supernatant", "control")) %>%
        mutate(sample = droplevels(sample)) %>%
        t_test(area ~ sample) %>%
        pull(p)}, error = function(e) {1})
    
    # Bacterial culture vs supernatant
    bacteria_supernatant_ttest <- tryCatch({tmp %>%
        filter(sample %in% c("bacterial culture", "supernatant")) %>%
        mutate(sample = droplevels(sample)) %>%
        t_test(area ~ sample) %>%
        pull(p)}, error = function(e) {1})
    
    df_anova <- df_pval %>% 
      mutate(PVAL_control_bacteria = control_bacteria_ttest,
             PVAL_control_supernatant = control_supernatant_ttest,
             PVAL_bacteria_supernatant = bacteria_supernatant_ttest)
    
    df_anova_list[[length(df_anova_list)+1]] <- df_anova
  }
}

# Bind data, insert NAs when pval is irrelevant
df_pval <- bind_rows(df_anova_list) %>% 
  mutate(PVAL_control_bacteria = ifelse(sample == "supernatant", NA, PVAL_control_bacteria),
         PVAL_control_supernatant = ifelse(sample == "bacterial culture", NA, PVAL_control_supernatant))

# Adjust p-values with fdr
df_adj <- df_pval %>% 
  select(-c(position, area)) %>% 
  distinct()

df_adj$FDR_control_bacteria <- p.adjust(df_adj$PVAL_control_bacteria, method = "fdr")
df_adj$FDR_control_supernatant <- p.adjust(df_adj$PVAL_control_supernatant, method = "fdr")
df_adj$FDR_bacteria_supernatant <- p.adjust(df_adj$PVAL_bacteria_supernatant, method = "fdr")

# Join data with controls
df_metabolites_full <- full_join(df_pval, df_adj) %>% 
  bind_rows(df_controls) %>% 
  write_csv(paste0(wd_path, "/output/2_metabolite_data_stats.csv"))

# Calculate mean
df_metabolites_norm <- df_metabolites_full %>%
  group_by(sample, condition, media, compound) %>% 
  mutate(mean_area = mean(area),
         sd_area = sd(area)) %>% 
  ungroup() %>% 
  group_by(media, compound) %>% 
  mutate(norm_area = area/max(mean_area)) %>% 
  select(1,2,3,4,6,8,5,7,14,15,16,9,10,11,12,13,14) %>% 
  write_csv(paste0(wd_path, "/output/2_metabolite_data_normalized.csv"))

# Clear environment
rm(list = ls())

### PLOT DATA ----

# Set paths to working directory and data folders
wd_path <- getwd()
setwd(wd_path)

df <- read_csv(paste0(wd_path, "/output/2_metabolite_data_normalized.csv"))

df$sample <- factor(df$sample, levels = c("control","bacterial culture","supernatant"))
df$condition <- factor(df$condition, levels = c("media","S. flexneri","F. nucleatum",
                                                "S. salivarius", "R. intestinalis",
                                                "C. comes", "G. elegans"))

species <- as.character(unique(df$condition))
species <- species[species != "media"]
compounds <- unique(df$compound)

#### PLOT NORM AREAS WITH PVALS PER SPECIES PER COMPOUND ----
for (s in species) {
  for (c in compounds) {
    bacteria <- gsub("\\. ", "", s)
    compound <- gsub(" ", "", c)
    
    pdf(paste0(wd_path, "/plots/2_metabolite_", bacteria, "_", compound, ".pdf"), width = 6, height = 5)
    
    # Find correct media for all species
    m <- df %>% 
      filter(condition%in%c(s)) %>% 
      pull(media) %>% 
      unique()
    
    # Check for and remove any NA or non-finite values in norm_area
    plot_data <- df %>%
      filter(compound == c & media == m) %>%
      filter(condition %in% c(s, "media")) %>%
      filter(!is.na(norm_area) & is.finite(norm_area))  # Remove rows with NA or non-finite values
    
    if (nrow(plot_data) > 0) {
      metabolite <- plot_data %>% pull(compound_analysed) %>% unique()
      # Build comparison data frame with positions and labels
      # Nothing passes 
      stat.test <- data.frame(
        group1 = c("control", "control", "bacterial culture"),
        group2 = c("bacterial culture", "supernatant", "supernatant"),
        p = c(
          plot_data %>% filter(!is.na(PVAL_control_bacteria) & is.finite(PVAL_control_bacteria)) %>% pull(PVAL_control_bacteria) %>% unique(),
          plot_data %>% filter(!is.na(PVAL_control_supernatant) & is.finite(PVAL_control_supernatant)) %>% pull(PVAL_control_supernatant) %>% unique(),
          plot_data %>% filter(!is.na(PVAL_bacteria_supernatant) & is.finite(PVAL_bacteria_supernatant)) %>% pull(PVAL_bacteria_supernatant) %>% unique()
        )
      ) %>%
        filter(!is.na(p)) %>%  # filter out missing p-values just in case
        mutate(
          p.signif = case_when(
            p <= 0.001 ~ "***",
            p <= 0.01 ~ "**",
            p <= 0.05 ~ "*",
            TRUE ~ "ns"
          ),
          y.position = max(plot_data$norm_area, na.rm = TRUE) * seq(1.05, 1.25, length.out = n())
        )
      
      p <- plot_data %>%
        group_by(condition, sample, compound) %>%
        ggplot(aes(x = sample, y = norm_area)) + 
        geom_boxplot() + 
        scale_y_continuous(limits = c(0, max(plot_data$norm_area, na.rm = TRUE) * 1.5)) +
        labs(x = "sample", y = paste(metabolite, "Relative Content")) +
        #scale_fill_viridis_d() + 
        theme_bw() + 
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14, vjust = 1, hjust = 0.5),
          axis.text.y = element_text(size = 14),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = "right",
          strip.text = element_text(size = 10)) + 
        ggtitle(paste0(s, " treated with ", c, " (", m, " media)")) + 
        stat_pvalue_manual(stat.test, label = "p.signif", 
                           tip.length = 0.01, 
                           bracket.size = 0.5, 
                           size = 5)
      
      print(p)
    }
    dev.off()
  }
}
