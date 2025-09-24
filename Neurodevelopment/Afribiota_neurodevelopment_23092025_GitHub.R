#################################
# AFRIBIOTA 16S fecal samples   #
# Jeanne TAMARELLE              #
# 19/10/2021                    #
# published 23/09/2025          #
#################################


# Overall structure
# I) Merging datasets and data management
  # 1) Sequence/count data
  # 2) Metadata
    # a) Correspondance file
    # b) Original metada
    # c) Neurodevelopment data
    # d) qPCR data
    # e) AAT, calpro
    # f) BCAA
  # 3) Formatting variables
  # 4) list of variables of interest

# II) Stunting and neurodevelopment
  # 1) Table sociodemo characteristics
  # 2) Regression model for stunting and development
  # 3) Univariate analysis for selection of variables for SEM
  # 4) Multivariate analysis for neurodevelopment

# III) Microbiome analyses
  # 1) Merge with microbiota and microbiota.asv datasets
  # 2) Description of reads and taxa
  # 3) alpha-diversity
    # a) univariate linear regression and plots
    # b) linear regression controlling for AgeASQ3 and haz
    # c) linear regression controlling for sequencing run
  # 4) Transformation of the data
    # a) at the lowest available rank
    # b) at the ASV level
    # c) at the genus level
    # d) at the family level
  # 5) Barplots
    # a) at the lowest available rank
    # b) at the genus level
    # c) at the family level
  # 6) Correlation plots
    # a) at the lowest available rank
    # b) at the genus level
    # c) at the family level
  # 7) DESEQ2
    # a) at the lowest available rank
    # b) at the genus level
    # c) at the family level
  # 8) ANCOM-BC
    # a) at the genus level
    # b) at the family level
  # 9) beta-diversity
  # 10) cluster analysis
    # a) hierarchical clustering, global heatmap
    # b) statistical testing

# IV) SEM
  # 1) Correlation between indicator variables for each latent variable
  # 2) CFA models
  # 3) SEM models
    # a) the simplest model
    # b) the simple model
    # c) the complex SEM
    # d) path analysis


rm(list=ls())

setwd("~/Documents/UNIL/Afribiota")
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(phyloseq)


# I) merging datasets and data management---------

    # 1) Sequence/count data------------

    count = read.csv(file = "Data/sequencetableAfribiotaMadafeces.csv", header=TRUE, sep=";")
    taxon = read.csv(file = "Data/taxtableAfribiotaMadafeces.csv", header=TRUE, sep=";") # What are the variables Accession & X.1?

    # Replace all "-" with "_" everywhere in the data frame
    taxon[] <- lapply(taxon, function(col) {
      if (is.character(col)) {
        gsub("-", "_", col)
      } else {
        col  # leave numeric/factor columns unchanged
      }
    })
    
    colnames(count)[which(colnames(count)=="X")] = "ASV"
    colnames(taxon)[which(colnames(taxon)=="X")] = "ASV"

    # Replace Rank7
    taxon$Rank7[which(taxon$Rank7 %in% c("Camelus_dromedarius_(Arabian_camel)","Echinogammarus_veneris"))]

    taxon$Accession[which(taxon$Rank7 == "Camelus_dromedarius_(Arabian_camel)")]
    taxon$Rank7[which(taxon$Rank7 == "Camelus_dromedarius_(Arabian_camel)")] <- "Escherichia_coli"

    taxon$Rank7[which(taxon$Rank7 == "Echinogammarus_veneris")]
    taxon$Rank7[which(taxon$Rank7 == "Echinogammarus_veneris")] <- NA

    # Replace unrelevant assignations
    taxon$Rank7[which(taxon$Rank7 %in% c("metagenome","mouse_gut_metagenome","wastewater_metagenome","groundwater_metagenome","soil_metagenome","gut_metagenome","human_gut_metagenome", "unidentified", "uncultured_bacterium","wallaby_gut_metagenome"))] = NA
    taxon$Rank6[which(taxon$Rank6 %in% c("uncultured"))] = NA
    taxon$Rank7[which(taxon$Rank7 %in% c("Trichuris_trichiura_(human_whipworm)"))] = NA
    
    # Create phyloseq object
    taxon2 = taxon
    rownames(taxon2) = taxon2$ASV
    taxon2 = as.matrix(taxon2)
    taxon2 <- tax_table(taxon2)
    
    count2 = count
    rownames(count2) = count2$ASV
    count2 = count2[,-which(colnames(count2) == "ASV")]
    count2 = otu_table(count2, taxa_are_rows = TRUE)
    
    ps = phyloseq(taxon2, count2)
    
    # Removing taxa with abundance 0
    get_taxa_unique(ps, "Rank1")  ## Bacteria, Archaea, unassigned
    # filter out Unassigned, Chloroplast, Mitochondria and Eukaryota and others
    ps_filtered <- ps %>%
      subset_taxa(Rank5 != "Mitochondria") %>%
      subset_taxa(Rank4 != "Chloroplast") %>%
      subset_taxa(Rank1 !="Unassigned") %>%
      subset_taxa(Rank1 != "Eukaryota")
    ps_filtered # new number of taxa: from 6178 to 5479
    ps_filtered <- prune_taxa(taxa_sums(ps_filtered)> 0, ps_filtered)
    ps_filtered # from 5479 to 1813 taxa only !
    
    # agglomerate phyloseq object at different levels
    # Relative abundance
    # ps_rel <- transform_sample_counts(ps_filtered, function(x) x / sum(x))
    ps_phy = tax_glom(ps_filtered, taxrank = "Rank2")
    # ps_phy_rel <- transform_sample_counts(ps_phy, function(x) x / sum(x))
      m_phy = as.data.frame(otu_table(ps_phy))
      m_phy$ASV = rownames(m_phy)
      m_phy2 = as.data.frame(tax_table(ps_phy))
      m_phy = inner_join(m_phy2[,which(colnames(m_phy2) %in% c("Rank2","ASV") )],m_phy, by = "ASV")
      m_phy = m_phy[,-which(colnames(m_phy) == "ASV")]
      m_phy %<>% 
        group_by(Rank2) %>% 
        summarise(across(everything(), sum))
      m_phy = tibble::column_to_rownames(m_phy, var = "Rank2")
      m_phy = as.data.frame(t(m_phy))
    ps_fam = tax_glom(ps_filtered, taxrank = "Rank5")
    # ps_fam_rel <- transform_sample_counts(ps_fam, function(x) x / sum(x))
      m_fam = as.data.frame(otu_table(ps_fam))
      m_fam$ASV = rownames(m_fam)
      m_fam2 = as.data.frame(tax_table(ps_fam))
      m_fam = inner_join(m_fam2[,which(colnames(m_fam2) %in% c("Rank5","ASV") )],m_fam, by = "ASV")
      m_fam = m_fam[,-which(colnames(m_fam) == "ASV")]
      m_fam %<>% 
        group_by(Rank5) %>% 
        summarise(across(everything(), sum))
      m_fam = tibble::column_to_rownames(m_fam, var = "Rank5")
      m_fam = as.data.frame(t(m_fam))
    ps_gen = tax_glom(ps_filtered, taxrank = "Rank6")
    # ps_gen_rel <- transform_sample_counts(ps_gen, function(x) x / sum(x))
      m_gen = as.data.frame(otu_table(ps_gen))
      m_gen$ASV = rownames(m_gen)
      m_gen2 = as.data.frame(tax_table(ps_gen))
      m_gen = inner_join(m_gen2[,which(colnames(m_gen2) %in% c("Rank6","ASV") )],m_gen, by = "ASV")
      m_gen = m_gen[,-which(colnames(m_gen) == "ASV")]
      m_gen %<>% 
        group_by(Rank6) %>% 
        summarise(across(everything(), sum))
      m_gen = tibble::column_to_rownames(m_gen, var = "Rank6")
      m_gen = as.data.frame(t(m_gen))
    
    rm(count2,taxon2,m_phy2,m_fam2,m_gen2)
    
    # Order by decreasing abundance
    m_phy = m_phy[,order(colSums(m_phy),decreasing = T)]
    m_fam = m_fam[,order(colSums(m_fam),decreasing = T)]
    m_gen = m_gen[,order(colSums(m_gen),decreasing = T)]
    
    # Create ID variable and move it at the beginning of table
    m_phy %<>%
      mutate(ID_metag = rownames(m_phy)) %>%
      dplyr::select("ID_metag", everything())
    rownames(m_phy) = NULL
    m_fam %<>%
      mutate(ID_metag = rownames(m_fam)) %>%
      dplyr::select("ID_metag", everything())
    rownames(m_fam) = NULL
    m_gen %<>%
      mutate(ID_metag = rownames(m_gen)) %>%
      dplyr::select("ID_metag", everything())
    rownames(m_gen) = NULL
    
    rm(ps, ps_phy, ps_fam, ps_gen)
    
    # Replace the highest level of taxonomy by the previous one
    taxon$Rank2[is.na(taxon$Rank2)] = paste("d_",taxon$Rank1[is.na(taxon$Rank2)], sep="")
    taxon$Rank3[is.na(taxon$Rank3)] = paste("p_",taxon$Rank2[is.na(taxon$Rank3)], sep="")
    taxon$Rank4[is.na(taxon$Rank4)] = paste("c_",taxon$Rank3[is.na(taxon$Rank4)], sep="")
    taxon$Rank5[is.na(taxon$Rank5)] = paste("o_",taxon$Rank4[is.na(taxon$Rank5)], sep="")
    taxon$Rank6[is.na(taxon$Rank6)] = paste("f_",taxon$Rank5[is.na(taxon$Rank6)], sep="")
    taxon$Rank7[is.na(taxon$Rank7)] = paste("g_",taxon$Rank6[is.na(taxon$Rank7)], sep="")

    # Remove the g_f_o_c_p prefixes for species level when unnecessary
    taxon$Rank7 = sub('g_f_o_c_p_d_', 'd_', taxon$Rank7)
    taxon$Rank7 = sub('g_f_o_c_p_', 'p_', taxon$Rank7)
    taxon$Rank7 = sub('g_f_o_c_', 'c_', taxon$Rank7)
    taxon$Rank7 = sub('g_f_o_', 'o_', taxon$Rank7)
    taxon$Rank7 = sub('g_f_', 'f_', taxon$Rank7)

    # Remove the f_o_c_p prefixes for genus level when unnecessary
    taxon$Rank6 = sub('f_o_c_p_d_', 'd_', taxon$Rank6)
    taxon$Rank6 = sub('f_o_c_p_', 'p_', taxon$Rank6)
    taxon$Rank6 = sub('f_o_c_', 'c_', taxon$Rank6)
    taxon$Rank6 = sub('f_o_', 'o_', taxon$Rank6)

    # Remove the f_o_c_p prefixes for family level when unnecessary
    taxon$Rank5 = sub('o_c_p_d_', 'd_', taxon$Rank5)
    taxon$Rank5 = sub('o_c_p_', 'p_', taxon$Rank5)
    taxon$Rank5 = sub('o_c_', 'c_', taxon$Rank5)
    taxon$Rank5 = sub('o_', 'o_', taxon$Rank5)

    # Remove special characters from Rank7's names (ex: [Bacteroides]_pectinophilus_ATCC_43243)
    taxon$Rank7 = sub("\\[", "", taxon$Rank7)
    taxon$Rank7 = sub("\\]", "", taxon$Rank7)
    taxon$Rank7 = sub("\\(", "", taxon$Rank7)
    taxon$Rank7 = sub("\\)", "", taxon$Rank7)

    # Remove special characters from Rank6's names
    taxon$Rank6 = sub("\\[", "", taxon$Rank6)
    taxon$Rank6 = sub("\\]", "", taxon$Rank6)
    taxon$Rank6 = sub("\\(", "", taxon$Rank6)
    taxon$Rank6 = sub("\\)", "", taxon$Rank6)

    # Remove special characters from Rank5's names
    taxon$Rank5 = sub("\\[", "", taxon$Rank5)
    taxon$Rank5 = sub("\\]", "", taxon$Rank5)
    taxon$Rank5 = sub("\\(", "", taxon$Rank5)
    taxon$Rank5 = sub("\\)", "", taxon$Rank5)

    # Remove dots. from names of Rank7
    grep("\\.", taxon$Rank7, value=TRUE) 
    taxon$Rank7 = gsub('\\.', '', taxon$Rank7)
    grep("=", taxon$Rank7, value=TRUE) 
    taxon$Rank7 = gsub('=', '', taxon$Rank7)
    grep(":", taxon$Rank7, value=TRUE) 
    taxon$Rank7 = gsub(':', '', taxon$Rank7)
    grep("-", taxon$Rank7, value=TRUE)
    taxon$Rank7 = gsub('-', '_', taxon$Rank7)
    
    # Check if all ASV of one table are included in the other table
    union(setdiff(taxon$ASV,count$ASV),setdiff(count$ASV,taxon$ASV))
    
    # For the ASV level, merge ASV + lowest rank for which information is available
    taxon$Rank8 = paste0(taxon$Rank7,"_",taxon$ASV)
    
    # remove ASVs that have a rowSums = 0
    count = count[which(rowSums(count[,-which(colnames(count) == "ASV")] ) != 0),] # we go down to 526 taxa instead of 1858
    
    # merge the two datasets
    microbiota = inner_join(taxon,count, by="ASV")

    # Keep a dataset at the ASV level
    microbiota.asv = microbiota[,-which(colnames(microbiota) %in% c("ASV","Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7","Accession"))]
    rownames(microbiota.asv) = microbiota.asv$Rank8
    microbiota.asv = microbiota.asv[,-which(colnames(microbiota.asv) == "Rank8")]
    microbiota.asv = as.data.frame(t(microbiota.asv))
    # Order microbiota.asv by abundance
    microbiota.asv = microbiota.asv[,order(colSums(microbiota.asv),decreasing = T)]
    # Create ID variable and move it at the beginning of table
    microbiota.asv %<>%
      mutate(ID_metag = rownames(microbiota.asv)) %>%
      dplyr::select("ID_metag",everything())
    rownames(microbiota.asv) = NULL
    
    # for microbiota at lowest available rank, sum the ASVs that correspond to the same species (taxa variable)
    microbiota2 = aggregate(x = microbiota[,11:ncol(microbiota)], by = list(microbiota$Rank7), FUN = sum)
    rownames(microbiota2) <- microbiota2$Group.1
    microbiota2 = microbiota2[,-which(colnames(microbiota2) == "Group.1")]
    microbiota = as.data.frame(t(microbiota2))
    # Order by decreasing abundance
    microbiota = microbiota[,order(colSums(microbiota),decreasing = T)]
    # Create ID variable and move it at the beginning of table
    microbiota %<>%
      mutate(ID_metag = rownames(microbiota)) %>%
      dplyr::select("ID_metag", everything())
    rownames(microbiota) = NULL

    rm(microbiota2)
    
    # 2) Metadata --------

      # a) Correspondance file------
      meta = read.csv(file = "Data/metadataAfribiotaAccessMadafeces.csv", header=TRUE, sep=";")
  
      # keep only certain variable
      meta = meta[,c("ID_afri","ID_metag","run")]
    
      # b) Original metadata----------
      meta2 = read.csv(file = "Data/AFRIBIOTA_parasites_data_jt.csv", header=TRUE, sep=",")
     
      # merge with correspondance file
      meta = left_join(meta2,meta, by=c("ID_afri")) # intersection = 445 patients
  
      rm(meta2)
      
      # c) Neurodevelopment data---------
      meta3 = read.csv(file = "Data/FINAL_DATA_2021/Afribiota_Develop+Clinic_FINAL_N372_UTF8.csv", header=TRUE, sep=",")
      
      # keep only certain variable
      meta3 = meta3[,c("ID_afri","Consentement","Date","Age_Dev","AgeASQ3","AGEASQ3_classes","Age_Clinic","sexe","haz_contin",
                       "C1","C2","C3","C4","C5","C6","Comm_tot","Comm_delay",
                       "GM1","GM2","GM3","GM4","GM5","GM6","GM_tot","GM_delay",
                       "FM1","FM2","FM3","FM4","FM5","FM6","FM_tot","FM_delay",
                       "PS1","PS2","PS3","PS4","PS5","PS6","PS_tot","PS_delay",
                       "PES1","PES2","PES3","PES4","PES5","PES6","PES_tot","PES_delay","Nº_domains_delays","Overal_score" )]
      
      mydata = inner_join(meta,meta3, by=c("ID_afri")) # intersection = 372 patients
  
      # check if sequencing data for all
      table(is.na(mydata$ID_metag)) # only 338 patients for which we have sequencing data

    
      # d) qPCR data---------
      qpcr = read.csv(file = "Data/MergeqPCRAfribiota_Munirv2.csv", header=TRUE, sep=";")
      
      # merge with mydata
      mydata = left_join(mydata, qpcr[,-which(colnames(qpcr) == "id")], by = "ID_afri")
      
      # e) AAT, calpro---------
      aatcalpro = read.csv(file = "Data/aatcalpro/AATdata.csv", header=TRUE, sep=",")
      
      # merge with mydata
      mydata = left_join(mydata, aatcalpro[,-which(colnames(aatcalpro) %in% c("id.x","id.y"))], by = "ID_afri")
      
      # f) BCAA--------
      bcaa = haven::read_dta("Data/bcaa/2018_07_13_Afribiota_Results_Envoi1et2_merge.dta") # 765 observations
      
      # merge with mydata
      mydata = left_join(mydata, bcaa[,-which(colnames(bcaa) == "id")], by = "ID_afri")
    
   
    rm(aatcalpro, bcaa, meta, meta3, qpcr)
     

    # 3) Formatting variables-------
    str(mydata)
    
    table(mydata$haz_cont,mydata$haz)
    mydata$haz = factor(mydata$haz, levels = unique(mydata$haz))
    mydata$haz = relevel(mydata$haz, ref = "1")
    mydata$haz = relevel(mydata$haz, ref = "0")
    
    
    mydata$Comm_tot = as.numeric(mydata$Comm_tot)
    mydata$Comm_tot[is.na(mydata$Comm_tot)==T] = sum(mydata[is.na(mydata$Comm_tot)==T,c("C1","C2","C3","C4","C5","C6")])
    mydata$Comm_delay = as.factor(mydata$Comm_delay)

    mydata$PS_tot = as.numeric(mydata$PS_tot)
    mydata$PS_tot[is.na(mydata$PS_tot)==T] = sum(mydata[is.na(mydata$PS_tot)==T,c("PS1","PS2","PS3","PS4","PS5","PS6")])
    mydata$PS_delay = as.factor(mydata$PS_delay)

    mydata$PES_tot = as.numeric(mydata$PES_tot)
    mydata$PES_delay = as.factor(mydata$PES_delay)
    
    mydata$FM_tot = as.numeric(mydata$FM_tot)
    mydata$FM_tot[is.na(mydata$FM_tot)==T] = sum(mydata[is.na(mydata$FM_tot)==T,c("FM1","FM2","FM3","FM4","FM5","FM6")])
    mydata$FM_delay = as.factor(mydata$FM_delay)
    
    mydata$GM_tot = as.numeric(mydata$GM_tot)
    mydata$GM_tot[is.na(mydata$GM_tot)==T] = sum(mydata[is.na(mydata$GM_tot)==T,c("GM1","GM2","GM3","GM4","GM5","GM6")])
    mydata$GM_delay = as.factor(mydata$GM_delay)
    
    mydata$Overal_score = as.numeric(mydata$Overal_score)
    mydata$Overal_score[is.na(mydata$Overal_score)] = rowSums(mydata[which(is.na(mydata$Overal_score)),c("Comm_tot","PS_tot","PES_tot","FM_tot","GM_tot")])
    # create overall delay variable
    mydata$Overal_delay = 0
    mydata$Overal_delay[which(mydata$Nº_domains_delays %in% c(1, 2, 3))]   = 1
    mydata$Overal_delay = as.factor(mydata$Overal_delay)
    table(mydata$Overal_delay, mydata$Nº_domains_delays, useNA = "always")
    
   
    # 4) List of variables of interest ---------

colnames(mydata)

household = c("score_socio_eco_bycountry","nb_pers","salaires_cum","eau","eau_traitee","eau_pot_uniq","nb_pieces")
  mydata$score_socio_eco_bycountry = factor(mydata$score_socio_eco_bycountry, levels = unique(mydata$score_socio_eco_bycountry))
  mydata$score_socio_eco_bycountry <- relevel(mydata$score_socio_eco_bycountry, ref = "1")
  mydata$eau = as.factor(mydata$eau)
  mydata$eau_traitee = as.factor(mydata$eau_traitee)
  mydata$eau_pot_uniq = as.factor(mydata$eau_pot_uniq)
  mydata$nb_pers = as.numeric(mydata$nb_pers)
  mydata$nb_pieces = as.numeric(mydata$nb_pieces)
labelled::var_label(mydata[,which(names(mydata) %in% household)]) = c("Number of people in household","Cumulated household weekly income","Running water","Number of rooms","Treated drinking water","Drinkable water exclusive consumption","Family economic status")
labelled::look_for(mydata[,which(names(mydata) %in% household)])
questionr::freq.na(mydata[,which(names(mydata) %in% household)]) 
    
family = c("age_prem_gross","age_prem_gross_15","statut_mere","mere_emploi","instruc_mere","instruc_mere_3","mere_sousnutris")
  mydata$age_prem_gross_15 = as.factor(mydata$age_prem_gross_15)
  mydata$age_prem_gross = as.numeric(mydata$age_prem_gross)
  mydata$statut_mere = factor(mydata$statut_mere, levels = unique(mydata$statut_mere))
  mydata$statut_mere <- relevel(mydata$statut_mere, ref = "2")
  mydata$statut_mere <- relevel(mydata$statut_mere, ref = "1")
  mydata$mere_emploi[which(mydata$mere_emploi == 9)] <- NA
  mydata$mere_emploi = as.factor(mydata$mere_emploi)
  mydata$instruc_mere = factor(mydata$instruc_mere, levels = unique(mydata$instruc_mere))
  mydata$instruc_mere = relevel(mydata$instruc_mere, ref = "0")
  mydata$instruc_mere_4 = factor(mydata$instruc_mere_4, levels = unique(mydata$instruc_mere_4))
  mydata$instruc_mere_4 = relevel(mydata$instruc_mere_4, ref = "1")
  mydata$instruc_mere_3 = factor(mydata$instruc_mere_4, levels=c("1","2","3","4"), labels = c("1","1","2","3"))
  mydata$mere_sousnutris = as.factor(mydata$mere_sousnutris)
labelled::var_label(mydata[,which(names(mydata) %in% family)]) = c("Maternal age at first pregnancy","Marital status","Maternal education level","Working mother","Undernourished mother","Maternal age at first pregnancy >15","Maternal education level cat")
labelled::look_for(mydata[,which(names(mydata) %in% family)])
questionr::freq.na(mydata[,which(names(mydata) %in% family)]) 

# age calculated from dt_naiss et dt_quest
child = c("sexe","age","haz_cont","haz","waz_cont","taille","poids","statut_nut","mere","atcd_malnut","poids_naiss","taille_naiss_p2","anemie2","ferritine2","ferritine2_seuil","allaite","age_allaite","allaite_exclu_6mois","voie_accou","hemoglobine2","hemoglobine2_seuil","sibo", "age_alim")
  mydata$sexe = as.factor(mydata$sexe)
  mydata$statut_nut = factor(mydata$statut_nut, levels = unique(mydata$statut_nut))
  mydata$statut_nut = relevel(mydata$statut_nut, ref = "2")
  mydata$statut_nut = relevel(mydata$statut_nut, ref = "1")
  mydata$mere = as.factor(mydata$mere)
  mydata$atcd_malnut = factor(mydata$atcd_malnut, levels = unique(mydata$atcd_malnut))
  mydata$anemie2 = as.factor(mydata$anemie2)
  mydata$ferritine2_seuil = as.factor(mydata$ferritine2_seuil)
  mydata$allaite = as.factor(mydata$allaite)
  mydata$allaite_exclu_6mois = as.factor(mydata$allaite_exclu_6mois)
  mydata$voie_accou = as.factor(mydata$voie_accou)
  mydata$hemoglobine2_seuil = as.factor(mydata$hemoglobine2_seuil)
  mydata$sibo = factor(mydata$sibo, levels = unique(mydata$sibo))
  mydata$sibo = relevel(mydata$sibo, ref = "0")
  mydata$taille_naiss_p2[which(mydata$taille_naiss_p2 == 9)] = NA
  anova(lm(mydata$poids_naiss ~ mydata$taille_naiss_p2))
  pdf(file = "Boxplot_birth_weight_reported.pdf")
  boxplot(mydata$poids_naiss ~ mydata$taille_naiss_p2)
  dev.off()
  mydata$taille_naiss_p2 = as.factor(mydata$taille_naiss_p2)
labelled::var_label(mydata[,which(names(mydata) %in% child)]) = c("Sex","Height","Weight","Nutritional status","Living with mother","Malnutrition history","Birth mode","Birth weight","Currently breastfed","Age of solid food introduction","Age of breastfeeding cessation","SIBO","Age","Weight-for-Age z-score","Height-for-Age z-score", "Height-for-Age z-score cat","Exclusive breastfeeding before 6 months","Haemoglobin","Anemia (adjusted)","Low haemoglobin","Ferritine","Low ferritine","Birth size (reported)")
labelled::look_for(mydata[,which(names(mydata) %in% child)])
questionr::freq.na(mydata[,which(names(mydata) %in% child)]) # problem with Birth weight (52% missing values)

development = c("Overal_score","Comm_tot","PS_tot","PES_tot","FM_tot","GM_tot","Nº_domains_delays")
labelled::var_label(mydata[,which(names(mydata) %in% development)]) = c("Communication score","Gross Motor score","Fine Motor score","Cognition-Problem Solving score","Personal-Social score","Number of domains with delay","Overall development score")
labelled::look_for(mydata[,which(names(mydata) %in% development)])
questionr::freq.na(mydata[,which(names(mydata) %in% development)]) 

# parasite = analyse parasitaire réalisée (only to make sure that you only include these ones if you include the parasites in the model)
parasite = c("parasite","ascaris","trichuris","giardiase","blastocystis","enterobius","microsporidies","hymenolepis")
labelled::var_label(mydata[,which(names(mydata) %in% parasite)]) = c("Parasite analysis done","Ascaris","Trichuris","Giardiasis","Blastocystis","Hymenolepis","Enterobius","Microsporidia")
labelled::look_for(mydata[,which(names(mydata) %in% parasite)])
questionr::freq.na(mydata[,which(names(mydata) %in% parasite)]) 

# qPCR
# Salmonella spp. and Shigella spp. were identified by amplification of the ompC and the invasion plasmid antigen H (ipaH) gene 
# For enterotoxigenic Escherichia coli (ETEC), heat-labile toxin (eltB) and heat-stable toxin (estA) coding regions were targeted. 
# For enteropathogenic E. coli (EPEC), the bundle-forming pilus (encoded by the bfpA gene) carried by the EPEC adherence factor (EAF) plasmid and the intimin (eae gene for EPEC attaching and effacing), an outer membrane adhesion essential for the intimate attachment of the EPEC or enterohemorrhagic E. coli (EHEC) to enterocytes were the targets. 
# For enteroaggregative E. coli (EAEC), aggR and aaiC genes were amplified. The aggR gene encodes a transcriptional activator of the aggregative adherence fimbriae expression in enteroaggregative E. coli and aaiC, is part of the aai gene cluster, encoding a type VI secretion system. When one of the two targets was detected for ETEC (estIa or eltB) or EAEC (aggR or aaiC), the pathogen was considered to be present. 
# Since eae can be present both in EPEC and EHEC, we considered only the presence of bfpA gene for EHEC detection.
# Fibronectin-binding protein (cadF) gene and the cholera toxin (CT) subunit A gene (ctxA) were the targeted genes for Campylobacter jejuni/coli and Vibrio cholerae, respectively
qpcr = c("ompC","ipaH","estla","eltB","AggR","aaiC","bfpA","eae","cadF","cta","No_genes")
mydata$No_genes = as.numeric(mydata$No_genes)
labelled::var_label(mydata[,which(names(mydata) %in% qpcr)]) = c("ompC","ipaH","estla","eltB","AggR","aaiC","bfpA","eae","cadF","cta","Number of enteropathogen genes")
labelled::look_for(mydata[,which(names(mydata) %in% qpcr)])
questionr::freq.na(mydata[,which(names(mydata) %in% qpcr)]) # 5% missing

# diet
diet = c("dds","ldd","asf")
  mydata$dds = as.numeric(mydata$dds)
  mydata$ldd = as.factor(mydata$ldd)
  mydata$asf = as.factor(mydata$asf)
labelled::var_label(mydata[,which(names(mydata) %in% diet)]) = c("Dietary diversity score","Low dietary diversity","Animal source foods")
labelled::look_for(mydata[,which(names(mydata) %in% diet)])
questionr::freq.na(mydata[,which(names(mydata) %in% diet)]) # 5% missing

# BCAA
bcaa = c("AlanineM","CitrulineM","ValineM","IsoleucineM","LeucineM","CitrulineM_cat")
labelled::var_label(mydata[,which(names(mydata) %in% bcaa)]) = c("Alanine","Citruline","Valine","Isoleucine","Leucine","Citruline cat")
labelled::look_for(mydata[,which(names(mydata) %in% bcaa)])
questionr::freq.na(mydata[,which(names(mydata) %in% bcaa)]) # 14% missing

# AAT, calpro
# create categories for aat and calpro
mydata <- mydata %>%
  mutate(AATlevel = case_when(
    AATmggPS > 2  ~ "elevated",
    AATmggPS < 1.25 | AATmggWET < 0.15   ~ "normal",
    AATmggPS < 2   ~ "grey zone", 
    is.na(AATmggPS) ~ "missing", )) 

mydata$AATwet_cat = cut(mydata$AATmggWET, breaks = c(0,0.149,100), labels = c("normal","elevated")) # values below 0.15 mg/g of fecal wet weight are considered normal
mydata$AATdry_cat = cut(mydata$AATmggPS, breaks = c(0,1.24,100), labels = c("normal","elevated")) # values below 1.25 mg/g of fecal dry weight are considered normal
table(mydata$AATmggWET,mydata$AATwet_cat,useNA="always") # 19 missing values
table(mydata$AATmggPS,mydata$AATdry_cat,useNA="always") # 34 missing values (including the 19 missing values according to wet weight)
table(mydata$AATdry_cat,mydata$AATwet_cat,useNA="always")
mydata$AATlevel2 = mydata$AATwet_cat
mydata$AATlevel2[mydata$AATdry_cat == "elevated"] = "elevated" # according to Pascale email: threshold AAT: We defined "elevated" if at least one of the two values was elevated. 
table(mydata$AATlevel2,useNA="always") # 19 missing values
# coherence between Munir's function and Pascale's definition
table(mydata$AATlevel, mydata$AATlevel2, useNA="always") # in pascale's definition, no grey zone and less missing values. Otherwise concordant --> use AATlevel2

mydata$CALwet_cat = ifelse( (mydata$age<36 & mydata$CALPROTECTINE..μg.g. < 150 ), "normal", ifelse((mydata$age>=36 & mydata$CALPROTECTINE..μg.g. < 100 ), "normal", "elevated") ) # values under 150 mg/g of fecal wet weight for children below 3 years old, and 100 mg/g for children of 3 years old and above are considered normal
table(mydata$CALwet_cat,useNA="always") # 12 missing values for calprotectine

mydata$crp2_cat = ifelse(mydata$crp2 < 10, "normal","elevated")

inflammation = c("crp2_cat","AATlevel2","AATwet_cat","AATdry_cat","CALwet_cat")
  mydata$CALwet_cat = factor(mydata$CALwet_cat, levels = unique(mydata$CALwet_cat))
  mydata$CALwet_cat = relevel(mydata$CALwet_cat, ref = "normal")
  mydata$crp2_cat = factor(mydata$crp2_cat, levels=unique(mydata$crp2_cat))
  mydata$crp2_cat = relevel(mydata$crp2_cat, ref = "normal")
labelled::var_label(mydata[,which(names(mydata) %in% inflammation)]) = c("AATwet_cat","AATdry_cat","AAT level","Calprotectine level cat (wet)","CRP level cat")
labelled::look_for(mydata[,which(names(mydata) %in% inflammation)])
questionr::freq.na(mydata[,which(names(mydata) %in% inflammation)]) 







# II) Stunting and neurodevelopment--------

    # 1) Table sociodemo characteristics ----------

    # suppl table 1
    family_simple = family[!family %in% c("instruc_mere")]
    child_simple = child[!child %in% c("statut_nut","atcd_malnut","poids_naiss","anemie2","sibo")]
    inflammation_simple = inflammation[!inflammation %in% c("AATwet_cat","AATdry_cat")]
    parasite_simple = parasite[!parasite %in% c("parasite")]
    
    mydata %>%
      gtsummary::tbl_summary(
        include = c(all_of(household),all_of(family_simple),all_of(child_simple),all_of(development),all_of(diet),all_of(parasite_simple),all_of(qpcr),all_of(inflammation_simple),all_of(bcaa) ),
        by = "haz", # my outcome variable in 3 modalities, which will be in columns 
        missing = "no",
        type = list(dds ~ "continuous", No_genes ~ "continuous", age_prem_gross ~ "continuous", nb_pers ~ "continuous", nb_pieces ~ "continuous"),
        statistic = list(all_continuous() ~ "{mean} ({min}-{max})"),
        digits = list(all_continuous() ~ 1, all_categorical() ~ c(0,1))) %>% 
      add_overall(last = TRUE) %>%
      add_p() %>% # makes the appropriate tests and indicates the corresponding p-value
      bold_labels()
    

    # 2) Regression model stunting and development---------
    # According to the literature: adjust for household wealth, rurality, size, access to improved sanitation and access to improved water source, maternal age, education and marital status, child age, sex, anemia and attendance of an early childhood education programme
    str(mydata[,which(colnames(mydata) %in% c("score_socio_eco_bycountry", "eau", "sexe", "age", "age_prem_gross", "instruc_mere_3", "anemie2"))])
    str(mydata[,which(colnames(mydata) %in% c(development))])
    dvt <- match(x=c(development),names(mydata))

      # Load the necessary library for the heteroskedasticity test
      library(lmtest)
    
      out_variable=rep(NA, length(dvt))
      out_MD1=rep(NA, length(dvt))
      out_lci1 = rep(NA, length(dvt))
      out_hci1 = rep(NA, length(dvt))
      out_pvalue1=rep(NA, length(dvt))
      out_MD2=rep(NA, length(dvt))
      out_lci2 = rep(NA, length(dvt))
      out_hci2 = rep(NA, length(dvt))
      out_pvalue2=rep(NA, length(dvt))
      out_MDa1=rep(NA, length(dvt))
      out_lcia1 = rep(NA, length(dvt))
      out_hcia1 = rep(NA, length(dvt))
      out_pvaluea1=rep(NA, length(dvt))
      out_MDa2=rep(NA, length(dvt))
      out_lcia2 = rep(NA, length(dvt))
      out_hcia2 = rep(NA, length(dvt))
      out_pvaluea2=rep(NA, length(dvt))
      hetero_pval_uni <- rep(NA, length(dvt))
      hetero_pval_multi <- rep(NA, length(dvt))
      number = 1
      
      # create pdf file to store lm diagnostics
      pdf(file = "lm_stunting_score_diagnostics.pdf")

      for (j in dvt) {
        print(j)
        
        outcome = colnames(mydata)[j]
        
        # Univariate model
        model <- lm(get(outcome) ~ haz, mydata)
        
        # Heteroskedasticity test for univariate model (Breusch-Pagan)
        bp_uni <- bptest(model)
        hetero_pval_uni[number] <- bp_uni$p.value  # Store p-value
        
        # extract results for "undernourished"
        MD1 <- summary(model)$coefficients[2,1]      # coefficient for count_ct
        lci1 <- confint(model)[2,1]  # lower bound for 95% CI, for count_ct
        hci1 <- confint(model)[2,2]  # higher bound for 95% CI, for count_ct
        pvalue1 <- summary(model)$coefficients[2,4]
        
        # extract results for severely undernourished
        MD2 <- summary(model)$coefficients[3,1]      # coefficient for count_ct
        lci2 <- confint(model)[3,1]  # lower bound for 95% CI, for count_ct
        hci2 <- confint(model)[3,2]  # higher bound for 95% CI, for count_ct
        pvalue2 <- summary(model)$coefficients[3,4]
        
        # store univariate results
        out_MD1[number] = as.numeric(MD1)
        out_lci1[number] = as.numeric(lci1)
        out_hci1[number] = as.numeric(hci1)
        out_pvalue1[number] = as.numeric(pvalue1)
        
        out_MD2[number] = as.numeric(MD2)
        out_lci2[number] = as.numeric(lci2)
        out_hci2[number] = as.numeric(hci2)
        out_pvalue2[number] = as.numeric(pvalue2)
        
        # Multivariate model
        modela <- lm(get(outcome) ~ haz + score_socio_eco_bycountry + eau_traitee + sexe + age +  age_prem_gross +  instruc_mere_3 , mydata)
        
        # Heteroskedasticity test for multivariate model (Breusch-Pagan)
        bp_multi <- bptest(modela)
        hetero_pval_multi[number] <- bp_multi$p.value  # Store p-value
        
        # extract results for "undernourished"
        MDa1 <- summary(modela)$coefficients[2,1]      # coefficient for count_ct
        lcia1 <- confint(modela)[2,1]  # lower bound for 95% CI, for count_ct
        hcia1 <- confint(modela)[2,2]  # higher bound for 95% CI, for count_ct
        pvaluea1 <- summary(modela)$coefficients[2,4]
        
        # extract results for "severely undernourished"
        MDa2 <- summary(modela)$coefficients[3,1]      # coefficient for count_ct
        lcia2 <- confint(modela)[3,1]  # lower bound for 95% CI, for count_ct
        hcia2 <- confint(modela)[3,2]  # higher bound for 95% CI, for count_ct
        pvaluea2 <- summary(modela)$coefficients[3,4]
        
        # store multivariate results
        out_MDa1[number] = as.numeric(MDa1)
        out_lcia1[number] = as.numeric(lcia1)
        out_hcia1[number] = as.numeric(hcia1)
        out_pvaluea1[number] = as.numeric(pvaluea1)
        
        out_MDa2[number] = as.numeric(MDa2)
        out_lcia2[number] = as.numeric(lcia2)
        out_hcia2[number] = as.numeric(hcia2)
        out_pvaluea2[number] = as.numeric(pvaluea2)
        
        out_variable[number] = outcome
        number = number + 1
        
        # Plot des diagnostics pour le modèle univarié
        par(mfrow = c(2, 2))  # 4 graphiques par page
        plot(model)  # Diagnostic plots for univariate model
        
        # Plot des diagnostics pour le modèle multivarié
        par(mfrow = c(2, 2))  # 4 graphiques par page
        plot(modela)  # Diagnostic plots for multivariate model
        
    }
     
      # close pdf
      dev.off()
      
      # Create final dataframe
      out = data.frame(out_variable, out_MD1, out_lci1, out_hci1, out_pvalue1, out_MDa1, out_lcia1, out_hcia1, out_pvaluea1, out_MD2, out_lci2, out_hci2, out_pvalue2, out_MDa2, out_lcia2, out_hcia2, out_pvaluea2)
      out$out_variable = as.character(out$out_variable)
      
      # Add the heteroskedasticity p-values to the data frame
      out$hetero_pval_uni <- hetero_pval_uni
      out$hetero_pval_multi <- hetero_pval_multi
      
      # round the results
      out %<>% mutate_at(vars(out_MD1, out_lci1, out_hci1,out_MDa1, out_lcia1, out_hcia1,out_MD2, out_lci2, out_hci2,out_MDa2, out_lcia2, out_hcia2), funs(round(., 1)))
      
    # table 1  
    write.csv(out, file = "Results_lm_stunting_score.csv")
      
    rm(number, out_variable, out_MD1, out_lci1, out_hci1, out_pvalue1, out_MD2, out_lci2, out_hci2, out_pvalue2
       ,out_MDa1, out_lcia1, out_hcia1, out_pvaluea1, out_MDa2, out_lcia2, out_hcia2, out_pvaluea2, MD1, lci1, hci1
       ,pvalue1, MD2,lci2,hci2,pvalue2,MDa1,lcia1,hcia1,pvaluea1,MDa2,lcia2,hcia2,pvaluea2, j,outcome
       ,hetero_pval_uni,hetero_pval_multi,bp_uni,bp_multi)
    
    rm(out)
    

    # 3) Univariate analyses for variable selection for SEM----------
    str(mydata[,which(colnames(mydata) %in% c(household,family_simple,child_simple,diet,parasite_simple,"No_genes",inflammation_simple,bcaa))]) # check that they are all num, int or two levels factors
    str(mydata[,which(colnames(mydata) %in% c(development))])
    # remove factors variables with more than two levels for the moment
    vec <- match(x=c(household[!household == 'score_socio_eco_bycountry'], family_simple[!family_simple %in% c('instruc_mere_3','statut_mere')], child_simple[ !child_simple %in% c('haz','taille_naiss_p2')], diet, parasite_simple, "No_genes",inflammation_simple, bcaa),names(mydata))
    dvt <- match(x=c(development[!development == 'Nº_domains_delays']),names(mydata))
    
    # First for loop for outcomes (dependent variables)
    for (j in dvt) {
      print(j)
      out_variable=rep(NA, length(vec))
      out_MD=rep(NA, length(vec))
      out_lci = rep(NA, length(vec))
      out_hci = rep(NA, length(vec))
      out_pvalue=rep(NA, length(vec))
      number = 1
      
        # second for loop for exposures (independent variables)
        for(i in vec){
          print(i)
          exposure = colnames(mydata)[i]
          outcome = colnames(mydata)[j]
          model <- lm(get(outcome) ~ get(exposure), mydata)
          
          # Check if the model has enough coefficients (i.e. valid model)
          if (length(summary(model)$coefficients) > 4) {
          
          # extract results for each binary covariate
          MD <- summary(model)$coefficients[2,1]      # coefficient for count_ct
          lci <- confint(model)[2,1]  # lower bound for 95% CI, for count_ct
          hci <- confint(model)[2,2]  # higher bound for 95% CI, for count_ct
          pvalue <- summary(model)$coefficients[2,4]
          
          # store results for each binary covariate
          out_MD[number] = as.numeric(MD)
          out_lci[number] = as.numeric(lci)
          out_hci[number] = as.numeric(hci)
          out_pvalue[number] = as.numeric(pvalue)
          out_variable[number] = exposure
          } else {
          out_MD[number] <- NA
          out_lci[number] <- NA
          out_hci[number] <- NA
          out_pvalue[number] <- NA
          out_variable[number] <- exposure
          }
          number = number + 1
        }
      
      # Dataframe for results of binary variables
      out = data.frame(out_variable, out_MD, out_lci, out_hci, out_pvalue)
      out$out_variable = as.character(out$out_variable)
      
      # Models for multilevel variables
      # taille_naiss_p2 (3 levels)
      model <- lm(get(outcome) ~ taille_naiss_p2, mydata)
      for(k in 2:3) {   # extractions for levels 2 and 3
        sm <- c(paste("Birth size (reported)",levels(mydata$taille_naiss_p2)[k])
                ,summary(model)$coefficients[k,1]
                ,confint(model)[k,1]
                ,confint(model)[k,2]
                ,summary(model)$coefficients[k,4])
        out <- rbind(out,sm)
      }
      
      # statut_mere (3 levels)
      model <- lm(get(outcome) ~ statut_mere, mydata)
      for(k in 2:3) {   # extractions for levels 2 and 3
        sm <- c(paste("Maternal marital status",levels(mydata$statut_mere)[k])
                ,summary(model)$coefficients[k,1]
                ,confint(model)[k,1]
                ,confint(model)[k,2]
                ,summary(model)$coefficients[k,4])
        out <- rbind(out,sm)
      }
      
      # score_socio_eco_bycountry (3 levels)
      model <- lm(get(outcome) ~ score_socio_eco_bycountry, mydata)
      for(k in 2:3) {   # extractions for levels 2 and 3
        sm <- c(paste("Family economic status",levels(mydata$score_socio_eco_bycountry)[k])
                ,summary(model)$coefficients[k,1]
                ,confint(model)[k,1]
                ,confint(model)[k,2]
                ,summary(model)$coefficients[k,4])
        out <- rbind(out,sm)
      }
      
      # instruc_mere_3 (3 levels)
      model <- lm(get(outcome) ~ instruc_mere_3, mydata)
      for(k in 2:3) {   # extractions for levels 2 and 3
        sm <- c(paste("Maternal education level cat",levels(mydata$instruc_mere_3)[k])
                ,summary(model)$coefficients[k,1]
                ,confint(model)[k,1]
                ,confint(model)[k,2]
                ,summary(model)$coefficients[k,4])
        out <- rbind(out,sm)
      }
      
      assign(paste0("out", j), out)
      
    }
    
    rm(out,out_MD,out_lci,out_hci,out_pvalue,out_variable,outcome)
    
    library(ggplot2)
    library(scales)
    
    # combine all results
    # df = rbind(out100, out103, out68, out76, out84, out92, )
    df = do.call("rbind",mget(ls(pattern = "^out*")))
    df$delay = factor(c(rep("PES",nrow(df)/6),rep("Overall",nrow(df)/6),rep("Comm",nrow(df)/6),rep("GM",nrow(df)/6),rep("FM",nrow(df)/6),rep("PS",nrow(df)/6)), levels = c("Overall","Comm","GM","FM","PS","PES"))
    df$out_MD = as.numeric(df$out_MD)
    df$out_lci = as.numeric(df$out_lci)
    df$out_hci = as.numeric(df$out_hci)
    df$out_pvalue = as.numeric(df$out_pvalue)
    write.csv(df, file="Results_lm_score_univariate.csv")
    # remove non-significant results
    df1 = df[which(df$out_pvalue <= 0.05),]
    
    pdf('Results_lm_score_univariate.pdf', height = 12, width = 10) 				# Write to PDF
    ggplot(df1, aes(x=out_variable, y=out_MD)) +
      geom_point(stat="identity", aes(colour=delay), position=position_dodge(width=.6)) +
      geom_errorbar(aes(ymin=out_lci, ymax=out_hci, colour=delay), width=.1, position=position_dodge(width=.6)) +
      facet_grid(delay ~  . , scales="free_y", switch = "y") +
      geom_hline(yintercept = 0) +
      ylab("Mean difference") +
      xlab("") +
      theme_bw() +						# remove grey background
      theme(legend.title=element_blank(), axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=.5)
            ,legend.position="none", panel.grid.major = element_blank()  # remove x and y major grid lines
            ,strip.text.x = element_text(size = 8), axis.text.y = element_text(angle = 0)
            ) +
      scale_y_continuous(position="right") +
      scale_x_discrete(limits = c("Family economic status 2","nb_pieces","salaires_cum","eau","eau_traitee"
                                  ,"age_prem_gross","age_prem_gross_15","mere_sousnutris","Maternal education level cat 2","Maternal education level cat 3"
                                  ,"sexe","age","haz_cont","waz_cont","taille","poids","Birth size (reported) 3","ferritine2","ferritine2_seuil","allaite","age_allaite","age_alim","hemoglobine2","hemoglobine2_seuil"
                                  ,"dds"
                                  ,"ascaris","trichuris","blastocystis","enterobius"
                                  ,"No_genes"
                                  ,"crp2_cat","AATlevel2"
                                  ,"AlanineM","CitrulineM","IsoleucineM","LeucineM","ValineM"),
                       labels=c( labelled::look_for(mydata[,which(names(mydata) %in% child_simple)])$label
                                   ,labelled::look_for(mydata[,which(names(mydata) %in% family_simple)])$label
                                   ,labelled::look_for(mydata[,which(names(mydata) %in% household)])$label
                                   ,labelled::look_for(mydata[,which(names(mydata) %in% inflammation_simple)])$label
                                   ,labelled::look_for(mydata[,which(names(mydata) %in% parasite_simple)])$label
                                   ,labelled::look_for(mydata[,which(names(mydata) %in% diet)])$label
                                   ,labelled::look_for(mydata[,which(names(mydata) %in% bcaa)])$label ) )
    dev.off()
    
    rm(sm,model,exposure,dvt, vec,i,j,MD,lci,hci,pvalue,number,df)
    rm(list = ls(, pattern = "out"))
  
    
    # 4) Multivariate model for neurodevelopment ----------
    library(GGally)
    df1$out_variable[df1$delay == "Overall"]
    ggpairs(mydata[,which(names(mydata) %in% df1$out_variable[df1$delay=="Overall"])]) # pairwise linear correlations between explanatory variables, for all significant variables in univariate analyses
    # remove poids (weight) and taille (height) because they are too correlated with haz and waz
    # remove valien because too correlated with leucine
    model1 <- lm(Overal_score ~ score_socio_eco_bycountry + salaires_cum + eau_traitee + nb_pieces + instruc_mere_3 + age_prem_gross + mere_sousnutris + haz_cont + waz_cont + hemoglobine2 + ascaris + LeucineM + taille_naiss_p2 , mydata)
    summary(model1)
    nobs(model1) # 279
    par(mfrow=c(2,2))
    plot(model1)
    bptest(model1) # residuals not normally distributed
    
    
    df1$out_variable[df1$delay == "Comm"]
    ggpairs(mydata[,which(names(mydata) %in% df1$out_variable[df1$delay=="Comm"])]) # pairwise linear correlations between explanatory variables, for all significant variables in univariate analyses
    model2 <- lm(Comm_tot ~ poids + eau_traitee + nb_pieces + age + ferritine2_seuil + ascaris + trichuris + CitrulineM + instruc_mere_3, data = mydata)
    summary(model2)
    nobs(model2) #291
    par(mfrow=c(2,2))
    plot(model2)
    bptest(model2) # residuals normally distributed
    
    df1$out_variable[df1$delay == "PS"]
    ggpairs(mydata[,which(names(mydata) %in% df1$out_variable[df1$delay=="PS"])]) # pairwise linear correlations between explanatory variables, for all significant variables in univariate analyses
    model3 <- lm(PS_tot ~ score_socio_eco_bycountry + salaires_cum + eau + eau_traitee + nb_pieces + instruc_mere_3 + age_prem_gross + mere_sousnutris + age + haz_cont + waz_cont + ascaris + trichuris + blastocystis + AlanineM , data = mydata)
    summary(model3)
    nobs(model3) #284
    par(mfrow=c(2,2))
    plot(model3)
    bptest(model3) # residuals normally distributed
    
    df1$out_variable[df1$delay == "PES"]
    ggpairs(mydata[,which(names(mydata) %in% df1$out_variable[df1$delay=="PES"])]) # # pairwise linear correlations between explanatory variables, for all significant variables in univariate analyses
    # remove poids (weight) and taille (height) because too correlated with haz and age
    # remove valine because too correlated with leucine
    model4 <- lm(PES_tot ~ nb_pieces + sexe + age + haz_cont + age_allaite + hemoglobine2 + dds + No_genes + CitrulineM + LeucineM, data = mydata)
    summary(model4)
    nobs(model4) #259
    par(mfrow=c(2,2))
    plot(model4)
    bptest(model4) # residuals not normally distributed
    
    df1$out_variable[df1$delay == "FM"]
    ggpairs(mydata[,which(names(mydata) %in% df1$out_variable[df1$delay=="FM"])]) # # pairwise linear correlations between explanatory variables, for all significant variables in univariate analyses
    # Remove poids (weight) and taille (height) because too correlated with haz and waz
    # remove valine because too correlated with leucine
    # remove age_prem_gross_15 because too correlated with age_prem_gross
    model5 <- lm(FM_tot ~ salaires_cum + nb_pieces + age_prem_gross + haz_cont + waz_cont + ascaris + trichuris + crp2_cat + LeucineM + score_socio_eco_bycountry + instruc_mere_3, data = mydata)
    summary(model5)
    nobs(model5) #280
    par(mfrow=c(2,2))
    plot(model5)
    bptest(model5) # residuals normally distributed
    
    df1$out_variable[df1$delay == "GM"]
    ggpairs(mydata[,which(names(mydata) %in% df1$out_variable[df1$delay=="GM"])]) # # pairwise linear correlations between explanatory variables, for all significant variables in univariate analyses
    # remove ValineM and leucineM because too correlated with isoleucineM
    # remove poids (weight) and taille (height) because too correlated with haz and waz
    # remove hemoglobine_seuil because too correlated with hemoglobine2
    # remove ferritine_seuil because too correlated with ferritine2
    model6 <- lm(GM_tot ~ age + haz_cont + waz_cont + ferritine2 + allaite + hemoglobine2 + No_genes +  AATlevel2 + IsoleucineM + taille_naiss_p2, data = mydata)
    summary(model6)
    nobs(model6) #271
    par(mfrow=c(2,2))
    plot(model6)
    bptest(model6) # residuals not normally distributed
    
    df2 = data.frame( rbind( cbind( rep("Overall",length(model1$coefficients)), model1$coefficients,confint(model1), summary(model1)$coefficients[, "Pr(>|t|)"])
                            ,cbind( rep("Comm",length(model2$coefficients)), model2$coefficients,confint(model2), summary(model2)$coefficients[, "Pr(>|t|)"])
                            ,cbind( rep("PS",length(model3$coefficients)), model3$coefficients,confint(model3), summary(model3)$coefficients[, "Pr(>|t|)"])
                            ,cbind( rep("PES",length(model4$coefficients)), model4$coefficients,confint(model4), summary(model4)$coefficients[, "Pr(>|t|)"])
                            ,cbind( rep("FM",length(model5$coefficients)), model5$coefficients,confint(model5), summary(model5)$coefficients[, "Pr(>|t|)"])
                            ,cbind( rep("GM",length(model6$coefficients)), model6$coefficients,confint(model6), summary(model6)$coefficients[, "Pr(>|t|)"]) ) )
    
    colnames(df2) <- c("delay","MD","lci","hci","pval")
    df2$Variable = rownames(df2)
    df2$MD = as.numeric(df2$MD)
    df2$lci = as.numeric(df2$lci)
    df2$hci = as.numeric(df2$hci)
    df2$pval = as.numeric(df2$pval)
    df2 = df2[-grep("Intercept",df2$Variable),]
    df2$Variable = gsub('\\.1', '', df2$Variable)
    df2$Variable = gsub('\\.2', '', df2$Variable)
    df2$Variable = gsub('\\.3', '', df2$Variable)
    df2$Variable = gsub('\\.4', '', df2$Variable)
    df2$signif = ifelse(df2$lci > 0, "8", ifelse(df2$hci < 0, "8" , "1"))
    df2$delay = factor(df2$delay, levels = c("Overall","Comm","GM","FM","PS","PES"))
    df2$Variable[which(df2$Variable %in% c("eau1", "eau_traitee1", "mere_sousnutris1", "ferritine2_seuil1", "allaite1"))] = gsub('1', '', df2$Variable[which(df2$Variable %in% c("eau1", "eau_traitee1", "mere_sousnutris1", "ferritine2_seuil1", "allaite1"))] ) 
    df2$Variable[df2$Variable == "sexe2"] = gsub('2', '', df2$Variable[df2$Variable == "sexe2"])
    df2$Variable[df2$Variable == "crp2_catelevated"] = gsub('elevated', '', df2$Variable[df2$Variable == "crp2_catelevated"])
    df2$Variable[df2$Variable == "AATlevel2elevated"] = gsub('elevated', '', df2$Variable[df2$Variable == "AATlevel2elevated"])
    df2$Variable[df2$Variable == "score_socio_eco_bycountry2"] <- "Family economic status 2"
    df2$Variable[df2$Variable == "score_socio_eco_bycountry3"] <- "Family economic status 3"
    df2$Variable[df2$Variable == "instruc_mere_32"] <- "Maternal education level cat 2"
    df2$Variable[df2$Variable == "instruc_mere_33"] <- "Maternal education level cat 3"
    df2$Variable[df2$Variable == "taille_naiss_p22"] <- "Birth size (reported) 2"
    df2$Variable[df2$Variable == "taille_naiss_p23"] <- "Birth size (reported) 3"
    write.csv(df2, file="Results_lm_score_multivariate.csv")

    pdf('Results_lm_score_multivariate.pdf', height = 12, width = 10) 				# Write to PDF
    ggplot(df2, aes(x=Variable, y=MD)) +
      geom_point(stat="identity", aes(colour=delay, shape = signif), position=position_dodge(width=.6), size = 2.5) +
      geom_errorbar(aes(ymin=lci, ymax=hci, colour=delay), width=.1, position=position_dodge(width=.6)) +
      facet_grid(delay ~  . , scales="free_y", switch = "y") +
      geom_hline(yintercept = 0) +
      ylab("Mean difference") +
      xlab("") +
      theme_bw() +						# remove grey background
      theme(legend.title=element_blank(), axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=.5)
            ,legend.position="bottom", panel.grid.major = element_blank()  # remove x and y major grid lines
            ,strip.text.x = element_text(size = 8), axis.text.y = element_text(angle = 0)
      ) +
      scale_shape_manual("", values = c("1" = 16, "8" = 8), labels = c("NS", "p-value < 0.05")) +
      scale_y_continuous(position="right") +
      scale_x_discrete(limits = c("Family economic status 2","Family economic status 3","nb_pieces","salaires_cum","eau","eau_traitee"
                                  ,"age_prem_gross","mere_sousnutris","Maternal education level cat 2","Maternal education level cat 3"
                                  ,"sexe","age","haz_cont","waz_cont","Birth size (reported) 2","Birth size (reported) 3","ferritine2","ferritine2_seuil","allaite","age_allaite","age_alim","hemoglobine2"
                                  ,"dds"
                                  ,"ascaris","trichuris","blastocystis"
                                  ,"No_genes"
                                  ,"crp2_cat","AATlevel2"
                                  ,"AlanineM","CitrulineM","IsoleucineM","LeucineM"),
                       labels=c( labelled::look_for(mydata[,which(names(mydata) %in% child_simple)])$label
                                 ,labelled::look_for(mydata[,which(names(mydata) %in% family_simple)])$label
                                 ,labelled::look_for(mydata[,which(names(mydata) %in% household)])$label
                                 ,labelled::look_for(mydata[,which(names(mydata) %in% inflammation_simple)])$label
                                 ,labelled::look_for(mydata[,which(names(mydata) %in% parasite_simple)])$label
                                 ,labelled::look_for(mydata[,which(names(mydata) %in% diet)])$label
                                 ,labelled::look_for(mydata[,which(names(mydata) %in% bcaa)])$label ) )
    dev.off()
    
    rm(model1,model2,model3,model4,model5,model6,df1, df2)
    
    # significant variables for at least one domain: haz (cont.), instruc_mere_3, age, ascaris, birth size (reported), eau (running water), eau_traitee (treated drinking water), blastocystis, sex, No_genes (nb of enteropathgen genes), 
    # other variables from literature: hemoglobine2, socioeconomic status and maternal factors: age_prem_gross, nb_pieces, score_socio_eco_bycountry
    

    
    
    
# III) Microbiome analyses --------

  # 1) merge with microbiota and microbiota.asv tables--------
  setdiff(mydata$ID_metag,microbiota$ID_metag)
  setdiff(microbiota$ID_metag,mydata$ID_metag)
  mytotal = inner_join(mydata,microbiota, by=c("ID_metag"))
  
  setdiff(mydata$ID_metag,microbiota.asv$ID_metag)
  setdiff(microbiota.asv$ID_metag,mydata$ID_metag) # the same as for microbiota dataset
  mytotal.asv = inner_join(mydata,microbiota.asv, by=c("ID_metag"))
  
  setdiff(mydata$ID_metag,m_fam$ID_metag)
  setdiff(m_fam$ID_metag,mydata$ID_metag)
  mytotal.fam = inner_join(mydata,m_fam, by=c("ID_metag"))
  
  setdiff(mydata$ID_metag,m_gen$ID_metag)
  setdiff(m_gen$ID_metag,mydata$ID_metag)
  mytotal.gen = inner_join(mydata,m_gen, by=c("ID_metag"))
  
  # check if samples with total abundance of 0
  rowSums(mytotal[,which(colnames(mytotal)=="g_Prevotella_9"):ncol(mytotal)])
  rowSums(mytotal.asv[,which(colnames(mytotal.asv)=="Prevotella_copri_DSM_18205_ASV3"):ncol(mytotal.asv)]) # Same as above, normal since for mytotal the ASVs were summed up, and nothing has been trimmed so far
  rowSums(mytotal[,which(colnames(mytotal)=="g_Prevotella_9"):ncol(mytotal)]) == rowSums(mytotal.asv[,which(colnames(mytotal.asv)=="Prevotella_copri_DSM_18205_ASV3"):ncol(mytotal.asv)]) # Same as above, normal since for mytotal the ASVs were summed up, and nothing has been trimmed so far
  
  rowSums(mytotal.fam[,which(colnames(mytotal.fam)=="Prevotellaceae"):ncol(mytotal.fam)])
  rowSums(mytotal.fam[,which(colnames(mytotal.fam)=="Prevotellaceae"):ncol(mytotal.fam)]) == rowSums(mytotal.asv[,which(colnames(mytotal.asv)=="Prevotella_copri_DSM_18205_ASV3"):ncol(mytotal.asv)]) # Same as above, normal since for mytotal the ASVs were summed up, and nothing has been trimmed so far
  rowSums(mytotal.fam[,which(colnames(mytotal.fam)=="Prevotellaceae"):ncol(mytotal.fam)]) - rowSums(mytotal.asv[,which(colnames(mytotal.asv)=="Prevotella_copri_DSM_18205_ASV3"):ncol(mytotal.asv)]) # Same as above, normal since for mytotal the ASVs were summed up, and nothing has been trimmed so far
  
  rowSums(mytotal.gen[,which(colnames(mytotal.gen)=="Prevotella_9"):ncol(mytotal.gen)])
  rowSums(mytotal.gen[,which(colnames(mytotal.gen)=="Prevotella_9"):ncol(mytotal.gen)]) == rowSums(mytotal.asv[,which(colnames(mytotal.asv)=="Prevotella_copri_DSM_18205_ASV3"):ncol(mytotal.asv)]) # Same as above, normal since for mytotal the ASVs were summed up, and nothing has been trimmed so far
  rowSums(mytotal.gen[,which(colnames(mytotal.gen)=="Prevotella_9"):ncol(mytotal.gen)]) - rowSums(mytotal.asv[,which(colnames(mytotal.asv)=="Prevotella_copri_DSM_18205_ASV3"):ncol(mytotal.asv)]) # Same as above, normal since for mytotal the ASVs were summed up, and nothing has been trimmed so far
  

  # 2) Description of reads and taxa----------
  # since we removed 349-317 = 32 participants with no cognitive data, some bacteria are not represented in the dataset anymore
  # order taxa by abundance - at the species level
  cs = colSums( mytotal[,which(colnames(mytotal)=="g_Prevotella_9"):ncol(mytotal)], na.rm=T)
  reads = mytotal[,which(colnames(mytotal)=="g_Prevotella_9"):ncol(mytotal)]
  reads = reads[,order(cs, decreasing=T)]
  rownames(reads) = mytotal$ID_afri
  mytotal = as.data.frame(cbind(mytotal[,1:(which(colnames(mytotal)=="g_Prevotella_9")-1)],reads))
  
  # remove bacteria with no occurrence in the dataset
  remove = colnames(reads)[which(colSums(reads, na.rm=T) == 0)]
  which(colnames(mytotal) %in% remove) # 40 columns
  mytotal = mytotal[,- (which(colnames(mytotal) %in% remove))]
  reads = reads[,- (which(colnames(reads) %in% remove))]
  
  # order taxa by abundance - at the ASV level
  cs = colSums( mytotal.asv[,which(colnames(mytotal.asv)=="Prevotella_copri_DSM_18205_ASV3"):ncol(mytotal.asv)],na.rm=T)
  reads.asv = mytotal.asv[,which(colnames(mytotal.asv)=="Prevotella_copri_DSM_18205_ASV3"):ncol(mytotal.asv)]
  reads.asv = reads.asv[,order(cs, decreasing=T)]
  rownames(reads.asv) = mytotal.asv$ID_afri
  mytotal.asv = as.data.frame(cbind(mytotal.asv[,1:(which(colnames(mytotal.asv)=="Prevotella_copri_DSM_18205_ASV3")-1)],reads.asv))
  
  # remove ASV with no occurrence in the dataset
  remove = colnames(reads.asv)[which(colSums(reads.asv,na.rm=T) == 0)]
  which(colnames(mytotal.asv) %in% remove)
  mytotal.asv = mytotal.asv[,- (which(colnames(mytotal.asv) %in% remove))]
  reads.asv = reads.asv[,- (which(colnames(reads.asv) %in% remove))]
  
  # order taxa by abundance - at the genus level
  cs = colSums( mytotal.gen[,which(colnames(mytotal.gen)=="Prevotella_9"):ncol(mytotal.gen)],na.rm=T)
  reads.gen = mytotal.gen[,which(colnames(mytotal.gen)=="Prevotella_9"):ncol(mytotal.gen)]
  reads.gen = reads.gen[,order(cs, decreasing=T)]
  rownames(reads.gen) = mytotal.gen$ID_afri
  mytotal.gen = as.data.frame(cbind(mytotal.gen[,1:(which(colnames(mytotal.gen)=="Prevotella_9")-1)],reads.gen))
  
  # remove genera with no occurrence in the dataset
  remove = colnames(reads.gen)[which(colSums(reads.gen,na.rm=T) == 0)]
  which(colnames(mytotal.gen) %in% remove)
  mytotal.gen = mytotal.gen[,- (which(colnames(mytotal.gen) %in% remove))]
  reads.gen = reads.gen[,- (which(colnames(reads.gen) %in% remove))]
  
  # order taxa by abundance - at the family level
  cs = colSums( mytotal.fam[,which(colnames(mytotal.fam)=="Prevotellaceae"):ncol(mytotal.fam)],na.rm=T)
  reads.fam = mytotal.fam[,which(colnames(mytotal.fam)=="Prevotellaceae"):ncol(mytotal.fam)]
  reads.fam = reads.fam[,order(cs, decreasing=T)]
  rownames(reads.fam) = mytotal.fam$ID_afri
  mytotal.fam = as.data.frame(cbind(mytotal.fam[,1:(which(colnames(mytotal.fam)=="Prevotellaceae")-1)],reads.fam))
  
  # remove families with no occurrence in the dataset
  remove = colnames(reads.fam)[which(colSums(reads.fam,na.rm=T) == 0)]
  which(colnames(mytotal.fam) %in% remove)
  mytotal.fam = mytotal.fam[,- (which(colnames(mytotal.fam) %in% remove))]
  reads.fam = reads.fam[,- (which(colnames(reads.fam) %in% remove))]
  
  rm(remove,cs)
  
  # save a vector of all the taxa names in the dataset
  bacteria = colnames(reads)
  asv = colnames(reads.asv)
  genera = colnames(reads.gen)
  families = colnames(reads.fam)
  # number of sequences, mean and sd
  sum(reads.asv) # 6 986 362 reads in the database
  mean = sum(reads.asv)/nrow(reads.asv) # 22 039 reads on average per sample
  row = rowSums(reads.asv)
  diff = (row-mean)^2
  sqrt(sum(diff)/(nrow(reads)-1)) # 18 867 reads as standard deviation
  
  rm(mean,row,diff)

  
  
  # 3) alpha diversity----------

  library(vegan)
  
  # At the ASV level only
  mytotal$ID_afri == mytotal.asv$ID_afri
  # inverse simpson index
  mytotal$invsimpson <- diversity(reads.asv,index = "invsimpson")
  hist(mytotal$invsimpson)
  # shannon index
  mytotal$shannon <- diversity(reads.asv,index = "shannon")
  hist(mytotal$shannon)    
  # richness
  # check if there are rare taxa in the dataset (with many samples that have counts of less than 10 and more than 0)
  k = 0
  rarespecies <- vector("numeric", 1630L)
  for (i in colnames(reads.asv)) {
    k = k+1
    rarespecies[k] = sum(rle(mytotal.asv[,i]>0 & mytotal.asv[,i]<=10)$values, na.rm=T)
    print(k)
  }
  rarespecies # gives the number of samples with non-null values inferior or equal to 10 for each ASV
  table(rarespecies) # 1626 - 883 = 743 ASVs that are rare in a sample
  max(rarespecies) # a max of 29 samples for which a species is rare, out of 317, so 9%
  
  mytotal$richness <- vegan::estimateR(reads.asv)[1,] # very weird because when NA values for microbiome data it calculates richness = 486 species
  hist(mytotal$richness) # normal distribution
  rm(i, k, rarespecies)
  
  # inverse simpson index
  mytotal.asv$invsimpson <- diversity(reads.asv,index = "invsimpson")
  hist(mytotal.asv$invsimpson)
  # shannon index
  mytotal.asv$shannon <- diversity(reads.asv,index = "shannon")
  hist(mytotal.asv$shannon)  
  #richness
  mytotal.asv$richness <- vegan::estimateR(reads.asv)[1,]
  hist(mytotal.asv$richness) # normal distribution
  

    # a) univariate linear regression and plots-----------
    ggplot( mytotal.asv, aes(y = PS_tot, x = invsimpson)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(label.x.npc = "left", label.y.npc = "bottom", size=2.5) +
      theme_minimal()
    
    mdl = lm(formula=PS_tot ~ invsimpson, data = mytotal.asv)
    summary(mdl)
    par(mfrow=c(2,2))
    plot(mdl)
  
    #-- residuals
    qqnorm(resid(mdl))
    hist(resid(mdl))
    
    #-- predict
    mytotal.asv$pred_PS = predict(mdl)
    ggplot(mytotal.asv) +
      geom_point(aes(y = pred_PS, x = PS_tot)) +
      geom_abline(slope=1, intercept=0) # très mauvais job --> add the AgeASQ3 variable to have a better prediction?
    
    plotStatsLin <- function(df, var2, var3)
    {
      var2 <- ensym(var2)
      var3 <- ensym(var3)
      
      #f <- as.formula(paste(as.character(var3), "~", as.character(var2)))
      
      gg = ggplot(df, aes(x = !!var2, y = !!var3)) +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "bottom", size=3.5) +
        theme_minimal()
      return(gg)
    }
    
    p1 = plotStatsLin(mytotal.asv, invsimpson, Comm_tot)
    p2 = plotStatsLin(mytotal.asv, shannon, Comm_tot)
    p3 = plotStatsLin(mytotal.asv, invsimpson, PS_tot)
    p4 = plotStatsLin(mytotal.asv, shannon, PS_tot)
    p5 = plotStatsLin(mytotal.asv, invsimpson, PES_tot)
    p6 = plotStatsLin(mytotal.asv, shannon, PES_tot)
    p7 = plotStatsLin(mytotal.asv, invsimpson, GM_tot)
    p8 = plotStatsLin(mytotal.asv, shannon, GM_tot)
    p9 = plotStatsLin(mytotal.asv, invsimpson, FM_tot)
    p10 = plotStatsLin(mytotal.asv, shannon, FM_tot)
    p11 = plotStatsLin(mytotal.asv, invsimpson, Overal_score)
    p12 = plotStatsLin(mytotal.asv, shannon, Overal_score)
    
    pdf(file="alpha_diversity_delay_cont.pdf", height=10, width=12)
    ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 + rremove("x.text"),
              labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
              ncol = 4, nrow = 3) # warnings() is because there is no value for lines that are NA
    dev.off() 
    
    rm(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12)
    
    p1 = plotStatsLin(mytotal.asv, richness, Comm_tot)
    p2 = plotStatsLin(mytotal.asv, richness, PS_tot)
    p3 = plotStatsLin(mytotal.asv, richness, PES_tot)
    p4 = plotStatsLin(mytotal.asv, richness, GM_tot)
    p5 = plotStatsLin(mytotal.asv, richness, FM_tot)
    p6 = plotStatsLin(mytotal.asv, richness, Overal_score)
    
    pdf(file="Richness_delay_cont.pdf", height=8, width=10)
    ggarrange(p1, p2, p3, p4, p5, p6 + rremove("x.text"),
              labels = c("A", "B", "C", "D", "E", "F"),
              ncol = 3, nrow = 2)
    dev.off() # warnings() is because there is no value for lines that are NA
    
    rm(p1, p2, p3, p4, p5, p6)
    rm(mdl)
    
    # b) Linear regression controlling for AgeASQ3 and haz------------
    summary(lm(formula = Comm_tot ~ invsimpson + AgeASQ3 + haz, data=mytotal.asv))
    summary(lm(formula = Comm_tot ~ shannon + AgeASQ3 + haz, data=mytotal.asv))
    summary(lm(formula = Comm_tot ~ richness + AgeASQ3 + haz, data=mytotal.asv))
    
    summary(lm(formula = PS_tot ~ invsimpson + AgeASQ3 + haz, data=mytotal.asv))
    summary(lm(formula = PS_tot ~ shannon + AgeASQ3 + haz, data=mytotal.asv))
    summary(lm(formula = PS_tot ~ richness + AgeASQ3 + haz, data=mytotal.asv))
    
    summary(lm(formula = PES_tot ~ invsimpson + AgeASQ3 + haz, data=mytotal.asv))
    summary(lm(formula = PES_tot ~ shannon + AgeASQ3 + haz, data=mytotal.asv))
    summary(lm(formula = PES_tot ~ richness + AgeASQ3 + haz, data=mytotal.asv))
    
    summary(lm(formula = GM_tot ~ invsimpson + AgeASQ3 + haz, data=mytotal.asv))
    summary(lm(formula = GM_tot ~ shannon + AgeASQ3 + haz, data=mytotal.asv))
    summary(lm(formula = GM_tot ~ richness + AgeASQ3 + haz, data=mytotal.asv))
    
    summary(lm(formula = FM_tot ~ invsimpson + AgeASQ3 + haz, data=mytotal.asv))
    summary(lm(formula = FM_tot ~ shannon + AgeASQ3 + haz, data=mytotal.asv))
    summary(lm(formula = FM_tot ~ richness + AgeASQ3 + haz, data=mytotal.asv))
    
    summary(lm(formula = Overal_score ~ invsimpson + AgeASQ3 + haz, data=mytotal.asv))
    summary(lm(formula = Overal_score ~ shannon + AgeASQ3 + haz, data=mytotal.asv))
    summary(lm(formula = Overal_score ~ richness + AgeASQ3 + haz, data=mytotal.asv))
    # the alpha-diversity indices are never significant
    
    # c) Linear regression controlling for sequencing run---------
    # make table
    library(broom)
    outcomes <- c("Comm_tot", "PS_tot", "PES_tot", "GM_tot", "FM_tot", "Overal_score")
    predictors <- c("invsimpson", "shannon", "richness")
    
    results_list <- list()
    number <- 1
    for (y in outcomes) {
      for (x in predictors) {
        # Fit model
        f <- as.formula(paste(y, "~", x, "+ run"))
        model <- lm(f, data = mytotal.asv)
        
        # Extract coefficient row for the predictor
        coef_table <- tidy(model)
        predictor_row <- coef_table[coef_table$term == x, ]
        
        # Extract model-level statistics
        model_stats <- glance(model)
        
        # Store results
        results_list[[number]] <- data.frame(
          outcome = y,
          predictor = x,
          beta = predictor_row$estimate,
          p_value = predictor_row$p.value,
          r_squared = model_stats$r.squared
        )
        number <- number + 1
      }
    }
    
    results <- do.call(rbind, results_list)
    results <- as.data.frame(results)
    write.csv(results, file = "Results_alpha_diversity_score_multivariate.csv")
    rm(results, results_list, predictors, outcomes, number, model_stats, model, f, coef_table, x, y)
    
 
  
  
  # 4) Transformation of the data ---------
      
    # a) at lowest available rank---------
  
        # 1) CLR transformation
        # three possibilities for taxa with 0 reads:
        # a) add 1 everywhere, substract 1 from resulting G(x) at the end
        # b) add 1 to zero values only
        # c) Impute zero counts with zCompositions package
        imp1 = zCompositions::cmultRepl(reads, label = 0, method = "BL", output = "p-counts")
        
        # "GBM" method doesn't work if only one positive value
        # Check how many columns have only zeros except for one non-zero value (singletons)
        checkNumZerosCol <- apply(reads,2,function(x) sum(x==0))
        cases <- which(checkNumZerosCol == (nrow(reads) - 1))
        length(cases) # 111 columns are all zeros but one positive value
        zCompositions::cmultRepl(reads[,-cases], label = 0, method = "GBM", output = "p-counts") # GBM imputation without them works but clr transformation can't be performed since GM vector not same length
        
        geoMeans = apply(imp1, 2, function(var) exp(mean(log(var))) ) # one GM for each taxon
        
        # création d'une fonction qui divise ma var1 par ma var2 puis renvoie le log de cette division
        divislog <- function(obs) {
          resuldiv <- obs/geoMeans
          resullog <- log(resuldiv)
          return(resullog)
        }
        reads_clr = apply(imp1, 1, divislog)
        reads_clr = as.data.frame(t(reads_clr))
        
        # d) leave zeros as blank. We use this last definition:
        # geoMeans = apply(reads, 2, function(var) if (all(var == 0)) 0 else exp(mean(log(var[var != 0]))))
        # clr = apply(reads, 2, function(x) if (all(x == 0)) 0 else log(x/geoMeans)) # pb with zeros
        # clr2 = as.data.frame(clr)
        
        rm(cases, checkNumZerosCol, imp1, geoMeans)
        
        # 2) Relative abundance transformation
        reads_ra = apply(reads,1,function(x) 100*x/sum(x))	# this function transposes the dataframe. We have to transpose it back
        # transpose matrix
        reads_ra = as.data.frame(t(reads_ra))
      
        
    # b) at ASV level---------
      
        # 1) CLR transformation
        imp1 = zCompositions::cmultRepl(reads.asv, label = 0, method = "BL", output = "p-counts")
        
        # "GBM" method doesn't work if only one positive value
        # Check how many columns have only zeros except for one non-zero value (singletons)
        checkNumZerosCol <- apply(reads.asv,2,function(x) sum(x==0))
        cases <- which(checkNumZerosCol == (nrow(reads.asv) - 1))
        length(cases) # 613 columns are all zeros but one positive value
        zCompositions::cmultRepl(reads.asv[,-cases], label = 0, method = "GBM", output = "p-counts") # GBM imputation without them works but clr transformation can't be performed since GM vector not same length
        
        geoMeans = apply(imp1, 2, function(var) exp(mean(log(var))) ) # one GM for each taxon
        
        reads.asv_clr = apply(imp1, 1, divislog)
        reads.asv_clr = as.data.frame(t(reads.asv_clr))
        
        rm(cases, checkNumZerosCol, imp1, geoMeans)
        
        # 2) Relative abundance transformation
        reads.asv_ra = apply(reads.asv,1,function(x) 100*x/sum(x))	# this function transposes the dataframe. We have to transpose it back
        # transpose matrix
        reads.asv_ra = as.data.frame(t(reads.asv_ra))
        
    # c) at the genus level---------
      
        # 1) CLR transformation
        imp1 = zCompositions::cmultRepl(reads.gen, label = 0, method = "BL", output = "p-counts")
        
        # "GBM" method doesn't work if only one positive value
        # Check how many columns have only zeros except for one non-zero value (singletons)
        checkNumZerosCol <- apply(reads.gen,2,function(x) sum(x==0))
        cases <- which(checkNumZerosCol == (nrow(reads.gen) - 1))
        length(cases) # 37 columns are all zeros but one positive value
        zCompositions::cmultRepl(reads.gen[,-cases], label = 0, method = "GBM", output = "p-counts") # GBM imputation without them works but clr transformation can't be performed since GM vector not same length
        
        geoMeans = apply(imp1, 2, function(var) exp(mean(log(var))) ) # one GM for each taxon
        
        reads.gen_clr = apply(imp1, 1, divislog)
        reads.gen_clr = as.data.frame(t(reads.gen_clr))
        
        rm(cases, checkNumZerosCol, imp1, geoMeans)
        
        # 2) Relative abundance transformation
        reads.gen_ra = apply(reads.gen,1,function(x) 100*x/sum(x))	# this function transposes the dataframe. We have to transpose it back
        # transpose matrix
        reads.gen_ra = as.data.frame(t(reads.gen_ra))
      
      
    # d) at the family level---------

        # 1) CLR transformation
        imp1 = zCompositions::cmultRepl(reads.fam, label = 0, method = "BL", output = "p-counts")
        
        # "GBM" method doesn't work if only one positive value
        # Check how many columns have only zeros except for one non-zero value (singletons)
        checkNumZerosCol <- apply(reads.fam,2,function(x) sum(x==0))
        cases <- which(checkNumZerosCol == (nrow(reads.fam) - 1))
        length(cases) # 9 columns are all zeros but one positive value
        zCompositions::cmultRepl(reads.fam[,-cases], label = 0, method = "GBM", output = "p-counts") # GBM imputation without them works but clr transformation can't be performed since GM vector not same length
        
        geoMeans = apply(imp1, 2, function(var) exp(mean(log(var))) ) # one GM for each taxon
        
        reads.fam_clr = apply(imp1, 1, divislog)
        reads.fam_clr = as.data.frame(t(reads.fam_clr))
        
        rm(cases, checkNumZerosCol, imp1, geoMeans)
        
        # 2) Relative abundance transformation
        reads.fam_ra = apply(reads.fam,1,function(x) 100*x/sum(x))	# this function transposes the dataframe. We have to transpose it back
        # transpose matrix
        reads.fam_ra = as.data.frame(t(reads.fam_ra))
      
        
  # 5) barplots ----------
  
    # a) at the highest available rank-----------
    # order samples by score values (from smaller to bigger and then plot the top species/bacteria to see if there is a trend)
    Top20 = colnames(mytotal[,which(colnames(mytotal) == "g_Prevotella_9"):(which(colnames(mytotal) == "g_Prevotella_9") + 19)]) 
    Top20
  
    mytotal1 = as.data.frame(cbind(mytotal[,c("ID_afri","Overal_score")], reads_ra[which(colnames(reads_ra) %in% Top20)]))
    mytotal1 = mytotal1[order(mytotal1$Overal_score),]
    mytotal1_long <- reshape2::melt(mytotal1, id.vars = c("ID_afri","Overal_score"))
    mytotal1_long$ID_afri <- factor(mytotal1_long$ID_afri, levels = unique(mytotal1_long$ID_afri))
    
    # Color palette creation from Brewer
    library(RColorBrewer)
    palbrew_info <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]  # Extract color info
    palbrew_all <- unlist(mapply(brewer.pal,                     # Create vector with all colors
                                 palbrew_info$maxcolors,
                                 rownames(palbrew_info)))
    pie(rep(1, 74), col = palbrew_all)
    # Set color palette (number can be changed)
    pal_20 <- palbrew_all[1:20]
    names(pal_20) = Top20
    
    phyloColorTbl = c("red1","lightseagreen",
                                  "yellow","chartreuse","gray75",
                                  "mediumblue","gray20","lightcoral",
                                  "firebrick","purple3","darkorange","darkseagreen1",
                                  "pink","sienna1","turquoise",
                                  "tan4","steelblue","thistle","aquamarine4","sandybrown")
                                  
    names(phyloColorTbl) = Top20
    
    # Barplot of relative abundances of top 20 species, ordered by Overal_score
    pdf(file="Barplot_Top20species_Overal_score.pdf", width=14)
    ggplot(mytotal1_long, aes(ID_afri)) +
      geom_bar(aes(weight=value, fill=variable)) +
      theme_bw() +
      theme(legend.position="bottom", legend.direction="horizontal", axis.text.x = element_blank()) +
      xlab("Individuals") + ylab("Phylotypes relative abundance") + labs(fill = "") +
      scale_fill_manual(values = phyloColorTbl) +
      geom_point(aes(y = Overal_score/3)) + 
      scale_y_continuous(sec.axis = sec_axis(~.*3, name="Overall score"))
    dev.off()  
  
    # b) at the genus level----------
    # order samples by score values (from smaller to bigger and then plot the top genera to see if there is a trend)
    Top20gen = colnames(mytotal.gen[,which(colnames(mytotal.gen) == "Prevotella_9"):(which(colnames(mytotal.gen) == "Prevotella_9") + 19)]) 
    Top20gen
    
    # Linear regression models with Top20gen and development outcomes, adjusted on run
    library(broom)
    outcomes <- c("Comm_tot", "PS_tot", "PES_tot", "GM_tot", "FM_tot", "Overal_score")
    results_list <- list()
    number <- 1
    
    for (y in outcomes) {
      for (x in Top20gen) {
        # Build formula safely using backticks
        f <- as.formula(sprintf("%s ~ `%s` + run", y, x))
        # Fit model
        model <- tryCatch(lm(f, data = mytotal.gen),
                          error = function(e) { 
                            warning(sprintf("Model failed for outcome=%s predictor=%s: %s", y, x, e$message))
                            return(NULL)
                          })
        if (is.null(model)) {
          results_list[[number]] <- data.frame(outcome = y, predictor = x,
                                                beta = NA_real_, p_value = NA_real_, r_squared = NA_real_,
                                                stringsAsFactors = FALSE)
          counter <- counter + 1
          next
        }
        # Extract coefficient row for the predictor
        coef_table <- tidy(model)
        # Use cleaned term names (remove backticks) to match
        coef_table$term_clean <- gsub("`", "", coef_table$term)
        predictor_row <- coef_table[coef_table$term_clean == x, , drop = FALSE]
        
        if (nrow(predictor_row) == 0) {
          beta <- NA_real_
          pval <- NA_real_
        } else {
          beta <- predictor_row$estimate[1]
          pval <- predictor_row$p.value[1]
        }
        # Extract model-level statistics
        model_stats <- glance(model)
        rsq <- model_stats$r.squared
        # Store results
        results_list[[number]] <- data.frame(outcome = y, predictor = x,
                                              beta = beta, p_value = pval, r_squared = rsq,
                                              stringsAsFactors = FALSE)
        number <- number + 1
      }
    }
    # Combine all results into a single data frame
    results <- do.call(rbind, results_list)
    write.csv(results, file = "Results_lm_Top20gen_scores_multivariate.csv")
    rm(results, results_list, outcomes, number, model_stats, model, f, coef_table, x, y, beta, pval, rsq, predictor_row)
    
    
    # dot plots + regression lines (univariate)
  
    mytotal.gen1 = as.data.frame(cbind(mytotal.gen[,c("ID_afri","Overal_score")], reads.gen_ra[which(colnames(reads.gen_ra) %in% Top20gen)]))
    mytotal.gen1 = mytotal.gen1[order(mytotal.gen1$Overal_score),]
    mytotal.gen1_long <- reshape2::melt(mytotal.gen1, id.vars = c("ID_afri","Overal_score"))
    mytotal.gen1_long$ID_afri <- factor(mytotal.gen1_long$ID_afri, levels = unique(mytotal.gen1_long$ID_afri))
    
    # Barplot of relative abundances of top 20 genera, ordered by Overal_score
    names(phyloColorTbl) = Top20gen
    
    pdf(file="Barplot_Top20genera_Overal_score.pdf", width=14)
    ggplot(mytotal.gen1_long, aes(ID_afri)) +
      geom_bar(aes(weight=value, fill=variable)) +
      theme_bw() +
      theme(legend.position="bottom", legend.direction="horizontal", axis.text.x = element_blank()) +
      xlab("Individuals") + ylab("Genera's relative abundance") + labs(fill = "") +
      scale_fill_manual(values = phyloColorTbl) +
      geom_point(aes(y = Overal_score/3)) + 
      scale_y_continuous(sec.axis = sec_axis(~.*3, name="Overall score"))
    dev.off()  
    
    
    # linear regression between top20 genera RA and each domain and overall
    # Overal_score
    p1 = ggplot( mytotal.gen1_long, aes(y = Overal_score, x = value)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
      facet_wrap(.~ variable) +
      theme_minimal() # Subdoligranulum (Oscillospiraceae) is significant
    
    # Comm_tot
    mytotal.gen2 = as.data.frame(cbind(mytotal.gen[,c("ID_afri","Comm_tot")], reads.gen_ra[which(colnames(reads.gen_ra) %in% Top20gen)]))
    mytotal.gen2_long <- reshape2::melt(mytotal.gen2, id.vars = c("ID_afri","Comm_tot"))
    mytotal.gen2_long$ID_afri <- factor(mytotal.gen2_long$ID_afri, levels = unique(mytotal.gen2_long$ID_afri))
    
    p2 = ggplot( mytotal.gen2_long, aes(y = Comm_tot, x = value)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
      facet_wrap(.~ variable) +
      theme_minimal() # Comm_tot increases when Bifidobacterium and Bacteroides increase
    
    # PS_tot
    mytotal.gen3 = as.data.frame(cbind(mytotal.gen[,c("ID_afri","PS_tot")], reads.gen_ra[which(colnames(reads.gen_ra) %in% Top20gen)]))
    mytotal.gen3_long <- reshape2::melt(mytotal.gen3, id.vars = c("ID_afri","PS_tot"))
    mytotal.gen3_long$ID_afri <- factor(mytotal.gen3_long$ID_afri, levels = unique(mytotal.gen3_long$ID_afri))
    
    p3 = ggplot( mytotal.gen3_long, aes(y = PS_tot, x = value)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
      facet_wrap(.~ variable) +
      theme_minimal() # PS_tot significantly increases when Faecalibacterium decreases, PS_tot significantly increases when Bifidobacterium and Subdoligranulum increase
    
    # PES_tot
    mytotal.gen4 = as.data.frame(cbind(mytotal.gen[,c("ID_afri","PES_tot")], reads.gen_ra[which(colnames(reads.gen_ra) %in% Top20gen)]))
    mytotal.gen4_long <- reshape2::melt(mytotal.gen4, id.vars = c("ID_afri","PES_tot"))
    mytotal.gen4_long$ID_afri <- factor(mytotal.gen4_long$ID_afri, levels = unique(mytotal.gen4_long$ID_afri))
    
    p4 = ggplot( mytotal.gen4_long, aes(y = PES_tot, x = value)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
      facet_wrap(.~ variable) +
      theme_minimal() # PES_tot significantly increases when Faecalibacterium decreases
    
    # FM_tot
    mytotal.gen5 = as.data.frame(cbind(mytotal.gen[,c("ID_afri","FM_tot")], reads.gen_ra[which(colnames(reads.gen_ra) %in% Top20gen)]))
    mytotal.gen5_long <- reshape2::melt(mytotal.gen5, id.vars = c("ID_afri","FM_tot"))
    mytotal.gen5_long$ID_afri <- factor(mytotal.gen5_long$ID_afri, levels = unique(mytotal.gen5_long$ID_afri))
    
    p5 = ggplot( mytotal.gen5_long, aes(y = FM_tot, x = value)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
      facet_wrap(.~ variable) +
      theme_minimal() # FM_tot significantly increases when Subdoligranulum increases
    
    # GM_tot
    mytotal.gen6 = as.data.frame(cbind(mytotal.gen[,c("ID_afri","GM_tot")], reads.gen_ra[which(colnames(reads.gen_ra) %in% Top20gen)]))
    mytotal.gen6_long <- reshape2::melt(mytotal.gen6, id.vars = c("ID_afri","GM_tot"))
    mytotal.gen6_long$ID_afri <- factor(mytotal.gen6_long$ID_afri, levels = unique(mytotal.gen6_long$ID_afri))
    
    p6 = ggplot( mytotal.gen6_long, aes(y = GM_tot, x = value)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
      facet_wrap(.~ variable) +
      theme_minimal() # GM_tot significantly increases when Escherichia−Shigella, Veillonella and Streptococcus decrease, and when Ruminococcaceae_UCG−002 increases
    
    pdf(file="lm_Top20genera_per_domain.pdf", height = 10, width =14 )
    p1
    p2
    p3
    p4
    p5
    p6
    dev.off()
    
    rm(p1, p2, p3, p4, p5, p6,mytotal.gen1,mytotal.gen1_long,mytotal.gen2,mytotal.gen2_long,mytotal.gen3,mytotal.gen3_long,mytotal.gen4,mytotal.gen4_long,mytotal.gen5,mytotal.gen5_long,mytotal.gen6,mytotal.gen6_long)
  
  
    # c) at the family level----------
    # order samples by score values (from smaller to bigger and then plot the top families to see if there is a trend)
    Top20fam = colnames(mytotal.fam[,which(colnames(mytotal.fam) == "Prevotellaceae"):(which(colnames(mytotal.fam) == "Prevotellaceae") + 19)]) 
    Top20fam
    
    # Linear regression models with Top20gen and development outcomes, adjusted on run 
    library(broom)
    outcomes <- c("Comm_tot", "PS_tot", "PES_tot", "GM_tot", "FM_tot", "Overal_score")
    results_list <- list()
    number <- 1
    
    for (y in outcomes) {
      for (x in Top20fam) {
        # Build formula safely using backticks
        f <- as.formula(sprintf("%s ~ `%s` + run", y, x))
        # Fit model
        model <- tryCatch(lm(f, data = mytotal.fam),
                          error = function(e) { 
                            warning(sprintf("Model failed for outcome=%s predictor=%s: %s", y, x, e$message))
                            return(NULL)
                          })
        if (is.null(model)) {
          results_list[[number]] <- data.frame(outcome = y, predictor = x,
                                               beta = NA_real_, p_value = NA_real_, r_squared = NA_real_,
                                               stringsAsFactors = FALSE)
          counter <- counter + 1
          next
        }
        # Extract coefficient row for the predictor
        coef_table <- tidy(model)
        # Use cleaned term names (remove backticks) to match
        coef_table$term_clean <- gsub("`", "", coef_table$term)
        predictor_row <- coef_table[coef_table$term_clean == x, , drop = FALSE]
        
        if (nrow(predictor_row) == 0) {
          beta <- NA_real_
          pval <- NA_real_
        } else {
          beta <- predictor_row$estimate[1]
          pval <- predictor_row$p.value[1]
        }
        # Extract model-level statistics
        model_stats <- glance(model)
        rsq <- model_stats$r.squared
        # Store results
        results_list[[number]] <- data.frame(outcome = y, predictor = x,
                                             beta = beta, p_value = pval, r_squared = rsq,
                                             stringsAsFactors = FALSE)
        number <- number + 1
      }
    }
    # Combine all results into a single data frame
    results <- do.call(rbind, results_list)
    write.csv(results, file = "Results_lm_Top20fam_scores_multivariate.csv")
    rm(results, results_list, outcomes, number, model_stats, model, f, coef_table, x, y, beta, pval, rsq, predictor_row)
    
    
    # dot plots + regression lines (univariate) 
    mytotal.fam1 = as.data.frame(cbind(mytotal.fam[,c("ID_afri","Overal_score")], reads.fam_ra[which(colnames(reads.fam_ra) %in% Top20fam)]))
    mytotal.fam1 = mytotal.fam1[order(mytotal.fam1$Overal_score),]
    mytotal.fam1_long <- reshape2::melt(mytotal.fam1, id.vars = c("ID_afri","Overal_score"))
    mytotal.fam1_long$ID_afri <- factor(mytotal.fam1_long$ID_afri, levels = unique(mytotal.fam1_long$ID_afri))
    
    # Barplot of relative abundances of top 20 families, ordered by Overal_score
    names(phyloColorTbl) = Top20fam
    
    pdf(file="Barplot_Top20families_Overal_score.pdf", width=14)
    ggplot(mytotal.fam1_long, aes(ID_afri)) +
      geom_bar(aes(weight=value, fill=variable)) +
      theme_bw() +
      theme(legend.position="bottom", legend.direction="horizontal", axis.text.x = element_blank()) +
      xlab("Individuals") + ylab("Families' relative abundance") + labs(fill = "") +
      scale_fill_manual(values = phyloColorTbl) +
      geom_point(aes(y = Overal_score/3)) + 
      scale_y_continuous(sec.axis = sec_axis(~.*3, name="Overall score"))
    dev.off()  
    
    
    
    # linear regression between top20 families RA and each domain and overall
          # Overal_score
          p1 = ggplot( mytotal.fam1_long, aes(y = Overal_score, x = value)) +
            geom_point() +
            geom_smooth(method='lm') +
            stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
            facet_wrap(.~ variable) +
            theme_minimal() # Verrucomicrobiaceae is significant
        
          # Comm_tot
          mytotal.fam2 = as.data.frame(cbind(mytotal.fam[,c("ID_afri","Comm_tot")], reads.fam_ra[which(colnames(reads.fam_ra) %in% Top20fam)]))
          mytotal.fam2_long <- reshape2::melt(mytotal.fam2, id.vars = c("ID_afri","Comm_tot"))
          mytotal.fam2_long$ID_afri <- factor(mytotal.fam2_long$ID_afri, levels = unique(mytotal.fam2_long$ID_afri))
          
          p2 = ggplot( mytotal.fam2_long, aes(y = Comm_tot, x = value)) +
            geom_point() +
            geom_smooth(method='lm') +
            stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
            facet_wrap(.~ variable) +
            theme_minimal() # none is significant
          
          # PS_tot
          mytotal.fam3 = as.data.frame(cbind(mytotal.fam[,c("ID_afri","PS_tot")], reads.fam_ra[which(colnames(reads.fam_ra) %in% Top20fam)]))
          mytotal.fam3_long <- reshape2::melt(mytotal.fam3, id.vars = c("ID_afri","PS_tot"))
          mytotal.fam3_long$ID_afri <- factor(mytotal.fam3_long$ID_afri, levels = unique(mytotal.fam3_long$ID_afri))
          
          p3 = ggplot( mytotal.fam3_long, aes(y = PS_tot, x = value)) +
            geom_point() +
            geom_smooth(method='lm') +
            stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
            facet_wrap(.~ variable) +
            theme_minimal() # PS_tot significantly increases when Rikenellaceae increases, PS_tot significantly increases when Verrucomicrobiaceae increases
     
          # PES_tot
          mytotal.fam4 = as.data.frame(cbind(mytotal.fam[,c("ID_afri","PES_tot")], reads.fam_ra[which(colnames(reads.fam_ra) %in% Top20fam)]))
          mytotal.fam4_long <- reshape2::melt(mytotal.fam4, id.vars = c("ID_afri","PES_tot"))
          mytotal.fam4_long$ID_afri <- factor(mytotal.fam4_long$ID_afri, levels = unique(mytotal.fam4_long$ID_afri))
          
          p4 = ggplot( mytotal.fam4_long, aes(y = PES_tot, x = value)) +
            geom_point() +
            geom_smooth(method='lm') +
            stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
            facet_wrap(.~ variable) +
            theme_minimal() # PES_tot significantly decreases when Streptococacceae increases  
          
          ggplot( mytotal.fam4, aes(y = PES_tot, x = Streptococcaceae)) +
            geom_point() +
            geom_smooth(method='lm') +
            stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
            theme_minimal() # PES_tot significantly decreases when Streptococacceae increases  
          
          # FM_tot
          mytotal.fam5 = as.data.frame(cbind(mytotal.fam[,c("ID_afri","FM_tot")], reads.fam_ra[which(colnames(reads.fam_ra) %in% Top20fam)]))
          mytotal.fam5_long <- reshape2::melt(mytotal.fam5, id.vars = c("ID_afri","FM_tot"))
          mytotal.fam5_long$ID_afri <- factor(mytotal.fam5_long$ID_afri, levels = unique(mytotal.fam5_long$ID_afri))
          
          p5 = ggplot( mytotal.fam5_long, aes(y = FM_tot, x = value)) +
            geom_point() +
            geom_smooth(method='lm') +
            stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
            facet_wrap(.~ variable) +
            theme_minimal() # FM_tot significantly increases when Bacteroidaceae increases
          
          # GM_tot
          mytotal.fam6 = as.data.frame(cbind(mytotal.fam[,c("ID_afri","GM_tot")], reads.fam_ra[which(colnames(reads.fam_ra) %in% Top20fam)]))
          mytotal.fam6_long <- reshape2::melt(mytotal.fam6, id.vars = c("ID_afri","GM_tot"))
          mytotal.fam6_long$ID_afri <- factor(mytotal.fam6_long$ID_afri, levels = unique(mytotal.fam6_long$ID_afri))
          
          p6 = ggplot( mytotal.fam6_long, aes(y = GM_tot, x = value)) +
            geom_point() +
            geom_smooth(method='lm') +
            stat_cor(method = "spearman", label.x.npc = "middle", label.y.npc = "bottom", size=2.5) +
            facet_wrap(.~ variable) +
            theme_minimal() # GM_tot significantly increases when Bacteroidales_S24-7_group increases
          
          pdf(file="lm_Top20families_per_domain.pdf", height = 10, width =14 )
          p1
          p2
          p3
          p4
          p5
          p6
          dev.off()
          
    rm(p1, p2, p3, p4, p5, p6,mytotal1, mytotal1_long,mytotal.fam1,mytotal.fam1_long,mytotal2, mytotal2_long,mytotal.fam2,mytotal.fam2_long,mytotal3, mytotal3_long,mytotal.fam3,mytotal.fam3_long,mytotal4, mytotal4_long,mytotal.fam4,mytotal.fam4_long,mytotal5, mytotal5_long,mytotal.fam5,mytotal.fam5_long,mytotal6, mytotal6_long,mytotal.fam6,mytotal.fam6_long)
    
  
  # 6) Correlation plots --------
  library(corrplot)
  
    # a) at the highest available rank----------
    CORR1 = as.data.frame(cbind(mytotal[,c("Overal_score","Comm_tot","PS_tot","PES_tot","FM_tot","GM_tot")], reads_ra[which(colnames(reads_ra) %in% Top20)]))
    CORR1 = as.matrix(CORR1)
    # calculate p-values
    testRes = cor.mtest(CORR1, conf.level = 0.95)
    
    pdf("Correlation_Neurodevelopment_Top20species.pdf")
    corrplot.mixed(cor(CORR1, use="pairwise.complete.obs") # work only on complete cases
                   ,lower = "number"                      # lower part of plot has correlation values
                   ,upper = "circle"                      # upper part of plot has circles
                   ,tl.col = "black"                      # color of text labels
                   ,tl.pos = "lt"                         # position of text labels (tl = top left)
                   ,tl.cex=0.5                            # size of text labels
                   #,order = "hclust"                     # clustering of variables
                   #,hclust.method = "ward.D"             # method for clustering
                   ,number.cex = 0.4                      # size of correlation values
                   ,p.mat = testRes$p                     # p-value matrix
                   ,insig = "blank")                      # erase corr values and circles that are not significant
    dev.off()
    
    # Correlation between all bacteria and each domain
    CORR2 = as.data.frame(cbind(mytotal[,c("Overal_score","Comm_tot","PS_tot","PES_tot","FM_tot","GM_tot")], reads_ra[which(colnames(reads_ra) %in% bacteria)]))
    CORR2res = cor(CORR2, use="pairwise.complete.obs")
    
    res <- outer(CORR2, CORR2, Vectorize(\(x, y) cor.test(x, y)$p.value)) |>
      as.table() |> as.data.frame() |> setNames(c("Var1", "Var2", "p")) |>
      {\(.) cbind(., p.adj=p.adjust(.$p, method='BH'))}()
    head(res)
    res.p = pivot_wider(res, id_cols = "Var1",names_from = "Var2",values_from = "p.adj")
    library(tidyverse)
    res.p %<>% 
      column_to_rownames(var="Var1")
    
    rownames(res.p)[which(res.p$Overal_score <= 0.05)] # Actinomyces_odontolyticus
    rownames(res.p)[which(res.p$PS_tot <= 0.05)] # Actinomyces_odontolyticus
    rownames(res.p)[which(res.p$Comm_tot <= 0.05)] # f_uncultured and g_Acidaminococcus
    rownames(res.p)[which(res.p$PES_tot <= 0.05)] # nothing
    rownames(res.p)[which(res.p$FM_tot <= 0.05)] # nothing
    rownames(res.p)[which(res.p$GM_tot <= 0.05)] # Streptococcus_parasanguinis and Prevotella_histicola_JCM_15637__DNF00424
    
    
    ggplot(mytotal, aes(y = Overal_score, x = Actinomyces_odontolyticus)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(label.x.npc = "left", label.y.npc = "bottom", size=2.5) +
      theme_minimal()
    
    ggplot(mytotal, aes(y = PS_tot, x = Actinomyces_odontolyticus)) +
      geom_point(color = mytotal$haz) +
      geom_smooth(method='lm') +
      stat_cor(label.x.npc = "left", label.y.npc = "bottom", size=2.5) +
      theme_minimal()
    
    ggplot(mytotal, aes(y = GM_tot, x = Streptococcus_parasanguinis)) +
      geom_point( color = mytotal$haz) +
      geom_smooth(method='lm') +
      stat_cor(label.x.npc = "left", label.y.npc = "bottom", size=2.5) +
      theme_minimal()
    
    
    rm(CORR1,testRes,CORR2,CORR2res,res,res.p)
  
  
    # b) at the genus level----------
    
    # Correlation plot between Top20 and each domain
    library(corrplot)
    CORR1 = as.data.frame(cbind(mytotal.gen[,c("Overal_score","Comm_tot","PS_tot","PES_tot","FM_tot","GM_tot")], reads.gen_ra[which(colnames(reads.gen_ra) %in% Top20gen)]))
    CORR1 = as.matrix(CORR1)
    # calculate p-values
    testRes = cor.mtest(CORR1, conf.level = 0.95)
    
    pdf("Correlation_Neurodevelopment_Top20genera.pdf")
    corrplot.mixed(cor(CORR1, use="pairwise.complete.obs") # work only on complete cases
                   ,lower = "number"                      # lower part of plot has correlation values
                   ,upper = "circle"                      # upper part of plot has circles
                   ,tl.col = "black"                      # color of text labels
                     ,tl.pos = "lt"                         # position of text labels (tl = top left)
                   ,tl.cex=0.5                            # size of text labels
                   #,order = "hclust"                     # clustering of variables
                   #,hclust.method = "ward.D"             # method for clustering
                   ,number.cex = 0.4                      # size of correlation values
                   ,p.mat = testRes$p                     # p-value matrix
                   ,insig = "blank")                      # erase corr values and circles that are not significant
    dev.off()
    
    # Correlation between all bacteria and each domain
    CORR2 = as.data.frame(cbind(mytotal.gen[,c("Overal_score","Comm_tot","PS_tot","PES_tot","FM_tot","GM_tot")], reads.gen_ra[which(colnames(reads.gen_ra) %in% genera)]))
    CORR2res = cor(CORR2, use="pairwise.complete.obs")
    
    res <- outer(CORR2, CORR2, Vectorize(\(x, y) cor.test(x, y)$p.value)) |>
      as.table() |> as.data.frame() |> setNames(c("Var1", "Var2", "p")) |>
      {\(.) cbind(., p.adj=p.adjust(.$p, method='BH'))}()
    head(res)
    res.p = pivot_wider(res, id_cols = "Var1",names_from = "Var2",values_from = "p.adj")
    library(tidyverse)
    res.p %<>% 
      column_to_rownames(var="Var1")
    
    rownames(res.p)[which(res.p$Overal_score <= 0.05)] # nothing
    rownames(res.p)[which(res.p$PS_tot <= 0.05)] # nothing
    rownames(res.p)[which(res.p$Comm_tot <= 0.05)] # nothing
    rownames(res.p)[which(res.p$PES_tot <= 0.05)] # nothing
    rownames(res.p)[which(res.p$FM_tot <= 0.05)] # nothing
    rownames(res.p)[which(res.p$GM_tot <= 0.05)] # Streptococcus and Anaerotruncus
    
    ggplot(mytotal.gen, aes(y = PES_tot, x = Streptococcus)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(label.x.npc = "left", label.y.npc = "bottom", size=2.5) +
      theme_minimal() # the more Streptococcus, the lower the PES score
    
    ggplot(mytotal.gen, aes(y = GM_tot, x = Streptococcus)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(label.x.npc = "left", label.y.npc = "bottom", size=2.5) +
      theme_minimal() # the more Streptococcus, the lower the GM score
    
    # keep only Top20 in CORR2
    CORR2 = CORR2[,1:26]
    res.p = res.p[1:26,1:26]
    
    pdf("Correlation_Neurodevelopment_Top20genera_q-value.pdf")
    corrplot.mixed(cor(CORR2, use="pairwise.complete.obs") # work only on complete cases
                   ,lower = "number"                      # lower part of plot has correlation values
                   ,upper = "circle"                      # upper part of plot has circles
                   ,tl.col = "black"                      # color of text labels
                     ,tl.pos = "lt"                         # position of text labels (tl = top left)
                   ,tl.cex=0.5                            # size of text labels
                   #,order = "hclust"                     # clustering of variables
                   #,hclust.method = "ward.D"             # method for clustering
                   ,number.cex = 0.4                      # size of correlation values
                   ,p.mat = as.matrix(res.p)      # p-value matrix
                   ,insig = "blank")                      # erase corr values and circles that are not significant
    dev.off()
    
    rm(CORR1,testRes,CORR2,CORR2res,res,res.p)
    
  
    # c) at the family level----------
    
    # Correlation plot between Top20 and each domain
    library(corrplot)
    CORR1 = as.data.frame(cbind(mytotal.fam[,c("Overal_score","Comm_tot","PS_tot","PES_tot","FM_tot","GM_tot")], reads.fam_ra[which(colnames(reads.fam_ra) %in% Top20fam)]))
    CORR1 = as.matrix(CORR1)
    # calculate p-values
    testRes = cor.mtest(CORR1, conf.level = 0.95)
    
    pdf("Correlation_Neurodevelopment_Top20families.pdf")
    corrplot.mixed(cor(CORR1, use="pairwise.complete.obs") # work only on complete cases
                   ,lower = "number"                      # lower part of plot has correlation values
                   ,upper = "circle"                      # upper part of plot has circles
                   ,tl.col = "black"                      # color of text labels
                   ,tl.pos = "lt"                         # position of text labels (tl = top left)
                   ,tl.cex=0.8                            # size of text labels
                   #,order = "hclust"                     # clustering of variables
                   #,hclust.method = "ward.D"             # method for clustering
                   ,number.cex = 0.5                      # size of correlation values
                   ,p.mat = testRes$p                     # p-value matrix
                   ,insig = "label_sig"                   # blank: erase corr values and circles that are not significant
                   )                      
    dev.off()
    
    pdf("correlation_plot.pdf", width = 12, height = 12)
    corrplot.mixed(
      cor(CORR1, use = "pairwise.complete.obs"),
      lower = "number",            # numbers in lower triangle
      upper = "circle",            # circles in upper triangle
      tl.col = "black",            # text labels in black
      tl.pos = "lt",               
      tl.cex = 0.8,
      number.cex = 0.7,            # slightly bigger numbers
      p.mat = testRes$p,
      insig = "label_sig",
      sig.level = c(0.05),   # significance cutoffs
      pch.cex = 0.8,               # shrink significance stars
      pch.col = "black"              # make them smaller + red so they don't hide circles
    )
    dev.off()
    
    # Correlation between all bacteria and each domain
    CORR2 = as.data.frame(cbind(mytotal.fam[,c("Overal_score","Comm_tot","PS_tot","PES_tot","FM_tot","GM_tot")], reads.fam_ra[which(colnames(reads.fam_ra) %in% families)]))
    CORR2res = cor(CORR2, use="pairwise.complete.obs")
    
    res <- outer(CORR2, CORR2, Vectorize(\(x, y) cor.test(x, y)$p.value)) |>
      as.table() |> as.data.frame() |> setNames(c("Var1", "Var2", "p")) |>
      {\(.) cbind(., p.adj=p.adjust(.$p, method='BH'))}()
    head(res)
    res.p = pivot_wider(res, id_cols = "Var1",names_from = "Var2",values_from = "p.adj")
    library(tidyverse)
    res.p %<>% 
      column_to_rownames(var="Var1")
    
    rownames(res.p)[which(res.p$Overal_score <= 0.05)] # nothing
    rownames(res.p)[which(res.p$PS_tot <= 0.05)] # uncultured
    rownames(res.p)[which(res.p$Comm_tot <= 0.05)] # only "uncultured"
    rownames(res.p)[which(res.p$PES_tot <= 0.05)] # only "Streptococcaceae"
    rownames(res.p)[which(res.p$FM_tot <= 0.05)] # nothing
    rownames(res.p)[which(res.p$GM_tot <= 0.05)] # nothing
    
    ggplot(mytotal.fam, aes(y = PES_tot, x = Streptococcaceae)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(label.x.npc = "left", label.y.npc = "bottom", size=2.5) +
      theme_minimal() # the more Streptococcaceae, the lower the PES score
    
    ggplot(mytotal.fam, aes(y = PES_tot, x = Pasteurellaceae)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(label.x.npc = "left", label.y.npc = "bottom", size=2.5) +
      theme_minimal() # the more Pasteurellaceae, the lower the PES score
    
    # keep only Top20 in CORR2
    CORR2 = CORR2[,1:26]
    res.p = res.p[1:26,1:26]
    
    pdf("Correlation_Neurodevelopment_Top20families_q-value.pdf")
    corrplot.mixed(cor(CORR2, use="pairwise.complete.obs") # work only on complete cases
                   ,lower = "number"                      # lower part of plot has correlation values
                   ,upper = "circle"                      # upper part of plot has circles
                   ,tl.col = "black"                      # color of text labels
                   ,tl.pos = "lt"                         # position of text labels (tl = top left)
                   ,tl.cex=0.5                            # size of text labels
                   #,order = "hclust"                     # clustering of variables
                   #,hclust.method = "ward.D"             # method for clustering
                   ,number.cex = 0.5                      # size of correlation values
                   ,p.mat = as.matrix(res.p)      # p-value matrix
                   ,insig = "blank")                      # erase corr values and circles that are not significant
    dev.off()
    
    rm(CORR1,testRes,CORR2,CORR2res,res,res.p)
  

    
  # 7) DESeq2 -----------
  library("DESeq2")
    
    # a) at the highest available rank -------------
    
    # keep only two extreme quartiles of the distribution in each domain
    quantile(mytotal$Comm_tot)
    mytotalComm = mytotal[which(mytotal$Comm_tot >= 60 | mytotal$Comm_tot <= 45),]
    mytotalComm$delay = ifelse(mytotalComm$Comm_tot>=60,"No",ifelse(mytotalComm$Comm_tot <= 45,"Yes","grey"))
    table(mytotalComm$delay,useNA="always")
    
    quantile(mytotal$PS_tot)
    mytotalPS = mytotal[which(mytotal$PS_tot >= 50 | mytotal$PS_tot <= 30),]
    mytotalPS$delay = ifelse(mytotalPS$PS_tot>=50,"No",ifelse(mytotalPS$PS_tot <= 30,"Yes","grey"))
    table(mytotalPS$delay,useNA="always")
    
    quantile(mytotal$PES_tot)
    mytotalPES = mytotal[which(mytotal$PES_tot >= 60 | mytotal$PES_tot <= 45),]
    mytotalPES$delay = ifelse(mytotalPES$PES_tot>=60,"No",ifelse(mytotalPES$PES_tot <= 45,"Yes","grey"))
    table(mytotalPES$delay,useNA="always")
    
    quantile(mytotal$FM_tot)
    mytotalFM = mytotal[which(mytotal$FM_tot >= 55 | mytotal$FM_tot <= 42),]
    mytotalFM$delay = ifelse(mytotalFM$FM_tot>=55,"No",ifelse(mytotalFM$FM_tot <= 42,"Yes","grey"))
    table(mytotalFM$delay,useNA="always")
    
    quantile(mytotal$GM_tot)
    mytotalGM = mytotal[which(mytotal$GM_tot >= 60 | mytotal$GM_tot <= 50),]
    mytotalGM$delay = ifelse(mytotalGM$GM_tot>=60,"No",ifelse(mytotalGM$GM_tot <= 50,"Yes","grey"))
    table(mytotalGM$delay,useNA="always")
      
    quantile(mytotal$Overal_score)
    mytotalOverall = mytotal[which(mytotal$Overal_score >= 265 | mytotal$Overal_score <= 225),]
    mytotalOverall$delay = ifelse(mytotalOverall$Overal_score>=265,"No",ifelse(mytotalOverall$Overal_score <= 225,"Yes","grey"))
    table(mytotalOverall$delay,useNA="always")
    
        # from mytotal, generate a count table with taxa in rows (cts), and a metadata table (coldata)
        colnames(mytotalOverall)
        
        cts =  t(mytotalOverall[,which(colnames(mytotalOverall) %in% bacteria)	])
        colnames(cts) = mytotalOverall$ID_afri
        coldata = mytotalOverall[,c("ID_afri","haz_cont","haz","delay","run")]
  
        rownames(coldata) = mytotalOverall$ID_afri
        
        all(rownames(coldata) %in% colnames(cts))
        all(rownames(coldata) == colnames(cts))   # if not : cts <- cts[, rownames(coldata)]
        
        dds <- DESeqDataSetFromMatrix(countData = cts,
                                      colData = coldata,
                                      design = ~ run + delay )   # design variable needs to be a factor
        # if we add a variable to adjust on, put the treatment condition at the end of formula
        dds
        
        all(apply(apply(cts, 2, function(x) x==0), 1, any)) # all columns have at least one zero --> DESeq2 cannot compute a log geometric means for estimating size factors
        
        # pre-filtering on taxa that have counts>10, to reduce memory size of dds object and increase speed
        # keep <- rowSums(counts(dds)) >= 10
        # dds <- dds[keep,] # stayed at 348
        # 
        # # change reference level, in case
        # dds$delay <- relevel(dds$delay, ref = "0")
        # 
        # # generate results
        # dds <- DESeq(dds) # OR the 3 steps, if impossible to get the geometric mean
        geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
        ddsLove<-estimateSizeFactors(dds, geoMeans=geoMeans)
        # or ddsLove = estimateSizeFactors(dds, type="iterate")
        dds <- estimateDispersions(ddsLove)
        dds <- nbinomWaldTest(dds)
        
        res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE) # or : res <- results(dds, contrast=c("haz_bin","normal","stunted"))
        
        # Log fold change shrinkage for visualization and ranking
        resultsNames(dds)
        resLFC <- lfcShrink(dds, coef="delay_Yes_vs_No", type="normal")
        resLFC
        resOrdered <- res[order(res$pvalue),]
        summary(res)
        sum(res$padj < 0.05, na.rm=TRUE) # how many species with adjusted pvalues below 0.05
        sum(res$pvalue < 0.05, na.rm=TRUE) # how many species with pvalues below 0.05
        
        # MA-plot
        plotMA(res, ylim=c(-25,10)) # species with padj < 0.1 in blue dots
        #idx <- identify(res$baseMean, res$log2FoldChange)
        #rownames(res)[idx]
        # terminate this command using esc/echap
        
        
        write.csv(as.data.frame(resOrdered), file="DESeq2_results_Overall_delay.csv")
        
        results = as.data.frame(resOrdered)
        results$species = rownames(results)
        results$species = as.factor(results$species)
        results$FC = 2^results$log2FoldChange
        results$hci = results$log2FoldChange + 1.96*results$lfcSE
        results$lci = results$log2FoldChange - 1.96*results$lfcSE
        results$species[which(results$padj <0.05)]
  
        #negative values of L2FC means lower expression in delay vs no_delay
        
        pdf("afribiota_DESEQ2_meancounts_Overall.pdf", width=12, height=8)
        ggplot() +
          geom_col(data=results[1:50,], aes(x=species, y=baseMean)) +
          theme_bw() +
          theme(axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=.5)) +
          scale_x_discrete(limits=results[1:50,]$species) +  # to order x-axis as in dataframe
          labs(title="Mean read counts", x = "Taxa", y = "baseMean")
        dev.off()
        
        pdf("afribiota_DESEQ2_log2FC_Overall.pdf", width=12, height=10)
        ggplot(results[1:50,], aes(x=species, y=log2FoldChange)) +
          geom_point(stat="identity", aes(colour=cut(padj,c(-Inf,0.05,Inf)) )) +
          geom_errorbar(aes(ymin=lci, ymax=hci, colour=cut(padj,c(-Inf,0.05,Inf)) ), width=.1) +
          geom_hline(yintercept = 0) +
          labs(title="Log2 Fold Change delayed vs not delayed", x = "Taxa", y = "log2 Fold Change") +
          theme_bw() +						# remove grey background
          theme(legend.position="bottom", axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=.5)) +
          scale_x_discrete(limits=results[1:50,]$species) + # to order x-axis as in dataframe
          scale_color_manual(name = "Adjusted P-value",
                             values = c("(-Inf,0.05]" = "red", "(0.05, Inf]" = "black"),
                             labels = c("<= 0.05", "> 0.05"),
                             drop = FALSE )
        dev.off()
        
        # plot counts for a specific species
        plotCounts(dds, gene=which.min(res$padj), intgroup="delay")
        plotCounts(dds, gene="Clostridiales_bacterium_VE202_07", intgroup="delay") # or d <- plotCounts(dds, gene=which.min(res$padj), intgroup="ROGRP", returnData=TRUE)
        mytotalOverall$Clostridiales_bacterium_VE202_07
        plotCounts(dds, gene="Paeniclostridium_sordellii_ATCC_9714", intgroup="delay")
        mytotalOverall$Paeniclostridium_sordellii_ATCC_9714
        plotCounts(dds, gene="g_Christensenellaceae_R_7_group", intgroup="delay")
        mytotalOverall$g_Christensenellaceae_R_7_group
  
        rm(mytotalComm,mytotalPS,mytotalPES,mytotalFM,mytotalGM,mytotalOverall)
        rm(res,resLFC,resOrdered,results,dds,ddsLove,cts,coldata,geoMeans)
      
      # b) at the genus level-----------
      # keep only two extreme quartiles of the distribution in each domain
      quantile(mytotal.gen$Comm_tot)
      mytotalComm = mytotal.gen[which(mytotal.gen$Comm_tot >= 60 | mytotal.gen$Comm_tot <= 45),]
      mytotalComm$delay = ifelse(mytotalComm$Comm_tot>=60,"No",ifelse(mytotalComm$Comm_tot <= 45,"Yes","grey"))
      table(mytotalComm$delay,useNA="always")
      
      quantile(mytotal.gen$PS_tot)
      mytotalPS = mytotal.gen[which(mytotal.gen$PS_tot >= 50 | mytotal.gen$PS_tot <= 30),]
      mytotalPS$delay = ifelse(mytotalPS$PS_tot>=50,"No",ifelse(mytotalPS$PS_tot <= 30,"Yes","grey"))
      table(mytotalPS$delay,useNA="always")
      
      quantile(mytotal.gen$PES_tot)
      mytotalPES = mytotal.gen[which(mytotal.gen$PES_tot >= 60 | mytotal.gen$PES_tot <= 45),]
      mytotalPES$delay = ifelse(mytotalPES$PES_tot>=60,"No",ifelse(mytotalPES$PES_tot <= 45,"Yes","grey"))
      table(mytotalPES$delay,useNA="always")
      
      quantile(mytotal.gen$FM_tot)
      mytotalFM = mytotal.gen[which(mytotal.gen$FM_tot >= 55 | mytotal.gen$FM_tot <= 42),]
      mytotalFM$delay = ifelse(mytotalFM$FM_tot>=55,"No",ifelse(mytotalFM$FM_tot <= 42,"Yes","grey"))
      table(mytotalFM$delay,useNA="always")
      
      quantile(mytotal.gen$GM_tot)
      mytotalGM = mytotal.gen[which(mytotal.gen$GM_tot >= 60 | mytotal.gen$GM_tot <= 50),]
      mytotalGM$delay = ifelse(mytotalGM$GM_tot>=60,"No",ifelse(mytotalGM$GM_tot <= 50,"Yes","grey"))
      table(mytotalGM$delay,useNA="always")
      
      quantile(mytotal.gen$Overal_score)
      mytotalOverall = mytotal.gen[which(mytotal.gen$Overal_score >= 265 | mytotal.gen$Overal_score <= 225),]
      mytotalOverall$delay = ifelse(mytotalOverall$Overal_score>=265,"No",ifelse(mytotalOverall$Overal_score <= 225,"Yes","grey"))
      table(mytotalOverall$delay,useNA="always")
      
      # from mytotal, generate a count table with taxa in rows (cts), and a metadata table (coldata)
      colnames(mytotalOverall)
      
      cts =  t(mytotalOverall[,which(colnames(mytotalOverall) %in% genera)	])
      colnames(cts) = mytotalOverall$ID_afri
      coldata = mytotalOverall[,c("ID_afri","haz_cont","haz","delay","run")]
      
      rownames(coldata) = mytotalOverall$ID_afri
      
      all(rownames(coldata) %in% colnames(cts))
      all(rownames(coldata) == colnames(cts))   # if not : cts <- cts[, rownames(coldata)]
      
      dds <- DESeqDataSetFromMatrix(countData = cts,
                                    colData = coldata,
                                    design = ~ run + delay)   # design variable needs to be a factor
      # if we add a variable to adjust on, put the treatment condition at the end of formula
      dds
      
      all(apply(apply(cts, 2, function(x) x==0), 1, any)) # all columns have at least one zero --> DESeq2 cannot compute a log geometric means for estimating size factors
      
      # pre-filtering on taxa that have counts>10, to reduce memory size of dds object and increase speed
      # keep <- rowSums(counts(dds)) >= 10
      # dds <- dds[keep,] # stayed at 348
      # 
      # # change reference level, in case
      # dds$delay <- relevel(dds$delay, ref = "0")
      # 
      # # generate results
      # dds <- DESeq(dds) # OR the 3 steps, if impossible to get the geometric mean
      geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
      ddsLove<-estimateSizeFactors(dds, geoMeans=geoMeans)
      # or ddsLove = estimateSizeFactors(dds, type="iterate")
      dds <- estimateDispersions(ddsLove)
      dds <- nbinomWaldTest(dds)
      
      res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE) # or : res <- results(dds, contrast=c("haz_bin","normal","stunted"))
      
      # Log fold change shrinkage for visualization and ranking
      resultsNames(dds)
      resLFC <- lfcShrink(dds, coef="delay_Yes_vs_No", type="normal")
      resLFC
      resOrdered <- res[order(res$pvalue),]
      summary(res)
      sum(res$padj < 0.05, na.rm=TRUE) # how many species with adjusted pvalues below 0.05
      sum(res$pvalue < 0.05, na.rm=TRUE) # how many species with pvalues below 0.05
      
      # MA-plot
      plotMA(res, ylim=c(-25,10)) # species with padj < 0.1 in blue dots
      #idx <- identify(res$baseMean, res$log2FoldChange)
      #rownames(res)[idx]
      # terminate this command using esc/echap
      
      
      write.csv(as.data.frame(resOrdered), file="DESeq2_results_Overall_delay_genera.csv")
      
      results = as.data.frame(resOrdered)
      results$species = rownames(results)
      results$species = as.factor(results$species)
      results$FC = 2^results$log2FoldChange
      results$hci = results$log2FoldChange + 1.96*results$lfcSE
      results$lci = results$log2FoldChange - 1.96*results$lfcSE
      results$species[which(results$padj <0.05)]
      
      #negative values of L2FC means lower expression in delay vs no_delay
      
      pdf("afribiota_DESEQ2_genera_meancounts_Overall.pdf", width=12, height=8)
      ggplot() +
        geom_col(data=results[1:50,], aes(x=species, y=baseMean)) +
        theme_bw() +
        theme(axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=.5)) +
        scale_x_discrete(limits=results[1:50,]$species) +  # to order x-axis as in dataframe
        labs(title="Mean read counts", x = "Taxa", y = "baseMean")
      dev.off()
      
      pdf("afribiota_DESEQ2_genera_log2FC_Overall.pdf", width=12, height=10)
      ggplot(results[1:50,], aes(x=species, y=log2FoldChange)) +
        geom_point(stat="identity", aes(colour=cut(padj,c(-Inf,0.05,Inf)) )) +
        geom_errorbar(aes(ymin=lci, ymax=hci, colour=cut(padj,c(-Inf,0.05,Inf)) ), width=.1) +
        geom_hline(yintercept = 0) +
        labs(title="Log2 Fold Change delayed vs not delayed", x = "Taxa", y = "log2 Fold Change") +
        theme_bw() +						# remove grey background
        theme(legend.position="bottom", axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=.5)) +
        scale_x_discrete(limits=results[1:50,]$species) + # to order x-axis as in dataframe
        scale_color_manual(name = "Adjusted P-value",
                           values = c("(-Inf,0.05]" = "red", "(0.05, Inf]" = "black"),
                           labels = c("<= 0.05", "> 0.05"),
                           drop = FALSE )
      dev.off()
      
      rm(res,resLFC,resOrdered,results,dds,ddsLove,cts,coldata,geoMeans)
      
      # DESeq2 at the genus level for GM domain
      # from mytotal, generate a count table with taxa in rows (cts), and a metadata table (coldata)
      colnames(mytotalGM)
      
      cts =  t(mytotalGM[,which(colnames(mytotalGM) %in% genera)	])
      colnames(cts) = mytotalGM$ID_afri
      coldata = mytotalGM[,c("ID_afri","haz_cont","haz","delay","run")]
      
      rownames(coldata) = mytotalGM$ID_afri
      
      all(rownames(coldata) %in% colnames(cts))
      all(rownames(coldata) == colnames(cts))   # if not : cts <- cts[, rownames(coldata)]
      
      dds <- DESeqDataSetFromMatrix(countData = cts,
                                    colData = coldata,
                                    design = ~ run + delay)   # design variable needs to be a factor
      # if we add a variable to adjust on, put the treatment condition at the end of formula
      dds
      
      all(apply(apply(cts, 2, function(x) x==0), 1, any)) # all columns have at least one zero --> DESeq2 cannot compute a log geometric means for estimating size factors
      
      geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
      ddsLove<-estimateSizeFactors(dds, geoMeans=geoMeans)
      # or ddsLove = estimateSizeFactors(dds, type="iterate")
      dds <- estimateDispersions(ddsLove)
      dds <- nbinomWaldTest(dds)
      
      res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE) # or : res <- results(dds, contrast=c("haz_bin","normal","stunted"))
      
      # Log fold change shrinkage for visualization and ranking
      resultsNames(dds)
      resLFC <- lfcShrink(dds, coef="delay_Yes_vs_No", type="normal")
      resLFC
      resOrdered <- res[order(res$pvalue),]
      summary(res)
      sum(res$padj < 0.05, na.rm=TRUE) # how many species with adjusted pvalues below 0.05
      sum(res$pvalue < 0.05, na.rm=TRUE) # how many species with pvalues below 0.05
      
      # MA-plot
      plotMA(res, ylim=c(-25,10)) # species with padj < 0.1 in blue dots
      #idx <- identify(res$baseMean, res$log2FoldChange)
      #rownames(res)[idx]
      # terminate this command using esc/echap
      
      
      write.csv(as.data.frame(resOrdered), file="DESeq2_results_GM_delay_genera.csv")
      
      results = as.data.frame(resOrdered)
      results$species = rownames(results)
      results$species = as.factor(results$species)
      results$FC = 2^results$log2FoldChange
      results$hci = results$log2FoldChange + 1.96*results$lfcSE
      results$lci = results$log2FoldChange - 1.96*results$lfcSE
      results$species[which(results$padj <0.05)]
      
      #negative values of L2FC means lower expression in delay vs no_delay
      
      pdf("afribiota_DESEQ2_genera_meancounts_GM.pdf", width=12, height=8)
      ggplot() +
        geom_col(data=results[1:50,], aes(x=species, y=baseMean)) +
        theme_bw() +
        theme(axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=.5)) +
        scale_x_discrete(limits=results[1:50,]$species) +  # to order x-axis as in dataframe
        labs(title="Mean read counts", x = "Genera", y = "baseMean")
      dev.off()
      
      pdf("afribiota_DESEQ2_genera_log2FC_GM.pdf", width=12, height=10)
      ggplot(results[1:50,], aes(x=species, y=log2FoldChange)) +
        geom_point(stat="identity", aes(colour=cut(padj,c(-Inf,0.05,Inf)) )) +
        geom_errorbar(aes(ymin=lci, ymax=hci, colour=cut(padj,c(-Inf,0.05,Inf)) ), width=.1) +
        geom_hline(yintercept = 0) +
        labs(title="Log2 Fold Change delayed vs not delayed", x = "Genera", y = "log2 Fold Change") +
        theme_bw() +						# remove grey background
        theme(legend.position="bottom", axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=.5)) +
        scale_x_discrete(limits=results[1:50,]$species) + # to order x-axis as in dataframe
        scale_color_manual(name = "Adjusted P-value",
                           values = c("(-Inf,0.05]" = "red", "(0.05, Inf]" = "black"),
                           labels = c("<= 0.05", "> 0.05"),
                           drop = FALSE )
      dev.off()
      
      # plot counts for a specific species
      plotCounts(dds, gene=which.min(res$padj), intgroup="delay")
      plotCounts(dds, gene="Streptococcus", intgroup="delay") # or d <- plotCounts(dds, gene=which.min(res$padj), intgroup="ROGRP", returnData=TRUE)
      plotCounts(dds, gene="Veillonella", intgroup="delay")
      
      rm(res,resLFC,resOrdered,results,dds,ddsLove,cts,coldata,geoMeans)
      rm(mytotalComm,mytotalPES,mytotalPS,mytotalFM,mytotalGM,mytotalOverall)
       
    # c) at the family level-----------
    # keep only two extreme quartiles of the distribution in each domain
    quantile(mytotal.fam$Comm_tot)
    mytotalComm = mytotal.fam[which(mytotal.fam$Comm_tot >= 60 | mytotal.fam$Comm_tot <= 45),]
    mytotalComm$delay = ifelse(mytotalComm$Comm_tot>=60,"No",ifelse(mytotalComm$Comm_tot <= 45,"Yes","grey"))
    table(mytotalComm$delay,useNA="always")
    
    quantile(mytotal.fam$PS_tot)
    mytotalPS = mytotal.fam[which(mytotal.fam$PS_tot >= 50 | mytotal.fam$PS_tot <= 30),]
    mytotalPS$delay = ifelse(mytotalPS$PS_tot>=50,"No",ifelse(mytotalPS$PS_tot <= 30,"Yes","grey"))
    table(mytotalPS$delay,useNA="always")
    
    quantile(mytotal.fam$PES_tot)
    mytotalPES = mytotal.fam[which(mytotal.fam$PES_tot >= 60 | mytotal.fam$PES_tot <= 45),]
    mytotalPES$delay = ifelse(mytotalPES$PES_tot>=60,"No",ifelse(mytotalPES$PES_tot <= 45,"Yes","grey"))
    table(mytotalPES$delay,useNA="always")
    
    quantile(mytotal.fam$FM_tot)
    mytotalFM = mytotal.fam[which(mytotal.fam$FM_tot >= 55 | mytotal.fam$FM_tot <= 42),]
    mytotalFM$delay = ifelse(mytotalFM$FM_tot>=55,"No",ifelse(mytotalFM$FM_tot <= 42,"Yes","grey"))
    table(mytotalFM$delay,useNA="always")
    
    quantile(mytotal.fam$GM_tot)
    mytotalGM = mytotal.fam[which(mytotal.fam$GM_tot >= 60 | mytotal.fam$GM_tot <= 50),]
    mytotalGM$delay = ifelse(mytotalGM$GM_tot>=60,"No",ifelse(mytotalGM$GM_tot <= 50,"Yes","grey"))
    table(mytotalGM$delay,useNA="always")
    
    quantile(mytotal.fam$Overal_score)
    mytotalOverall = mytotal.fam[which(mytotal.fam$Overal_score >= 265 | mytotal.fam$Overal_score <= 225),]
    mytotalOverall$delay = ifelse(mytotalOverall$Overal_score>=265,"No",ifelse(mytotalOverall$Overal_score <= 225,"Yes","grey"))
    table(mytotalOverall$delay,useNA="always")
        
        # from mytotal, generate a count table with taxa in rows (cts), and a metadata table (coldata)
        colnames(mytotalOverall)
        
        cts =  t(mytotalOverall[,which(colnames(mytotalOverall) %in% families)	])
        colnames(cts) = mytotalOverall$ID_afri
        coldata = mytotalOverall[,c("ID_afri","haz_cont","haz","delay","run")]
        
        rownames(coldata) = mytotalOverall$ID_afri
        
        all(rownames(coldata) %in% colnames(cts))
        all(rownames(coldata) == colnames(cts))   # if not : cts <- cts[, rownames(coldata)]
        
        dds <- DESeqDataSetFromMatrix(countData = cts,
                                      colData = coldata,
                                      design = ~ run + delay)   # design variable needs to be a factor
        # if we add a variable to adjust on, put the treatment condition at the end of formula
        dds
        
        all(apply(apply(cts, 2, function(x) x==0), 1, any)) # all columns have at least one zero --> DESeq2 cannot compute a log geometric means for estimating size factors
        
        # pre-filtering on taxa that have counts>10, to reduce memory size of dds object and increase speed
        # keep <- rowSums(counts(dds)) >= 10
        # dds <- dds[keep,] # stayed at 348
        # 
        # # change reference level, in case
        # dds$delay <- relevel(dds$delay, ref = "0")
        # 
        # # generate results
        # dds <- DESeq(dds) # OR the 3 steps, if impossible to get the geometric mean
        geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
        ddsLove<-estimateSizeFactors(dds, geoMeans=geoMeans)
        # or ddsLove = estimateSizeFactors(dds, type="iterate")
        dds <- estimateDispersions(ddsLove)
        dds <- nbinomWaldTest(dds)
        
        res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE) # or : res <- results(dds, contrast=c("haz_bin","normal","stunted"))
        
        # Log fold change shrinkage for visualization and ranking
        resultsNames(dds)
        resLFC <- lfcShrink(dds, coef="delay_Yes_vs_No", type="normal")
        resLFC
        resOrdered <- res[order(res$pvalue),]
        summary(res)
        sum(res$padj < 0.05, na.rm=TRUE) # how many species with adjusted pvalues below 0.05
        sum(res$pvalue < 0.05, na.rm=TRUE) # how many species with pvalues below 0.05
        
        # MA-plot
        plotMA(res, ylim=c(-25,10)) # species with padj < 0.1 in blue dots
        #idx <- identify(res$baseMean, res$log2FoldChange)
        #rownames(res)[idx]
        # terminate this command using esc/echap
        
        
        write.csv(as.data.frame(resOrdered), file="DESeq2_results_Overall_delay_families.csv")
        
        results = as.data.frame(resOrdered)
        results$species = rownames(results)
        results$species = as.factor(results$species)
        results$FC = 2^results$log2FoldChange
        results$hci = results$log2FoldChange + 1.96*results$lfcSE
        results$lci = results$log2FoldChange - 1.96*results$lfcSE
        results$species[which(results$padj <0.05)]
        
        #negative values of L2FC means lower expression in delay vs no_delay
        
        pdf("afribiota_DESEQ2_families_meancounts_Overall.pdf", width=12, height=8)
        ggplot() +
          geom_col(data=results[1:50,], aes(x=species, y=baseMean)) +
          theme_bw() +
          theme(axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=.5)) +
          scale_x_discrete(limits=results[1:50,]$species) +  # to order x-axis as in dataframe
          labs(title="Mean read counts", x = "Taxa", y = "baseMean")
        dev.off()
        
        pdf("afribiota_DESEQ2_families_log2FC_Overall.pdf", width=12, height=10)
        ggplot(results[1:50,], aes(x=species, y=log2FoldChange)) +
          geom_point(stat="identity", aes(colour=cut(padj,c(-Inf,0.05,Inf)) )) +
          geom_errorbar(aes(ymin=lci, ymax=hci, colour=cut(padj,c(-Inf,0.05,Inf)) ), width=.1) +
          geom_hline(yintercept = 0) +
          labs(title="Log2 Fold Change delayed vs not delayed", x = "Taxa", y = "log2 Fold Change") +
          theme_bw() +						# remove grey background
          theme(legend.position="bottom", axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=.5)) +
          scale_x_discrete(limits=results[1:50,]$species) + # to order x-axis as in dataframe
          scale_color_manual(name = "Adjusted P-value",
                             values = c("(-Inf,0.05]" = "red", "(0.05, Inf]" = "black"),
                             labels = c("<= 0.05", "> 0.05"),
                             drop = FALSE )
        dev.off()
        
        rm(res,resLFC,resOrdered,results,dds,ddsLove,cts,coldata,geoMeans)
        
        # for PES domain
        # from mytotal, generate a count table with taxa in rows (cts), and a metadata table (coldata)
        colnames(mytotalPES)
        
        cts =  t(mytotalPES[,which(colnames(mytotalPES) %in% families)	])
        colnames(cts) = mytotalPES$ID_afri
        coldata = mytotalPES[,c("ID_afri","haz_cont","haz","delay","run")]
        
        rownames(coldata) = mytotalPES$ID_afri
        
        all(rownames(coldata) %in% colnames(cts))
        all(rownames(coldata) == colnames(cts))   # if not : cts <- cts[, rownames(coldata)]
        
        dds <- DESeqDataSetFromMatrix(countData = cts,
                                      colData = coldata,
                                      design = ~ run + delay)   # design variable needs to be a factor
        # if we add a variable to adjust on, put the treatment condition at the end of formula
        dds
        
        all(apply(apply(cts, 2, function(x) x==0), 1, any)) # all columns have at least one zero --> DESeq2 cannot compute a log geometric means for estimating size factors
        
        geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
        ddsLove<-estimateSizeFactors(dds, geoMeans=geoMeans)
        # or ddsLove = estimateSizeFactors(dds, type="iterate")
        dds <- estimateDispersions(ddsLove)
        dds <- nbinomWaldTest(dds)
        
        res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE) # or : res <- results(dds, contrast=c("haz_bin","normal","stunted"))
        
        # Log fold change shrinkage for visualization and ranking
        resultsNames(dds)
        resLFC <- lfcShrink(dds, coef="delay_Yes_vs_No", type="normal")
        resLFC
        resOrdered <- res[order(res$pvalue),]
        summary(res)
        sum(res$padj < 0.05, na.rm=TRUE) # how many species with adjusted pvalues below 0.05
        sum(res$pvalue < 0.05, na.rm=TRUE) # how many species with pvalues below 0.05
        
        # MA-plot
        plotMA(res, ylim=c(-25,10)) # species with padj < 0.1 in blue dots
        #idx <- identify(res$baseMean, res$log2FoldChange)
        #rownames(res)[idx]
        # terminate this command using esc/echap
        
        
        write.csv(as.data.frame(resOrdered), file="DESeq2_results_PES_delay_families.csv")
        
        results = as.data.frame(resOrdered)
        results$species = rownames(results)
        results$species = as.factor(results$species)
        results$FC = 2^results$log2FoldChange
        results$hci = results$log2FoldChange + 1.96*results$lfcSE
        results$lci = results$log2FoldChange - 1.96*results$lfcSE
        results$species[which(results$padj <0.05)]
        
        #negative values of L2FC means lower expression in delay vs no_delay
        
        pdf("afribiota_DESEQ2_families_meancounts_PES.pdf", width=12, height=8)
        ggplot() +
          geom_col(data=results[1:50,], aes(x=species, y=baseMean)) +
          theme_bw() +
          theme(axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=.5)) +
          scale_x_discrete(limits=results[1:50,]$species) +  # to order x-axis as in dataframe
          labs(title="Mean read counts", x = "Taxa", y = "baseMean")
        dev.off()
        
        pdf("afribiota_DESEQ2_families_log2FC_PES.pdf", width=12, height=10)
        ggplot(results[1:50,], aes(x=species, y=log2FoldChange)) +
          geom_point(stat="identity", aes(colour=cut(padj,c(-Inf,0.05,Inf)) )) +
          geom_errorbar(aes(ymin=lci, ymax=hci, colour=cut(padj,c(-Inf,0.05,Inf)) ), width=.1) +
          geom_hline(yintercept = 0) +
          labs(title="Log2 Fold Change delayed vs not delayed", x = "Taxa", y = "log2 Fold Change") +
          theme_bw() +						# remove grey background
          theme(legend.position="bottom", axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust=.5)) +
          scale_x_discrete(limits=results[1:50,]$species) + # to order x-axis as in dataframe
          scale_color_manual(name = "Adjusted P-value",
                             values = c("(-Inf,0.05]" = "red", "(0.05, Inf]" = "black"),
                             labels = c("<= 0.05", "> 0.05"),
                             drop = FALSE )
        dev.off()
        
        # plot counts for a specific species
        plotCounts(dds, gene=which.min(res$padj), intgroup="delay")
        plotCounts(dds, gene="Streptococcaceae", intgroup="delay") # or d <- plotCounts(dds, gene=which.min(res$padj), intgroup="ROGRP", returnData=TRUE)
        plotCounts(dds, gene="Pasteurellaceae", intgroup="delay")
  
        rm(res,resLFC,resOrdered,results,dds,ddsLove,cts,coldata,geoMeans)
        rm(mytotalComm,mytotalPS,mytotalPES,mytotalFM,mytotalGM,mytotalOverall)
      

  # 8) ANCOM-BC -----------
    library(ANCOMBC)
      
      # a) at the genus level ----
      # Overall -----
        
      # keep only two extreme quartiles of the distribution in each domain
      quantile(mytotal.gen$Overal_score)
      mytotalOverall = mytotal.gen[which(mytotal.gen$Overal_score >= 265 | mytotal.gen$Overal_score <= 225),]
      mytotalOverall$delay = ifelse(mytotalOverall$Overal_score>=265,"No",ifelse(mytotalOverall$Overal_score <= 225,"Yes","grey"))
      table(mytotalOverall$delay,useNA="always")
      
      # from mytotal, generate a count table with taxa in rows (cts), and a metadata table (coldata)
      colnames(mytotalOverall)
      
      cts =  t(mytotalOverall[,which(colnames(mytotalOverall) %in% genera)	])
      colnames(cts) = mytotalOverall$ID_afri
      coldata = mytotalOverall[,c("ID_afri","haz_cont","haz","delay","run")]
      
      # Check dimensions / names
      all(rownames(coldata) == colnames(cts))  # if FALSE, reorder columns of cts
      if(!all(rownames(coldata) == colnames(cts))){
        cts <- cts[, rownames(coldata)]
      }
      
      # Create a minimal phyloseq object expected by ANCOM-BC (OTU table + sample_data).
      # If you have taxonomic table, add it as tax_table(); otherwise create a dummy tax_table
      otu <- otu_table(as.matrix(cts), taxa_are_rows = TRUE)
      sdata <- sample_data(coldata)
      
      ps <- phyloseq(otu, sdata)
      
      # --- filtering similar to your DESeq2 pre-filtering ---
      # remove taxa with extremely low total counts (e.g., sum < 10)
      keep_taxa <- rowSums(otu_table(ps)) >= 10
      psf <- prune_taxa(keep_taxa, ps) # from 226 genera to 196
      
      # Also optionally remove samples with very low library size
      lib_sizes <- sample_sums(psf)
      psf <- prune_samples(lib_sizes > 0, psf) # none was removed
      
      # --- Run ANCOM-BC ---
      # formula: include covariates if you need; treatment/group should be included in formula.
      # If you want to adjust on e.g., haz or age, use formula = "delay + haz" (put confounders before main)
      res_ancombc <- ancombc(
        phyloseq = psf,
        formula = "delay + run",           # or "delay + haz + other_covariates"
        p_adj_method = "BH",         # or "BH"
        zero_cut = 0.90,             # taxa with >90% zeros removed (adjust as needed)
        lib_cut = 1000,              # remove samples with library < 1000 reads (tune to your data)
        group = "delay",             # group variable; helpful for some outputs
        struc_zero = TRUE,           # detect structural zeros
        neg_lb = TRUE,               # use negative lower bounds for CI (recommended)
        alpha = 0.05
      )
      
      # res_ancombc is a list;
      # Extract the main components
      beta_mat <- res_ancombc$res$beta       # taxa x covariates
      se_mat   <- res_ancombc$res$se
      pval_mat <- res_ancombc$res$p_val
      qval_mat <- res_ancombc$res$q_val
      
      # Extract only the "delay" column
      delay_results <- data.frame(
        taxa = rownames(beta_mat),
        beta  = beta_mat[, "delayYes"],
        se    = se_mat[, "delayYes"],
        p_val = pval_mat[, "delayYes"],
        q_val = qval_mat[, "delayYes"],
        stringsAsFactors = FALSE
      )
      
      # ANCOM-BC 'beta' is (log) fold change for the coefficient. Confirm whether it's log or log2:
      # The package returns "beta" as log fold-change on natural log scale — convert to base-2 if you want:
      delay_results$log2FC <- delay_results$beta / log(2)   # if beta is ln FC
      delay_results$FC <- exp(delay_results$beta)          # fold change on natural scale
      
      # approximate 95% CI on log(ln) scale
      delay_results$hci <- delay_results$beta + 1.96 * delay_results$se
      delay_results$lci <- delay_results$beta - 1.96 * delay_results$se
      # convert CIs to log2 if desired:
      delay_results$hci_log2 <- delay_results$hci / log(2)
      delay_results$lci_log2 <- delay_results$lci / log(2)
      
      # Order by q-value (adjusted p-value)
      delay_results_ordered <- delay_results %>% arrange(q_val)
      
      # Count significant taxa
      sum(delay_results_ordered$q_val < 0.05, na.rm = TRUE) # [Ruminococcus]_torques_group, Ruminococcus, Weissella at the genus level
      
      # write csv
      write.csv(delay_results_ordered, file = "ANCOMBC_results_Overall_delay_genera.csv", row.names = FALSE)
      
      # --- plotting ---
      # 1) mean counts (baseMean equivalent)
      # Get counts as a matrix
      counts_mat <- as.matrix(otu_table(psf))
      # Compute mean counts only for taxa present in ANCOM-BC results
      mean_counts <- rowMeans(counts_mat[rownames(counts_mat) %in% delay_results_ordered$taxa, , drop = FALSE])
      # Merge mean counts into results table
      delay_results_ordered$mean_count <- mean_counts[delay_results_ordered$taxa]
      
      # plot top 50 by q-value or mean count
      topn <- 50
      top_taxa <- delay_results_ordered$taxa[1:topn]
      
      pdf("afribiota_ANCOMBC_genera_meancounts_Overall.pdf", width=12, height=8)
      ggplot(delay_results_ordered[1:topn, ], aes(x = factor(taxa, levels = top_taxa), y = mean_count)) +
        geom_col() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = "Mean read counts (top 50 Genera)", x = "Genera", y = "Mean count")
      dev.off()    
      
      # 2) log2 fold changes with CI (top 50)
      pdf("afribiota_ANCOMBC_genera_log2FC_Overall.pdf", width=12, height=8)
      ggplot(delay_results_ordered[1:topn, ], aes(x = factor(taxa, levels = top_taxa), y = log2FC)) +
        geom_point(aes(color = q_val <= 0.05)) +
        geom_errorbar(aes(ymin = lci_log2, ymax = hci_log2, color = q_val <= 0.05), width = .1) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = "ANCOM-BC log2 Fold Change (delay vs no delay)", x = "Genera", y = "log2 Fold Change")
      dev.off()
      
      rm(coldata, cts, delay_results, delay_results_ordered,mytotalOverall,ps,psf,pval_mat,qval_mat,beta_mat,se_mat,res_ancombc,sdata,counts_mat,keep_taxa,lib_taxa,mean_counts,otu,top_taxa,topn)
      
      
      
      # b) at the family level ----
      # Overall -----
      # keep only two extreme quartiles of the distribution in each domain
      quantile(mytotal.fam$Overal_score)
      mytotalOverall = mytotal.fam[which(mytotal.fam$Overal_score >= 265 | mytotal.fam$Overal_score <= 225),]
      mytotalOverall$delay = ifelse(mytotalOverall$Overal_score>=265,"No",ifelse(mytotalOverall$Overal_score <= 225,"Yes","grey"))
      table(mytotalOverall$delay,useNA="always")
      
      # from mytotal, generate a count table with taxa in rows (cts), and a metadata table (coldata)
      colnames(mytotalOverall)
      
      cts =  t(mytotalOverall[,which(colnames(mytotalOverall) %in% families)	])
      colnames(cts) = mytotalOverall$ID_afri
      coldata = mytotalOverall[,c("ID_afri","haz_cont","haz","delay","run")]
      
      # Check dimensions / names
      all(rownames(coldata) == colnames(cts))  # if FALSE, reorder columns of cts
      if(!all(rownames(coldata) == colnames(cts))){
        cts <- cts[, rownames(coldata)]
      }
      
      # Create a minimal phyloseq object expected by ANCOM-BC (OTU table + sample_data).
      # If you have taxonomic table, add it as tax_table(); otherwise create a dummy tax_table
      otu <- otu_table(as.matrix(cts), taxa_are_rows = TRUE)
      sdata <- sample_data(coldata)
    
      ps <- phyloseq(otu, sdata)
      
      # --- filtering similar to your DESeq2 pre-filtering ---
      # remove taxa with extremely low total counts (e.g., sum < 10)
      keep_taxa <- rowSums(otu_table(ps)) >= 10
      psf <- prune_taxa(keep_taxa, ps) # from 74 families to 65
      
      # Also optionally remove samples with very low library size
      lib_sizes <- sample_sums(psf)
      psf <- prune_samples(lib_sizes > 0, psf) # none was removed
      
      # --- Run ANCOM-BC ---
      # formula: include covariates if you need; treatment/group should be included in formula.
      # If you want to adjust on e.g., haz or age, use formula = "delay + haz" (put confounders before main)
      res_ancombc <- ancombc(
        phyloseq = psf,
        formula = "delay + run",           # or "delay + haz + other_covariates"
        p_adj_method = "BH",         # or "BH"
        zero_cut = 0.90,             # taxa with >90% zeros removed (adjust as needed)
        lib_cut = 1000,              # remove samples with library < 1000 reads (tune to your data)
        group = "delay",             # group variable; helpful for some outputs
        struc_zero = TRUE,           # detect structural zeros
        neg_lb = TRUE,               # use negative lower bounds for CI (recommended)
        alpha = 0.05
      )
      
      # res_ancombc is a list;
      # Extract the main components
      beta_mat <- res_ancombc$res$beta       # taxa x covariates
      se_mat   <- res_ancombc$res$se
      pval_mat <- res_ancombc$res$p_val
      qval_mat <- res_ancombc$res$q_val
      
      # Extract only the "delay" column
      delay_results <- data.frame(
        taxa = rownames(beta_mat),
        beta  = beta_mat[, "delayYes"],
        se    = se_mat[, "delayYes"],
        p_val = pval_mat[, "delayYes"],
        q_val = qval_mat[, "delayYes"],
        stringsAsFactors = FALSE
      )
      
      # ANCOM-BC 'beta' is (log) fold change for the coefficient. Confirm whether it's log or log2:
      # The package returns "beta" as log fold-change on natural log scale — convert to base-2 if you want:
      delay_results$log2FC <- delay_results$beta / log(2)   # if beta is ln FC
      delay_results$FC <- exp(delay_results$beta)          # fold change on natural scale
      
      # approximate 95% CI on log(ln) scale
      delay_results$hci <- delay_results$beta + 1.96 * delay_results$se
      delay_results$lci <- delay_results$beta - 1.96 * delay_results$se
      # convert CIs to log2 if desired:
      delay_results$hci_log2 <- delay_results$hci / log(2)
      delay_results$lci_log2 <- delay_results$lci / log(2)
      
      # Order by q-value (adjusted p-value)
      delay_results_ordered <- delay_results %>% arrange(q_val)
      
      # Count significant taxa
      sum(delay_results_ordered$q_val < 0.05, na.rm = TRUE) # [Ruminococcus]_torques_group, Ruminococcus, Weissella at the genus level
      
      # write csv
      write.csv(delay_results_ordered, file = "ANCOMBC_results_Overall_delay_families.csv", row.names = FALSE)
      
      # --- plotting ---
      # 1) mean counts (baseMean equivalent)
      # Get counts as a matrix
      counts_mat <- as.matrix(otu_table(psf))
      # Compute mean counts only for taxa present in ANCOM-BC results
      mean_counts <- rowMeans(counts_mat[rownames(counts_mat) %in% delay_results_ordered$taxa, , drop = FALSE])
      # Merge mean counts into results table
      delay_results_ordered$mean_count <- mean_counts[delay_results_ordered$taxa]
      
      # plot top 50 by q-value or mean count
      topn <- 50
      top_taxa <- delay_results_ordered$taxa[1:topn]
      
      pdf("afribiota_ANCOMBC_families_meancounts_Overall.pdf", width=12, height=8)
      ggplot(delay_results_ordered[1:topn, ], aes(x = factor(taxa, levels = top_taxa), y = mean_count)) +
        geom_col() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = "Mean read counts (top 50 Genera)", x = "Families", y = "Mean count")
      dev.off()    
      
      # 2) log2 fold changes with CI (top 50)
      pdf("afribiota_ANCOMBC_families_log2FC_Overall.pdf", width=12, height=8)
      ggplot(delay_results_ordered[1:topn, ], aes(x = factor(taxa, levels = top_taxa), y = log2FC)) +
        geom_point(aes(color = q_val <= 0.05)) +
        geom_errorbar(aes(ymin = lci_log2, ymax = hci_log2, color = q_val <= 0.05), width = .1) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = "ANCOM-BC log2 Fold Change (delay vs no delay)", x = "Families", y = "log2 Fold Change")
      dev.off()
      
      rm(coldata, cts, delay_results, delay_results_ordered,mytotalOverall,ps,psf,pval_mat,qval_mat,beta_mat,se_mat,res_ancombc,sdata,counts_mat,keep_taxa,lib_taxa,mean_counts,otu,top_taxa,topn)
      
      
      # PES -----
      # keep only two extreme quartiles of the distribution in each domain
      quantile(mytotal.fam$PES_tot)
      mytotalPES = mytotal.fam[which(mytotal.fam$PES_tot >= 60 | mytotal.fam$PES_tot <= 45),]
      mytotalPES$delay = ifelse(mytotalPES$PES_tot>=60,"No",ifelse(mytotalPES$PES_tot <= 45,"Yes","grey"))
      table(mytotalPES$delay,useNA="always")
      
      # from mytotal, generate a count table with taxa in rows (cts), and a metadata table (coldata)
      colnames(mytotalPES)
      
      cts =  t(mytotalPES[,which(colnames(mytotalPES) %in% families)	])
      colnames(cts) = mytotalPES$ID_afri
      coldata = mytotalPES[,c("ID_afri","haz_cont","haz","delay","run")]
      
      # Check dimensions / names
      all(rownames(coldata) == colnames(cts))  # if FALSE, reorder columns of cts
      if(!all(rownames(coldata) == colnames(cts))){
        cts <- cts[, rownames(coldata)]
      }
      
      # Create a minimal phyloseq object expected by ANCOM-BC (OTU table + sample_data).
      # If you have taxonomic table, add it as tax_table(); otherwise create a dummy tax_table
      otu <- otu_table(as.matrix(cts), taxa_are_rows = TRUE)
      sdata <- sample_data(coldata)
      
      ps <- phyloseq(otu, sdata)
      
      # --- filtering similar to your DESeq2 pre-filtering ---
      # remove taxa with extremely low total counts (e.g., sum < 10)
      keep_taxa <- rowSums(otu_table(ps)) >= 10
      psf <- prune_taxa(keep_taxa, ps) # from 74 to 66 families
      
      # Also optionally remove samples with very low library size
      lib_sizes <- sample_sums(psf)
      psf <- prune_samples(lib_sizes > 0, psf) # none was removed
      
      # --- Run ANCOM-BC ---
      # formula: include covariates if you need; treatment/group should be included in formula.
      # If you want to adjust on e.g., haz or age, use formula = "delay + haz" (put confounders before main)
      res_ancombc <- ancombc(
        phyloseq = psf,
        formula = "delay + run",           # or "delay + haz + other_covariates"
        p_adj_method = "BH",         # or "BH"
        zero_cut = 0.90,             # taxa with >90% zeros removed (adjust as needed)
        lib_cut = 1000,              # remove samples with library < 1000 reads (tune to your data)
        group = "delay",             # group variable; helpful for some outputs
        struc_zero = TRUE,           # detect structural zeros
        neg_lb = TRUE,               # use negative lower bounds for CI (recommended)
        alpha = 0.05
      )
      
      # res_ancombc is a list;
      # Extract the main components
      beta_mat <- res_ancombc$res$beta       # taxa x covariates
      se_mat   <- res_ancombc$res$se
      pval_mat <- res_ancombc$res$p_val
      qval_mat <- res_ancombc$res$q_val
      
      # Extract only the "delay" column
      delay_results <- data.frame(
        taxa = rownames(beta_mat),
        beta  = beta_mat[, "delayYes"],
        se    = se_mat[, "delayYes"],
        p_val = pval_mat[, "delayYes"],
        q_val = qval_mat[, "delayYes"],
        stringsAsFactors = FALSE
      )
      
      # ANCOM-BC 'beta' is (log) fold change for the coefficient. Confirm whether it's log or log2:
      # The package returns "beta" as log fold-change on natural log scale — convert to base-2 if you want:
      delay_results$log2FC <- delay_results$beta / log(2)   # if beta is ln FC
      delay_results$FC <- exp(delay_results$beta)          # fold change on natural scale
      
      # approximate 95% CI on log(ln) scale
      delay_results$hci <- delay_results$beta + 1.96 * delay_results$se
      delay_results$lci <- delay_results$beta - 1.96 * delay_results$se
      # convert CIs to log2 if desired:
      delay_results$hci_log2 <- delay_results$hci / log(2)
      delay_results$lci_log2 <- delay_results$lci / log(2)
      
      # Order by q-value (adjusted p-value)
      delay_results_ordered <- delay_results %>% arrange(q_val)
      
      # Count significant taxa
      sum(delay_results_ordered$q_val < 0.05, na.rm = TRUE) # [Ruminococcus]_torques_group, Ruminococcus, Weissella at the genus level
      
      # write csv
      write.csv(delay_results_ordered, file = "ANCOMBC_results_PES_delay_families.csv", row.names = FALSE)
      
      # --- plotting ---
      # 1) mean counts (baseMean equivalent)
      # Get counts as a matrix
      counts_mat <- as.matrix(otu_table(psf))
      # Compute mean counts only for taxa present in ANCOM-BC results
      mean_counts <- rowMeans(counts_mat[rownames(counts_mat) %in% delay_results_ordered$taxa, , drop = FALSE])
      # Merge mean counts into results table
      delay_results_ordered$mean_count <- mean_counts[delay_results_ordered$taxa]
      
      # plot top 50 by q-value or mean count
      topn <- 50
      top_taxa <- delay_results_ordered$taxa[1:topn]
      
      pdf("afribiota_ANCOMBC_families_meancounts_PES.pdf", width=12, height=8)
      ggplot(delay_results_ordered[1:topn, ], aes(x = factor(taxa, levels = top_taxa), y = mean_count)) +
        geom_col() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = "Mean read counts (top 50 Genera)", x = "Families", y = "Mean count")
      dev.off()    
      
      # 2) log2 fold changes with CI (top 50)
      pdf("afribiota_ANCOMBC_families_log2FC_PES.pdf", width=12, height=8)
      ggplot(delay_results_ordered[1:topn, ], aes(x = factor(taxa, levels = top_taxa), y = log2FC)) +
        geom_point(aes(color = q_val <= 0.05)) +
        geom_errorbar(aes(ymin = lci_log2, ymax = hci_log2, color = q_val <= 0.05), width = .1) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = "ANCOM-BC log2 Fold Change (delay vs no delay)", x = "Families", y = "log2 Fold Change")
      dev.off()
      
      rm(coldata, cts, delay_results, delay_results_ordered,mytotalPES,ps,psf,pval_mat,qval_mat,beta_mat,se_mat,res_ancombc,sdata,counts_mat,keep_taxa,lib_sizes,mean_counts,otu,top_taxa,topn)
    
      
        
  # 9) Beta-diversity -----------
  # at ASV level only
      
  # dissimilarity matrix
  distance_matrix  = vegan::vegdist(reads.asv_ra, "euclidean" )
  class(distance_matrix)
  hist(distance_matrix)
  
  # PCoA with stats package
  my_pcoa_ra = stats::cmdscale(distance_matrix, k=3, eig=TRUE)
  # or with vegan package: my_pcoa2 = vegan::wcmdscale(distance_matrix, k=2, eig = TRUE)
  print(my_pcoa_ra)
  my_pcoa_ra <- as.data.frame(my_pcoa_ra$points)
  names(my_pcoa_ra)[1:3] <- c('PC1','PC2','PC3')
  which(rownames(my_pcoa_ra) != mytotal$ID_afri) # check that order of IDs is the same in both datasets
  which(rownames(my_pcoa_ra) != mytotal.asv$ID_afri) # check that order of IDs is the same in both datasets
  my_pcoa_ra$Comm_delay <- mytotal$Comm_delay
  
  # plot the Principal Coordinates Analysis Results
  ggplot(my_pcoa_ra) +
    geom_point(aes(x = PC1, y = PC2, colour = Comm_delay),size = 1) +
    stat_ellipse(geom = "polygon",
                 aes(x = PC1, y = PC2, fill = Comm_delay),
                 alpha = 0.25) +
    theme_minimal()
  
  # PERMANOVA test
  set.seed(452)
  perm <- adonis2(reads.asv_ra ~ Comm_delay, data = mytotal.asv, permutations = 999, method="euclidean")
  perm[1,5]
  
  
  
  # on clr transformed data: Euclidean distance applied to clr transformed data =
  # dissimilarity matrix
  distance_matrix  = vegan::vegdist(reads.asv_clr, "euclidean" )
  hist(distance_matrix)
  
  # PCoA with stats package
  my_pcoa = stats::cmdscale(distance_matrix, k=3, eig=TRUE)
  my_pcoa <- as.data.frame(my_pcoa$points)
  names(my_pcoa)[1:3] <- c('PC1','PC2','PC3')
  which(rownames(my_pcoa) != mytotal$ID_afri) # check that order of IDs is the same in both datasets
  which(rownames(my_pcoa) != mytotal.asv$ID_afri) # check that order of IDs is the same in both datasets
  my_pcoa$Comm_delay <- mytotal$Comm_delay
  
  # plot the Principal Coordinates Analysis Results
  ggplot(my_pcoa) +
    geom_point(aes(x = PC1, y = PC2, colour = Comm_delay),size = 1) +
    stat_ellipse(geom = "polygon",
                 aes(x = PC1, y = PC2, fill = Comm_delay),
                 alpha = 0.25) +
    theme_minimal()
  
  # PERMANOVA test
  set.seed(452)
  perm <- adonis2(reads.asv_clr ~ Comm_delay, data = mytotal.asv, permutations = 999, method="euclidean")
  perm[1,5]
  
  
  
  # Create a function to do it on all domains
  plotPcoa = function (df, distmeth, dfvar, var) {
    
    # dissimilarity matrix
    distance_matrix  = vegan::vegdist(df, distmeth )
    
    # PCoA with stats package
    my_pcoa = stats::cmdscale(distance_matrix, k=2, eig=T)
    my_pcoa1 = stats::cmdscale(distance_matrix, k=1, eig=T)
    GOF1 = round(my_pcoa1$GOF[1],2)*100
    GOF2 = round((my_pcoa$GOF[1] - my_pcoa1$GOF[1]),2)*100
    col1 = paste0('PC1',' ','[', GOF1, '%]')
    col2 = paste0('PC2',' ','[', GOF2, '%]')
    
    my_pcoa <- as.data.frame(my_pcoa$points)
    my_pcoa = cbind(my_pcoa, dfvar[,var])
    names(my_pcoa)[1:3] <- c(col1,col2, var)
    
    set.seed(452)
    form = as.formula(paste("df", "~", var))
    perm <- adonis2(form, data = dfvar, permutations = 999, method= distmeth)
    perm[1,5]
    
    # plot the Principal Coordinates Analysis Results
    var <- ensym(var)
    
    gg = ggplot(my_pcoa) +
      geom_point(aes(x = my_pcoa[,1], y = my_pcoa[,2], colour = !!var),size = 1) +
      stat_ellipse(geom = "polygon",
                   aes(x = my_pcoa[,1], y = my_pcoa[,2], fill = !!var), alpha = 0.25) +
      ggtitle(paste0(distmeth,","," ","p=", perm[1,5])) +
      theme_minimal() +
      xlab(col1) +
      ylab(col2) +
      theme(plot.title = element_text(vjust = - 5, size=10))
    
    return(gg)
    
  }
  
  p1 = plotPcoa(reads.asv_ra, "euclidean", mytotal.asv, "Comm_delay")
  p2 = plotPcoa(reads.asv_ra, "bray", mytotal.asv, "Comm_delay")
  p3 = plotPcoa(reads.asv_ra, "euclidean", mytotal.asv, "PS_delay")
  p4 = plotPcoa(reads.asv_ra, "bray", mytotal.asv, "PS_delay")
  p5 = plotPcoa(reads.asv_ra, "euclidean", mytotal.asv, "PES_delay")
  p6 = plotPcoa(reads.asv_ra, "bray", mytotal.asv, "PES_delay")
  p7 = plotPcoa(reads.asv_ra, "euclidean", mytotal.asv, "GM_delay")
  p8 = plotPcoa(reads.asv_ra, "bray", mytotal.asv, "GM_delay")
  p9 = plotPcoa(reads.asv_ra, "euclidean", mytotal.asv, "FM_delay")
  p10 = plotPcoa(reads.asv_ra, "bray", mytotal.asv, "FM_delay")
  p11 = plotPcoa(reads.asv_ra, "euclidean", mytotal.asv, "Overal_delay")
  p12 = plotPcoa(reads.asv_ra, "bray", mytotal.asv, "Overal_delay")
  
  
  pdf(file="PCoA_ASV_ra.pdf", height=10, width=12)
  ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 + rremove("x.text"),
            labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
            ncol = 4, nrow = 3)
  dev.off()
  
  
  rm(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, distance_matrix)
  
  
  # with clr transformed data
  p1 = plotPcoa(reads.asv_clr, "euclidean", mytotal.asv, "Comm_delay")
  p3 = plotPcoa(reads.asv_clr, "euclidean", mytotal.asv, "PS_delay")
  p5 = plotPcoa(reads.asv_clr, "euclidean", mytotal.asv, "PES_delay")
  p7 = plotPcoa(reads.asv_clr, "euclidean", mytotal.asv, "GM_delay")
  p9 = plotPcoa(reads.asv_clr, "euclidean", mytotal.asv, "FM_delay")
  p11 = plotPcoa(reads.asv_clr, "euclidean", mytotal.asv, "Overal_delay")
  
  
  pdf(file="PCoA_ASV_clr.pdf", height=8, width=10)
  ggarrange(p1, p3, p5, p7, p9, p11 + rremove("x.text"),
            labels = c("A", "B", "C", "D", "E", "F"),
            ncol = 3, nrow = 2)
  dev.off()
  
  
  rm(p1, p3, p5, p7, p9, p11)
  rm(perm,my_pcoa,data.sites,example_NMDS)

  # PCoA with euclidean distance on RA data generates a PC1 than explains 18% of the variance, more than other distance metrics
  # Add PC1 (euclidean, RA data) to mytotal
  my_pcoa_ra$ID_afri <- rownames(my_pcoa_ra)
  my_pcoa_ra <- my_pcoa_ra[,-which(colnames(my_pcoa_ra) == "Comm_delay")]
  mytotal = left_join(mytotal,my_pcoa_ra, by = "ID_afri")
  mytotal.asv = left_join(mytotal,my_pcoa_ra, by = "ID_afri")
  mytotal.fam = left_join(mytotal.fam,my_pcoa_ra, by = "ID_afri")
  
  
  # 10) Cluster analysis--------

    # a) hierarchical clustering, global heatmap--------
    source("heatmapLib2.R")
    library(gplots)
    library("stats")
    library("vegan")
    
    # check if taxa with total abundance of 0 or constant value
    names(reads.asv_ra[,apply(reads.asv_ra, 2, sum, na.rm=TRUE) == 0])
    names(reads.asv_ra[,apply(reads.asv_ra, 2, var, na.rm=TRUE) == 0]) # or names(mydata[, sapply(mydata, function(v) var(v, na.rm=TRUE)==0)])
    
    # hierarchical clustering
    is.na(reads.asv_ra) %>% table()
    is.matrix(reads.asv_ra)
    # JSdist=vegdist(mydata3,method="jensen-shannon",useShrinkage=TRUE)
    hc=hclust(dist(reads.asv_ra),method="ward.D2")
    #plot(hc, hang=-1)
    plot(hc)
    dend=as.dendrogram(hc)
    
    
    # determine optimal number of clusters
    library(clusterCrit)
    
    CH <- vector("numeric", 9L)
    for (k in 2:10){
    clrs=cutree(hc,k=k)
    CH[k-1] = intCriteria(as.matrix(reads.asv_ra),clrs,"all")$calinski_harabasz
    # wss <- WSS( data = reads_ra, groups =  clrs )
    # tss <- Distance( cluster = reads_ra )
    # bss <- tss - wss
    # CH[k-1] = bss*(333-k)/(wss*(k-1))
    }
    
    pdf(file="CHplot.pdf")
    ggplot() +
    geom_point(aes(x=c(2:10),y=CH)) +
    geom_line(aes(x=c(2:10),y=CH)) +
    xlab("Number of clusters (k)") +
    ylab("Calinski-Harabasz index") +
    theme_minimal()
    dev.off()
    
    nClrs=2
    clrs=cutree(hc,k=nClrs)
    table(clrs)
    
    # Sidebars
    # number of domains with delays
    delayf = factor(mytotal$Overal_delay)
    levels(delayf) # "0" "1"
    delaysTbl=c()
    delaysTbl[levels(delayf)[1]]="green"
    delaysTbl[levels(delayf)[2]]="red"
      
    
    # haz status
    hazf = factor(mytotal$haz)
    levels(hazf) # "0" "1" "2"
    hazTbl=c()
    hazTbl[levels(hazf)[1]] = "green"
    hazTbl[levels(hazf)[2]] = "yellow"
    hazTbl[levels(hazf)[3]] = "red"
          
    # clusters
    clrsf=factor(clrs)
    levels(clrsf)
    clrsTbl=c()
    clrsTbl[levels(clrsf)[1]]="red"
    clrsTbl[levels(clrsf)[2]]="blue"
        
          
    sideBars.ward=cbind(clrsTbl[clrsf], delaysTbl[delayf], hazTbl[hazf])
    colnames(sideBars.ward)=c("Clusters", "Overall delay", "Haz status")
    
    nCols=20
    pdf("heatmap_afribiota.pdf")	#
    heatmap2(as.matrix(reads_ra)[,1:nCols],   # show only nCols first columns of the table tbl
             col=rainbow(50,start=1/6,end=0), # set heatmap color palette
             Colv=NA,                       	# comment this out if you want to have columns clustered as well
             Rowv=dend,			                  # use our h. clustering
             RowSideColors=sideBars.ward,     # this sets three color columns: clustering, pH and Nugent score
             RowSideTitleCex=0.5,             # scaling factor for the titles of the color bars
             RowSideTitleLine=0.5,            # the distance of the titles for the color bars
             margins=c(13,15),                # =c(bottom margin, right margin)
             #labRow=rownames(mydata),        # add row labels
             labRow=NA,                       # suppress row labels
             cexCol=0.9,                      # magnification factor for column labels
             xlas=2,                          # column labels are to be perpendicular to the x-axis
             main=""                         # main title
             )                         
    legend(0.75,0.75, legend=c("1","2"), bg="white", fill=clrsTbl, title="clusters", cex=0.7)
    legend(0.75,0.5, legend=c("0","1"), bg="white", fill=delaysTbl, title="Overall delay", cex=0.7)
    legend(0.75,0.25, legend=c("Normal","Moderate","Severe"), bg="white", fill=hazTbl, title="Haz status", cex=0.7)
    dev.off()
    
    mytotal$clusters = as.factor(clrs)
    mytotal.asv$clusters = as.factor(clrs)
    table(mytotal$clusters)
    table(clrs)
    
    rm(nClrs,hc, dend, k, nCols, hazf, hazTbl, delayf, delaysTbl,sideBars.ward, CH)
    
          
    # b) statistical testing-------

    # library(gtsummary)
    mytotal %>%
      gtsummary::tbl_summary(
        include = c("Overal_score","Comm_tot", "PS_tot","PES_tot","FM_tot","GM_tot"),
        by = "clusters", # ma variable en 5 modalités qui sera en colonnes
        missing_text ="NA", # choisir l'intitulé des données manquantes
        statistic = list(all_continuous() ~ "{mean} ({min}-{max})")) %>% # choisir lle format des statistiques de variables continues
      add_overall(last = TRUE) %>%
      add_p() %>% # fait les tests stat appropriés et indique la pvalue correspondante
      bold_labels()
    # none is significant
    
    # Using delays as continuous variables
    boxplot(Comm_tot ~ clusters, data=mytotal.asv)
    boxplot(PS_tot ~ clusters, data=mytotal.asv)
    boxplot(PES_tot ~ clusters, data=mytotal.asv)
    boxplot(GM_tot ~ clusters, data=mytotal.asv)
    boxplot(FM_tot ~ clusters, data=mytotal.asv)
    boxplot(Overal_score ~ clusters, data=mytotal.asv)
    
    summary(aov(Comm_tot ~ clusters + run, data=mytotal.asv))
    summary(aov(PS_tot ~ clusters + run, data=mytotal.asv))
    summary(aov(PES_tot ~ clusters + run, data=mytotal.asv))
    summary(aov(GM_tot ~ clusters + run, data=mytotal.asv))
    summary(aov(FM_tot ~ clusters + run, data=mytotal.asv))
    summary(aov(Overal_score ~ clusters + run, data=mytotal.asv))
    t.test(FM_tot ~ clusters, data=mytotal.asv)
    
    # graphical representation
    library(rstatix)
    library(ggpubr)
    
    stat.test <- mytotal.asv %>%
      rstatix::t_test(FM_tot ~ clusters) %>%
      rstatix::add_xy_position( x = "clusters", fun = "max", step.increase=0.5)
    
    ggboxplot(mytotal.asv, x = "clusters", y = "FM_tot", fill = "clusters") +
      stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01) +
      ylim(min=0, max=80)
    
    # using student t test
    plotStatsStud = function(df, var2, var3)
    {
      var2 <- ensym(var2)
      var3 <- ensym(var3)
      
      f <- as.formula(paste(as.character(var2), "~", as.character(var3)))
      stat.test <- df %>%
        rstatix::t_test(f) %>%
        rstatix::add_xy_position(step.increase = 0.5)
      
      ggplot(df, aes(x = !!var3, y = !!var2, fill = factor(!!var3))) +
        geom_boxplot() +
        scale_x_discrete(expand = c(0.5, 0.5)) +
        geom_bracket(data = stat.test, inherit.aes = FALSE,
                     aes(xmin = group1, xmax = group2, label = p), size=0.3, label.size = 3, tip.length = 0) +
        labs(fill = as.character(var3) ) +
        theme(panel.spacing = unit(0, "points"),
              strip.placement = "outside",
              strip.background = element_blank()) +
        theme_minimal()
    }
    
    plotStatsStud(mytotal.asv, FM_tot, clusters)
    p1 = plotStatsStud(mytotal.asv, Comm_tot, clusters)
    p2 = plotStatsStud(mytotal.asv, PS_tot, clusters)
    p3 = plotStatsStud(mytotal.asv, PES_tot, clusters)
    p4 = plotStatsStud(mytotal.asv, GM_tot, clusters)
    p5 = plotStatsStud(mytotal.asv, FM_tot, clusters)
    p6 = plotStatsStud(mytotal.asv, Overal_score, clusters)
    
    pdf(file="Clusters_cont.pdf", height=8, width=10)
    ggarrange(p1, p2, p3, p4, p5, p6 + rremove("x.text"),
              labels = c("A", "B", "C", "D", "E", "F"),
              ncol = 3, nrow = 2)
    dev.off()
    
    rm(p1, p2, p3, p4, p5, p6, stat.test)
    
    # using wilcoxon test
    plotStats <- function(df, var2, var3)
    {
      var2 <- ensym(var2)
      var3 <- ensym(var3)

      f <- as.formula(paste(as.character(var2), "~", as.character(var3)))
      stat.test <- df %>%
        rstatix::wilcox_test(f) %>%
        rstatix::add_xy_position(step.increase = 0.5)

      ggplot(df, aes(x = !!var3, y = !!var2, fill = factor(!!var3))) +
        geom_boxplot() +
        scale_x_discrete(expand = c(0.5, 0.5)) +
        geom_bracket(data = stat.test, inherit.aes = FALSE,
                     aes(xmin = group1, xmax = group2, label = p), size=0.3, label.size = 3, tip.length = 0) +
        labs(fill = as.character(var3) ) +
        theme(panel.spacing = unit(0, "points"),
              strip.placement = "outside",
              strip.background = element_blank()) +
        theme_minimal()
    }
    
    p7 = plotStats(mytotal.asv, Comm_tot, clusters)
    p8 = plotStats(mytotal.asv, PS_tot, clusters)
    p9 = plotStats(mytotal.asv, PES_tot, clusters)
    p10= plotStats(mytotal.asv, GM_tot, clusters)
    p11= plotStats(mytotal.asv, FM_tot, clusters)
    p12= plotStats(mytotal.asv, Overal_score, clusters)
    
    
    pdf(file="Clusters_cont_wilc.pdf", height=8, width=10)
    ggarrange(p7, p8, p9, p10, p11, p12 + rremove("x.text"),
              labels = c("A", "B", "C", "D", "E", "F"),
              ncol = 3, nrow = 2)
    dev.off()
    
    rm(p7, p8, p9, p10, p11, p12)
    


    
    
    
    
# IV) SEM ----------

   
    # 1) Correlation between indicator variables for each latent variable --------

    # mixed variable type or factors only
    require(tidyverse)
    require(lsr)
    # Calculate a pairwise association between all variables in a data-frame. In particular nominal vs nominal with Chi-square, numeric vs numeric with Pearson correlation, and nominal vs numeric with ANOVA.
    # Adopted from https://stackoverflow.com/a/52557631/590437
    mixed_assoc = function(df, cor_method="spearman"){
      df_comb = expand.grid(names(df), names(df),  stringsAsFactors = F) 
      colnames(df_comb) <- c("X1","X2")
      
      is_nominal = function(x) class(x) %in% c("factor", "character")
      # https://community.rstudio.com/t/why-is-purr-is-numeric-deprecated/3559
      # https://github.com/r-lib/rlang/issues/781
      is_numeric <- function(x) { is.integer(x) || is_double(x)}
      
      f = function(xName,yName) {
        x =  pull(df, xName)
        y =  pull(df, yName)
        
        result = if(is_nominal(x) && is_nominal(y)){
          # use bias corrected cramersV as described in https://rdrr.io/cran/rcompanion/man/cramerV.html
          cv = cramersV(as.character(x), as.character(y))
          data.frame(xName, yName, assoc=cv, type="cramersV")
          
        }else if(is_numeric(x) && is_numeric(y)){
          correlation = cor(x, y, method=cor_method, use="complete.obs")
          data.frame(xName, yName, assoc=correlation, type="correlation")
          
        }else if(is_numeric(x) && is_nominal(y)){
          # from https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable/124618#124618
          r_squared = summary(lm(x ~ y))$r.squared
          data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
          
        }else if(is_nominal(x) && is_numeric(y)){
          r_squared = summary(lm(y ~x))$r.squared
          data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
          
        }else {
          warning(paste("unmatched column type combination: ", class(x), class(y)))
        }
        
        # finally add complete obs number and ratio to table
        result %>% 
          #mutate(complete_obs_pairs=sum(!is.na(x) & !is.na(y)), complete_obs_ratio=complete_obs_pairs/length(x)) %>%
          rename("xName"="x", "yName"="y")
      }
      
      # apply function to each variable combination
      map2_df(df_comb$X1, df_comb$X2, f)
    }
    
    cor1 = mixed_assoc(mytotal[,c("Overal_score","PS_tot","PES_tot","FM_tot","GM_tot","Comm_tot")])
    cor2 = mixed_assoc(mytotal[,c("score_socio_eco_bycountry","instruc_mere_3","nb_pieces","eau_traitee","age_prem_gross")])
    cor3 = mixed_assoc(mytotal[,c("asf","dds")]) # 0.47
    cor4 = mixed_assoc(mytotal[,c("AATlevel2","CALwet_cat","crp2_cat","CitrulineM_cat")])
    cor5 = mixed_assoc(mytotal[,c("ascaris","trichuris","giardiase","No_genes")])
    cor6 = mixed_assoc(mytotal[,c("ompC","ipaH","estla","eltB","AggR","aaiC","bfpA","eae","cadF","cta")])
    cor7 = mixed_assoc(mytotal[,c("AlanineM","CitrulineM","ValineM","IsoleucineM","LeucineM")])
    cor8 = mixed_assoc(mytotal.fam[,c("Lachnospiraceae","Enterobacteriaceae","Streptococcaceae","Lactobacillaceae")])
    
    custom_labels <- c("Overal_score"="Overall score", "PS_tot"="PS score", "PES_tot"="PES score", "FM_tot" = "FM score", "GM_tot" = "GM score", "Comm_tot" = "Comm score")
    # Replace x and y in cor1
    cor1$x <- custom_labels[as.character(cor1$x)]
    cor1$y <- custom_labels[as.character(cor1$y)]
    # reorder levels
    cor1$x = factor(cor1$x, levels = c("Overall score","PS score","PES score","FM score","GM score","Comm score"), labels = c("Overall score","PS score","PES score","FM score","GM score","Comm score"))
    cor1$y = factor(cor1$y, levels = c("Overall score","PS score","PES score","FM score","GM score","Comm score"), labels = c("Overall score","PS score","PES score","FM score","GM score","Comm score"))
    
    custom_labels <- c("score_socio_eco_bycountry"="Family economic status", "instruc_mere_3"="Maternal education level cat", "nb_pieces"="Number of rooms", "eau_traitee" = "Treated drinking water", "age_prem_gross" = "Maternal age at first pregnancy") 
    # Replace x and y in cor2
    cor2$x <- custom_labels[as.character(cor2$x)]
    cor2$y <- custom_labels[as.character(cor2$y)]
    
    custom_labels <- c("AATlevel2" = "AAT level","CALwet_cat" = "Calprotectine level cat (wet)","crp2_cat" = "CRP level cat","CitrulineM_cat" = "Citruline cat") 
    # Replace x and y in cor4
    cor4$x <- custom_labels[as.character(cor4$x)]
    cor4$y <- custom_labels[as.character(cor4$y)]
    
    custom_labels <- c("ascaris"="Ascaris","trichuris"="Trichuris","giardiase"="Giardiasis", "No_genes"="Number of enteropathogens genes") 
    # Replace x and y in cor5
    cor5$x <- custom_labels[as.character(cor5$x)]
    cor5$y <- custom_labels[as.character(cor5$y)]
    # reorder levels
    cor5$x = factor(cor5$x, levels = c("Ascaris","Trichuris","Giardiasis","Number of enteropathogens genes"), labels = c("Ascaris","Trichuris","Giardiasis","Number of enteropathogens genes"))
    cor5$y = factor(cor5$y, levels = c("Ascaris","Trichuris","Giardiasis","Number of enteropathogens genes"), labels = c("Ascaris","Trichuris","Giardiasis","Number of enteropathogens genes"))
    
    custom_labels <- c("AlanineM" = "Alanine","CitrulineM"="Citrulline","ValineM"="Valine","IsoleucineM"="Isoleucine","LeucineM"="Leucine") 
    # Replace x and y in cor7
    cor7$x <- custom_labels[as.character(cor7$x)]
    cor7$y <- custom_labels[as.character(cor7$y)]
    
    library(ggcorrplot)
    pdf("correlation_plots_indicator_variables.pdf", width = 9)
    #ggcorrplot(cor1,show.diag=TRUE, type="full", lab=TRUE, lab_size=4,colors = c("red", "white", "blue"))
    cor1 %>%
      ggplot(aes(x,y,fill=assoc)) +
      geom_tile() +
      geom_text(aes(x,y,label=round(assoc,2))) +
      scale_fill_gradientn(values=c(1, 0, -1), colours=c("blue", "white", "red"), name = "Correlation") +
      theme_classic()  +
      labs(x = "", y = "") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12), axis.text.y = element_text( size=12 ) )
    cor2 %>%
      ggplot(aes(x,y,fill=assoc)) +
      geom_tile() +
      geom_text(aes(x,y,label=round(assoc,2))) +
      scale_fill_gradientn(values=c(1, 0, -1), colours=c("blue", "white", "red"), name = "Correlation") +
      theme_classic()  +
      labs(x = "", y = "") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12), axis.text.y = element_text( size=12 ) )
    cor4 %>%
      ggplot(aes(x,y,fill=assoc)) +
      geom_tile() +
      geom_text(aes(x,y,label=round(assoc,2))) +
      scale_fill_gradientn(values=c(1, 0, -1), colours=c("blue", "white", "red"), name = "Correlation") +
      theme_classic() +
      labs(x = "", y = "") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12), axis.text.y = element_text( size=12 ) )
    cor5 %>%
      ggplot(aes(x,y,fill=assoc)) +
      geom_tile() +
      geom_text(aes(x,y,label=round(assoc,2))) +
      scale_fill_gradientn(values=c(1, 0, -1), colours=c("blue", "white", "red"), name = "Correlation") +
      theme_classic() +
      labs(x = "", y = "") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12), axis.text.y = element_text( size=12 ) )
    cor6 %>%
      ggplot(aes(x,y,fill=assoc)) +
      geom_tile() +
      geom_text(aes(x,y,label=round(assoc,2))) +
      scale_fill_gradientn(values=c(1, 0, -1), colours=c("blue", "white", "red"), name = "Correlation") +
      theme_classic() +
      labs(x = "", y = "") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12), axis.text.y = element_text( size=12 ) )
    cor7 %>%
      ggplot(aes(x,y,fill=assoc)) +
      geom_tile() +
      geom_text(aes(x,y,label=round(assoc,2))) +
      scale_fill_gradientn(values=c(1, 0, -1), colours=c("blue", "white", "red"), name = "Correlation") +
      theme_classic() +
      labs(x = "", y = "") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12), axis.text.y = element_text( size=12 ) )
    dev.off()
    
    rm(cor1, cor2, cor3, cor4, cor5, cor6,cor7,cor8)
    
    
    
    
    # 2) CFA model ------
    library(lavaan)
    # order factor variables that are unordered: haz, AATlevel2 CALwet_cat asf score_socio_eco_bycountry instruc_mere_4 eau_traitee hemoglobine2_seuil
    mytotal[,c("haz", "AATlevel2", "CALwet_cat","crp2_cat","CitrulineM_cat", "score_socio_eco_bycountry", "instruc_mere_3", "eau_traitee", "hemoglobine2_seuil")] <- lapply(mytotal[,c("haz", "AATlevel2", "CALwet_cat","crp2_cat","CitrulineM_cat", "score_socio_eco_bycountry", "instruc_mere_3", "eau_traitee", "hemoglobine2_seuil")], ordered)
    
    cfa.model1 <- 'Neurodevelopment =~ Comm_tot + PS_tot + PES_tot + FM_tot + GM_tot'
    fitcfa1 <- lavaan::cfa(cfa.model1, data = mytotal, estimator = "MLR") # MLR is an robust maximum likelihood estimator, better when the indicators are not super well normally distributed
    lavaan::summary(fitcfa1, fit.measures = TRUE, standardized=TRUE)
    
    cfa.model2 <- 'Enteric  =~ ascaris + trichuris + giardiase + No_genes'
    fitcfa2 <- lavaan::cfa(cfa.model2, data = mytotal)
    lavaan::summary(fitcfa2, fit.measures = TRUE, standardized=TRUE) # really bad
    
    cfa.model3 <- 'qpcr  =~  ompC + ipaH  + estla  + eltB  + AggR  + aaiC  + bfpA  + eae  + cadF  + cta'
    fitcfa3 <- lavaan::cfa(cfa.model3, data = mytotal)
    lavaan::summary(fitcfa3, fit.measures = TRUE, standardized=TRUE)
    
    cfa.model4 <- 'Inflammation  =~  AATlevel2 + CALwet_cat + crp2_cat + CitrulineM_cat'
    fitcfa4 <- lavaan::cfa(cfa.model4, data = mytotal)
    lavaan::summary(fitcfa4, fit.measures = TRUE, standardized=TRUE) # really bad. Better to use just AATlevel2 or CALwet_cat
    questionr::freq.na(mytotal[,which(names(mytotal) %in% c("AATlevel2","CALwet_cat","crp2_cat","CitrulineM_cat"))]) # CitrulineM_cat has 12% missing variables
    
    cfa.model5 <- 'SES  =~   score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross'
    fitcfa5 <- lavaan::cfa(cfa.model5, data = mytotal)
    lavaan::summary(fitcfa5, fit.measures = TRUE, standardized=TRUE)  # Good
    
    cfa.model6 <- 'BCAA  =~ CitrulineM + ValineM + IsoleucineM + LeucineM'
    fitcfa6 <- lavaan::cfa(cfa.model6, data = mytotal, estimator = "MLR")
    lavaan::summary(fitcfa6, fit.measures = TRUE, standardized=TRUE)  # Good except RMSEA
  
    
    # Try with relative abundance data
    mytotal.fam_ra = mytotal.fam[,-c(which(colnames(mytotal.fam) %in% colnames(reads.fam_ra)))]
    reads.fam_ra$ID_afri = rownames(reads.fam_ra)
    mytotal.fam_ra = inner_join(mytotal.fam_ra, reads.fam_ra, by="ID_afri")
    reads.fam_ra = reads.fam_ra[,-which(colnames(reads.fam_ra) == "ID_afri")]
    
    mytotal.fam_ra[,c("haz", "AATlevel2", "CALwet_cat","crp2_cat","CitrulineM_cat", "score_socio_eco_bycountry", "instruc_mere_3", "eau_traitee", "hemoglobine2_seuil")] <- lapply(mytotal.fam_ra[,c("haz", "AATlevel2", "CALwet_cat","crp2_cat","CitrulineM_cat", "score_socio_eco_bycountry", "instruc_mere_3", "eau_traitee", "hemoglobine2_seuil")], ordered)
    
    cfa.model7 <- '
    Microbiome =~ Lachnospiraceae + Enterobacteriaceae + Streptococcaceae + Lactobacillaceae ' 
    fitcfa7 <- lavaan::cfa(cfa.model7, data = mytotal.fam_ra)
    lavaan::summary(fitcfa7, fit.measures = TRUE, standardized=TRUE)
    
    # Try with relative abundance data
    mytotal_ra = mytotal[,-c(which(colnames(mytotal) %in% colnames(reads_ra)))]
    reads_ra$ID_afri = rownames(reads_ra)
    mytotal_ra = inner_join(mytotal_ra, reads_ra, by="ID_afri")
    reads_ra = reads_ra[,-which(colnames(reads_ra) == "ID_afri")]
    
    cfa.model8 <- '
    Microbiome =~ g_Prevotella_9  +  Prevotella_copri_DSM_18205 + g_Faecalibacterium + Escherichia_coli + g_Alloprevotella + Bacteroides_xylanisolvens + g_Mycoplasmataceae + g_Eubacterium_rectale_group + g_Prevotella_2 + Prevotellaceae_bacterium_DJF_LS10 + g_Succinivibrio + g_Subdoligranulum + g_Phascolarctobacterium + g_Bifidobacterium + g_Megasphaera + g_Roseburia + Veillonella_sp_ICM51a + Prevotella_sp_Marseille_P2439 + Streptococcus_salivarius_subsp_thermophilus + g_Ruminococcaceae_UCG_002
    SES  =~   score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
    Neurodevelopment =~ Comm_tot + PS_tot + PES_tot + FM_tot + GM_tot ' 
    fitcfa8 <- lavaan::cfa(cfa.model8, data = mytotal_ra)
    lavaan::summary(fitcfa8, fit.measures = TRUE, standardized=TRUE)
    
    cfa.model9 <- '
    SES  =~   score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
    Neurodevelopment =~ Comm_tot + PS_tot + PES_tot + FM_tot + GM_tot ' 
    fitcfa9 <- lavaan::cfa(cfa.model9, data = mytotal)
    lavaan::summary(fitcfa9, fit.measures = TRUE, standardized=TRUE) # very good fit
    
    cfa.model10 <- '
    SES  =~   score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
    Neurodevelopment =~ Comm_tot + PS_tot + PES_tot + FM_tot + GM_tot
    BCAA  =~ CitrulineM + ValineM + IsoleucineM + LeucineM' 
    fitcfa10 <- lavaan::cfa(cfa.model10, data = mytotal)
    lavaan::summary(fitcfa10, fit.measures = TRUE, standardized=TRUE) # very good fit
    
    
    fitcfa1gof = c("Neurodevelopment"
                  ,lavaan::summary(fitcfa1, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                  ,as.vector(fitMeasures(fitcfa1, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
    
    fitcfa5gof = c("SES"
                  ,lavaan::summary(fitcfa5, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                  ,as.vector(fitMeasures(fitcfa5, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
    
    fitcfa6gof = c("BCAA"
                   ,lavaan::summary(fitcfa6, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                   ,as.vector(fitMeasures(fitcfa6, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
    
    fitcfa9gof = c("Neurodevelopment + SES"
                  ,lavaan::summary(fitcfa9, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                  ,as.vector(fitMeasures(fitcfa9, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
    
    fitcfa10gof = c("Neurodevelopment + SES + BCAA"
                   ,lavaan::summary(fitcfa10, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                   ,as.vector(fitMeasures(fitcfa10, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
    
    cfagof = as.data.frame(rbind(fitcfa1gof,fitcfa5gof,fitcfa6gof,fitcfa9gof,fitcfa10gof))
    colnames(cfagof) = c('model','estimator','npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')
    rm(fitcfa1gof,fitcfa5gof,fitcfa6gof,fitcfa9gof,fitcfa10gof)
    write.csv(cfagof, file = "CFAs_Goodness_of_fit.csv")
    
    
    rm(cfa.model1,fitcfa1,cfa.model2,fitcfa2,cfa.model3,fitcfa3,cfa.model4,fitcfa4,cfa.model5,fitcfa5,cfa.model6,fitcfa6,cfa.model7,fitcfa7,cfa.model8,fitcfa8,cfa.model9,fitcfa9,fitcfa10,cfa.model10,cfagof)
    
    
    
    # 3) SEM model---------
    
        # a) the simplest model -------
        # Model specification
        mytotal_ra$haz <- factor( mytotal_ra$haz , ordered = TRUE ) #endogenous variables should be ordered
    
        model1.0 <- '
        # measurement model
        Neurodevelopment =~ Comm_tot + PS_tot + PES_tot + FM_tot + GM_tot 
        SES =~ score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
        # regressions
        Neurodevelopment ~ a*haz_cont + b*Streptococcaceae + c*SES
        '
        fit1.0 <- lavaan::sem(model1.0, data=mytotal.fam_ra)
        lavaan::summary(fit1.0, fit.measures = TRUE, standardized=TRUE) # was also tried with enterobacteriaceae, lactobacillaceae and Lachnospiraceae but didn't give any significant results
        
        # Look at model fit
        fitMeasures(fit1.0, c('chisq','df','pvalue','cfi','rmsea','srmr','AIC'))
        # chi square (Model Test User Model) measures the discrepancy between the observed correlations and those implied by the model. For this test we want a NS result.
        # the df is equal to the number of variances + covariances minus the number of parameters estimated. In a just-identified situation, you get no fit statistics because df = 0. It is better to have an overidentified model.
        # lavaan also provides a built-in chisquare test which compares the current model to a model with all paths = 0, similar to the likelihood ratio test where the comparison model is an intercept only model. For that test we want a significant result.
        # Comparative Fit Index (CFI) compares the fitted model to a null model that assumes no relationship among the measured items (complete independence). It ranges from 0 to 1 and CFI values larger than .9 are typically desired.
        # Tucker-Lewis fit index. Both CFI and TLI of these are sensitive to the average correlations, to the point that adding meaningful covariates could actually worsen fit if they aren’t notably correlated with other variables in the model. 
        # They also may improve simply by having a more complex model, leading rise to ‘parsimony-adjusted’ versions. You only need to report one, as they are generally notably correlated. Given that TLI can actually fall outside the 0-1 range, you might as well prefer CFI,
        # RMSEA is a measure that also centers on the model-implied vs. sample covariance matrix.look for values less than .05
        # Lavaan also provides a one-sided test that the RMSEA is ≤.05, or ‘test of close fit’, for which the corresponding p-value ideally would be high. The confidence interval is enough for reporting purposes, and an upper bound beyond .10 may indicate a problem.
        # the SRMR is the mean absolute correlation residual, i.e. the difference between the observed and model-implied correlations.
        # Historical suggestions are to also look for values less than .05, but it is better to simply inspect the residuals and note where there are large discrepancies
        residuals(fit1.0, type = 'cor')$res.cov
        # Beware: there is little correlation between fit indices and the type of misspecification, and predictive ability and fit are essentially independent. with SEM the primary focus should be on the parameter estimates and their corresponding uncertainty.
        
        # Model comparison
        # AIC is a good way to compare models, where a penalty is imposed on the likelihood based on the number of parameters estimated, i.e. model complexity. The value by itself is not useful, but the ‘better’ model will have a lower value.
        # Even if not nested, models must have the same variables for them to be compared, which may require you to constrain certain paths to be zero. Also, AIC can be used to compare just-identified models to over-identified ones.
        # BIC is similar but has a Bayesian interpretation, then AIC is usually preferable. BIC has a different penalty than AIC and is not a measure based on predictive performance, which is what we typically want in a model selection
        # semTools package for easy side by side comparison via the compareFit function
        
        fit1.0gof = c("Simplest model","Streptococcaceae"
                      ,lavaan::summary(fit1.0, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit1.0, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit1.0res = summary(fit1.0, fit.measures = TRUE, standardized=TRUE)$pe
        fit1.0res$model = "Simplest model"
        fit1.0res$microbiota = "Streptococcaceae" 
        
        model1.1 <- '
        # measurement model
        Neurodevelopment =~ Comm_tot + PS_tot + PES_tot + FM_tot + GM_tot 
        SES =~ score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
        # regressions
        Neurodevelopment ~ a*haz_cont + b*PC1 + c*SES
        '
        fit1.1 <- lavaan::sem(model1.1, data=mytotal)
        lavaan::summary(fit1.1, fit.measures = TRUE, standardized=TRUE) 
        # summary(fit1.1, fit.measures = TRUE, standardized=TRUE)$optim$estimator
        fitMeasures(fit1.1, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        fit1.1gof = c("Simplest model","PC1"
                      ,lavaan::summary(fit1.1, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit1.1, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit1.1res = summary(fit1.1, fit.measures = TRUE, standardized=TRUE)$pe
        fit1.1res$model = "Simplest model"
        fit1.1res$microbiota = "PC1"
        
        model1.2 <- '
        # measurement model
        Neurodevelopment =~ Comm_tot + PS_tot + PES_tot + FM_tot + GM_tot 
        SES =~ score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
        # regressions
        Neurodevelopment ~ a*haz_cont + b*shannon + c*SES
        '
        fit1.2 <- lavaan::sem(model1.2, data=mytotal)
        lavaan::summary(fit1.2, fit.measures = TRUE, standardized=TRUE) 
        fitMeasures(fit1.2, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        fit1.2gof = c("Simplest model","shannon"
                      ,lavaan::summary(fit1.2, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit1.2, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit1.2res = summary(fit1.2, fit.measures = TRUE, standardized=TRUE)$pe
        fit1.2res$model = "Simplest model"
        fit1.2res$microbiota = "shannon"
        
        
        mytotal$clusters <- factor( mytotal$clusters , ordered = FALSE ) # exogenous variables should be unordered
        model1.3 <- '
        # measurement model
        Neurodevelopment =~ Comm_tot + PS_tot + PES_tot + FM_tot + GM_tot 
        SES =~ score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
        # regressions
        Neurodevelopment ~ a*haz_cont + b*clusters + c*SES
        '
        fit1.3 <- lavaan::sem(model1.3, data=mytotal)
        lavaan::summary(fit1.3, fit.measures = TRUE, standardized=TRUE) 
        fitMeasures(fit1.3, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        fit1.3gof = c("Simplest model","clusters"
                      ,lavaan::summary(fit1.3, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit1.3, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit1.3res = summary(fit1.3, fit.measures = TRUE, standardized=TRUE)$pe
        fit1.3res$model = "Simplest model"
        fit1.3res$microbiota = "clusters"
        
        fit1res = as.data.frame(rbind(fit1.1res,fit1.2res,fit1.3res,fit1.0res))
        #fit1res = fit1res[which(fit1res$op == "~"),]
        rm(fit1.1res,fit1.2res,fit1.3res,fit1.0res)
        
        
        
        # b) the simple model -------
        # Model specification
        model2.0 <- '
        # measurement model
        Neurodevelopment =~ PS_tot + PES_tot + FM_tot + GM_tot + Comm_tot
        SES =~ score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
        # regressions
        Neurodevelopment ~  a*haz_cont + b*Streptococcaceae + e*SES
        haz_cont ~ c*SES
        Streptococcaceae ~ d*SES
        # defined parameters
        Indirect_H := a*c
        Indirect_M := b*d
        Total := a*c + b*d + e
        '
        fit2.0 <- lavaan::sem(model2.0, data=mytotal.fam_ra)
        lavaan::summary(fit2.0, fit.measures = TRUE, standardized=TRUE)
        fitMeasures(fit2.0, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        # many parameters to estimate
        
        fit2.0gof = c("Simple model","Streptococcaceae"
                      ,lavaan::summary(fit2.0, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit2.0, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit2.0res = summary(fit2.0, fit.measures = TRUE, standardized=TRUE)$pe
        fit2.0res$model = "Simple model"
        fit2.0res$microbiota = "Streptococcaceae" 
        
        model2.1 <- '
        # measurement model
        Neurodevelopment =~ Comm_tot + PS_tot + PES_tot + FM_tot + GM_tot 
        SES =~ score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
        # regressions
        Neurodevelopment ~  a*haz_cont + b*PC1 + e*SES
        haz_cont ~ c*SES
        PC1 ~ d*SES
        # defined parameters
        Indirect_H := a*c
        Indirect_M := b*d
        Total := a*c + b*d + e
        '
        fit2.1 <- lavaan::sem(model2.1, data=mytotal)
        lavaan::summary(fit2.1, fit.measures = TRUE, standardized=TRUE) 
        fitMeasures(fit2.1, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        fit2.1gof = c("Simple model","PC1"
                      ,lavaan::summary(fit2.1, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit2.1, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit2.1res = summary(fit2.1, fit.measures = TRUE, standardized=TRUE)$pe
        fit2.1res$model = "Simple model"
        fit2.1res$microbiota = "PC1"
        
        model2.2 <- '
        # measurement model
        Neurodevelopment =~ Comm_tot + PS_tot + PES_tot + FM_tot + GM_tot 
        SES =~ score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
        # regressions
        Neurodevelopment ~  a*haz_cont + b*shannon +e*SES
        haz_cont ~ c*SES
        shannon ~ d*SES
        # defined parameters
        Indirect_H := a*c
        Indirect_M := b*d
        Total := a*c + b*d + e
        '
        fit2.2 <- lavaan::sem(model2.2, data=mytotal)
        lavaan::summary(fit2.2, fit.measures = TRUE, standardized=TRUE) 
        fitMeasures(fit2.2, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        fit2.2gof = c("Simple model","shannon"
                      ,lavaan::summary(fit2.2, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit2.2, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit2.2res = summary(fit2.2, fit.measures = TRUE, standardized=TRUE)$pe
        fit2.2res$model = "Simple model"
        fit2.2res$microbiota = "shannon"
        
        
        mytotal$clusters <- factor( mytotal$clusters , ordered = TRUE ) #endogenous variables should be ordered
        model2.3 <- '
        # measurement model
        Neurodevelopment =~ Comm_tot + PS_tot + PES_tot + FM_tot + GM_tot 
        SES =~ score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
        # regressions
        Neurodevelopment ~  a*haz_cont + b*clusters + e*SES
        haz_cont ~ c*SES
        clusters ~ d*SES
        # defined parameters
        Indirect_H := a*c
        Indirect_M := b*d
        Total := a*c + b*d + e
        '
        fit2.3 <- lavaan::sem(model2.3, data=mytotal)
        lavaan::summary(fit2.3, fit.measures = TRUE, standardized=TRUE) 
        fitMeasures(fit2.3, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        fit2.3gof = c("Simple model","clusters"
                      ,lavaan::summary(fit2.3, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit2.3, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit2.3res = summary(fit2.3, fit.measures = TRUE, standardized=TRUE)$pe
        fit2.3res$model = "Simple model"
        fit2.3res$microbiota = "clusters"
        
        fit2res = as.data.frame(rbind(fit2.1res,fit2.2res,fit2.3res,fit2.0res))
        #fit2res = fit2res[which(fit2res$op == "~"),]
        rm(fit2.1res,fit2.2res,fit2.3res,fit2.0res)
        
        
        
        # c) The complex SEM model ----------
        mytotal_ra$AATlevel2 <- factor( mytotal_ra$AATlevel2 , ordered = FALSE ) #exogenous variables should be unordered
        mytotal$AATlevel2 <- factor( mytotal$AATlevel2 , ordered = FALSE ) #exogenous variables should be unordered
        mytotal$taille_naiss_p3 <- factor( mytotal$taille_naiss_p2 , levels = c("1","2","3"), labels = c("1","2","2"), ordered = FALSE ) #exogenous variables should be unordered
        mytotal_ra$taille_naiss_p3 <- factor( mytotal_ra$taille_naiss_p2 , levels = c("1","2","3"), labels = c("1","2","2"), ordered = FALSE ) #exogenous variables should be unordered
        mytotal.fam_ra$taille_naiss_p3 <- factor( mytotal.fam_ra$taille_naiss_p2 , levels = c("1","2","3"), labels = c("1","2","2"), ordered = FALSE ) #exogenous variables should be unordered
        
        #Model Specification
        model3.0 <- '
        #measurement model
        Neurodevelopment =~ PS_tot + PES_tot + FM_tot + GM_tot + Comm_tot
        SES =~ score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
        BCAA  =~ CitrulineM + ValineM + IsoleucineM + LeucineM
        # regressions
        haz_cont ~ hemoglobine2 + d*Streptococcaceae + c*SES + BCAA
        Neurodevelopment ~ e*SES + a*haz_cont + b*Streptococcaceae + AATlevel2 + age + taille_naiss_p3
        Streptococcaceae ~ BCAA
        #residual correlations
        # defined parameters
        Indirect_M := a*d
        Total_M := b + a*d
        Indirect_SES := a*c
        Total_SES := a*c + e
        '
        fit3.0 <- lavaan::sem(model3.0, data=mytotal.fam_ra)
        summary(fit3.0, fit.measures = TRUE, standardized=TRUE)
        fitMeasures(fit3.0, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        
        fit3.0gof = c("Complex SEM","Streptococcaceae"
                      ,lavaan::summary(fit3.0, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit3.0, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit3.0res = summary(fit3.0, fit.measures = TRUE, standardized=TRUE)$pe
        fit3.0res$model = "Complex SEM"
        fit3.0res$microbiota = "Streptococcaceae"
        
        model3.1 <- '
        #measurement model
        Neurodevelopment =~ PS_tot + PES_tot + FM_tot + GM_tot + Comm_tot
        SES =~ score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
        BCAA  =~ CitrulineM + ValineM + IsoleucineM + LeucineM
        # regressions
        haz_cont ~ hemoglobine2 + d*PC1 + c*SES + BCAA
        Neurodevelopment ~ e*SES + a*haz_cont + b*PC1 + AATlevel2 + age + taille_naiss_p3
        PC1 ~ BCAA
        #residual correlations
        # defined parameters
        Indirect_M := a*d
        Total_M := b + a*d
        Indirect_SES := a*c
        Total_SES := a*c + e
        '
        fit3.1 <- lavaan::sem(model3.1, data=mytotal)
        summary(fit3.1, fit.measures = TRUE, standardized=TRUE)
        fitMeasures(fit3.1, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        fit3.1gof = c("Complex SEM","PC1"
                      ,lavaan::summary(fit3.1, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit3.1, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit3.1res = summary(fit3.1, fit.measures = TRUE, standardized=TRUE)$pe
        fit3.1res$model = "Complex SEM"
        fit3.1res$microbiota = "PC1"
        
        model3.2 <- '
        #measurement model
        Neurodevelopment =~ PS_tot + PES_tot + FM_tot + GM_tot + Comm_tot
        SES =~ score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
        BCAA  =~ CitrulineM + ValineM + IsoleucineM + LeucineM
        # regressions
        haz_cont ~ hemoglobine2 + d*shannon + c*SES + BCAA
        Neurodevelopment ~ e*SES + a*haz_cont + b*shannon + AATlevel2 + age + taille_naiss_p3
        shannon ~ BCAA
        #residual correlations
        # defined parameters
        Indirect_M := a*d
        Total_M := b + a*d
        Indirect_SES := a*c
        Total_SES := a*c + e
        '
        fit3.2 <- lavaan::sem(model3.2, data=mytotal)
        summary(fit3.2, fit.measures = TRUE, standardized=TRUE)
        fitMeasures(fit3.2, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        fit3.2gof = c("Complex SEM","shannon"
                      ,lavaan::summary(fit3.2, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit3.2, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit3.2res = summary(fit3.2, fit.measures = TRUE, standardized=TRUE)$pe
        fit3.2res$model = "Complex SEM"
        fit3.2res$microbiota = "shannon"
        
        model3.3 <- '
        #measurement model
        Neurodevelopment =~ PS_tot + PES_tot + FM_tot + GM_tot + Comm_tot
        SES =~ score_socio_eco_bycountry + instruc_mere_3 + nb_pieces + eau_traitee + age_prem_gross
        BCAA  =~ CitrulineM + ValineM + IsoleucineM + LeucineM
        # regressions
        haz_cont ~ hemoglobine2 + d*clusters + c*SES + BCAA
        Neurodevelopment ~ e*SES + a*haz_cont + b*clusters + AATlevel2 + age + taille_naiss_p3
        clusters ~ BCAA
        #residual correlations
        # defined parameters
        Indirect_M := a*d
        Total_M := b + a*d
        Indirect_SES := a*c
        Total_SES := a*c + e
       '
        fit3.3 <- lavaan::sem(model3.3, data=mytotal)
        summary(fit3.3, fit.measures = TRUE, standardized=TRUE)
        fitMeasures(fit3.3, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        fit3.3gof = c("Complex SEM","clusters"
                      ,lavaan::summary(fit3.3, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit3.3, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit3.3res = summary(fit3.3, fit.measures = TRUE, standardized=TRUE)$pe
        fit3.3res$model = "Complex SEM"
        fit3.3res$microbiota = "clusters"
        
        fit3res = as.data.frame(rbind(fit3.1res,fit3.2res,fit3.3res,fit3.0res))
        rm(fit3.1res,fit3.2res,fit3.3res,fit3.0res)
        
        
        
        # d) The complex model, alternative 2: path analysis -----------
        
        # Replace latent variables with observed variables
          # SES
          summary(lm(Overal_score ~ score_socio_eco_bycountry, data = mytotal))
          summary(lm(Overal_score ~ instruc_mere_3, data = mytotal)) # S for lower category
          summary(lm(Overal_score ~ nb_pieces, data = mytotal)) # very significant
          summary(lm(Overal_score ~ eau_traitee, data = mytotal))
          summary(lm(Overal_score ~ age_prem_gross, data = mytotal)) # S
          #-- plot
          ggplot( mydata, aes(y = Overal_score, x = nb_pieces , col = haz )) +
            geom_point() +
            geom_smooth(method='lm')
          ggplot( mydata, aes(y = Overal_score, x = age_prem_gross , col = haz )) +
            geom_point() +
            geom_smooth(method='lm')
          # keep nb_pieces and age_prem_gross

        
        # Model Specification
        mytotal$AATlevel2 <- factor( mytotal$AATlevel2 , ordered = FALSE ) #exogenous variables should be unordered
        mytotal.fam_ra$AATlevel2 <- factor( mytotal.fam_ra$AATlevel2 , ordered = FALSE ) #exogenous variables should be unordered
        mytotal$hemoglobine2_seuil <- factor(mytotal$hemoglobine2_seuil, ordered = FALSE)
        
        model4.0 <- '
        # measurement model
        # regressions
        haz_cont ~ e*age_prem_gross + c*nb_pieces + hemoglobine2 + d*Streptococcaceae + LeucineM
        Streptococcaceae ~ LeucineM
        Overal_score ~ f*nb_pieces + g*age_prem_gross + a*haz_cont + b*Streptococcaceae + AATlevel2 + age + taille_naiss_p3
        # residual correlations
        # defined parameters
        Indirect_M := a*d
        Total_M := b + a*d
        Indirect_SES := a*c
        Total_SES := f + a*c
        Indirect_Mat := a*e
        Total_Mat := g + a*e
        '
        fit4.0 <- lavaan::sem(model4.0, data=mytotal.fam_ra, estimator = "MLR")
        summary(fit4.0, fit.measures = TRUE, standardized=TRUE)
        
        fit4.0gof = c("Complex path analysis","Streptococcaceae"
                      ,lavaan::summary(fit4.0, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit4.0, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit4.0res = summary(fit4.0, fit.measures = TRUE, standardized=TRUE)$pe
        fit4.0res$model = "Complex path analysis"
        fit4.0res$microbiota = "Streptococcaceae"
        
        model4.1 <- '
        # measurement model
        # regressions
        haz_cont ~ e*age_prem_gross + c*nb_pieces + hemoglobine2 + d*PC1 + LeucineM
        PC1 ~ LeucineM
        Overal_score ~ f*nb_pieces + g*age_prem_gross + a*haz_cont + b*PC1 + AATlevel2 + age + taille_naiss_p3
        # residual correlations
        # defined parameters
        Indirect_M := a*d
        Total_M := b + a*d
        Indirect_SES := a*c
        Total_SES := f + a*c
        Indirect_Mat := a*e
        Total_Mat := g + a*e
        '
        fit4.1 <- lavaan::sem(model4.1, data=mytotal, estimator = "MLR")
        summary(fit4.1, fit.measures = TRUE, standardized=TRUE)
        fitMeasures(fit4.1, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        # horrible fit !!!!
        fit4.1gof = c("Complex path analysis","PC1"
                      ,lavaan::summary(fit4.1, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit4.1, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit4.1res = summary(fit4.1, fit.measures = TRUE, standardized=TRUE)$pe
        fit4.1res$model = "Complex path analysis"
        fit4.1res$microbiota = "PC1"
        
        
        model4.2 <- '
        # measurement model
        # regressions
        haz_cont ~ e*age_prem_gross + c*nb_pieces + hemoglobine2 + d*shannon + LeucineM
        shannon ~ LeucineM
        Overal_score ~ f*nb_pieces + g*age_prem_gross + a*haz_cont + b*shannon + AATlevel2 + age + taille_naiss_p3
        # residual correlations
        # defined parameters
        Indirect_M := a*d
        Total_M := b + a*d
        Indirect_SES := a*c
        Total_SES := f + a*c
        Indirect_Mat := a*e
        Total_Mat := g + a*e
        '
        fit4.2 <- lavaan::sem(model4.2, data=mytotal, estimator = "MLR")
        summary(fit4.2, fit.measures = TRUE, standardized=TRUE)
        fitMeasures(fit4.2, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        # horrible fit !!!!
        fit4.2gof = c("Complex path analysis","shannon"
                      ,lavaan::summary(fit4.2, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit4.2, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit4.2res = summary(fit4.2, fit.measures = TRUE, standardized=TRUE)$pe
        fit4.2res$model = "Complex path analysis"
        fit4.2res$microbiota = "shannon"
        
        model4.3 <- '
        # measurement model
        # regressions
        haz_cont ~ e*age_prem_gross + c*nb_pieces + hemoglobine2 + d*clusters + LeucineM
        clusters ~ LeucineM
        Overal_score ~ f*nb_pieces + g*age_prem_gross + a*haz_cont + b*clusters + AATlevel2 + age + taille_naiss_p3
        # residual correlations
        # defined parameters
        Indirect_M := a*d
        Total_M := b + a*d
        Indirect_SES := a*c
        Total_SES := f + a*c
        Indirect_Mat := a*e
        Total_Mat := g + a*e
        '
        fit4.3 <- lavaan::sem(model4.3, data=mytotal)
        summary(fit4.3, fit.measures = TRUE, standardized=TRUE)
        fitMeasures(fit4.3, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled'))
        # horrible fit !!!!
        fit4.3gof = c("Complex path analysis","clusters"
                      ,lavaan::summary(fit4.3, fit.measures = TRUE, standardized=TRUE)$optim$estimator
                      ,as.vector(fitMeasures(fit4.3, c('npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')) ))
        
        fit4.3res = summary(fit4.3, fit.measures = TRUE, standardized=TRUE)$pe
        fit4.3res$model = "Complex path analysis"
        fit4.3res$microbiota = "clusters"
        
        fit4res = as.data.frame(rbind(fit4.1res,fit4.2res,fit4.3res,fit4.0res))
        # fit4res = fit4res[which(fit4res$op == "~"),]
        rm(fit4.1res,fit4.2res,fit4.3res,fit4.0res)
        
        rm(fit1.0, fit1.1,fit1.2,fit1.3,fit2.0,fit2.1,fit2.2,fit2.3,fit3.0,fit3.1,fit3.2,fit3.3,fit4.0,fit4.1,fit4.2,fit4.3)
        rm(model1.0,model1.1,model1.2,model1.3,model2.0,model2.1,model2.2,model2.3,model3.0,model3.1,model3.2,model3.3,model4.0,model4.1,model4.2,model4.3)
        
        gof = as.data.frame(rbind(fit1.1gof,fit1.2gof,fit1.3gof,fit1.0gof,fit2.1gof,fit2.2gof,fit2.3gof,fit2.0gof,fit3.1gof,fit3.2gof,fit3.3gof,fit3.0gof,fit4.1gof,fit4.2gof,fit4.3gof,fit4.0gof))
        colnames(gof) = c('model','microbiota','estimator','npar','chisq','df','pvalue','baseline.pvalue','cfi','tli','rmsea','rmsea.pvalue','srmr','chisq.scaled','pvalue.scaled','baseline.pvalue.scaled','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.pvalue.scaled')
        rm(fit1.1gof,fit1.2gof,fit1.3gof,fit2.1gof,fit2.2gof,fit2.3gof,fit3.1gof,fit3.2gof,fit3.3gof,fit4.1gof,fit4.2gof,fit4.3gof)
        write.csv(gof, file = "SEMs_Goodness_of_fit.csv")
        
        fitres = as.data.frame(rbind(fit1res,fit2res,fit3res,fit4res))
        fitres = fitres[,c((ncol(fitres)-1),ncol(fitres),1:(ncol(fitres)-2))] # reorder columns
        fitres.indirect = fitres[which(fitres$op == "~" | fitres$op == ":="),]
        fitres = fitres[which(fitres$op == "~"),]
        rm(fit1res,fit2res,fit3res,fit4res)
        write.csv(fitres, file = "SEMs_parameters_estimates.csv")
        
        # graphical representation of results
        is.numeric(fitres$est)
        is.numeric(fitres$se)
        is.numeric(fitres$pvalue)
        is.numeric(fitres$std.all)
        # remove non-significant results
        #df1 = df[which(df$out_pvalue <= 0.05),]
        # remove lhs != Neurodevelopment
        fitres1 = fitres[which(fitres$lhs == "Neurodevelopment" | fitres$lhs == "Overal_score"),]
        fitres1 = fitres1[which(fitres1$microbiota == "shannon"),]
        fitres1$model = factor(fitres1$model, levels = c("Simplest model","Simple model","Complex SEM","Complex path analysis"))
        
        pdf('Results_sem_parameters_neurodevelopment_shannon.pdf', height = 12, width = 10) 				# Write to PDF
        ggplot(fitres1, aes(x=rhs, y=est, colour=ifelse(pvalue<0.05,"red","black"))) +
          #geom_point(stat="identity", aes(color = ifelse(pvalue<0.05,"red","black")), position=position_dodge(width =1)) +
          #geom_errorbar(aes(ymin=est - 1.96*se, ymax=est + 1.96*se, colour=ifelse(pvalue<0.05,"red","black")), width=.1, position=position_dodge(width=1)) +
          geom_pointrange(aes(ymin=est - 1.96*se, ymax=est + 1.96*se ), position=position_dodge(0.9), size = 0.4) + 
          coord_flip() + 
          facet_grid(. ~  model , scales="free_x", switch = "x") + 
          geom_hline(yintercept = 0) +
          ylab("Parameter estimate") +
          xlab("") +
          theme_bw() +						# remove grey background
          theme(legend.title=element_blank(), axis.text.x = element_text(size=10, angle = 0)
                ,legend.position="bottom", panel.grid.major = element_blank()  # remove x and y major grid lines
                ,strip.text.x = element_text(size = 8), axis.text.y = element_text(angle = 0)
          ) +
          scale_color_identity() +
          scale_y_continuous(position="right") +
          scale_x_discrete(limits = c("haz_cont","SES","nb_pieces","age_prem_gross","age","taille_naiss_p3","AATlevel2","shannon"),
                           labels=c( "Stunting (cont.)","Socioeconomic status","Number of rooms","Maternal age at first pregnancy","Age","Reported birth size","AAT level (cat.)","Microbiome shannon diversity") )
        dev.off()
        
        fitres2 = fitres[which(fitres$lhs == "Neurodevelopment" | fitres$lhs == "Overal_score"),]
        fitres2 = fitres2[which(fitres2$microbiota == "shannon" | fitres2$microbiota == "PC1" | fitres2$microbiota == "clusters" | fitres2$microbiota == "Streptococcaceae"),]
        fitres2$model = factor(fitres2$model, levels = c("Simplest model","Simple model","Complex SEM","Complex path analysis"))
        
        
        pdf('Results_sem_parameters_neurodevelopment.pdf', height = 12, width = 10) 				# Write to PDF
        ggplot(fitres2, aes(x=rhs, y=est, colour=ifelse(pvalue<0.05,"red","black"), group = microbiota, shape = microbiota)) +
          geom_pointrange(aes(ymin=est - 1.96*se, ymax=est + 1.96*se ), position=position_dodge(0.9), size = 0.4) +
          coord_flip() +
          facet_grid(. ~  model , scales="free_x", switch = "x") +
          geom_hline(yintercept = 0) +
          ylab("Parameter estimate") +
          xlab("") +
          theme_bw() +						# remove grey background
          theme(legend.title=element_blank(), axis.text.x = element_text(size=10, angle = 0)
                ,legend.position="bottom", panel.grid.major = element_blank()  # remove x and y major grid lines
                ,strip.text.x = element_text(size = 8), axis.text.y = element_text(angle = 0)
          ) +
          scale_color_identity() +
          scale_y_continuous(position="right") +
          scale_x_discrete(limits = c("haz_cont","SES","nb_pieces","age_prem_gross","age","taille_naiss_p3","AATlevel2","shannon","PC1","clusters","Streptococcaceae"),
                           labels=c( "Stunting (cont.)","Socioeconomic status","Number of rooms","Maternal age at first pregnancy","Age","Reported birth size","AAT level (cat.)","Microbiota shannon diversity", "Microbiota first principal component", "Microbiota clusters", "Streptococcaceae") ) 
       dev.off()

       
        write.csv(fitres.indirect, file = "SEM_indirect_effects.csv")
        
        