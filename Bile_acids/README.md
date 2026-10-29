READ ME
-
R scripts used to analyse the data and produce the figures from: Lithocholic acid modulates the growth of butyrate-producing bacteria and is decreased in the feces of stunted children

DOI: https://doi.org/10.64898/2026.05.04.722639

Datasets:
-
16s rRNA gene amplicon sequencing data of the Afribtiota are available on ENA under project number: PRJEB48119 

16S rRNA gene amplicon sequencing data of the stool-derived in vitro communities are available on ENA under project number: PRJNA1468128

Workflow:
- 
## Data cleaning and DADA2

Demultiplexed raw 16S rRNA gene amplicon sequencing data were processed using the DADA2 pipeline. Parameters and database versions are described in the method section.

The fecal 16S rRNA gene amplicon sequencing raw data from the Afribiota project was previously processed using the script: *Afribiota/Dada2_16SPNASVonaeschetal2022.R* (DOI: https://doi.org/10.1073/pnas.220958911 ).

The SICs raw data were processed in DADA2 using the script: *DADA2_SICs.Rmd*

The resulting ASVs and taxonomy tables, as well as the metadata associated to these datasets are available as supplementary tables:

  > Afribiota data: Supplementary Data S1-S3

  > Stool-derived in vitro community: Supplementary Data S11-S13
  
  > Stool samples used to produce SICs: Supplementary Data S8-S9 
  
The bile acid abundance per samples is included in the metadata table in Supplementary Data S3.

## Analysis

Analysis of the bile acid profiles :

  > Analysis script: *BileAcidBatchEffects.Rmd*, *BileAcidCode_CleanedCombined.Rmd*

Analysis of the microbiome profiles from the Afribiota study in correlation with the bile acid profiles:

  > Analysis script: *A_BA_Afribiota.Rmd*

Analysis of the stool-derived in vitro communities (SICs) with additional bile acids: 

  > Analysis script: *BA_SIC_ms.Rmd*

Analysis of the bacterial isolates exposed to bile acids:

  > Analysis script: *BA_Isolates_ms.Rmd*

Metabolomic analysis:

  > Analysis scripts: *1_metabolomic_analysis_parent_data.R*, *2_metabolomic_analysis_metabolit_data.R*, *3_main_figures.R*
  




 
