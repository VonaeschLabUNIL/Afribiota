### Analysis of 16S study Afribiota June 2021 #########
#### Analysis performed by Pascale Vonaesch and code simplified by Munir Winkel ##########
#### General data import and cleaning ####
##### install all necessary packages ####

library(ggplot2)
library(magrittr)
library(biomformat)
library(lme4)
library(microbiome)
library(factoextra)
library(vegan)

##### set the work directory and set working environment ####
rm(list=ls(all=TRUE)) # go get a clean workspace
setwd("/Users/admin/Desktop/Papers and manuscripts/under revision/Lipid: 16S paper/testlipidpaper copy")
source("useful_functions.R")


#### import data ####
tax = read.delim("/Users/admin/Desktop/Papers and manuscripts/under revision/Lipid: 16S paper/Submission to PNAS/Supplementary Files/Supplementary File S3.csv", header = TRUE, sep = ",")
colnames(tax)[1] = "ASV"
which(duplicated(row.names(tax))==TRUE) # there are no duplicates in term of all variables
dim(tax) #6178 ASVs and 10 taxonomic levels overall
row.names(tax) = tax[, 1]
View(head(tax))


taxonomy = phyloseq::tax_table(tax)
# #View(taxonomy)
row.names(taxonomy) = tax[, 1]
colnames(taxonomy) = colnames(tax)

ASV = read.delim("/Users/admin/Desktop/Papers and manuscripts/under revision/Lipid: 16S paper/Submission to PNAS/Supplementary Files/Supplementary File S2.csv", header = TRUE, sep = ";")
# #View(ASV)
colnames(ASV)[1] = "ASV"
which(duplicated(row.names(ASV))==TRUE) # there are no duplicates
which(duplicated(colnames(ASV))==TRUE) # there are no duplicates
dim(ASV) #6178 ASVs and 1173 samples (including the ASV column)
row.names(ASV) = ASV$ASV
ASV = ASV[, -1] %>% as.matrix()
View(head(ASV))

taxa   = phyloseq::otu_table(ASV, taxa_are_rows = TRUE) %>% as.matrix()
physeq = phyloseq::phyloseq(taxa, taxonomy)

meta   = read.delim("/Users/admin/Desktop/Supplementary File S2.csv", header = TRUE, sep = ",")

dim(meta) #1172 samples and 68 metadata categories
row.names(meta) = meta$id_metag

meta = phyloseq::sample_data(meta)
df   = phyloseq::merge_phyloseq(meta, physeq)

df # check if all makes sense
phyloseq::tax_table(df) = phyloseq::tax_table(df)[, -1]

#this produces a maximum possible prevalence count per species per ASV
df2 = df %>% phyloseq::transform_sample_counts(fun = prune_singlets)

length(which(phyloseq::sample_data(df2)$sampletype=="gastric")) # we start with 248 samples
length(which(phyloseq::sample_data(df2)$sampletype=="duodenal")) # we start with 170 samples
length(which(phyloseq::sample_data(df2)$sampletype=="feces")) # we start with 737 samples
length(which(phyloseq::sample_data(df2)$project=="Control")) # we start with 14 samples
dim( phyloseq::otu_table(df))

df_control = df2 %>% phyloseq::subset_samples(project== "Control")

df_mock  = phyloseq::subset_samples(df_control, sampletype=="control")
df_water = phyloseq::subset_samples(df_control, sampletype=="water")

average = mean(phyloseq::sample_sums(df_water)) # especially the Mada water samples have a lot of sequences!
average ## 35078
average2 = average * 7/9 ## this is the actual average, was we do have 9 samples in analysis
average2 ##27282.89

median = median( phyloseq::sample_sums(df_water))
median ## 23199

average = mean( phyloseq::sample_sums(df_mock))
average ## 52963 more than double the average water sample

#### clean up your metadata and keep only samples with at least 5000 seq ####
#create a variable age in years
phyloseq::sample_data(df2)$ageyears = cut(
  phyloseq::sample_data(df2)$age, c(24,36,48,61)
  , include.lowest = TRUE, right=TRUE, dig.lab = 5, ordered_result = TRUE)

which(is.na(phyloseq::sample_data(df2)$ageyears)) # these are the controls
levels(phyloseq::sample_data(df2)$ageyears)
levels(phyloseq::sample_data(df2)$haz)
levels(phyloseq::sample_data(df2)$pays)

# make labels on variables to better interpret results
phyloseq::sample_data(df2)$ageyears = factor(phyloseq::sample_data(df2)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))
phyloseq::sample_data(df2)$haz  = as.factor(phyloseq::sample_data(df2)$haz)
phyloseq::sample_data(df2)$haz  = factor(phyloseq::sample_data(df2)$haz, labels =c("non-stunted", "moderately stunted", "severely stunted"))

levels(phyloseq::sample_data(df2)$haz)

phyloseq::sample_data(df2)$stunted = ""
phyloseq::sample_data(df2)$stunted[which(phyloseq::sample_data(df2)$haz == "moderately stunted")] <- "stunted"
phyloseq::sample_data(df2)$stunted[which(phyloseq::sample_data(df2)$haz == "severely stunted")] <- "stunted"
phyloseq::sample_data(df2)$stunted[which(phyloseq::sample_data(df2)$haz == "non-stunted")] <- "non-stunted"
phyloseq::sample_data(df2)$stunted<-as.factor(phyloseq::sample_data(df2)$stunted)
levels(phyloseq::sample_data(df2)$stunted)

phyloseq::sample_data(df2)$anemie = factor(phyloseq::sample_data(df2)$anemie, labels =c("Non", "Oui"))

phyloseq::sample_data(df2)$ph_estomac  = as.numeric(phyloseq::sample_data(df2)$ph_estomac)
phyloseq::sample_data(df2)$ph_intestin = as.numeric(phyloseq::sample_data(df2)$ph_intestin)

sampledata = as.data.frame(phyloseq::sample_data(df2))
mean(phyloseq::sample_sums(df2))

plot(sort(phyloseq::sample_sums(df2)), ylim=c(0, 200000))

df_controls = phyloseq::subset_samples(df2, project == "Control")

dff = phyloseq::prune_samples(phyloseq::sample_sums(df2)>=5000, df2)
dim(phyloseq::sample_data(dff)) ## we are left with 1089 samples

length(which(phyloseq::sample_data(dff)$sampletype=="gastric")) # we start with 241 samples
length(which(phyloseq::sample_data(dff)$sampletype=="duodenal")) # we start with 153 samples
length(which(phyloseq::sample_data(dff)$sampletype=="feces")) # we start with 683 samples
length(which(phyloseq::sample_data(dff)$project=="Control")) # we start with 11 samples, water samples are still here

# kick-out the control samples
dff = phyloseq::subset_samples(dff, project != "Control")

#filter your dataset to keep only relevant data
phyloseq::get_taxa_unique(dff, "Rank1") ## only bacteria, archae, unassigned in dataset

# filter out Unassigned, Chloroplast etc.
dff = dff %>%
 phyloseq::subset_taxa(Rank5 != "Mitochondria") %>%
 phyloseq::subset_taxa(Rank4 != "Chloroplast") %>%
 phyloseq::subset_taxa(Rank1 != "Unassigned")
dff # 5479 taxa

# kick out samples that do not match the criteria
phyloseq::sample_data(dff)$ph_estomac  = as.numeric(as.character(phyloseq::sample_data(dff)$ph_estomac))
phyloseq::sample_data(dff)$ph_intestin = as.numeric(as.character(phyloseq::sample_data(dff)$ph_intestin))

dffiltered2 = dff %>%
 phyloseq::subset_samples(
   # to filter out samples which are not recruited in the community
   raison_hospi=="Recrutement communautaire")  %>%
   phyloseq::subset_samples(whz_cont <=2) %>%
  # to keep only samples with no acute undernutrition or obesity
   phyloseq::subset_samples(whz_cont >=-2)

sample_data(dff)


dfa = dffiltered2 %>%
 phyloseq::subset_samples(
          sampletype=="feces" | (
            sampletype == "gastric" & ph_estomac <= 4 & ph_estomac!="")
          )

dfb = phyloseq::subset_samples(dff
          , sampletype == "duodenal" & ph_intestin > 4 & ph_intestin != "")

table(phyloseq::sample_data(dfb)$sampletype)
table(phyloseq::sample_data(dfb)$sampletype,
      phyloseq::sample_data(dfb)$ph_intestin)
table(phyloseq::sample_data(dfa)$sampletype,
      phyloseq::sample_data(dfa)$ph_estomac)

dffiltered2 = phyloseq::merge_phyloseq(dfa, dfb) %>%
  # to keep only children below age 5
 phyloseq::subset_samples(age <= 60) %>%
  # to keep only children above two years
 phyloseq::subset_samples(age >= 24)

dffiltered2 # 924 samples

#track the evolution of your sample numbers
levels( phyloseq::sample_data(dffiltered2)$sampletype)

length( which( phyloseq::sample_data(dffiltered2)$ sampletype == "gastric")) # 168
length( which( phyloseq::sample_data(dffiltered2)$ sampletype == "duodenal")) # 129
length( which( phyloseq::sample_data(dffiltered2)$ sampletype == "feces")) # 627

colnames(
  phyloseq::tax_table(dffiltered2)
         ) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Accession", "ASV_sequence")

# kick out "empty" taxa from taxonomy table
dffiltered2 = phyloseq::filter_taxa(dffiltered2, function(x) mean(x) > 0, TRUE)

#### who is there in the dataset? Analyze and filter accordingly ####
phyloseq::get_taxa_unique(dffiltered2, "Kingdom")
dfArcheae = phyloseq::subset_taxa(dffiltered2, Kingdom=="Archaea")

phyloseq::get_taxa_unique(dfArcheae, "Phylum") # to see how many phylums we have in Euryarchaeota
phyloseq::get_taxa_unique(dfArcheae, "Class") # to see how many classes: Methanobacteria
phyloseq::get_taxa_unique(dfArcheae, "Order") # "Methanobacteriales"
phyloseq::get_taxa_unique(dfArcheae, "Family") # "Methanobacteriaceae"
phyloseq::get_taxa_unique(dfArcheae, "Genus") # "Methanobrevibacter"
phyloseq::get_taxa_unique(dfArcheae, "Species") # "Methanobrevibacter_smithii_DSM_2375"

dfBacteriae = phyloseq::subset_taxa(dffiltered2, Kingdom=="Bacteria")
length( phyloseq::get_taxa_unique(dfBacteriae, "Phylum")) # to see how many phylums: 28 phyla
length( phyloseq::get_taxa_unique(dfBacteriae, "Class")) # to see how many classes: we have 56 classes
length( phyloseq::get_taxa_unique(dfBacteriae, "Order")) # to see how many orders, we have 102 orders
length( phyloseq::get_taxa_unique(dfBacteriae, "Family")) # to see how many families, 175 different families
length( phyloseq::get_taxa_unique(dfBacteriae, "Genus")) # to see how many genera, 471 different genera
length( phyloseq::get_taxa_unique(dfBacteriae, "Species")) # 841 different species

dffiltered3 = dffiltered2

# how frequent are the different taxa? Assess to filter more, if needed
# Compute prevalence of each feature
prevdf = apply(X = phyloseq::otu_table(dffiltered3),
        MARGIN = ifelse(phyloseq::taxa_are_rows(dffiltered3), yes = 1, no = 2),
        FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
          TotalAbundance = phyloseq::taxa_sums(dffiltered3),
          phyloseq::tax_table(dffiltered3))
# #View(prevdf)
# write.csv(prevdf, "prevalencetable.csv")

# Subset to the remaining phyla

prevdf1 = subset(prevdf, Phylum %in% phyloseq::get_taxa_unique(dffiltered3, "Phylum"))

ggplot(prevdf1, aes(TotalAbundance, Prevalence / phyloseq::nsamples(dffiltered3),color=Phylum)) +
 # Include a guess for parameter
 geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
 scale_x_log10() +
  xlab("Total Abundance") +
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) +
  theme(legend.position = "none")
# comment: there are quite a few phyla that are only present in a given dataset. However, the data is also heterogenous, as the samples are from different sources and countries. Apply only a very permissive filtering step here and maybe make more stringent filters on more restricted datasets

## Define prevalence threshold as 1% of total samples
prevalenceThreshold = 0.01 * phyloseq::nsamples(dffiltered3)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa    = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
dffiltered4 = phyloseq::prune_taxa(keepTaxa, dffiltered3)

dffiltered4 # we are now down to 1021 taxa!!

prevdf2 = subset(prevdf, Phylum %in% phyloseq::get_taxa_unique(dffiltered4, "Phylum"))

ggplot(prevdf1, aes(TotalAbundance, Prevalence / phyloseq::nsamples(dffiltered4),color=Phylum)) +
 # Include a guess for parameter
 geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
 scale_x_log10() +
  xlab("Total Abundance") +
  ylab("Prevalence [Frac. Samples]") +
 facet_wrap(~Phylum) +
  theme(legend.position="none")

# comment: there are quite a few phyla that are only present in a given dataset. However, the data is also heterogenous, as the samples are from different sources and countries. Apply only a very permissive filtering step here and maybe make more stringent filters on more restricted datasets


#### make subsets of samples for analysis ####

phyloseq::sample_data(dffiltered4)$read_count = rowSums(t(
  phyloseq::otu_table(dffiltered4)))

table(phyloseq::sample_data(dffiltered4)$pays)
table(phyloseq::sample_data(dffiltered4)$sampletype)
table(phyloseq::sample_data(dffiltered4)$sampletype, phyloseq::sample_data(dffiltered4)$pays)

# make labels on variables to better interpret results
phyloseq::sample_data(dffiltered4)$ageyears = factor(
  phyloseq::sample_data(dffiltered4)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))

phyloseq::sample_data(dffiltered4)$stunted = factor(
  phyloseq::sample_data(dffiltered4)$stunted, labels =c("non-stunted", "stunted"))

phyloseq::sample_data(dffiltered4)$haz = factor(
  phyloseq::sample_data(dffiltered4)$haz, labels =c("non-stunted", "moderately stunted", "severely stunted"))

phyloseq::sample_data(dffiltered4)$ph_estomac = as.numeric(
  phyloseq::sample_data(dffiltered4)$ph_estomac )

phyloseq::sample_data(dffiltered4)$ph_intestin = as.numeric(
  phyloseq::sample_data(dffiltered4)$ph_intestin )

##### saving dffiltered4 #####
saveRDS( object = dffiltered4, file = "dffiltered4.RDS")

dffilteredfeces = dffiltered4 %>% phyloseq::subset_samples(sampletype == "feces")

dffilteredfeces #627 samples

dffilteredfecesBangui = dffilteredfeces %>% phyloseq::subset_samples(pays == "RCA")
dffilteredfecesBangui # 248 feces

dffilteredfecesMada = dffilteredfeces %>% phyloseq::subset_samples(pays == "Madagascar")
dffilteredfecesMada # 379 samples

dffilteredduodenal = dffiltered4 %>%
 phyloseq::subset_samples(sampletype == "duodenal" & stunted!="non-stunted") %>%
 phyloseq::subset_samples(haz!="non-stunted")

dffilteredduodenal # 128 samples

table(phyloseq::sample_data(dffilteredduodenal)$sibo) # 109 samples with data
table(phyloseq::sample_data(dff)$sibo, phyloseq::sample_data(dff)$sampletype) # we lost number of samples due to low sequence reads
table(phyloseq::sample_data(df)$sibo, phyloseq::sample_data(df)$sampletype) # we lost number of samples due to low sequence reads

dffilteredduodenalBangui <- dffiltered4 %>%
 phyloseq::subset_samples(
   sampletype == "duodenal" & pays == "RCA" & stunted != "non-stunted") #

dffilteredduodenalBangui
table(phyloseq::sample_data(dffilteredduodenalBangui)$sibo) # 70 samples with data


dffilteredduodenalMada <- dffiltered4 %>%
 phyloseq::subset_samples(
   sampletype == "duodenal" & pays=="Madagascar" & stunted!="non-stunted") #

dffilteredduodenalMada
table(phyloseq::sample_data(dffilteredduodenalMada)$sibo) # 39 samples with data

dffilteredgastric <- dffiltered4 %>%
 phyloseq::subset_samples(sampletype == "gastric" & stunted!="non-stunted") #

dffilteredgastric # 167 samples

dffilteredgastricBangui <- dffiltered4 %>%
 phyloseq::subset_samples(sampletype == "gastric" & pays=="RCA" & stunted!="non-stunted") #
dffilteredgastricBangui # 67 samples

dffilteredgastricMada <- dffiltered4 %>%
 phyloseq::subset_samples(sampletype == "gastric" & pays=="Madagascar" & stunted!="non-stunted") #
dffilteredgastricMada #100 samples


#### make a substset of subjects with samples for all three compartments ####
test         = phyloseq::sample_data(dffiltered4)
test_gastric = data.frame(phyloseq::sample_data(dffilteredgastric))
test_duodenal= data.frame(phyloseq::sample_data(dffilteredduodenal))
test_feces   = data.frame(phyloseq::sample_data(dffilteredfeces))

test_gastric2 = test_gastric[  test_gastric$ id_individual %in% test_duodenal$id_individual, ]
test_gastric3 = test_gastric2[ test_gastric2$id_individual %in% test_feces$id_individual, ]
test_gastric3 = test_gastric2[ test_gastric2$id_individual %in% test_feces$id_individual, ]
dim(test_gastric3) # 50 samples remaining!

test_duodenal2 = test_duodenal[test_duodenal$id_individual %in% test_gastric3$id_individual, ]
dim(test_duodenal2) # 50 samples remaining!

test_feces2<-test_feces[test_feces$id_individual %in% test_gastric3$id_individual, ]
dim(test_feces2) # 50 samples remaining!

test_feces3<-test_feces2[test_feces2$id_individual %in% test_duodenal2$id_individual, ]
dim(test_feces3) # 50 samples remaining!

dffiltered4_red = phyloseq::subset_samples(dffiltered4, phyloseq::sample_data(dffiltered4)$id_individual %in% test_feces3$id_individual)
dffiltered4_red # 150 samples are remaining!

table(phyloseq::sample_data(dffiltered4_red)$sampletype) # ok!
table(phyloseq::sample_data(dffiltered4_red)$sampletype,
      phyloseq::sample_data(dffiltered4_red)$pays) # ok!



#### get metadata of the final datasets ####
## fecal dataset
phyloseq::sample_data(dffilteredfeces)$ageyears = cut(
  phyloseq::sample_data(dffilteredfeces)$age, c(24,36,48,61)
  , include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)

which(is.na(phyloseq::sample_data(dffilteredfeces)$ageyears)) # these are the controls

levels(phyloseq::sample_data(dffilteredfeces)$ageyears)

# make labels on variables to better interpret results
phyloseq::sample_data(dffilteredfeces)$ageyears = factor(
  phyloseq::sample_data(dffilteredfeces)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))

table(phyloseq::sample_data(dffilteredfeces)$pays)
chi_square_table( dffilteredfeces, "pays", "sexe")
chi_square_table( dffilteredfeces, "pays", "ageyears")
chi_square_table( dffilteredfeces, "pays", "haz")
chi_square_table( dffilteredfeces, "pays", "sibo")
chi_square_table( dffilteredfeces, "pays", "anemie")
chi_square_table( dffilteredfeces, "pays", "crp_seuil")
chi_square_table( dffilteredfeces, "pays", "calprotectinelevel")
chi_square_table( dffilteredfeces, "pays", "alphaantitrypsinlevel")
chi_square_table( dffilteredfeces, "pays", "current_breastfeeding")



table(phyloseq::sample_data(dffilteredfeces)$alpha_calpro_connection) ##
prop.table(table(phyloseq::sample_data(dffilteredfeces)$alpha_calpro_connection))*100



cor(phyloseq::sample_data(dffilteredfeces)$aatmggdeps
    , phyloseq::sample_data(dffilteredfeces)$calprotectineggdeps, method = c("spearman"))
cor.test(phyloseq::sample_data(dffilteredfeces)$aatmggdeps
         , phyloseq::sample_data(dffilteredfeces)$calprotectineggdeps, method=c( "spearman"))

shapiro.test(phyloseq::sample_data(dffilteredfeces)$aatmggdeps) # => p-value < 2.2e-16
shapiro.test(phyloseq::sample_data(dffilteredfeces)$calprotectineggdeps) # => p-value < 2.2e-16

ggpubr::ggscatter(phyloseq::sample_data(dffilteredfeces), x = "calprotectineggdeps", y = "aatmggdeps",
     add = "reg.line", conf.int = TRUE,
     cor.coef = TRUE, cor.method = "spearmann",
     xlab = "calprotectin levels", ylab = "AAT levels")

dffilteredfeces_test = phyloseq::subset_samples(dffilteredfeces, phyloseq::sample_data(dffilteredfeces)$calprotectineggdeps<3000)

cor(phyloseq::sample_data(dffilteredfeces_test)$aatmggdeps, phyloseq::sample_data(dffilteredfeces_test)$calprotectineggdeps, method = c("spearman"))
cor.test(phyloseq::sample_data(dffilteredfeces_test)$aatmggdeps, phyloseq::sample_data(dffilteredfeces_test)$calprotectineggdeps, method=c( "spearman"))

shapiro.test(phyloseq::sample_data(dffilteredfeces_test)$aatmggdeps) # => p-value < 2.2e-16
shapiro.test(phyloseq::sample_data(dffilteredfeces_test)$calprotectineggdeps) # => p-value < 2.2e-16

ggpubr::ggscatter(phyloseq::sample_data(dffilteredfeces_test), x = "calprotectineggdeps", y = "aatmggdeps",
     add = "reg.line", conf.int = TRUE,
     cor.coef = TRUE, cor.method = "spearmann",
     xlab = "calprotectin levels", ylab = "AAT levels")

dffilteredfeces_Bangui = phyloseq::subset_samples(dffilteredfeces, pays=="RCA")

cor(phyloseq::sample_data(dffilteredfeces_Bangui)$aatmggdeps, phyloseq::sample_data(dffilteredfeces_Bangui)$calprotectineggdeps, method = c("spearman"))
cor.test(phyloseq::sample_data(dffilteredfeces_Bangui)$aatmggdeps, phyloseq::sample_data(dffilteredfeces_Bangui)$calprotectineggdeps, method=c( "spearman"))

shapiro.test(phyloseq::sample_data(dffilteredfeces_Bangui)$aatmggdeps) # => p-value = 4.845e-08
shapiro.test(phyloseq::sample_data(dffilteredfeces_Bangui)$calprotectineggdeps) # => p-value < 2.2e-16

ggpubr::ggscatter(phyloseq::sample_data(dffilteredfeces_Bangui), x = "calprotectineggdeps", y = "aatmggdeps",
     add = "reg.line", conf.int = TRUE,
     cor.coef = TRUE, cor.method = "spearmann",
     xlab = "calprotectin levels", ylab = "AAT levels")

dffilteredfeces_Tana = phyloseq::subset_samples(dffilteredfeces, pays=="Madagascar")

cor(phyloseq::sample_data(dffilteredfeces_Tana)$aatmggdeps, phyloseq::sample_data(dffilteredfeces_Tana)$calprotectineggdeps, method = c("spearman"))
cor.test(phyloseq::sample_data(dffilteredfeces_Tana)$aatmggdeps, phyloseq::sample_data(dffilteredfeces_Tana)$calprotectineggdeps, method=c( "spearman"))

shapiro.test(phyloseq::sample_data(dffilteredfeces_Tana)$aatmggdeps) # => p-value < 2.2e-16
shapiro.test(phyloseq::sample_data(dffilteredfeces_Tana)$calprotectineggdeps) # => p-value < 2.2e-16

ggpubr::ggscatter(phyloseq::sample_data(dffilteredfeces_Tana), x = "calprotectineggdeps", y = "aatmggdeps",
     add = "reg.line", conf.int = TRUE,
     cor.coef = TRUE, cor.method = "spearmann",
     xlab = "calprotectin levels", ylab = "AAT levels")


## duodenal dataset
phyloseq::sample_data(dffilteredduodenal)$ageyears<-cut(phyloseq::sample_data(dffilteredduodenal)$age, c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(phyloseq::sample_data(dffilteredduodenal)$ageyears)) # these are the controls
levels(phyloseq::sample_data(dffilteredduodenal)$ageyears)

# make labels on variables to better interpret results
phyloseq::sample_data(dffilteredduodenal)$ageyears<-factor(phyloseq::sample_data(dffilteredduodenal)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))
phyloseq::sample_data(dffilteredduodenal)$ageyears

table(phyloseq::sample_data(dffilteredduodenal)$pays)

chi_square_table( dffilteredduodenal, "pays", "sexe")
chi_square_table( dffilteredduodenal, "pays", "ageyears")
chi_square_table( dffilteredduodenal, "pays", "haz")
chi_square_table( dffilteredduodenal, "pays", "sibo")
chi_square_table( dffilteredduodenal, "pays", "anemie")
chi_square_table( dffilteredduodenal, "pays", "crp_seuil")


cor(phyloseq::sample_data(dffilteredduodenal)$aatsurldmgml, phyloseq::sample_data(dffilteredduodenal)$calprotectinesurliquideduodn, method = c("spearman"))
cor.test(phyloseq::sample_data(dffilteredduodenal)$aatsurldmgml, phyloseq::sample_data(dffilteredduodenal)$calprotectinesurliquideduodn, method=c( "spearman"))

shapiro.test(phyloseq::sample_data(dffilteredduodenal)$aatsurldmgml) # => p-value = 1.247e-14
shapiro.test(phyloseq::sample_data(dffilteredduodenal)$calprotectinesurliquideduodn) # => p-value = 2.247e-15

ggpubr::ggscatter(phyloseq::sample_data(dffilteredduodenal), x = "calprotectineggdeps", y = "aatmggdeps",
     add = "reg.line", conf.int = TRUE,
     cor.coef = TRUE, cor.method = "spearmann",
     xlab = "calprotectin levels in duodenum", ylab = "AAT levels in duodenum")

calprolevelsBangui = as.numeric(phyloseq::sample_data(dffilteredduodenalBangui)$calprotectinesurliquideduodn)
mean(na.omit(calprolevelsBangui))
calprolevelsTana<-as.numeric(phyloseq::sample_data(dffilteredduodenalMada)$calprotectinesurliquideduodn)
mean(na.omit(calprolevelsTana))
wilcox.test(na.omit(calprolevelsBangui), na.omit(calprolevelsTana))

AATlevelsBangui = as.numeric(phyloseq::sample_data(dffilteredduodenalBangui)$aatsurldmgml)
mean(na.omit(AATlevelsBangui))
AATlevelsTana = as.numeric(phyloseq::sample_data(dffilteredduodenalMada)$aatsurldmgml)
mean(na.omit(AATlevelsTana))
wilcox.test(na.omit(AATlevelsBangui), na.omit(AATlevelsTana))

cfumlBangui = as.numeric(phyloseq::sample_data(dffilteredduodenalBangui)$ufcml)
mean(na.omit(cfumlBangui))
min(na.omit(cfumlBangui))
max(na.omit(cfumlBangui))
IQR(na.omit(cfumlBangui))

cfumlTana = as.numeric(phyloseq::sample_data(dffilteredduodenalMada)$ufcml)
mean(na.omit(cfumlTana))
min(na.omit(cfumlTana))
max(na.omit(cfumlTana))
IQR(na.omit(cfumlTana))

wilcox.test(na.omit(cfumlBangui), na.omit(cfumlTana))

#' line
ph_intestinBangui<-as.numeric(phyloseq::sample_data(dffilteredduodenalBangui)$ph_intestin)
mean(na.omit(ph_intestinBangui))
ph_intestinTana<-as.numeric(phyloseq::sample_data(dffilteredduodenalMada)$ph_intestin)
mean(na.omit(ph_intestinTana))
wilcox.test(na.omit(ph_intestinBangui), na.omit(ph_intestinTana))


## gastric dataset
phyloseq::sample_data(dffilteredgastric)$ageyears = cut(phyloseq::sample_data(dffilteredgastric)$age, c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(phyloseq::sample_data(dffilteredgastric)$ageyears)) # these are the controls
levels(phyloseq::sample_data(dffilteredgastric)$ageyears)

# make labels on variables to better interpret results
phyloseq::sample_data(dffilteredgastric)$ageyears = factor(phyloseq::sample_data(dffilteredgastric)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))
phyloseq::sample_data(dffilteredgastric)$ageyears

chi_square_table( dffilteredgastric, "pays", "sexe")
chi_square_table( dffilteredgastric, "pays", "ageyears")
chi_square_table( dffilteredgastric, "pays", "haz")
chi_square_table( dffilteredgastric, "pays", "sibo")
chi_square_table( dffilteredgastric, "pays", "anemie")
chi_square_table( dffilteredgastric, "pays", "crp_seuil")


ph_estomacBangui = as.numeric(phyloseq::sample_data(dffilteredgastricBangui)$ph_estomac)
mean(na.omit(ph_estomacBangui))
ph_estomacTana = as.numeric(phyloseq::sample_data(dffilteredgastricMada)$ph_estomac)
mean(na.omit(ph_estomacTana))
wilcox.test(na.omit(ph_estomacBangui), na.omit(ph_estomacTana))

#### Get information on sequencing depths ####

dffiltered4 = phyloseq::prune_samples( phyloseq::sample_sums(dffiltered4) > 0, dffiltered4 )
minimum= min((phyloseq::sample_sums(dffiltered4)))
minimum ## minimum sample depth is 4581 sequences (remember: we kicked out everyone with less sequences!)

average= mean(phyloseq::sample_sums(dffiltered4))
average #25104.86

.median= median(phyloseq::sample_sums(dffiltered4))
.median #22239.5

dffilteredfeces = phyloseq::prune_samples( phyloseq::sample_sums(dffilteredfeces) > 0, dffilteredfeces )
minimumfeces= min(phyloseq::sample_sums(dffilteredfeces))
minimumfeces ## minimum sample depth is 5032 sequences

averagefeces= mean(phyloseq::sample_sums(dffilteredfeces))
averagefeces #20645.36

medianfeces= median(phyloseq::sample_sums(dffilteredfeces))
medianfeces ##15530

minimumduodenal= min(phyloseq::sample_sums(dffilteredduodenal))
minimumduodenal ## minimum sample depth is 4971 sequences

averageduodenal= mean(phyloseq::sample_sums(dffilteredduodenal))
averageduodenal #32625.01

medianduodenal= median(phyloseq::sample_sums(dffilteredduodenal))
medianduodenal ##30421 sequences

minimumgastric= min(phyloseq::sample_sums(dffilteredgastric))
minimumgastric ## minimum sample depth is 4581 sequences

averagegastric= mean(phyloseq::sample_sums(dffilteredgastric))
averagegastric #36037.12

mediangastric= median(phyloseq::sample_sums(dffilteredgastric))
mediangastric ##31270 sequences

pdf("samplesequeningdepthfulldataset.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
sample_sum_dffiltered <- data.frame(sum = phyloseq::sample_sums(dffiltered4))
ggplot(sample_sum_dffiltered, aes(x = sum)) +
 geom_histogram(color = "black", fill = "blue", binwidth = 500) +
 ggtitle("Distribution of sample sequencing depth before rarefaction") +
 xlab("Read counts") + theme(axis.text = element_text(size = 16), axis.title = element_text(size=18, face="bold"), plot.title=element_text(size=22, face="bold")) +
 theme(axis.title.y = element_blank())
dev.off()


#### Generate rarefaction curves ####
raremax = min(rowSums(t(data.frame(phyloseq::otu_table(dffiltered4)))))
raremax
dffiltered4

col  = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty  = c("solid", "dashed", "longdash", "dotdash")
pars = expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)


### UNCOMMENT:
# out = with(pars[1:57, ],
#       vegan::rarecurve(t(data.frame(phyloseq::otu_table(dffiltered4))), step = 20, sample = 10000, col = "red",
#            lty = "solid", label = FALSE))
### we can also make rarefaction curves colored by Sample Type
# vegan::rarecurve(t(data.frame(phyloseq::otu_table(dffiltered4))), step=100
#           , col=phyloseq::sample_data(dffiltered4)$sampletype, lwd=2, ylab="ASVs", label=F)
# abline(v=(min(rowSums(t(data.frame(phyloseq::otu_table(dffiltered4)))))))

#### Alpha diversity measures ####
#### Prepare untrimmed dataset but rarefied to 5000 sequences ####
df_red = phyloseq::subset_samples(df,
        sampletype=="feces" | sampletype=="duodenal" | sampletype=="gastric")

phyloseq::sample_data(df_red)$ageyears = cut(
  phyloseq::sample_data(df_red)$age, c(24,36,48,61)
  , include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(phyloseq::sample_data(df_red)$ageyears)) # these are the controls
levels(phyloseq::sample_data(df_red)$ageyears)

phyloseq::sample_data(df_red)$haz<-as.factor(phyloseq::sample_data(df_red)$haz)
levels(phyloseq::sample_data(df_red)$haz)
levels(phyloseq::sample_data(df_red)$pays) # we need to transform this to factor and give this labels

# make labels on variables to better interpret results
phyloseq::sample_data(df_red)$ageyears<-factor(phyloseq::sample_data(df_red)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))
phyloseq::sample_data(df_red)$haz<-factor(phyloseq::sample_data(df_red)$haz, labels =c("non-stunted", "moderately stunted", "severely stunted"))
phyloseq::sample_data(df_red)$stunted<-""
phyloseq::sample_data(df_red)$stunted[which(phyloseq::sample_data(df_red)$haz == "moderately stunted")] <- "stunted"
phyloseq::sample_data(df_red)$stunted[which(phyloseq::sample_data(df_red)$haz == "severely stunted")] <- "stunted"
phyloseq::sample_data(df_red)$stunted[which(phyloseq::sample_data(df_red)$haz == "non-stunted")] <- "non-stunted"
phyloseq::sample_data(df_red)$stunted<-as.factor(phyloseq::sample_data(df_red)$stunted)
levels(phyloseq::sample_data(df_red)$stunted)

phyloseq::sample_data(df_red)$ph_estomac<-as.numeric(phyloseq::sample_data(df_red)$ph_estomac)
phyloseq::sample_data(df_red)$ph_intestin<-as.numeric(phyloseq::sample_data(df_red)$ph_intestin)
phyloseq::sample_data(df_red)$anemie<-factor(phyloseq::sample_data(df_red)$anemie, labels =c("Non", "Oui"))


#filter your dataset to keep only relevant data
phyloseq::get_taxa_unique(df_red, "Rank1") ## only bacteria, archae in dataset


# filter out Unassigned, Chloroplast etc.
df_red = df_red %>%
 phyloseq::subset_taxa(Rank5 != "Mitochondria") %>%
 phyloseq::subset_taxa(Rank4 != "Chloroplast") %>%
 phyloseq::subset_taxa(Rank1 != "Unassigned") # was "Kingdom" before
df_red # 5479 taxa

df_redfiltered = df_red %>%
   # to filter out samples which are not recruited in the community
 phyloseq::subset_samples( raison_hospi == "Recrutement communautaire") %>%
  # to keep only samples with good pH range (<=4 for gastric, >=5 for duodenal))
 phyloseq::subset_samples(
   (sampletype == "feces") | ph_intestin >= 5 | ph_estomac <= 4 )  %>%
# kick-out all overweight children
 phyloseq::subset_samples( whz_cont <= 2 )  %>%
#kick-out all children above age 5 years and below two years
 phyloseq::subset_samples(age <= 60) %>%
 phyloseq::subset_samples(age >= 24) # to keep only children above two years

table(phyloseq::sample_data(df_redfiltered)$pays)

#- 670
# Now make a loop and rarefy to even depth in repetitive manner to avoid inducing bias. It is always important to set a seed when you subsample so your result is replicable
set.seed(3)
for (i in 1:1) {
 # Subsample
 df_red_rar  = phyloseq::rarefy_even_depth(
   df_redfiltered, sample.size = 5000, verbose = FALSE, replace = TRUE)
 }

length(which(phyloseq::sample_data(df_red_rar)$pays=="RCA")) #436 are from RCA
length(which(phyloseq::sample_data(df_red_rar)$pays=="Madagascar")) #556 are from Mada

df_red_rar_gastric = phyloseq::subset_samples(df_red_rar, sampletype=="gastric")
df_red_rar_gastric # 214 samples ###  212 

df_red_rar_duodenal = phyloseq::subset_samples(df_red_rar, sampletype=="duodenal")
df_red_rar_duodenal # 139 samples ### 137 

df_red_rar_feces = phyloseq::subset_samples(df_red_rar, sampletype=="feces")
df_red_rar_feces #634 

#
df_red_rar_feces_Mada = phyloseq::subset_samples(df_red_rar_feces, pays=="Madagascar")
df_red_rar_feces_Mada # 380 

df_red_rar_feces_CAR = phyloseq::subset_samples(df_red_rar_feces, pays=="RCA")
df_red_rar_feces_CAR # 254 

df_red_rar_gastric_Mada = phyloseq::subset_samples(df_red_rar_gastric, pays=="Madagascar")
df_red_rar_gastric_Mada #  113 

df_red_rar_gastric_CAR = phyloseq::subset_samples(df_red_rar_gastric, pays=="RCA")
df_red_rar_gastric_CAR # 99 

df_red_rar_duodenal_Mada = phyloseq::subset_samples(df_red_rar_duodenal, pays=="Madagascar")
df_red_rar_duodenal_Mada # 60 samples

df_red_rar_duodenal_CAR = phyloseq::subset_samples(df_red_rar_duodenal, pays=="RCA")
df_red_rar_duodenal_CAR #  77 

#### Make a subset of rarefied, untrimmed samples which do have all three compartments covered for untrimmed dataset  ####
test          = phyloseq::sample_data(df_red_rar)
test_gastric  = phyloseq::sample_data(df_red_rar_gastric)
test_duodenal = phyloseq::sample_data(df_red_rar_duodenal)
test_feces    = phyloseq::sample_data(df_red_rar_feces)

test_gastric2 = test_gastric[test_gastric$id_individual %in% test_duodenal$id_individual, ]
test_gastric3 = test_gastric2[test_gastric2$id_individual %in% test_feces$id_individual, ]

test_gastric3 = test_gastric2[test_gastric2$id_individual %in% test_feces$id_individual, ]
dim(test_gastric3) # 86 samples remaining!

test_duodenal2 = test_duodenal[test_duodenal$id_individual %in% test_gastric3$id_individual, ]
dim(test_duodenal2) # 86 samples remaining!

test_feces2 = test_feces[test_feces$id_individual %in% test_gastric3$id_individual, ]
dim(test_feces2) # 86 samples remaining!

df_red_rar_allthree = phyloseq::subset_samples(df_red_rar
        , phyloseq::sample_data(df_red_rar)$id_individual %in% test_feces2$id_individual)
df_red_rar_allthree # 258 samples are remaining!

table(phyloseq::sample_data(df_red_rar_allthree)$sampletype) # ok!

#### now assess the alpha diversity in the full dataset for sampletype ####

phyloseq::sample_data(df_red_rar)$sampletype =  factor(
phyloseq::sample_data(df_red_rar)$sampletype,
                      levels = c("gastric", "duodenal", "feces"))

table(phyloseq::sample_data(df_red_rar)$sampletype)


out = richness_functions( df_red_rar, variable = "sampletype"
                    , title =  "Alpha Diversity according to Sample Type"
                    , write_table = TRUE, save_pdf = TRUE
                    , pdfname = "alpha_diversity_sampletype.pdf"
                  )
out


#### now assess the alpha diversity in the reduced dataset with only subjects with samples in all three compartments for sampletype  ####
#### Changed this so it works
phyloseq::sample_data(df_red_rar_allthree)$sampletype =  factor(
  phyloseq::sample_data(df_red_rar_allthree)$sampletype,
                      levels = c("gastric", "duodenal", "feces"))

table(phyloseq::sample_data(df_red_rar_allthree)$sampletype)


out = richness_functions( df_red_rar_allthree, variable = "sampletype"
                    , title =  "Alpha Diversity according to Sample Type"
                    , write_table = TRUE, save_pdf = TRUE
                    , pdfname = "alphadiversityASVlevelsampletypereduced.pdf"
                    , filename = "DataAlphasampletypereduced.txt"
                  )
out


#### now assess the alpha diversity in the fecal dataset ####

out = richness_functions( df_red_rar_allthree, variable = "pays"
                    , title =  "Alpha Diversity according to Country of Origin"
                    , write_table = TRUE, save_pdf = TRUE
                    , pdfname = "alphadiversityASVlevelpays.pdf"
                    , filename = "DataAlphapays.txt"
                  )
out


pdf("alphadiversityASVlevelgutinflammation.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 17, # define plot width and height. completely up to user.
  height = 8)
p = plot_richness(df_red_rar_feces
                  , x = "alpha_calpro_connection"
                  , measures = c("Observed", "Observed", "Chao1", "Shannon", "InvSimpson", "InvSimpson")
                  , title = "Alpha Diversity according to Country of origin")
p + geom_boxplot(data = p$data
              , aes(x = alpha_calpro_connection, y = value, color = NULL), alpha = 0.1) +
 theme(axis.text=element_text(size=18)
       , axis.title = element_text(size=18, face="bold")
       , plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0))
       , strip.text.x = element_text(size=18, face="bold" )) +
 xlab("Gut inflammation")

dev.off()

# dataset = df_red_rar_feces; variable = "alphaantitrypsinlevel"

dd = exclude_missing( df_red_rar_feces, "alphaantitrypsinlevel")
outs = richness_loop( df_red_rar_feces
                      , dataset_name = "feces"
                      , variables = c( "stunted" ,"haz", "ageyears",
                                        "calprotectinelevel", "alphaantitrypsinlevel", "current_breastfeeding"))
outs


# df_red_rar_fecesaat = phyloseq::subset_samples(df_red_rar_feces, alphaantitrypsinlevel!="")
# phyloseq::sample_data(df_red_rar_fecesaat)$alphaantitrypsinlevel = as.character(phyloseq::sample_data(df_red_rar_fecesaat)$alphaantitrypsinlevel)
# phyloseq::sample_data(df_red_rar_fecesaat)$alphaantitrypsinlevel[phyloseq::sample_data(df_red_rar_fecesaat)$alphaantitrypsinlevel=="normal"]<-"normal"
# phyloseq::sample_data(df_red_rar_fecesaat)$alphaantitrypsinlevel[phyloseq::sample_data(df_red_rar_fecesaat)$alphaantitrypsinlevel=="elevated"]<-"too_high"
# phyloseq::sample_data(df_red_rar_fecesaat)$alphaantitrypsinlevel<-as.factor(phyloseq::sample_data(df_red_rar_fecesaat)$alphaantitrypsinlevel)

out = richness_functions(  phyloseq::subset_samples(df_red_rar_feces, sibo!="" & stunted=="stunted")
                          , variable = 'sibo'
                          , pdfname = "alphadiversitySIBOASV.pdf"
                          # , filename = "DataAlphagutinflammation.txt"
                          , title ="Alpha Diversity according to Sibo status, fecal samples"
                         )
out

#
df_red_rar_fecesanemie = phyloseq::subset_samples(df_red_rar_feces, anemie!="")
phyloseq::sample_data(df_red_rar_fecesanemie)$anemie =
  as.factor(phyloseq::sample_data(df_red_rar_fecesanemie)$anemie)

out = richness_functions(df_red_rar_fecesanemie, variable = 'anemie'
                          , pdfname = "alphadiversityanemieASV.pdf"
                          # , filename = "DataAlphagutinflammation.txt"
                          , title ="Alpha Diversity according to anemia status, fecal samples"
                         )
out

df_red_rar_fecesanemiestunted = phyloseq::subset_samples(df_red_rar_fecesanemie, haz!="non-stunted")

out = richness_functions(df_red_rar_fecesanemiestunted, variable = 'anemie'
                          , pdfname = "alphadiversityanemieASVstuntedonly.pdf"
                          # , filename = "DataAlphagutinflammation.txt"
                          , title ="Alpha Diversity according to anemia status, fecal samples"
                         )
out

df_red_rar_fecesanemienonstunted =  phyloseq::subset_samples(df_red_rar_fecesanemie, haz=="non-stunted")

#' line
out = richness_functions(df_red_rar_fecesanemienonstunted, variable = 'anemie'
                          , pdfname = "alphadiversityanemieASVnonstuntedonly.pdf"
                          # , filename = "DataAlphagutinflammation.txt"
                          , title ="Alpha Diversity according to anemia status, fecal samples"
                         )
out

df_red_rar_fecesaat =  phyloseq::subset_samples(df_red_rar_feces,
                                                              alphaantitrypsinlevel != "" )

phyloseq::sample_data(df_red_rar_fecesaat)$aatstunted = ""
phyloseq::sample_data(df_red_rar_fecesaat)$aatstunted[which((phyloseq::sample_data(df_red_rar_fecesaat)$alphaantitrypsinlevel != "normal") & (phyloseq::sample_data(df_red_rar_fecesaat)$haz != "non-stunted"))] <- "highAATandstunted"
phyloseq::sample_data(df_red_rar_fecesaat)$aatstunted[which((phyloseq::sample_data(df_red_rar_fecesaat)$alphaantitrypsinlevel != "normal") & (phyloseq::sample_data(df_red_rar_fecesaat)$haz == "non-stunted"))] <- "highAATnonstunted"
phyloseq::sample_data(df_red_rar_fecesaat)$aatstunted[which((phyloseq::sample_data(df_red_rar_fecesaat)$alphaantitrypsinlevel == "normal") & (phyloseq::sample_data(df_red_rar_fecesaat)$haz != "non-stunted"))] <- "normalAATandstunted"
phyloseq::sample_data(df_red_rar_fecesaat)$aatstunted[which((phyloseq::sample_data(df_red_rar_fecesaat)$alphaantitrypsinlevel == "normal") & (phyloseq::sample_data(df_red_rar_fecesaat)$haz == "non-stunted"))] <- "normalAATnonstunted"
table(phyloseq::sample_data(df_red_rar_fecesaat)$aatstunted, phyloseq::sample_data(df_red_rar_fecesaat)$pays)


### NOTE: this has more than 3 levels of a variable (need more code to handle it)
out = richness_functions(df_red_rar_fecesaat, variable = "aatstunted"
                          , pdfname = "alphadiversityaatstuntedASV.pdf"
                          # , filename = "DataAlphagutinflammation.txt"
                          , title ="Alpha Diversity according to AAT levels and stunting, fecal samples"
                         )
out

pdf("alphadiversityaatstuntedASV.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 17, # define plot width and height. completely up to user.
  height = 8)
p = phyloseq::plot_richness(df_red_rar_fecesaat, x= "aatstunted", measures=c("Observed", "Chao1", "Shannon", "InvSimpson")
                            , title = "Alpha Diversity according to AAT levels and stunting, fecal samples")
p + geom_boxplot(data = p$data, aes(x = aatstunted, y = value, color = NULL), alpha = 0.1) +
 theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
 xlab("AAT levels & stunting")
dev.off()

#- 884
results = phyloseq::estimate_richness(df_red_rar_fecesaat)
d       = phyloseq::sample_data(df_red_rar_fecesaat)
highAATandstunted = results[d[,'aatstunted'] == 'highAATandstunted',]
normalAATnonstunted = results[d[,'aatstunted'] == 'normalAATnonstunted',]

wilcox.test(highAATandstunted$Shannon, normalAATnonstunted$Shannon) # p-value = 0.001672
wilcox.test(highAATandstunted$Chao1, normalAATnonstunted$Chao1) # p-value = 0.01766
wilcox.test(highAATandstunted$Observed, normalAATnonstunted$Observed) # p-value = 0.006739
wilcox.test(highAATandstunted$InvSimpson, normalAATnonstunted$InvSimpson) # p-value = 0.05312

highAATnonstunted = results[d[,'aatstunted'] == 'highAATnonstunted',]

wilcox.test(highAATandstunted$Shannon, highAATnonstunted$Shannon) # ns
wilcox.test(highAATandstunted$Chao1, highAATnonstunted$Chao1) # ns
wilcox.test(highAATandstunted$Observed, highAATnonstunted$Observed) # ns
wilcox.test(highAATandstunted$InvSimpson, highAATnonstunted$InvSimpson) # ns

#### now assess the alpha diversity in the fecal dataset, Mada only ####
df_red_rar_feces_Mada = phyloseq::subset_samples(df_red_rar_feces, pays=="Madagascar")

#
outs = richness_loop( df_red_rar_feces_Mada
                      , dataset_name = "feces_Mada"
                      , variables = c("stunted","haz", "ageyears"
                                      , "calprotectinelevel", "alphaantitrypsinlevel", "current_breastfeeding"))
outs


df_red_rar_feces_Madasibo = phyloseq::subset_samples(df_red_rar_feces_Mada, sibo!="" & stunted=="stunted")
out = richness_functions( df_red_rar_feces_Madasibo, variable = "sibo"
                          , pdfname = "alphadiversity_sibo_ASVMada.pdf"
                          # , filename = "DataAlphagutinflammation.txt"
                          , title = "Alpha Diversity according to sibo levels, fecal samples"
                         )
out

#' line
df_red_rar_feces_Madaanemie = phyloseq::subset_samples(df_red_rar_feces_Mada, anemie!="")
phyloseq::sample_data(df_red_rar_feces_Madaanemie)$anemie<-as.factor(phyloseq::sample_data(df_red_rar_feces_Madaanemie)$anemie)

out = richness_functions( df_red_rar_feces_Madaanemie, variable = "anemie"
                          , pdfname = "alphadiversityASV_anemie_fecesMada.pdf"
                          # , filename = "DataAlphagutinflammation.txt"
                          , title = "Alpha Diversity according to anemie status, fecal samples"
                         )
out


#### now assess the alpha diversity in the fecal dataset, CAR only ####
outs = richness_loop(
  phyloseq::subset_samples(
    df_red_rar_feces, pays=="RCA" )
                      , dataset_name = "feces_RCA"
                      , variables = c("stunted","haz", "ageyears"
                                      , "calprotectinelevel", "alphaantitrypsinlevel"))
outs

outs = richness_loop(
  phyloseq::subset_samples(
                       df_red_rar_duodenal, pays == "RCA" & stunted == "stunted"
                        )
                      , dataset_name = "duodenal_stunted_RCA"
                      , variables = c("sibo", "anemie"))

outs

#### now assess the alpha diversity in the duodenal dataset ####
outs = richness_loop( df_red_rar_duodenal
                      , dataset_name = "duodenal"
                      , variables = c( "pays", "stunted"
                                       ,"haz", "ageyears"
                                      , "calprotectinelevel", "alphaantitrypsinlevel"
                                      ))
outs

#- 966
outs = richness_loop(
  phyloseq::subset_samples(
                        df_red_rar_duodenal, stunted == "stunted"
                        )
                      , dataset_name = "duodenal_stunted"
                      , variables = c("sibo", "anemie"))
outs

#### now assess the alpha diversity in the duodenal dataset, Mada only ####
outs = richness_loop(
  phyloseq::subset_samples(
                        df_red_rar_duodenal, pays == "Madagascar"
                        )
                      , dataset_name = "duodenal_mada"
                      , variables = c( "stunted","haz"
                                       , "ageyears"
                                      , "calprotectinelevel", "alphaantitrypsinlevel"
                                      ))
outs

outs = richness_loop(
  phyloseq::subset_samples(
                        df_red_rar_duodenal, pays == "Madagascar" & stunted == "stunted"
                        )
                      , dataset_name = "duodenal_stunted_mada"
                      , variables = c("sibo", "anemie"))
outs

#### now assess the alpha diversity in the duodenal dataset, CAR only ####
outs = richness_loop(
  phyloseq::subset_samples(
                        df_red_rar_duodenal, pays == "RCA"
                        )
                      , dataset_name = "duodenal_RCA"
                      , variables = c(
                        # only stunted children in RCA
                        # "stunted",
                                       "haz", "ageyears"
                                      , "calprotectinelevel", "alphaantitrypsinlevel"
                                      ))
outs

outs = richness_loop(
  phyloseq::subset_samples(
                        df_red_rar_duodenal, pays == "RCA" & stunted == "stunted"
                        )
                      , dataset_name = "duodenal_stunted_RCA"
                      , variables = c("sibo", "anemie"))

outs

#### Who is there? Descriptive tables and plots for the overall sample composition ####
#### calculate prevalence and abundance for species, phylum, genus for all samples ####
dffiltered_s = phyloseq::tax_glom(dffiltered4, "Species")

# fecal samples subset to 129 samples (to match with duodenal samples)
length(which(phyloseq::sample_data(dffiltered_s)$ sampletype == "duodenal")) #129 samples
length(which(phyloseq::sample_data(dffiltered_s)$ sampletype == "gastric")) #168 samples
length(which(phyloseq::sample_data(dffiltered_s)$ sampletype == "feces")) #627 samples, reduce and restrict on malnourished samples

# dffiltered_g = phyloseq::tax_glom(dffiltered4, "Genus")
# dffiltered_p = phyloseq::tax_glom(dffiltered4, "Phylum")
# dffiltered_s_red    = phyloseq::tax_glom(dffiltered4_red, "Species")
# dffiltered4_Bangui = phyloseq::subset_samples(dffiltered4, pays == "RCA")
# dffiltered4_Tana   = phyloseq::subset_samples(dffiltered4, pays=="Madagascar")

colnames( phyloseq::tax_table(dffiltered4)) <- c(
              "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Accession", "ASV_sequence")


for (level in c("Species", "Phylum", "Genus")){
res = prevalence_abundance( dataset = dffiltered4
                      , level = level
                      , name = ""
                      , group = "sampletype"
                      , write_csv = TRUE )
}

#### calculate prevalence for Species for subjects with samples from all three compartments ####
colnames(phyloseq::tax_table(dffiltered4_red)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Accession", "ASV_sequence")

##### prevalence for Species, Genus, Phylum #####
for( level in c("Species","Genus","Phylum")){
res = prevalence_abundance( dataset = dffiltered4_red
                      , level = level
                      , name = ""
                      , group = "sampletype"
                      , write_csv = TRUE )
}
#- 1050
#### calculate prevalence for phylum for all samples CAR, Madagascar ####
for( pays in c("RCA" , "Madagascar")){
res = prevalence_abundance( dataset =  phyloseq::subset_samples(dffiltered4, pays == pays)
                      , level = "Phylum"
                      , name = ""
                      , group = "sampletype"
                      , write_csv = TRUE )
}

#### calculate prevalence for Phylum, Genus, Species for fecal samples and country ####

for( level in  c( "Phylum" ,"Genus", "Species" )) {
res = prevalence_abundance( dataset =  dffilteredfeces
                      , level = level
                      , name = "feces"
                      , group = "pays"
                      , write_csv = TRUE )

}
#### calculate prevalence for Phylum, Genus, Species for duodenal samples and pays ####

colnames(phyloseq::tax_table(dffilteredduodenal)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Accession", "ASV_sequence")

for( level in  c( "Phylum" ,"Genus", "Species" )) {
res = prevalence_abundance( dataset =  dffilteredduodenal
                      , level = level
                      , name = "duodenal"
                      , group = "pays"
                      , write_csv = TRUE )
}

#### calculate prevalence for Phylum, Genus, Species for gastric samples and pays ####
colnames(phyloseq::tax_table(dffilteredfeces)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Accession", "ASV_sequence")

for( level in c( "Phylum" ,"Genus", "Species" )) {
res = prevalence_abundance( dataset =  dffilteredgastric
                      , level = level
                      , name = "gastric"
                      , group = "pays"
                      , write_csv = TRUE )
}


#### Abundance plots: make an abundance plot at the phylum level for different categories (country, stunting status, Sample Type) ####
phyla_counts_tab  = phyloseq::otu_table(phyloseq::tax_glom(dffiltered4, taxrank="Phylum"))
phyla_tax_vec <- as.vector(phyloseq::tax_table(phyloseq::tax_glom(dffiltered4, taxrank="Phylum"))[,2]) # set phyla names as row name
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)
unclassified_tax_counts <- colSums(phyloseq::otu_table(dffiltered4)) - colSums(phyla_counts_tab) # these are the unclassified phyla
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)
dim(phyla_and_unidentified_counts_tab) ## we have 18 phyla. This is ok for a graph!
taxa_proportions_tab <- apply(phyla_and_unidentified_counts_tab, 2, function(x) x/sum(x)*100) # transform in rel. abundance

Afribiota16Sft  = phyloseq::transform_sample_counts(dffiltered4, function(x) x/sum(x))

phylum_colors <- c(
 "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
 "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
 "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "#8A7C64", "blue"
)
#- 1112


Afribiota_phylum <- dffiltered4 %>%
 phyloseq::tax_glom(taxrank = "Phylum") %>%           # agglomerate at phylum level
 phyloseq::transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

Afribiota_phylum$sampletype<-factor(Afribiota_phylum$sampletype, levels=c("gastric", "duodenal", "feces"))

# Plot according to Sample Type and Country
pdf("Phylumcomposition of Afribiota Samples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 30, # define plot width and height. completely up to user.
  height = 10)
ggplot(Afribiota_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) +
 facet_wrap(pays~sampletype, scales = "free_x") +
 theme(strip.text = element_text(colour = 'black', size=18), axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.title.y = element_text(size=18), title =element_text(size=25, face='bold')) +
 geom_bar(stat = "identity") +
 scale_fill_manual(values = phylum_colors) +
 scale_x_discrete(
  drop = FALSE
 ) +
 # Remove x axis title
 theme(axis.title.x = element_blank()) +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
 ylab("Relative Abundance") +
 ggtitle("Phylum Composition of Afribiota samples")
dev.off()

# now look at average, using malnourished samples only for comparison
Afribiota16Sfmalonly = phyloseq::subset_samples(dffiltered4, stunted="stunted")
Afribiota16Sfmmalonlytyp = phyloseq::merge_samples(Afribiota16Sfmalonly, "sampletype")
phyloseq::sample_data(Afribiota16Sfmmalonlytyp)$sampletype <- levels(phyloseq::sample_data(Afribiota16Sfmmalonlytyp)$sampletype)

# Phylum level and Sample Type in stunted children

#- 1150
Afribiota16Sftmalonlytyp_p  = phyloseq::tax_glom(Afribiota16Sfmmalonlytyp, "Phylum")

Afribiota16Sftmalonlytyp_p = phyloseq::transform_sample_counts(Afribiota16Sftmalonlytyp_p, function(x) 100 * x/sum(x))

Afribiota16Sftmalonlytyp_p_s= Afribiota16Sftmalonlytyp_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

Afribiota16Sftmalonlytyp_p_s$Sample<-factor(Afribiota16Sftmalonlytyp_p_s$Sample, levels=c("gastric", "duodenal", "feces"))

pdf("AbundancetablePhylumsampletype.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 20, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftmalonlytyp_p_s, aes(x =Sample, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Sample Type") + ggtitle ("Relative Abundance of phyla according to Sample Type")
dev.off()

#- 1175

# now look at average, using only subjects which have samples for all three compartments
# phyloseq::sample_data(dffiltered4_red)$sampletype <- levels(phyloseq::sample_data(dffiltered4_red)$sampletype)
Afribiota16Sshared = phyloseq::merge_samples(dffiltered4_red, "sampletype")

# Phylum level and Sample Type in stunted children

Afribiota16Ssharedtyp_p  = phyloseq::tax_glom(Afribiota16Sshared, "Phylum")

Afribiota16Ssharedtyp_p = phyloseq::transform_sample_counts(Afribiota16Ssharedtyp_p, function(x) 100 * x/sum(x))

Afribiota16Ssharedtyp_p_s= Afribiota16Ssharedtyp_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

Afribiota16Ssharedtyp_p_s$Sample<-factor(Afribiota16Ssharedtyp_p_s$Sample, levels=c("gastric", "duodenal", "feces"))

pdf("AbundancetablePhylumsampletypesharedonly.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 20, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Ssharedtyp_p_s, aes(x =Sample, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Sample Type") + ggtitle ("Relative Abundance of phyla according to Sample Type")
dev.off()


# Phylum level in feces and country
Afribiota16Sftfeces = phyloseq::subset_samples(dffilteredfeces, sampletype=="feces")
phyloseq::sample_data(Afribiota16Sftfeces)$pays<-as.factor(phyloseq::sample_data(Afribiota16Sftfeces)$pays)
Afribiota16Sftfecespays = phyloseq::merge_samples(Afribiota16Sftfeces, "pays")

Afribiota16Sftfeces_p  = phyloseq::tax_glom(Afribiota16Sftfecespays, "Phylum")

Afribiota16Sftfeces_p = phyloseq::transform_sample_counts(Afribiota16Sftfeces_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftfeces_p)$pays[(phyloseq::sample_data(Afribiota16Sftfeces_p)$pays=="1")] <- "Madagascar"
phyloseq::sample_data(Afribiota16Sftfeces_p)$pays[(phyloseq::sample_data(Afribiota16Sftfeces_p)$pays=="2")] <- "CAR"

Afribiota16Sftfeces_p_s= Afribiota16Sftfeces_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

pdf("AbundancetablePhylumfecescountryoforigin.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 20, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftfeces_p_s, aes(x =pays, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Country of origin") + ggtitle ("Relative Abundance of phyla in feces according to country of origin")
dev.off()

# Phylum level in duodenal samples and country
phyloseq::sample_data(dffilteredduodenal)$pays<-as.factor(phyloseq::sample_data(dffilteredduodenal)$pays)

Afribiota16Sftduodenalpays = phyloseq::merge_samples(dffilteredduodenal, "pays")

Afribiota16Sftduodenal_p  = phyloseq::tax_glom(Afribiota16Sftduodenalpays, "Phylum")

Afribiota16Sftduodenal_p = phyloseq::transform_sample_counts(Afribiota16Sftduodenal_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftduodenal_p)$pays[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$pays=="1")] <- "Madagascar"
phyloseq::sample_data(Afribiota16Sftduodenal_p)$pays[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$pays=="2")] <- "CAR"

Afribiota16Sftduodenal_p_s= Afribiota16Sftduodenal_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

pdf("AbundancetablePhylumduodenalcountryoforigin.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 6, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftduodenal_p_s, aes(x =pays, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Country of origin") + ggtitle ("Relative Abundance of phyla in duodenal aspirates according to country of origin")
dev.off()

#- 1270
dffilteredduodenal2 = phyloseq::subset_samples(dffilteredduodenal, sibo!="")
Afribiota16Sftduodenalsibo = phyloseq::merge_samples(dffilteredduodenal, "sibo")

Afribiota16Sftduodenal_p  = phyloseq::tax_glom(Afribiota16Sftduodenalsibo, "Phylum")

Afribiota16Sftduodenal_p = phyloseq::transform_sample_counts(Afribiota16Sftduodenal_p, function(x) 100 * x/sum(x))

Afribiota16Sftduodenal_p_s= Afribiota16Sftduodenal_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

pdf("AbundancetablePhylumduodenalsibo.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 6, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftduodenal_p_s, aes(x =Sample, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Sibo yes/no") + ggtitle ("Relative Abundance of phyla in duodenal aspirates according to sibo")
dev.off()


# Phylum level in gastric samples and country

phyloseq::sample_data(dffilteredgastric)$pays<-as.factor(phyloseq::sample_data(dffilteredgastric)$pays)

Afribiota16Sftgastricpays = phyloseq::merge_samples(dffilteredgastric, "pays")

Afribiota16Sftgastric_p  = phyloseq::tax_glom(Afribiota16Sftgastricpays, "Phylum")

Afribiota16Sftgastric_p = phyloseq::transform_sample_counts(Afribiota16Sftgastric_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftgastric_p)$pays[(phyloseq::sample_data(Afribiota16Sftgastric_p)$pays=="1")] <- "Madagascar"
phyloseq::sample_data(Afribiota16Sftgastric_p)$pays[(phyloseq::sample_data(Afribiota16Sftgastric_p)$pays=="2")] <- "CAR"

Afribiota16Sftgastric_p_s= Afribiota16Sftgastric_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

pdf("AbundancetablePhylumgastriccountryoforigin.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 20, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftgastric_p_s, aes(x =pays, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Country of origin") + ggtitle ("Relative Abundance of phyla in gastric aspirates according to country of origin")
dev.off()

# Phylum level in gastric samples and stunting status, both countries merged
dffilteredgastric = phyloseq::subset_samples(dffilteredgastric, haz!="0")
Afribiota16Sftgastrichaz = phyloseq::merge_samples(dffilteredgastric, "haz")

Afribiota16Sftgastric_p  = phyloseq::tax_glom(Afribiota16Sftgastrichaz, "Phylum")

Afribiota16Sftgastric_p = phyloseq::transform_sample_counts(Afribiota16Sftgastric_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftgastric_p)$haz[(phyloseq::sample_data(Afribiota16Sftgastric_p)$haz=="1")] <- "moderately stunted"
phyloseq::sample_data(Afribiota16Sftgastric_p)$haz[(phyloseq::sample_data(Afribiota16Sftgastric_p)$haz=="2")] <- "severely stunted"

Afribiota16Sftgastric_p_s= Afribiota16Sftgastric_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

pdf("AbundancetablePhylumgastricstuntingstatus.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 11, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftgastric_p_s, aes(x =haz, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Stunting status") + ggtitle ("Relative Abundance of phyla in gastric aspirates according to stunting status")
dev.off()

#Phylum level in gastric samples and stunting status, CAR only
Afribiota16Sftgastrichaz = phyloseq::merge_samples(dffilteredgastricBangui, "haz")

Afribiota16Sftgastric_p  = phyloseq::tax_glom(Afribiota16Sftgastrichaz, "Phylum")

Afribiota16Sftgastric_p = phyloseq::transform_sample_counts(Afribiota16Sftgastric_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftgastric_p)$haz[(phyloseq::sample_data(Afribiota16Sftgastric_p)$haz=="1")] <- "moderately stunted"
phyloseq::sample_data(Afribiota16Sftgastric_p)$haz[(phyloseq::sample_data(Afribiota16Sftgastric_p)$haz=="2")] <- "severely stunted"

Afribiota16Sftgastric_p_s= Afribiota16Sftgastric_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

pdf("AbundancetablePhylumgastricstuntingstatusBanguionly.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 11, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftgastric_p_s, aes(x =haz, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Stunting status") + ggtitle ("Relative Abundance of phyla in gastric aspirates from Bangui, CAR according to stunting status")
dev.off()

# Phylum level in gastric samples and stunting status, Mada only
Afribiota16Sftgastrichaz = phyloseq::merge_samples(dffilteredgastricMada, "haz")

Afribiota16Sftgastric_p  = phyloseq::tax_glom(Afribiota16Sftgastrichaz, "Phylum")

Afribiota16Sftgastric_p = phyloseq::transform_sample_counts(Afribiota16Sftgastric_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftgastric_p)$haz[(phyloseq::sample_data(Afribiota16Sftgastric_p)$haz=="0")] <- "non-stunted"
phyloseq::sample_data(Afribiota16Sftgastric_p)$haz[(phyloseq::sample_data(Afribiota16Sftgastric_p)$haz=="1")] <- "moderately stunted"
phyloseq::sample_data(Afribiota16Sftgastric_p)$haz[(phyloseq::sample_data(Afribiota16Sftgastric_p)$haz=="2")] <- "severely stunted"

Afribiota16Sftgastric_p_s= Afribiota16Sftgastric_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

pdf("AbundancetablePhylumgastricstuntingstatusTanaonly.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 15, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftgastric_p_s, aes(x =haz, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Stunting status") + ggtitle ("Relative Abundance of phyla in gastric aspirates from Antananarivo, Madagascar according to stunting status")
dev.off()


#-1420

# Phylum level in duodenal samples and stunting status, both countries merged

Afribiota16Sftduodenalhaz = phyloseq::merge_samples(dffilteredduodenal, "haz")

Afribiota16Sftduodenal_p  = phyloseq::tax_glom(Afribiota16Sftduodenalhaz, "Phylum")

Afribiota16Sftduodenal_p = phyloseq::transform_sample_counts(Afribiota16Sftduodenal_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftduodenal_p)$haz[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$haz=="1")] <- "moderately stunted"
phyloseq::sample_data(Afribiota16Sftduodenal_p)$haz[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$haz=="2")] <- "severely stunted"

Afribiota16Sftduodenal_p_s= Afribiota16Sftduodenal_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

pdf("AbundancetablePhylumduodenalstuntingstatus.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 11, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftduodenal_p_s, aes(x =haz, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Stunting status") + ggtitle ("Relative Abundance of phyla in duodenal aspirates according to stunting status")
dev.off()

# Phylum level in duodenal samples and stunting status, CAR only
Afribiota16Sftduodenalhaz = phyloseq::merge_samples(dffilteredduodenalBangui, "haz")

Afribiota16Sftduodenal_p  = phyloseq::tax_glom(Afribiota16Sftduodenalhaz, "Phylum")

Afribiota16Sftduodenal_p = phyloseq::transform_sample_counts(Afribiota16Sftduodenal_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftduodenal_p)$haz[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$haz=="1")] <- "moderately stunted"
phyloseq::sample_data(Afribiota16Sftduodenal_p)$haz[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$haz=="2")] <- "severely stunted"

Afribiota16Sftduodenal_p_s= Afribiota16Sftduodenal_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

Afribiota16Sftduodenal_p_s_f = dplyr::filter(Afribiota16Sftduodenal_p_s, Abundance>0.1)

#There are more than 18 phyla here (31!). Restrict on taxa with at least an abundance of 1%

pdf("AbundancetablePhylumduodenalstuntingstatusBanguionly.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 11, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftduodenal_p_s_f, aes(x =haz, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Stunting status") + ggtitle ("Relative Abundance of phyla in duodenal aspirates from Bangui, CAR according to stunting status")
dev.off()

#- 1487
# Phylum level in duodenal samples and stunting status, Mada only
Afribiota16Sftduodenalhaz = phyloseq::merge_samples(dffilteredduodenalMada, "haz")

Afribiota16Sftduodenal_p  = phyloseq::tax_glom(Afribiota16Sftduodenalhaz, "Phylum")

Afribiota16Sftduodenal_p = phyloseq::transform_sample_counts(Afribiota16Sftduodenal_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftduodenal_p)$haz[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$haz=="1")] <- "moderately stunted"
phyloseq::sample_data(Afribiota16Sftduodenal_p)$haz[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$haz=="2")] <- "severely stunted"

Afribiota16Sftduodenal_p_s= Afribiota16Sftduodenal_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

Afribiota16Sftduodenal_p_s_f = dplyr::filter(Afribiota16Sftduodenal_p_s, Abundance>0.1)


pdf("AbundancetablePhylumduodenalstuntingstatusTanaonly.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 11, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftduodenal_p_s_f, aes(x =haz, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Stunting status") + ggtitle ("Relative Abundance of phyla in duodenal aspirates from Antananarivo, Madagascar according to stunting status")
dev.off()

# Phylum level in feces samples and stunting status, both countries merged

table((phyloseq::sample_data(dffilteredfeces)$haz))
Afribiota16Sftfeceshaz = phyloseq::merge_samples(dffilteredfeces, "haz")
table((phyloseq::sample_data(Afribiota16Sftfeceshaz)$haz))

Afribiota16Sftfeces_p  = phyloseq::tax_glom(Afribiota16Sftfeceshaz, "Phylum")

Afribiota16Sftfeces_p = phyloseq::transform_sample_counts(Afribiota16Sftfeces_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftfeces_p)$haz[(phyloseq::sample_data(Afribiota16Sftfeces_p)$haz=="2")] <- "moderately stunted"
phyloseq::sample_data(Afribiota16Sftfeces_p)$haz[(phyloseq::sample_data(Afribiota16Sftfeces_p)$haz=="3")] <- "severely stunted"
phyloseq::sample_data(Afribiota16Sftfeces_p)$haz[(phyloseq::sample_data(Afribiota16Sftfeces_p)$haz=="1")] <- "non stunted"

#' line
Afribiota16Sftfeces_p_s= Afribiota16Sftfeces_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

Afribiota16Sftfeces_p_s$haz<-factor(Afribiota16Sftfeces_p_s$haz, levels=c("non stunted", "moderately stunted", "severely stunted"))


pdf("AbundancetablePhylumfecesstuntingstatus.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 15, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftfeces_p_s, aes(x =haz, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Stunting status") + ggtitle ("Relative Abundance of phyla in feces according to stunting status")
dev.off()

# Phylum level in feces samples and stunting status, CAR only
Afribiota16Sftfeceshaz = phyloseq::merge_samples(dffilteredfecesBangui, "haz")

Afribiota16Sftfeces_p  = phyloseq::tax_glom(Afribiota16Sftfeceshaz, "Phylum")

Afribiota16Sftfeces_p = phyloseq::transform_sample_counts(Afribiota16Sftfeces_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftfeces_p)$haz[(phyloseq::sample_data(Afribiota16Sftfeces_p)$haz=="2")] <- "moderately stunted"
phyloseq::sample_data(Afribiota16Sftfeces_p)$haz[(phyloseq::sample_data(Afribiota16Sftfeces_p)$haz=="3")] <- "severely stunted"
phyloseq::sample_data(Afribiota16Sftfeces_p)$haz[(phyloseq::sample_data(Afribiota16Sftfeces_p)$haz=="1")] <- "non stunted"

Afribiota16Sftfeces_p_s= Afribiota16Sftfeces_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

Afribiota16Sftfeces_p_s_f = dplyr::filter(Afribiota16Sftfeces_p_s, Abundance>0.1)

Afribiota16Sftfeces_p_s_f$haz<-factor(Afribiota16Sftfeces_p_s_f$haz, levels=c("non stunted", "moderately stunted", "severely stunted"))


#There are more than 18 phyla here (31!). Restrict on taxa with at least an abundance of 1%

pdf("AbundancetablePhylumfecesstuntingstatusBanguionly.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 15, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftfeces_p_s_f, aes(x =haz, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Stunting status") + ggtitle ("Relative Abundance of phyla in feces aspirates from Bangui, CAR according to stunting status")
dev.off()

# Phylum level in feces samples and stunting status, Mada only
Afribiota16Sftfeceshaz = phyloseq::merge_samples(dffilteredfecesMada, "haz")

Afribiota16Sftfeces_p  = phyloseq::tax_glom(Afribiota16Sftfeceshaz, "Phylum")

Afribiota16Sftfeces_p = phyloseq::transform_sample_counts(Afribiota16Sftfeces_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftfeces_p)$haz[(phyloseq::sample_data(Afribiota16Sftfeces_p)$haz=="2")] <- "moderately stunted"
phyloseq::sample_data(Afribiota16Sftfeces_p)$haz[(phyloseq::sample_data(Afribiota16Sftfeces_p)$haz=="3")] <- "severely stunted"
phyloseq::sample_data(Afribiota16Sftfeces_p)$haz[(phyloseq::sample_data(Afribiota16Sftfeces_p)$haz=="1")] <- "non stunted"

Afribiota16Sftfeces_p_s= Afribiota16Sftfeces_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Phylum)

Afribiota16Sftfeces_p_s_f = dplyr::filter(Afribiota16Sftfeces_p_s, Abundance>0.1)

Afribiota16Sftfeces_p_s_f$haz<-factor(Afribiota16Sftfeces_p_s_f$haz, levels=c("non stunted", "moderately stunted", "severely stunted"))


pdf("AbundancetablePhylumfecesstuntingstatusTanaonly.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 15, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftfeces_p_s_f, aes(x =haz, y = Abundance, fill = Phylum)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 scale_fill_manual(values = phylum_colors) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Stunting status") + ggtitle ("Relative Abundance of phyla in feces aspirates from Antananarivo, Madagascar according to stunting status")
dev.off()

# Genus level in duodenal samples and sibo
phyloseq::sample_data(dffilteredduodenal)$sibo<-as.factor(phyloseq::sample_data(dffilteredduodenal)$sibo)

Afribiota16Sftduodenalsibo = phyloseq::merge_samples(dffilteredduodenal, "sibo")

Afribiota16Sftduodenal_p  = phyloseq::tax_glom(Afribiota16Sftduodenalsibo, "Genus")

Afribiota16Sftduodenal_p = phyloseq::transform_sample_counts(Afribiota16Sftduodenal_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftduodenal_p)$sibo[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$sibo=="1")] <- "No SIBO"
phyloseq::sample_data(Afribiota16Sftduodenal_p)$sibo[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$sibo=="2")] <- "SIBO"

Afribiota16Sftduodenal_p_s= Afribiota16Sftduodenal_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Genus)


pdf("AbundancetableGenusduodenalsibo.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 50, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftduodenal_p_s, aes(x =sibo, y = Abundance, fill = Genus)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Sibo yes/no") + ggtitle ("Relative Abundance of phyla in duodenal aspirates according to sibo")
dev.off()

#' 2220 line


# now restricted to top 10
TopNOTUs <- names(sort(taxa_sums(Afribiota16Sftduodenal_p), TRUE)[1:10])
taxtable<-as.data.frame(phyloseq::tax_table(Afribiota16Sftduodenal_p))
TopNOTUannotated = dplyr::filter(taxtable, (row.names(taxtable) %in% TopNOTUs))

duodenal10 <- prune_taxa(TopNOTUs, Afribiota16Sftduodenal_p)
print(duodenal10)

duodenal10= duodenal10 %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Genus)

pdf("AbundancetableGenusduodenalsibotop10.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 6, # define plot width and height. completely up to user.
  height = 8)
ggplot(duodenal10, aes(x =sibo, y = Abundance, fill = Genus)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_fill_manual(values = c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "blue"
 )) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Sibo yes/no") + ggtitle ("Relative Abundance of phyla in duodenal aspirates according to sibo")
dev.off()

# Genus level in duodenal samples and sibo CAR only
phyloseq::sample_data(dffilteredduodenalBangui)$sibo<-as.factor(phyloseq::sample_data(dffilteredduodenalBangui)$sibo)

Afribiota16Sftduodenalsibo = phyloseq::merge_samples(dffilteredduodenalBangui, "sibo")

Afribiota16Sftduodenal_p  = phyloseq::tax_glom(Afribiota16Sftduodenalsibo, "Genus")

Afribiota16Sftduodenal_p = phyloseq::transform_sample_counts(Afribiota16Sftduodenal_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftduodenal_p)$sibo[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$sibo=="1")] <- "No SIBO"
phyloseq::sample_data(Afribiota16Sftduodenal_p)$sibo[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$sibo=="2")] <- "SIBO"

Afribiota16Sftduodenal_p_s= Afribiota16Sftduodenal_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Genus)


pdf("AbundancetableGenusduodenalsiboBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 50, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftduodenal_p_s, aes(x =sibo, y = Abundance, fill = Genus)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Sibo yes/no") + ggtitle ("Relative Abundance of phyla in duodenal aspirates according to sibo")
dev.off()


# now restricted to top 10
TopNOTUs <- names(sort(taxa_sums(Afribiota16Sftduodenal_p), TRUE)[1:10])
taxtable<-as.data.frame(phyloseq::tax_table(Afribiota16Sftduodenal_p))
TopNOTUannotated = dplyr::filter(taxtable, (row.names(taxtable) %in% TopNOTUs))

duodenal10 <- prune_taxa(TopNOTUs, Afribiota16Sftduodenal_p)
print(duodenal10)

duodenal10= duodenal10 %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Genus)

pdf("AbundancetableGenusduodenalsibotop10Bangui.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 12, # define plot width and height. completely up to user.
  height = 8)
ggplot(duodenal10, aes(x =sibo, y = Abundance, fill = Genus)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_fill_manual(values = c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "blue"
 )) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Sibo yes/no") + ggtitle ("Relative Abundance of phyla in duodenal aspirates according to sibo")
dev.off()

# Genus level in duodenal samples and sibo CAR only
phyloseq::sample_data(dffilteredduodenalMada)$sibo<-as.factor(phyloseq::sample_data(dffilteredduodenalMada)$sibo)

Afribiota16Sftduodenalsibo = phyloseq::merge_samples(dffilteredduodenalMada, "sibo")

Afribiota16Sftduodenal_p  = phyloseq::tax_glom(Afribiota16Sftduodenalsibo, "Genus")

Afribiota16Sftduodenal_p = phyloseq::transform_sample_counts(Afribiota16Sftduodenal_p, function(x) 100 * x/sum(x))

phyloseq::sample_data(Afribiota16Sftduodenal_p)$sibo[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$sibo=="1")] <- "No SIBO"
phyloseq::sample_data(Afribiota16Sftduodenal_p)$sibo[(phyloseq::sample_data(Afribiota16Sftduodenal_p)$sibo=="2")] <- "SIBO"

Afribiota16Sftduodenal_p_s= Afribiota16Sftduodenal_p %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Genus)


pdf("AbundancetableGenusduodenalsiboTana.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 50, # define plot width and height. completely up to user.
  height = 8)
ggplot(Afribiota16Sftduodenal_p_s, aes(x =sibo, y = Abundance, fill = Genus)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Sibo yes/no") + ggtitle ("Relative Abundance of phyla in duodenal aspirates according to sibo")
dev.off()


# now restricted to top 10
TopNOTUs <- names(sort(taxa_sums(Afribiota16Sftduodenal_p), TRUE)[1:10])
taxtable<-as.data.frame(phyloseq::tax_table(Afribiota16Sftduodenal_p))
TopNOTUannotated = dplyr::filter(taxtable, (row.names(taxtable) %in% TopNOTUs))

duodenal10 <- prune_taxa(TopNOTUs, Afribiota16Sftduodenal_p)
print(duodenal10)

#' line
duodenal10= duodenal10 %>%
 phyloseq::psmelt() %>%                     # Melt to long format
 dplyr::arrange(Genus)

pdf("AbundancetableGenusduodenalsibotop10Tana.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 12, # define plot width and height. completely up to user.
  height = 8)
ggplot(duodenal10, aes(x =sibo, y = Abundance, fill = Genus)) +
 theme_bw() +
 theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
 theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
 geom_bar(stat = "identity", width = 1.0) +
 scale_fill_manual(values = c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "blue"
 )) +
 scale_y_continuous(expand = c(0.01,0.01)) +
 theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) +
 ylab("Relative Abundance") +
 xlab("Sibo yes/no") + ggtitle ("Relative Abundance of phyla in duodenal aspirates according to sibo")
dev.off()




#### Ordination plots ####
#### Make a PERMANOVA to look for factors influencing dispersion of whole sample set ####

dffiltered5 = phyloseq::subset_samples(dffiltered4, (sampletype=="feces" | sampletype=="duodenal" | sampletype=="gastric"))

phyloseq::sample_data(dffiltered5)$ageyears = cut(phyloseq::sample_data(dffiltered5)$age, c(24,36,48,61)
                                        , include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(phyloseq::sample_data(dffiltered5)$ageyears)) # these are the controls
levels(phyloseq::sample_data(dffiltered5)$ageyears)

#' line
set.seed(3)
for (i in 1:1) {
# Subsample
dffiltered5_rar = phyloseq::rarefy_even_depth(dffiltered5
                                              , sample.size = 5000
                                              , verbose = FALSE, replace = TRUE)
}


#- 1839
controls = c("run", "haz", "sampletype", "anemie", "calprotectinelevel"
  , "ageyears", "read_count", "pays", "sexe", "alphaantitrypsinlevel")
out = beta_diversity_setup( dataset = dffiltered5_rar, controls = controls
                            , pdfname = "Contributionfulldataset.pdf"
                            , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

#- 1846

##### Plot Ordination with significant taxa as explanatory arrows on Genus, Family, Order level- sampletype ####
for ( level in c("Genus","Family","Order")){
# level = "Family"
out = setup_ordination( dataset = dffiltered5_rar, level = level )
out2= setup_pcoa_scores( dist = out$dist, test4 = out$test4, meta = out$meta, sig_level = 0.01 )

pdf(file="Genussampletype16S_explantory_arrowsrelabun.pdf",
  width = 8, height = 5)


make_biplot(       vectorPCOA = out2$vectorPCOA
              , arrow_dataset = out2$df_envfit
             , color_vector = out$meta$sampletype
             , shape_vector = out$meta$pays
             , label_var = "Order", mult = .18
             , title = paste("PCoA bi-plot of", level)
             ) + labs( subtitle = "All samples")

dev.off()

}

#### Make a PERMANOVA to look for factors influencing dispersion of reduced sample set ####
phyloseq::sample_data(dffiltered4_red)$read_count<-rowSums(t(phyloseq::otu_table(dffiltered4_red)))
dffiltered5_red = phyloseq::subset_samples(dffiltered4_red, (sampletype=="feces" | sampletype=="duodenal" | sampletype=="gastric"))

phyloseq::sample_data(dffiltered5_red)$ageyears<-cut(phyloseq::sample_data(dffiltered5_red)$age, c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(phyloseq::sample_data(dffiltered5_red)$ageyears)) # these are the controls

# make labels on variables to better interpret results
phyloseq::sample_data(dffiltered5_red)$ageyears<-factor(phyloseq::sample_data(dffiltered5_red)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))

set.seed(3)
for (i in 1:1) {
 # Subsample
 dffiltered5_red_rar = phyloseq::rarefy_even_depth(dffiltered5_red, sample.size = 5000, verbose = FALSE, replace = TRUE)
}

controls   = c("run","haz","elas_status","anemie","calprotectinelevel","ageyears","read_count","sexe","pays","alphaantitrypsinlevel")
res.adonis = beta_diversity_setup( dffiltered5_red_rar, controls = controls
                                   , pdfname = "Contribution_full_reduced.pdf"
                                   , csvname = "Disperions_full_reduced.csv"
                                   , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

####' line stuff

#### Make a PERMANOVA to look for factors influencing dispersion of CAR sample set ####
dffiltered4_CAR = phyloseq::subset_samples(dffiltered4, pays=="RCA")
dffiltered5_CAR = phyloseq::subset_samples(dffiltered4_CAR, (sampletype=="feces" | sampletype=="duodenal" | sampletype=="gastric"))

phyloseq::sample_data(dffiltered5_CAR)$ageyears<-cut(phyloseq::sample_data(dffiltered5_CAR)$age, c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(phyloseq::sample_data(dffiltered5_CAR)$ageyears)) # these are the controls
levels(phyloseq::sample_data(dffiltered5_CAR)$ageyears)

# make labels on variables to better interpret results
phyloseq::sample_data(dffiltered5_CAR)$ageyears<-factor(phyloseq::sample_data(dffiltered5_CAR)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))
#- 1897
set.seed(3)
for (i in 1:1) {
 # Subsample
 dffiltered5_CAR_rar  = phyloseq::rarefy_even_depth(dffiltered5_CAR, sample.size = 5000, verbose = FALSE, replace = TRUE)
 }


controls   = c("run","haz","sampletype","anemie","calprotectinelevel","ageyears","read_count","sexe","alphaantitrypsinlevel")
res.adonis = beta_diversity_setup( dffiltered5_CAR_rar, controls = controls
                                   , pdfname = "ContributionfulldatasetCAR.pdf"
                                   , csvname = "DisperionstestallsamplesCAR.csv"
                                   , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

##### Plot Ordination with significant taxa as explanatory arrows on Order level- sampletype, CAR only ####
for ( level in c( "Order")){

out = setup_ordination( dataset = dffiltered5_CAR_rar, level = level )
out2= setup_pcoa_scores( dist = out$dist, test4 = out$test4, meta = out$meta, sig_level = 0.01 )

pdf(file="Genussampletype16S_explantory_arrowsrelabun.pdf",
  width = 8, height = 5)

make_biplot(       vectorPCOA = out2$vectorPCOA
              , arrow_dataset = out2$df_envfit
             , color_vector = out$meta$sampletype
             , shape_vector = out$meta$stunted
             , label_var = "Order", mult = .2
             , title = paste("PCoA bi-plot of", level)
             ) + labs( subtitle = "in CAR")

dev.off()

}

#### Make a PERMANOVA to look for factors influencing dispersion of Mada sample set ####
dffiltered4_Mada = phyloseq::subset_samples(dffiltered4, pays == "Madagascar")
dffiltered5_Mada = phyloseq::subset_samples(dffiltered4_Mada, (sampletype=="feces" | sampletype=="duodenal" | sampletype=="gastric"))

phyloseq::sample_data(dffiltered5_Mada)$ageyears<-cut(phyloseq::sample_data(dffiltered5_Mada)$age, c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(phyloseq::sample_data(dffiltered5_Mada)$ageyears)) # these are the controls
levels(phyloseq::sample_data(dffiltered5_Mada)$ageyears)

# make labels on variables to better interpret results
phyloseq::sample_data(dffiltered5_Mada)$ageyears<-factor(phyloseq::sample_data(dffiltered5_Mada)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))

set.seed(3)
for (i in 1:1) {
 # Subsample
 dffiltered5_Mada_rar  = phyloseq::rarefy_even_depth(dffiltered5_Mada, sample.size = 5000, verbose = FALSE, replace = TRUE)
}

#- 1966

controls   = c("run","haz","sampletype","anemie","calprotectinelevel","ageyears"
               ,"read_count","sexe","alphaantitrypsinlevel")

res.adonis = beta_diversity_setup( dffiltered5_Mada_rar, controls = controls
                                   , pdfname = "Contribution_fulldata_mada.pdf"
                                   , csvname = "Disperions_testall_mada.csv"
                                   , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

##### Plot Ordination with significant taxa as explanatory arrows on Order level- sampletype, Mada only ####
for ( level in c( "Order")){

out = setup_ordination( dataset = dffiltered5_Mada_rar, level = level )
out2= setup_pcoa_scores( dist = out$dist, test4 = out$test4, meta = out$meta , sig_level = 0.01)

pdf(file= "Ordersampletype16S_explantory_arrowsrelabun_Mada2.pdf",
  width = 8, height = 5)

make_biplot(       vectorPCOA = out2$vectorPCOA
              , arrow_dataset = out2$df_envfit
             , color_vector = out$meta$sampletype
             , shape_vector = out$meta$stunted
             , label_var = "Order", mult = .15
             , title = paste("PCoA bi-plot of", level)
             ) + labs( subtitle = "in Madagascar")

dev.off()

}



#### Make a PERMANOVA to look for factors influencing dispersion of fecal sample set ####
phyloseq::sample_data(dffilteredfeces)$read_count<-rowSums(t(phyloseq::otu_table(dffilteredfeces)))
phyloseq::sample_data(dffilteredfeces)$ageyears<-cut(phyloseq::sample_data(dffilteredfeces)$age, c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(phyloseq::sample_data(dffilteredfeces)$ageyears)) # these are the controls
levels(phyloseq::sample_data(dffilteredfeces)$ageyears)

# make labels on variables to better interpret results
phyloseq::sample_data(dffilteredfeces)$ageyears<-factor(phyloseq::sample_data(dffilteredfeces)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))

#- 2012
set.seed(3)
for (i in 1:1) {
 # Subsample
 dffilteredfeces_rar  = phyloseq::rarefy_even_depth(dffilteredfeces, sample.size = 5000, verbose = FALSE, replace = TRUE)
}

controls   = c("run","stunted","elas_status","anemie","calprotectinelevel","ageyears","read_count","sexe","pays","alphaantitrypsinlevel")
res.adonis = beta_diversity_setup( dffilteredfeces_rar, controls = controls
                                   , pdfname = "Contribution_feces2.pdf"
                                   , csvname = "Disperions_test_feces2.csv"
                                   , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

##### Plot Ordination with significant taxa as explanatory arrows on Order level- pays ####

out = setup_ordination( dataset = dffilteredfeces_rar, level = "Order")
out2= setup_pcoa_scores( dist = out$dist, test4 = out$test4, meta = out$meta, sig_level = 0.01 )

pdf(file="Orderpays_RCA_vs_Madagascar6S_fecal.pdf",
  width = 5,
  height = 5)

make_biplot(       vectorPCOA = out2$vectorPCOA
              , arrow_dataset = out2$df_envfit
             , color_vector = out$meta$pays
             , shape_vector = out$meta$stunted
             , label_var = "Order", mult = .15
             , title = "PCoA bi-plot of Order"
             ) + labs( subtitle = "Fecal Samples")

dev.off()


#### Make a PERMANOVA to look for factors influencing dispersion of fecal sample set CAR ####
phyloseq::sample_data(dffilteredfecesBangui)$read_count<-rowSums(t(phyloseq::otu_table(dffilteredfecesBangui)))
phyloseq::sample_data(dffilteredfecesBangui)$ageyears<-cut(phyloseq::sample_data(dffilteredfecesBangui)$age, c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(phyloseq::sample_data(dffilteredfecesBangui)$ageyears)) # these are the controls
levels(phyloseq::sample_data(dffilteredfecesBangui)$ageyears)

# make labels on variables to better interpret results
phyloseq::sample_data(dffilteredfecesBangui)$ageyears<-factor(phyloseq::sample_data(dffilteredfecesBangui)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))

set.seed(3)
for (i in 1:1) {
 # Subsample
 dffilteredfecesBangui_rar  = phyloseq::rarefy_even_depth(dffilteredfecesBangui, sample.size = 5000, verbose = FALSE, replace = TRUE)
}

#- 3198
controls   = c("run","stunted","elas_status","anemie","calprotectinelevel","ageyears","read_count","sexe","alphaantitrypsinlevel")
res.adonis = beta_diversity_setup( dffilteredfecesBangui_rar, controls = controls
                                   , pdfname = "Contribution_feces_CAR.pdf"
                                   , csvname = "Disperionstestall_feces_car.csv"
                                   , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

#### Make a PERMANOVA to look for factors influencing dispersion of fecal sample set Mada ####
phyloseq::sample_data(dffilteredfecesMada)$read_count<-rowSums(t(phyloseq::otu_table(dffilteredfecesMada)))
phyloseq::sample_data(dffilteredfecesMada)$ageyears<-cut(phyloseq::sample_data(dffilteredfecesMada)$age, c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(phyloseq::sample_data(dffilteredfecesMada)$ageyears)) # these are the controls
levels(phyloseq::sample_data(dffilteredfecesMada)$ageyears)

# make labels on variables to better interpret results
phyloseq::sample_data(dffilteredfecesMada)$ageyears<-factor(phyloseq::sample_data(dffilteredfecesMada)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))
dffilteredfecesMadabf<-subset_samples(dffilteredfecesMada, current_breastfeeding!="")

set.seed(3)
for (i in 1:1) {
 # Subsample
 dffilteredfecesMada_rar  = phyloseq::rarefy_even_depth(dffilteredfecesMadabf, sample.size = 5000, verbose = FALSE, replace = TRUE)
}


#- 2089
controls   = c("run","stunted","elas_status","anemie","calprotectinelevel","ageyears","read_count","sexe","alphaantitrypsinlevel", "current_breastfeeding")
res.adonis = beta_diversity_setup( dffilteredfecesMada_rar, controls = controls
                                   , pdfname = "Contribution_feces_mada.pdf"
                                   , csvname = "Disperions_feces_mada.csv"
                                   , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

#### Make a PERMANOVA to look for factors influencing dispersion of duodenal sample set ####
phyloseq::sample_data(dffilteredduodenal)$read_count<-rowSums(t(phyloseq::otu_table(dffilteredduodenal)))
phyloseq::sample_data(dffilteredduodenal)$ageyears<-cut(phyloseq::sample_data(dffilteredduodenal)$age, c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(phyloseq::sample_data(dffilteredduodenal)$ageyears)) # these are the controls
levels(phyloseq::sample_data(dffilteredduodenal)$ageyears)

# make labels on variables to better interpret results
phyloseq::sample_data(dffilteredduodenal)$ageyears<-factor(phyloseq::sample_data(dffilteredduodenal)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))

set.seed(3)
for (i in 1:1) {
 # Subsample
 dffilteredduodenal_rar  = phyloseq::rarefy_even_depth(dffilteredduodenal, sample.size = 5000, verbose = FALSE, replace = TRUE)
 }

controls   = c("run","sibo","haz","elas_status","anemie","ageyears","read_count","sexe","pays")
res.adonis = beta_diversity_setup( dffilteredduodenal_rar, controls = controls
                                   , pdfname = "Contribution_duodenal.pdf"
                                   , csvname = "Disperions_duodenal.csv"
                                   , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

##### Plot Ordination with significant taxa as explanatory arrows on Order level- pays ####

out = setup_ordination( dataset = dffilteredduodenal_rar, level = "Order")

out2= setup_pcoa_scores( dist = out$dist, test4 = out$test4, meta = out$meta, sig_level = 0.01 )

pdf(file="Orderpays_RCA_vs_Madagascar6S_explantory_arrowsrelabun_duodenal.pdf",
  width = 8, height = 5)

make_biplot(       vectorPCOA = out2$vectorPCOA
              , arrow_dataset = out2$df_envfit
             , color_vector = out$meta$pays
             , shape_vector = out$meta$stunted
             , label_var = "Order", mult = 0.1
             , title = "PCoA bi-plot of Order"
             ) + labs( subtitle = "Duodenal Samples")

dev.off()


##### Plot Ordination with significant taxa as explanatory arrows on Family level- pays ####

out = setup_ordination( dataset = dffilteredduodenal_rar, level = "Family")
out2= setup_pcoa_scores( dist = out$dist, test4 = out$test4, meta = out$meta , sig_level = 0.01)

pdf(file="Familypays_RCA_vs_Madagascar6S_explantory_arrowsrelabun_duodenal.pdf",
  width = 5,
  height = 5)

make_biplot(       vectorPCOA = out2$vectorPCOA
              , arrow_dataset = out2$df_envfit
             , color_vector = out$meta$pays
             , shape_vector = out$meta$sampletype
             , label_var = "Order", mult = .2
             , title = "PCoA bi-plot of Family"
             ) + labs( subtitle = "Duodenal Samples")

dev.off()


#### Make a PERMANOVA to look for factors influencing dispersion of duodenal sample set without sibo as co-factor ####
# set.seed(3)
# for (i in 1:1) {
#  # Subsample
#  dffilteredduodenal_rar  = phyloseq::rarefy_even_depth(dffilteredduodenal, sample.size = 5000, verbose = FALSE, replace = TRUE)
# }

controls   = c("run","haz","elas_status","anemie","ageyears","read_count","sexe","pays")
res.adonis = beta_diversity_setup( dffilteredduodenal_rar, controls = controls
                                   , pdfname = "Contribution_duodenal_nosibo.pdf"
                                   , csvname = "Disperions_duodenal_nosibo.csv"
                                   , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

controls   = c("run","haz","elas_status","anemie","ageyears","read_count","sexe")

#- 3310
res.adonis = beta_diversity_setup(
  # RCA subset
  dffilteredduodenal_rar %>% phyloseq::subset_samples( pays == "RCA")
                                   , controls = controls
                                   , pdfname = "Contribution_duodenal_nosibo_CAR.pdf"
                                   , csvname = "Disperions_duodenal_nosibo_CAR.csv"
                                   , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

res.adonis = beta_diversity_setup(
  # Madagascar subset
  dffilteredduodenal_rar %>% phyloseq::subset_samples( pays == "Madagascar")
                                   , controls = controls
                                   , pdfname = "Contribution_duodenal_nosibo_Mada.pdf"
                                   , csvname = "Disperions_duodenal_nosibo_Mada.csv"
                                   , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

#### Make a PERMANOVA to look for factors influencing dispersion of gastric sample set ####
phyloseq::sample_data(dffilteredgastric)$read_count<-rowSums(t(phyloseq::otu_table(dffilteredgastric)))
phyloseq::sample_data(dffilteredgastric)$ageyears<-cut(phyloseq::sample_data(dffilteredgastric)$age, c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(phyloseq::sample_data(dffilteredgastric)$ageyears)) # these are the controls
levels(phyloseq::sample_data(dffilteredgastric)$ageyears)

# make labels on variables to better interpret results
phyloseq::sample_data(dffilteredgastric)$ageyears<-factor(phyloseq::sample_data(dffilteredgastric)$ageyears, labels =c("2-3 years", "3-4 years", "4-5 years"))

set.seed(3)
for (i in 1:1) {
 # Subsample
 dffilteredgastric_rar  = phyloseq::rarefy_even_depth(dffilteredgastric, sample.size = 5000, verbose = FALSE, replace = TRUE)
}

#- 3343
controls   = c("run","haz","elas_status","anemie","ageyears","read_count","sexe","pays")
res.adonis = beta_diversity_setup( dffilteredgastric_rar, controls = controls
                                   , pdfname = "Contribution_gastric.pdf"
                                   , csvname = "Disperions_gastric.csv"
                                   , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

#now the adonis test to see if there is a signficant difference according to different variables with sibo
controls   = c("run","haz","elas_status","anemie","ageyears","read_count","sexe","pays","sibo")
res.adonis = beta_diversity_setup( dffilteredgastric_rar, controls = controls
                                   , pdfname = "Contribution_gastric_sibo.pdf"
                                   , csvname = "Disperions_gastric_sibo.csv"
                                   , save_pdf = TRUE, write_csv = TRUE, perms = 499 )

##### Plot Ordination with significant taxa as explanatory arrows on Order level- pays, gastric ####
out = setup_ordination( dataset = dffilteredgastric_rar, level = "Order")
out2= setup_pcoa_scores( dist = out$dist, test4 = out$test4, meta = out$meta )

pdf(file="Orderpays_RCA_vs_Madagascar6S_explantory_order_gastric.pdf",
  width = 5,
  height = 5)

make_biplot(       vectorPCOA = out2$vectorPCOA
              , arrow_dataset = out2$df_envfit
             , color_vector = out$meta$pays
             , shape_vector = out$meta$sampletype
             , label_var = "Order", mult = .15
             , title = "PCoA bi-plot of Order"
             ) + labs( subtitle = "Gastric Samples")

dev.off()


##### Plot Ordination with significant taxa as explanatory arrows on Family level- pays, gastric ####
out = setup_ordination( dataset = dffilteredgastric_rar, level = "Family")
out2= setup_pcoa_scores( dist = out$dist, test4 = out$test4, meta = out$meta, sig_level = 0.01 )

pdf(file="Orderpays_RCA_vs_Madagascar6S_explantory_family_gastric.pdf",
  width = 5,
  height = 5)

make_biplot(       vectorPCOA = out2$vectorPCOA
              , arrow_dataset = out2$df_envfit
             , color_vector = out$meta$pays
             , shape_vector = out$meta$sampletype
             , label_var = "Order", mult = .15
             , title = "PCoA bi-plot of Family"
             ) + labs( subtitle = "Gastric Samples")

dev.off()



##### Analyse the core microbiota in the duodenum using the microbiome package #####
#### first define the core microbiota of each individual compartment: ASV level ####
## gastric samples
# find taxa which are at least shared by 90% of all samples and have an aundance of at least 0.01%
dffilteredgastric_rel = phyloseq::transform_sample_counts(dffilteredgastric, function(x) x/sum(x))
core.taxa.gastric = microbiome::core(dffilteredgastric_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.gastric # 9 taxa

core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))
#View(core.taxa.gastric.mat)
write.csv(core.taxa.gastric.mat, "Coretaxagastricmerged.csv")


####' line stuff

# Bangui only gastric
dffilteredgastricBangui_rel = phyloseq::transform_sample_counts(dffilteredgastricBangui, function(x) x/sum(x))
core.taxa.gastric = microbiome::core(dffilteredgastricBangui_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.gastric # 4 taxa

core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))
#View(core.taxa.gastric.mat)
write.csv(core.taxa.gastric.mat, "CoretaxagastricmergedBangui.csv")

# Tana only gastric
dffilteredgastricTana_rel = phyloseq::transform_sample_counts(dffilteredgastricMada, function(x) x/sum(x))
core.taxa.gastric = microbiome::core(dffilteredgastricTana_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.gastric # 13 taxa

core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))
#View(core.taxa.gastric.mat)
write.csv(core.taxa.gastric.mat, "CoretaxagastricmergedTana.csv")

# find taxa which are at least shared by 75% of all samples and have an aundance of at least 0.01%

dffilteredgastric_rel = phyloseq::transform_sample_counts(dffilteredgastric, function(x) x/sum(x))
core.taxa.gastric = microbiome::core(dffilteredgastric_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.gastric # 23 taxa

core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))
#View(core.taxa.gastric.mat)
write.csv(core.taxa.gastric.mat, "Coretaxagastricmerged75.csv")

# Bangui only gastric
dffilteredgastricBangui_rel = phyloseq::transform_sample_counts(dffilteredgastricBangui, function(x) x/sum(x))
core.taxa.gastric = microbiome::core(dffilteredgastricBangui_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.gastric # 19 taxa

core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))
#View(core.taxa.gastric.mat)
write.csv(core.taxa.gastric.mat, "CoretaxagastricmergedBangui75.csv")

# Tana only gastric
dffilteredgastricTana_rel = phyloseq::transform_sample_counts(dffilteredgastricMada, function(x) x/sum(x))
core.taxa.gastric = microbiome::core(dffilteredgastricTana_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.gastric # 29 taxa

core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))
#View(core.taxa.gastric.mat)
write.csv(core.taxa.gastric.mat, "CoretaxagastricmergedTana75.csv")

## duodenal samples
# find taxa which are at least shared by 90% of all samples and have an aundance of at least 0.01%
dffilteredduodenal_rel = phyloseq::transform_sample_counts(dffilteredduodenal, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredduodenal_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.duodenal # 9 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "Coretaxaduodenalmerged.csv")

# Bangui only duodenal
dffilteredduodenalBangui_rel = phyloseq::transform_sample_counts(dffilteredduodenalBangui, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredduodenalBangui_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.duodenal # 8 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "CoretaxaduodenalmergedBangui.csv")

# Tana only duodenal
dffilteredduodenalTana_rel = phyloseq::transform_sample_counts(dffilteredduodenalMada, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredduodenalTana_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.duodenal # 13 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "CoretaxaduodenalmergedTana.csv")

# find taxa which are at least shared by 75% of all samples and have an aundance of at least 0.01%

dffilteredduodenal_rel = phyloseq::transform_sample_counts(dffilteredduodenal, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredduodenal_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal # 26 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
write.csv(core.taxa.duodenal.mat, "Coretaxaduodenalmerged75.csv")

#' line
# Bangui only duodenal
dffilteredduodenalBangui_rel = phyloseq::transform_sample_counts(dffilteredduodenalBangui, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredduodenalBangui_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal # 23 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "CoretaxaduodenalmergedBangui75.csv")

# Tana only duodenal
dffilteredduodenalTana_rel = phyloseq::transform_sample_counts(dffilteredduodenalMada, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredduodenalTana_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal # 30 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "CoretaxaduodenalmergedTana75.csv")


## fecal samples
# find taxa which are at least shared by 90% of all samples and have an aundance of at least 0.01%
dffilteredfeces_rel = phyloseq::transform_sample_counts(dffilteredfeces, function(x) x/sum(x))
core.taxa.feces = microbiome::core(dffilteredfeces_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.feces # 3 taxa

core.taxa.feces.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces))
#View(core.taxa.feces.mat)
write.csv(core.taxa.feces.mat, "Coretaxafecesmerged.csv")

# Bangui only feces
dffilteredfecesBangui_rel = phyloseq::transform_sample_counts(dffilteredfecesBangui, function(x) x/sum(x))
core.taxa.feces = microbiome::core(dffilteredfecesBangui_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.feces # 3 taxa

core.taxa.feces.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces))
write.csv(core.taxa.feces.mat, "CoretaxafecesmergedBangui.csv")

# Tana only feces
dffilteredfecesTana_rel = phyloseq::transform_sample_counts(dffilteredfecesMada, function(x) x/sum(x))
core.taxa.feces = microbiome::core(dffilteredfecesTana_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.feces # 4 taxa

core.taxa.feces.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces))
write.csv(core.taxa.feces.mat, "CoretaxafecesmergedTana.csv")

# find taxa which are at least shared by 75% of all samples and have an aundance of at least 0.01%

dffilteredfeces_rel = phyloseq::transform_sample_counts(dffilteredfeces, function(x) x/sum(x))
core.taxa.feces = microbiome::core(dffilteredfeces_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.feces # 11 taxa

core.taxa.feces.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces))
#View(core.taxa.feces.mat)
write.csv(core.taxa.feces.mat, "Coretaxafecesmerged75.csv")

# Bangui only feces
dffilteredfecesBangui_rel = phyloseq::transform_sample_counts(dffilteredfecesBangui, function(x) x/sum(x))
core.taxa.feces = microbiome::core(dffilteredfecesBangui_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.feces # 10 taxa

core.taxa.feces.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces))
#View(core.taxa.feces.mat)
write.csv(core.taxa.feces.mat, "CoretaxafecesmergedBangui75.csv")

# Tana only feces
dffilteredfecesTana_rel = phyloseq::transform_sample_counts(dffilteredfecesMada, function(x) x/sum(x))
core.taxa.feces = microbiome::core(dffilteredfecesTana_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.feces # 11 taxa

core.taxa.feces.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces))
#View(core.taxa.feces.mat)
write.csv(core.taxa.feces.mat, "CoretaxafecesmergedTana75.csv")


# fecal samples subset to 137 samples (to match with duodenal samples)
length(which(phyloseq::sample_data(dffiltered4)$sampletype=="duodenal")) #129 samples
length(which(phyloseq::sample_data(dffiltered4)$sampletype=="gastric")) #168 samples
length(which(phyloseq::sample_data(dffiltered4)$sampletype=="feces")) #627 samples, reduce and restrict on malnourished samples


dffilteredfecesmaln_rel = phyloseq::subset_samples(dffilteredfeces_rel, stunted="stunted")
dataframe=dffilteredfecesmaln_rel

dataframe=as.data.frame(phyloseq::sample_data(dataframe))

subset = dataframe[sample(nrow(dataframe), 129), ]
vector = as.vector(subset$ID_metag)


dffilteredfecesmaln129only_rel=prune_samples(vector, dffilteredfecesmaln_rel)
dffilteredfecesmaln129only_rel

core.taxa.feces129 = microbiome::core(dffilteredfecesmaln129only_rel, detection = 0.0001, prevalence = 90/100)
core.taxa.feces129 # 3 taxa

core.taxa.feces.129.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces129))
#View(core.taxa.feces.129.mat)
write.csv(core.taxa.feces.129.mat, "Coretaxafecesmerged129.csv")

core.taxa.feces129 = microbiome::core(dffilteredfecesmaln129only_rel, detection = 0.0001, prevalence = 75/100)
core.taxa.feces129 # 13 taxa

core.taxa.feces.129.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces129))
#View(core.taxa.feces.129.mat)
write.csv(core.taxa.feces.129.mat, "Coretaxafecesmerged129_75.csv")

## now see who is conserved in between different compartments with 75% threshold, merged dataset
core.taxa.gastric = microbiome::core(dffilteredgastric_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.gastric
core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))

core.taxa.duodenal = microbiome::core(dffilteredduodenal_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal
core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))


commongastduod <- cbind(core.taxa.gastric.mat[ intersect(rownames(core.taxa.gastric.mat), rownames(core.taxa.duodenal.mat)), ])
#View(commongastduod)
nrow(commongastduod) # 21 taxa are shared between gastric and duodenal

commongastfeces <- cbind( core.taxa.gastric.mat[ intersect(rownames(core.taxa.gastric.mat), rownames(core.taxa.feces.129.mat)), ])
#View(commongastfeces)
nrow(commongastfeces) # 2 taxa are shared between gastric and feces

commonduodfeces <- cbind( core.taxa.duodenal.mat[ intersect(rownames(core.taxa.duodenal.mat), rownames(core.taxa.feces.129.mat)), ])
#View(commonduodfeces)
nrow(commonduodfeces) # 2 taxa are shared between duodenal and feces


#### first define the core microbiota of each individual compartment: Species level ####
## phyloseq::tax_glom on Species level

dffilteredSpeciesgastric = phyloseq::tax_glom(dffilteredgastric, "Species")
dffilteredSpeciesduodenal = phyloseq::tax_glom(dffilteredduodenal, "Species")
dffilteredSpeciesfeces = phyloseq::tax_glom(dffilteredfeces, "Species")

dffilteredSpeciesgastricBangui = phyloseq::tax_glom(dffilteredgastricBangui, "Species")
dffilteredSpeciesduodenalBangui = phyloseq::tax_glom(dffilteredduodenalBangui, "Species")
dffilteredSpeciesfecesBangui = phyloseq::tax_glom(dffilteredfecesBangui, "Species")

dffilteredSpeciesgastricMada = phyloseq::tax_glom(dffilteredgastricMada, "Species")
dffilteredSpeciesduodenalMada = phyloseq::tax_glom(dffilteredduodenalMada, "Species")
dffilteredSpeciesfecesMada = phyloseq::tax_glom(dffilteredfecesMada, "Species")

## gastric samples

# find taxa which are at least shared by 90% of all samples and have an aundance of at least 0.01%
dffilteredSpeciesgastric_rel = phyloseq::transform_sample_counts(dffilteredSpeciesgastric, function(x) x/sum(x))
core.taxa.gastric = microbiome::core(dffilteredSpeciesgastric_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.gastric # 9 taxa

core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))
#View(core.taxa.gastric.mat)
write.csv(core.taxa.gastric.mat, "CoretaxagastricmergedSpecies.csv")

# Bangui only gastric
dffilteredSpeciesgastricBangui_rel = phyloseq::transform_sample_counts(dffilteredSpeciesgastricBangui, function(x) x/sum(x))
core.taxa.gastric = microbiome::core(dffilteredSpeciesgastricBangui_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.gastric # 3 taxa

core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))
#View(core.taxa.gastric.mat)
write.csv(core.taxa.gastric.mat, "CoretaxagastricmergedBanguiSpecies.csv")

# Tana only gastric
dffilteredSpeciesgastricTana_rel = phyloseq::transform_sample_counts(dffilteredSpeciesgastricMada, function(x) x/sum(x))
core.taxa.gastric = microbiome::core(dffilteredSpeciesgastricTana_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.gastric # 15 taxa

core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))
#View(core.taxa.gastric.mat)
write.csv(core.taxa.gastric.mat, "CoretaxagastricmergedTanaSpecies.csv")

# find taxa which are at least shared by 75% of all samples and have an aundance of at least 0.01%

dffilteredSpeciesgastric_rel = phyloseq::transform_sample_counts(dffilteredSpeciesgastric, function(x) x/sum(x))
core.taxa.gastric = microbiome::core(dffilteredSpeciesgastric_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.gastric # 23 taxa

core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))
#View(core.taxa.gastric.mat)
write.csv(core.taxa.gastric.mat, "Coretaxagastricmerged75Species.csv")

# Bangui only gastric
dffilteredSpeciesgastricBangui_rel = phyloseq::transform_sample_counts(dffilteredSpeciesgastricBangui, function(x) x/sum(x))
core.taxa.gastric = microbiome::core(dffilteredSpeciesgastricBangui_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.gastric # 18 taxa

core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))
#View(core.taxa.gastric.mat)
write.csv(core.taxa.gastric.mat, "CoretaxagastricmergedBangui75Species.csv")

# Tana only gastric
dffilteredSpeciesgastricTana_rel = phyloseq::transform_sample_counts(dffilteredSpeciesgastricMada, function(x) x/sum(x))
core.taxa.gastric = microbiome::core(dffilteredSpeciesgastricTana_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.gastric # 26 taxa

core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))
#View(core.taxa.gastric.mat)
write.csv(core.taxa.gastric.mat, "CoretaxagastricmergedTana75Species.csv")


## duodenal samples
# find taxa which are at least shared by 90% of all samples and have an abundance of at least 0.01%
dffilteredSpeciesduodenal_rel = phyloseq::transform_sample_counts(dffilteredSpeciesduodenal, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenal_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.duodenal # 9 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "CoretaxaduodenalmergedSpecies.csv")

# Bangui only duodenal
dffilteredSpeciesduodenalBangui_rel = phyloseq::transform_sample_counts(dffilteredSpeciesduodenalBangui, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenalBangui_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.duodenal # 7 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "CoretaxaduodenalmergedBanguiSpecies.csv")

# Tana only duodenal
dffilteredSpeciesduodenalTana_rel = phyloseq::transform_sample_counts(dffilteredSpeciesduodenalMada, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenalTana_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.duodenal # 13 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "CoretaxaduodenalmergedTanaSpecies.csv")

# find taxa which are at least shared by 75% of all samples and have an aundance of at least 0.01%

dffilteredSpeciesduodenal_rel = phyloseq::transform_sample_counts(dffilteredSpeciesduodenal, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenal_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal # 24 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "Coretaxaduodenalmerged75Species.csv")

# Bangui only duodenal
dffilteredSpeciesduodenalBangui_rel = phyloseq::transform_sample_counts(dffilteredSpeciesduodenalBangui, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenalBangui_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal # 25 taxa
#- 2606 (3752)
core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "CoretaxaduodenalmergedBangui75Species.csv")

# Tana only duodenal
dffilteredSpeciesduodenalTana_rel = phyloseq::transform_sample_counts(dffilteredSpeciesduodenalMada, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenalTana_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal # 28 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "CoretaxaduodenalmergedTana75Species.csv")

## now by sibo status
# find taxa which are at least shared by 75% of all samples and have an aundance of at least 0.01%, with sibo
dffilteredSpeciesduodenalsibo = phyloseq::subset_samples(dffilteredSpeciesduodenal, sibo=="yes")
dffilteredSpeciesduodenal_rel = phyloseq::transform_sample_counts(dffilteredSpeciesduodenalsibo, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenal_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal # 24 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(dffilteredSpeciesduodenalsibo))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "Coretaxaduodenalmerged75Speciessibo.csv")

# find taxa which are at least shared by 75% of all samples and have an aundance of at least 0.01%, no sibo
dffilteredSpeciesduodenalnosibo = phyloseq::subset_samples(dffilteredSpeciesduodenal, sibo=="no")
dffilteredSpeciesduodenal_rel = phyloseq::transform_sample_counts(dffilteredSpeciesduodenalnosibo, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenal_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal # 24 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "Coretaxaduodenalmerged75Speciesnosibo.csv")

# Bangui only duodenal, sibo
dffilteredSpeciesduodenalBanguisibo = phyloseq::subset_samples(dffilteredSpeciesduodenalBangui, sibo=="yes")
dffilteredSpeciesduodenalBangui_rel = phyloseq::transform_sample_counts(dffilteredSpeciesduodenalBanguisibo, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenalBangui_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal # 25 taxa


#- 2662

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "CoretaxaduodenalmergedBangui75Speciessibo.csv")

# Bangui only duodenal, no sibo
dffilteredSpeciesduodenalBanguinosibo = phyloseq::subset_samples(dffilteredSpeciesduodenalBangui, sibo=="no")
dffilteredSpeciesduodenalBangui_rel = phyloseq::transform_sample_counts(dffilteredSpeciesduodenalBanguinosibo, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenalBangui_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal # 23 taxa

core.taxa.duodenal.mat2 <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat2)
write.csv(core.taxa.duodenal.mat2, "CoretaxaduodenalmergedBangui75Speciesnosibo.csv")

#- 2678

# Tana only duodenal, sibo
dffilteredSpeciesduodenalMadasibo = phyloseq::subset_samples(dffilteredSpeciesduodenalMada, sibo=="yes")
dffilteredSpeciesduodenalTana_rel = phyloseq::transform_sample_counts(dffilteredSpeciesduodenalMadasibo, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenalTana_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal # 28 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "CoretaxaduodenalmergedTana75Speciessibo.csv")


# Tana only duodenal, no sibo
dffilteredSpeciesduodenalMadanosibo = phyloseq::subset_samples(dffilteredSpeciesduodenalMada, sibo=="no")
dffilteredSpeciesduodenalTana_rel = phyloseq::transform_sample_counts(dffilteredSpeciesduodenalMadanosibo, function(x) x/sum(x))
core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenalTana_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal # 33 taxa

core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))
#View(core.taxa.duodenal.mat)
write.csv(core.taxa.duodenal.mat, "CoretaxaduodenalmergedTana75Speciesnosibo.csv")


## fecal samples
# find taxa which are at least shared by 90% of all samples and have an aundance of at least 0.01%
dffilteredSpeciesfeces_rel = phyloseq::transform_sample_counts(dffilteredSpeciesfeces, function(x) x/sum(x))
core.taxa.feces = microbiome::core(dffilteredSpeciesfeces_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.feces # 3 taxa

core.taxa.feces.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces))
write.csv(core.taxa.feces.mat, "CoretaxafecesmergedSpecies.csv")

# Bangui only feces
dffilteredSpeciesfecesBangui_rel = phyloseq::transform_sample_counts(dffilteredSpeciesfecesBangui, function(x) x/sum(x))
core.taxa.feces = microbiome::core(dffilteredSpeciesfecesBangui_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.feces # 5 taxa

core.taxa.feces.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces))
write.csv(core.taxa.feces.mat, "CoretaxafecesmergedBanguiSpecies.csv")

# Tana only feces
dffilteredSpeciesfecesTana_rel = phyloseq::transform_sample_counts(dffilteredSpeciesfecesMada, function(x) x/sum(x))
core.taxa.feces = microbiome::core(dffilteredSpeciesfecesTana_rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE)
core.taxa.feces # 4 taxa

core.taxa.feces.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces))
write.csv(core.taxa.feces.mat, "CoretaxafecesmergedTanaSpecies.csv")

# find taxa which are at least shared by 75% of all samples and have an aundance of at least 0.01%

dffilteredSpeciesfeces_rel = phyloseq::transform_sample_counts(dffilteredSpeciesfeces, function(x) x/sum(x))
core.taxa.feces = microbiome::core(dffilteredSpeciesfeces_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.feces # 11 taxa

core.taxa.feces.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces))
#View(core.taxa.feces.mat)
write.csv(core.taxa.feces.mat, "Coretaxafecesmerged75Species.csv")

# Bangui only feces
dffilteredSpeciesfecesBangui_rel = phyloseq::transform_sample_counts(dffilteredSpeciesfecesBangui, function(x) x/sum(x))
core.taxa.feces = microbiome::core(dffilteredSpeciesfecesBangui_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.feces # 11 taxa

core.taxa.feces.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces))
#View(core.taxa.feces.mat)
write.csv(core.taxa.feces.mat, "CoretaxafecesmergedBangui75Species.csv")

# Tana only feces
dffilteredSpeciesfecesTana_rel = phyloseq::transform_sample_counts(dffilteredSpeciesfecesMada, function(x) x/sum(x))
core.taxa.feces = microbiome::core(dffilteredSpeciesfecesTana_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.feces # 12 taxa

core.taxa.feces.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces))
#View(core.taxa.feces.mat)
write.csv(core.taxa.feces.mat, "CoretaxafecesmergedTana75Species.csv")

dffilteredSpeciesfecesmaln_rel = phyloseq::subset_samples(dffiltered_s, stunted="stunted")
dataframe=dffilteredSpeciesfecesmaln_rel

dataframe=as.data.frame(phyloseq::sample_data(dataframe))

subset= dataframe[sample(nrow(dataframe), 129), ]
vector = as.vector(subset$ID_metag)

dffilteredSpeciesfecesmaln129only_rel=prune_samples(vector, dffilteredSpeciesfecesmaln_rel)
dffilteredSpeciesfecesmaln129only_rel

core.taxa.feces129 = microbiome::core(dffilteredSpeciesfecesmaln129only_rel, detection = 0.0001, prevalence = 90/100)
core.taxa.feces129 # 11 taxa

core.taxa.feces.129.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces129))
#View(core.taxa.feces.129.mat)
write.csv(core.taxa.feces.129.mat, "Coretaxafecesmerged129Species.csv")

core.taxa.feces129 = microbiome::core(dffilteredSpeciesfecesmaln129only_rel, detection = 0.0001, prevalence = 75/100)
core.taxa.feces129 # 7 taxa

core.taxa.feces.129.mat <- as.data.frame(phyloseq::tax_table(core.taxa.feces129))
#View(core.taxa.feces.129.mat)
write.csv(core.taxa.feces.129.mat, "Coretaxafecesmerged129Species75.csv")

## now see who is conserved in between different compartments with 75% threshold, merged dataset
core.taxa.gastric = microbiome::core(dffilteredSpeciesgastric_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.gastric
core.taxa.gastric.mat <- as.data.frame(phyloseq::tax_table(core.taxa.gastric))

core.taxa.duodenal = microbiome::core(dffilteredSpeciesduodenal_rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
core.taxa.duodenal
core.taxa.duodenal.mat <- as.data.frame(phyloseq::tax_table(core.taxa.duodenal))

commongastduod <- cbind(core.taxa.gastric.mat[ intersect(rownames(core.taxa.gastric.mat), rownames(core.taxa.duodenal.mat)), ])
#View(commongastduod)
nrow(commongastduod) # 19 taxa are shared between gastric and duodenal

commongastfeces <- cbind( core.taxa.gastric.mat[ intersect(rownames(core.taxa.gastric.mat), rownames(core.taxa.feces.129.mat)), ])
#View(commongastfeces)
nrow(commongastfeces) # 3

commonduodfeces <- cbind( core.taxa.duodenal.mat[ intersect(rownames(core.taxa.duodenal.mat), rownames(core.taxa.feces.129.mat)), ])
#View(commonduodfeces)
nrow(commonduodfeces) # 3

#- 2811 deseq
#### DESEQ2 ANALYSIS ####
#### make the Deseq analysis on all gastric samples and pays of origin on Genus level all gastric samples ####
level = "Genus"
name  = "gastric"
variable_of_interest = "pays"

res = deseq_cooks( dataset = dffilteredgastric, level = level
             , full_model = c("age","sexe","run", variable_of_interest )
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )

res$results_filtered

merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 1
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

payscorr


#make the actual plot

pdf("differentially_present_taxa_Genus_bycountryoforigin_corr_gender_age_gastric_LRTallsamples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot( payscorr, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) +
  ggtitle("Differences by country of origin (CAR vs. Madagascar), controlled gender, age, LRT model")
dev.off()

# make graph plot
pdf("differentially_present_genera by country of origin_gastric samples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 8, # define plot width and height. completely up to user.
  height = 6)
ggplot(data= payscorr , aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
 geom_bar(position="dodge",stat="identity", color="black") +
 coord_flip() +
 theme(axis.text.x = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 xlab("")+
 theme(axis.title.y = element_text(size=14))+
 theme(title = element_text(size=16, face="bold"))+
 ylab("log2 fold change")+
 ggtitle("Significantly different genera by country of origin in gastric samples")
dev.off()

#- 2858 (4002)
#### make the Deseq analysis on all gastric samples and sibo on Species level all gastric samples, no filter for foldchange ####
level = "Species"
name  = "gastric"
variable_of_interest = "sibo"

res = dffilteredgastric  %>%
         phyloseq::subset_samples( sibo == "no" | sibo == "yes") %>%
         filter_prune( ) %>%
         deseq_cooks( level = level
             , full_model = c("age","sexe","run","haz","pays", variable_of_interest )
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )

res$results_filtered

merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

payscorr


#make the actual plot
#- 2885
pdf("differentially_present_taxa_Species_bysibo_corr_gender_age_gastric_LRTnofilterfoldchangeallsamples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(payscorr, aes(x=Species, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by sibo status (yes vs. no) in gastric samples, controlled gender, age, haz LRT model")
dev.off()

pdf("differentially_present_species by SIBO status_gastric samples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 8, # define plot width and height. completely up to user.
  height = 6)
ggplot(data=payscorr, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
 geom_bar(position="dodge",stat="identity", color="black") +
 coord_flip() +
 theme(axis.text.x = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 xlab("")+
 theme(axis.title.y = element_text(size=14))+
 theme(title = element_text(size=16, face="bold"))+
 ylab("log2 fold change")+
 ggtitle("Significantly different species by SIBO status in gastric samples")
dev.off()


#### make the Deseq analysis on gastric from Bangui samples and sibo on Species level gastric from Bangui samples, no filter for foldchange ####
level = "Species"
name  = "gastric"
variable_of_interest = "sibo"

res = dffilteredgastricBangui %>%
  phyloseq::subset_samples( sibo == "no" | sibo == "yes") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run","haz", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
res$results_filtered

merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrsibogenderage_reduced_prevabdeseq <- payscorr

#- 2938
#make the actual plot





pdf("differentially_present_taxa_Species_bysibo_corr_gender_age_gastric_LRTnofilterfoldchangeallsamplesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrsibogenderage_reduced_prevabdeseq, aes(x=Species, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by sibo status (yes vs. no) in gastric samples, controlled gender, age, haz LRT model")
dev.off()

corrsibogenderage_reduced_prevabdeseq$taxonomy<- paste0(corrsibogenderage_reduced_prevabdeseq$Species
                                                        , "/", corrsibogenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_species by SIBO status_gastric samplesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 8, # define plot width and height. completely up to user.
  height = 6)

differential_plot(corrsibogenderage_reduced_prevabdeseq, "Significantly different species by SIBO status in gastric samples, Bangui" )
dev.off()

#- 2967 (4111)
#### make the Deseq analysis on gastric from Mada samples and sibo on Species level gastric from Mada samples, no filter for foldchange ####
level = "Species"
variable_of_interest = "sibo"
name = "gastric"

res = dffilteredgastricMada %>%
  phyloseq::subset_samples( sibo == "no" | sibo == "yes") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run",variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
res$results_filtered

merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrsibogenderage_reduced_prevabdeseq <- payscorr


#make the actual plot

pdf("differentially_present_taxa_Species_bysibo_corr_gender_age_gastric_LRTnofilterfoldchangeallsamplesMada.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrsibogenderage_reduced_prevabdeseq, aes(x=Species, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by sibo status (yes vs. no) in gastric samples, controlled gender, age, haz LRT model")
dev.off()

corrsibogenderage_reduced_prevabdeseq$taxonomy<- paste0(corrsibogenderage_reduced_prevabdeseq$Species, "/", corrsibogenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_species by SIBO status_gastric samplesTana.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 8, # define plot width and height. completely up to user.
  height = 6)
differential_plot(corrsibogenderage_reduced_prevabdeseq, "Significantly different species by SIBO status in gastric samples, Tana" )

dev.off()

#- 3021
#### duodenal samples ####
#### make the Deseq analysis on all duodenal samples and pays of origin on Genus level all duodenal samples ####
level = "Genus"
variable_of_interest = "pays"
name = "duodenal"

res = deseq_cooks( dataset = dffilteredduodenal
             , level = level
             , full_model = c("age", "sexe", "run",variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
res$results_filtered

merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

sibocorrcountrygenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Genus_bycountryoforigin_corr_gender_age_duodenal_LRTallsamples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(sibocorrcountrygenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by country of origin (CAR vs. Madagascar), controlled gender, age, LRT model")
dev.off()

# make graph plot
sibocorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(sibocorrcountrygenderage_reduced_prevabdeseq$Genus, "/", sibocorrcountrygenderage_reduced_prevabdeseq$Family)
pdf("differentially_present_genera by country of origin_duodenal samples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 6, # define plot width and height. completely up to user.
  height = 1.3)

differential_plot(sibocorrcountrygenderage_reduced_prevabdeseq, "Significantly different genera by country of origin in duodenal samples" )

dev.off()


#### make the Deseq analysis on all duodenal samples and sibo on Species level all duodenal samples, no filter for foldchange ####
level = "Species"
name = "duodenal"
variable_of_interest = "sibo"

res = dffilteredduodenal %>%
  phyloseq::subset_samples( sibo == "no" | sibo == "yes") %>%
  filter_prune( ) %>%
  deseq_cooks( level = "Species"
             , full_model = c("age","sexe","run","haz","pays",  variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
res$results_filtered

merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrsibogenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Species_bysibo_corr_gender_age_duodenal_LRTnofilterfoldchangeallsamples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrsibogenderage_reduced_prevabdeseq, aes(x=Species, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by sibo status (yes vs. no) in duodenal samples, controlled gender, age, haz LRT model")
dev.off()

corrsibogenderage_reduced_prevabdeseq$taxonomy<- paste0(corrsibogenderage_reduced_prevabdeseq$Species, "/", corrsibogenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_species by SIBO status_duodenal samples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 8, # define plot width and height. completely up to user.
  height = 6)


differential_plot(corrsibogenderage_reduced_prevabdeseq,"Significantly different species by SIBO status in duodenal samples" )
dev.off()

#- 3127
#### make the Deseq analysis on all duodenal samples and sibo on Genus level all duodenal samples, no filter for foldchange ####
level = "Genus"
variable_of_interest = "sibo"
name = "duodenal"

res = dffilteredduodenal %>%
  phyloseq::subset_samples( sibo == "no" | sibo == "yes") %>%
  filter_prune( ) %>%
  deseq_cooks( level = "Genus"
             , full_model = c("age","sexe","run", variable_of_interest )
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )

res$results_filtered

merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrsibogenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Genus_bysibo_corr_gender_age_duodenal_LRTnofilterfoldchangeallsamples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrsibogenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by sibo status (yes vs. no) in duodenal samples, controlled gender, age, haz LRT model")
dev.off()

corrsibogenderage_reduced_prevabdeseq$taxonomy<- paste0(corrsibogenderage_reduced_prevabdeseq$Genus, "/", corrsibogenderage_reduced_prevabdeseq$Family)
pdf("differentially_present_genera by SIBO status_duodenal samples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 8, # define plot width and height. completely up to user.
  height = 6)


differential_plot(corrsibogenderage_reduced_prevabdeseq,"Significantly different genera by SIBO status in duodenal samples" )
dev.off()



####make the Deseq analysis on duodenal from Bangui samples and sibo on Species level duodenal from Bangui samples, no filter for foldchange ####
level = "Species"
name= "duodenal_mada"
variable_of_interest = "sibo"

res = dffilteredduodenalMada %>%
  phyloseq::subset_samples( sibo == "no" | sibo == "yes") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )

res$results_filtered

merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrsibogenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Species_bysibo_corr_gender_age_duodenal_LRTnofilterfoldchangeallsamplesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrsibogenderage_reduced_prevabdeseq, aes(x=Species, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by sibo status (yes vs. no) in duodenal samples, controlled gender, age, haz LRT model")
dev.off()

corrsibogenderage_reduced_prevabdeseq$taxonomy<- paste0(corrsibogenderage_reduced_prevabdeseq$Species, "/", corrsibogenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_species by SIBO status_duodenal samplesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 8, # define plot width and height. completely up to user.
  height = 6)

differential_plot(corrsibogenderage_reduced_prevabdeseq,"Significantly different species by SIBO status in duodenal samples, Bangui" )
dev.off()



#### make the Deseq analysis on duodenal from Bangui samples and sibo on Genus level duodenal from Bangui samples, no filter for foldchange ####
level = "Genus"
name= "duodenal_mada"
variable_of_interest = "sibo"

res = dffilteredduodenalBangui %>%
  phyloseq::subset_samples( sibo == "no" | sibo == "yes") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
res$results_filtered

merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrsibogenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Genus_bysibo_corr_gender_age_duodenal_LRTnofilterfoldchangeallsamplesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrsibogenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by sibo status (yes vs. no) in duodenal samples, controlled gender, age, haz LRT model")
dev.off()

corrsibogenderage_reduced_prevabdeseq$taxonomy<- paste0(corrsibogenderage_reduced_prevabdeseq$Genus, "/", corrsibogenderage_reduced_prevabdeseq$Family)
pdf("differentially_present_genera by SIBO status_duodenal samplesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 8, # define plot width and height. completely up to user.
  height = 6)


differential_plot(corrsibogenderage_reduced_prevabdeseq,"Significantly different genera by SIBO status in duodenal samples, Bangui" )
dev.off()



#### make the Deseq analysis on duodenal from Mada samples and sibo on Species level duodenal from Mada samples, no filter for foldchange ####
level = "Species"
name = "duodnenal_mada"
variable_of_interest = "sibo"
res = dffilteredduodenalMada %>%
  phyloseq::subset_samples( sibo == "no" | sibo == "yes") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             , fold_change = 0
             )
res$results_filtered

merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrsibogenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Species_bysibo_corr_gender_age_duodenal_LRTnofilterfoldchangeallsamplesMada.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrsibogenderage_reduced_prevabdeseq, aes(x=Species, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by sibo status (yes vs. no) in duodenal samples, controlled gender, age, haz LRT model")
dev.off()

corrsibogenderage_reduced_prevabdeseq$taxonomy<- paste0(corrsibogenderage_reduced_prevabdeseq$Species, "/", corrsibogenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_species by SIBO status_duodenal samplesTana.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 8, # define plot width and height. completely up to user.
  height = 6)


differential_plot(corrsibogenderage_reduced_prevabdeseq, "Significantly different species by SIBO status in duodenal samples, Tana" )
dev.off()



#### make the Deseq analysis on duodenal from Mada samples and sibo on Genus level duodenal from Mada samples, no filter for foldchange ####
 level = "Genus"
name = "duodenal_mada"
variable_of_interest = "sibo"

res = dffilteredduodenalMada %>%
  phyloseq::subset_samples( sibo == "no" | sibo == "yes") %>%
  filter_prune( ) %>%
  deseq_cooks( level = "Genus"
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
res$results_filtered

merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrsibogenderage_reduced_prevabdeseq <- payscorr
#make the actual plot

pdf("differentially_present_taxa_Genus_bysibo_corr_gender_age_duodenal_LRTnofilterfoldchangeallsamplesMada.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrsibogenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by sibo status (yes vs. no) in duodenal samples, controlled gender, age, haz LRT model")
dev.off()

corrsibogenderage_reduced_prevabdeseq$taxonomy<- paste0(corrsibogenderage_reduced_prevabdeseq$Genus, "/", corrsibogenderage_reduced_prevabdeseq$Family)
pdf("differentially_present_genera by SIBO status_duodenal samples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 8, # define plot width and height. completely up to user.
  height = 6)


differential_plot(corrsibogenderage_reduced_prevabdeseq,"Significantly different genera by SIBO status in duodenal samples, Bangui" )
dev.off()



#- 3399
#### fecal samples ####
#### make the Deseq analysis on all feces samples and pays of origin on Genus level all feces samples ####
level = "Genus"
name = "feces"
variable_of_interest = "pays"
res = dffilteredfeces %>%
  phyloseq::subset_samples( sibo == "no" | sibo == "yes") %>%
  filter_prune( ) %>%
  deseq_cooks( level = "Genus"
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

sibocorrcountrygenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Genus_bycountryoforigin_corr_gender_age_feces_LRTallsamples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(sibocorrcountrygenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) +
  ggtitle("Differences by country of origin (CAR vs. Madagascar), controlled gender, age, LRT model")
dev.off()

#' line
# make graph plot
sibocorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(sibocorrcountrygenderage_reduced_prevabdeseq$Genus, "/", sibocorrcountrygenderage_reduced_prevabdeseq$Family)
pdf("differentially_present_genera by country of origin_fecal samples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 8, # define plot width and height. completely up to user.
  height =12
)
differential_plot(sibocorrcountrygenderage_reduced_prevabdeseq, "Significantly different  genera by country of origin in fecal samples" )
dev.off()

#- 3335

#### make the Deseq analysis on all feces samples and stunted_stunted_vs_non.stunted on Genus level all feces samples, no filter for foldchange ####
phyloseq::sample_data(dffilteredfeces)$haz<-as.factor(phyloseq::sample_data(dffilteredfeces)$haz)
phyloseq::sample_data(dffilteredfeces)$haz<-factor(phyloseq::sample_data(dffilteredfeces)$haz, labels =c("non-stunted", "moderately stunted", "severely stunted"))
levels(phyloseq::sample_data(dffilteredfeces)$haz)

phyloseq::sample_data(dffilteredfeces)$stunted<-""
phyloseq::sample_data(dffilteredfeces)$stunted[which(phyloseq::sample_data(dffilteredfeces)$haz == "moderately stunted")] <- "stunted"
phyloseq::sample_data(dffilteredfeces)$stunted[which(phyloseq::sample_data(dffilteredfeces)$haz == "severely stunted")] <- "stunted"
phyloseq::sample_data(dffilteredfeces)$stunted[which(phyloseq::sample_data(dffilteredfeces)$haz == "non-stunted")] <- "non-stunted"
phyloseq::sample_data(dffilteredfeces)$stunted<-as.factor(phyloseq::sample_data(dffilteredfeces)$stunted)
levels(phyloseq::sample_data(dffilteredfeces)$stunted)

 level = "Genus"
 variable_of_interest = "stunted"
 name = "feces"

res = dffilteredfeces %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run","pays",variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             , fold_change = 0
             )
merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Genus_bystunted_stunted_vs_non.stunted_corr_gender_age_feces_LRTnofilterfoldchangeallsamples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by stunting status (stunted_stunted_vs_non.stunted vs. non-stunted_stunted_vs_non.stunted), controlled gender, age, LRT model")
dev.off()

# make graph plot
corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq$taxonomy<- paste0(corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq$Family, "/", corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_genera by stunting status-merged dataset.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 12, # define plot width and height. completely up to user.
  height = 8)
ggplot(data=corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
 geom_bar(position="dodge",stat="identity", color="black") +
 coord_flip() +
 theme(axis.text.x = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 xlab("")+
 theme(axis.title.y = element_text(size=14))+
 theme(title = element_text(size=16, face="bold"))+
 ylab("log2 fold change")+
 ggtitle("Significantly different genera by stunting status, both countries merged")
dev.off()


#- 4663
#### make the Deseq analysis on feces from Bangui samples and stunted_stunted_vs_non.stunted on Genus level feces from Bangui samples, no filter for foldchange ####
phyloseq::sample_data(dffilteredfecesBangui)$haz<-as.factor(phyloseq::sample_data(dffilteredfecesBangui)$haz)
phyloseq::sample_data(dffilteredfecesBangui)$haz<-factor(phyloseq::sample_data(dffilteredfecesBangui)$haz, labels =c("non-stunted", "moderately stunted", "severely stunted"))
levels(phyloseq::sample_data(dffilteredfecesBangui)$haz)

phyloseq::sample_data(dffilteredfecesBangui)$stunted<-""
phyloseq::sample_data(dffilteredfecesBangui)$stunted[which(phyloseq::sample_data(dffilteredfecesBangui)$haz == "moderately stunted")] <- "stunted"
phyloseq::sample_data(dffilteredfecesBangui)$stunted[which(phyloseq::sample_data(dffilteredfecesBangui)$haz == "severely stunted")] <- "stunted"
phyloseq::sample_data(dffilteredfecesBangui)$stunted[which(phyloseq::sample_data(dffilteredfecesBangui)$haz == "non-stunted")] <- "non-stunted"
phyloseq::sample_data(dffilteredfecesBangui)$stunted<-as.factor(phyloseq::sample_data(dffilteredfecesBangui)$stunted)
levels(phyloseq::sample_data(dffilteredfecesBangui)$stunted)

level = "Genus"
 variable_of_interest = "stunted"
 name = "feces_Mada"

res = dffilteredfecesBangui %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             , fold_change = 0
             )
merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq <- payscorr
#- 3553
#make the actual plot

pdf("differentially_present_taxa_Genus_bystunted_stunted_vs_non.stunted_corr_gender_age_feces_LRTnofilterfoldchangeallsamplesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by stunting status (stunted_stunted_vs_non.stunted vs. non-stunted_stunted_vs_non.stunted), controlled gender, age, LRT model")
dev.off()

# make graph plot
corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq$taxonomy<- paste0(corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq$Family, "/", corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_genera by stunting status CAR.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 12, # define plot width and height. completely up to user.
  height = 8)
ggplot(data=corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
 geom_bar(position="dodge",stat="identity", color="black") +
 coord_flip() +
 theme(axis.text.x = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 xlab("")+
 theme(axis.title.y = element_text(size=14))+
 theme(title = element_text(size=16, face="bold"))+
 ylab("log2 fold change")+
 ggtitle("Significantly different genera by stunting status, CAR")
dev.off()



#### make the Deseq analysis on feces from Mada samples and stunted_stunted_vs_non.stunted on Genus level feces from Mada samples, no filter for foldchange ####
phyloseq::sample_data(dffilteredfecesMada)$haz<-as.factor(phyloseq::sample_data(dffilteredfecesMada)$haz)
phyloseq::sample_data(dffilteredfecesMada)$haz<-factor(phyloseq::sample_data(dffilteredfecesMada)$haz, labels =c("non-stunted", "moderately stunted", "severely stunted"))
levels(phyloseq::sample_data(dffilteredfecesMada)$haz)

phyloseq::sample_data(dffilteredfecesMada)$stunted<-""
phyloseq::sample_data(dffilteredfecesMada)$stunted[which(phyloseq::sample_data(dffilteredfecesMada)$haz == "moderately stunted")] <- "stunted"
phyloseq::sample_data(dffilteredfecesMada)$stunted[which(phyloseq::sample_data(dffilteredfecesMada)$haz == "severely stunted")] <- "stunted"
phyloseq::sample_data(dffilteredfecesMada)$stunted[which(phyloseq::sample_data(dffilteredfecesMada)$haz == "non-stunted")] <- "non-stunted"
phyloseq::sample_data(dffilteredfecesMada)$stunted<-as.factor(phyloseq::sample_data(dffilteredfecesMada)$stunted)
levels(phyloseq::sample_data(dffilteredfecesMada)$stunted)

variable_of_interest = "stunted"
name = "feces_Mada"
level = "Genus"

res = dffilteredfecesMada %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             , fold_change = 0
             )
merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Genus_bystunted_stunted_vs_non.stunted_corr_gender_age_feces_LRTnofilterfoldchangeallsamplesMada.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by stunting status (stunted_stunted_vs_non.stunted vs. non-stunted_stunted_vs_non.stunted), controlled gender, age, LRT model")
dev.off()

# make graph plot
corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq$taxonomy<- paste0(corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq$Family, "/", corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_genera by stunting status Madagascar.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 12, # define plot width and height. completely up to user.
  height = 8)
ggplot(data=corrstunted_stunted_vs_non.stuntedgenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
 geom_bar(position="dodge",stat="identity", color="black") +
 coord_flip() +
 theme(axis.text.x = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 xlab("")+
 theme(axis.title.y = element_text(size=14))+
 theme(title = element_text(size=16, face="bold"))+
 ylab("log2 fold change")+
 ggtitle("Significantly different genera by stunting status, Madagascar")
dev.off()



#### make the Deseq analysis on all feces samples and alphaantitrypsinlevel on Genus level all feces samples, no filter for foldchange ####
 variable_of_interest = "alphaantitrypsinlevel"
name = "feces"
level = "Genus"

res = dffilteredfeces %>%
  phyloseq::subset_samples( alphaantitrypsinlevel != "") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run","stunted","pays", variable_of_interest )
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )

merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corralphaantitrypsinlevelgenderage_reduced_prevabdeseq <- payscorr


#make the actual plot

pdf("differentially_present_taxa_Genus_byalphaantitrypsinlevel_corr_gender_age_feces_LRTnofilterfoldchangeallsamples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corralphaantitrypsinlevelgenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by alphaantitrypsin level (high alphaantitrypsinlevel vs. low alphaantitrypsinlevel), controlled gender, age, LRT model")
dev.off()


# make graph plot
corralphaantitrypsinlevelgenderage_reduced_prevabdeseq$taxonomy<- paste0(corralphaantitrypsinlevelgenderage_reduced_prevabdeseq$Family, "/", corralphaantitrypsinlevelgenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_genera by AAT status-merged dataset.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 12, # define plot width and height. completely up to user.
  height = 8)
ggplot(data=corralphaantitrypsinlevelgenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
 geom_bar(position="dodge",stat="identity", color="black") +
 coord_flip() +
 theme(axis.text.x = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 xlab("")+
 theme(axis.title.y = element_text(size=14))+
 theme(title = element_text(size=16, face="bold"))+
 ylab("log2 fold change")+
 ggtitle("Significantly different genera by fecal alpha-antitrypsin levels, merged dataset")
dev.off()


#### make the Deseq analysis on feces from Bangui samples and alphaantitrypsinlevel on Genus level feces from Bangui samples, no filter for foldchange ####
dffilteredfecesBanguiG = phyloseq::tax_glom(dffilteredfecesBangui, "Genus")
dffilteredfecesBanguiG <- dffilteredfecesBanguiG %>%
 phyloseq::subset_samples(alphaantitrypsinlevel == "normal" | alphaantitrypsinlevel=="elevated")
phyloseq::sample_data(dffilteredfecesBanguiG)$alphaantitrypsinlevel<-as.character(phyloseq::sample_data(dffilteredfecesBanguiG)$alphaantitrypsinlevel)
phyloseq::sample_data(dffilteredfecesBanguiG)$alphaantitrypsinlevel[phyloseq::sample_data(dffilteredfecesBanguiG)$alphaantitrypsinlevel=="normal"]<-"normal"
phyloseq::sample_data(dffilteredfecesBanguiG)$alphaantitrypsinlevel[phyloseq::sample_data(dffilteredfecesBanguiG)$alphaantitrypsinlevel=="elevated"]<-"too_high"
phyloseq::sample_data(dffilteredfecesBanguiG)$alphaantitrypsinlevel<-as.factor(phyloseq::sample_data(dffilteredfecesBanguiG)$alphaantitrypsinlevel)


variable_of_interest = "alphaantitrypsinlevel"
name = "feces_Bangui"
 level = "Genus"

res = dffilteredfecesBangui %>%
  phyloseq::subset_samples( alphaantitrypsinlevel != "") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corralphaantitrypsinlevelgenderage_reduced_prevabdeseq <- payscorr
#- 3734
#make the actual plot

pdf("differentially_present_taxa_Genus_byalphaantitrypsinlevel_corr_gender_age_feces_LRTnofilterfoldchangeallsamplesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corralphaantitrypsinlevelgenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by alphaantitrypsin level (high alphaantitrypsinlevel vs. low alphaantitrypsinlevel), controlled gender, age, LRT model")
dev.off()

# make graph plot
corralphaantitrypsinlevelgenderage_reduced_prevabdeseq$taxonomy<- paste0(corralphaantitrypsinlevelgenderage_reduced_prevabdeseq$Family, "/", corralphaantitrypsinlevelgenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_genera by AAT status-CAR.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 6, # define plot width and height. completely up to user.
  height = 1)
ggplot(data=corralphaantitrypsinlevelgenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
 geom_bar(position="dodge",stat="identity", color="black") +
 coord_flip() +
 theme(axis.text.x = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 xlab("")+
 theme(axis.title.y = element_text(size=14))+
 theme(title = element_text(size=16, face="bold"))+
 ylab("log2 fold change")+
 ggtitle("Significantly different genera by fecal alpha-antitrypsin levels, CAR")
dev.off()

#- 4907
#### make the Deseq analysis on feces from Mada samples and alphaantitrypsinlevel on Genus level feces from Mada samples, no filter for foldchange ####
variable_of_interest = "alphaantitrypsinlevel"
name = "feces_Mada"
 level = "Genus"

res = dffilteredfecesMada %>%
  phyloseq::subset_samples( alphaantitrypsinlevel != "") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corralphaantitrypsinlevelgenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Genus_byalphaantitrypsinlevel_corr_gender_age_feces_LRTnofilterfoldchangeallsamplesMada.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corralphaantitrypsinlevelgenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by alphaantitrypsin level (high alphaantitrypsinlevel vs. low alphaantitrypsinlevel), controlled gender, age, LRT model")
dev.off()

# make graph plot
corralphaantitrypsinlevelgenderage_reduced_prevabdeseq$taxonomy<- paste0(corralphaantitrypsinlevelgenderage_reduced_prevabdeseq$Family, "/", corralphaantitrypsinlevelgenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_genera by AAT status-Mada.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 12, # define plot width and height. completely up to user.
  height = 8)
ggplot(data=corralphaantitrypsinlevelgenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
 geom_bar(position="dodge",stat="identity", color="black") +
 coord_flip() +
 theme(axis.text.x = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 xlab("")+
 theme(axis.title.y = element_text(size=14))+
 theme(title = element_text(size=16, face="bold"))+
 ylab("log2 fold change")+
 ggtitle("Significantly different genera by fecal alpha-antitrypsin levels, Madagascar")
dev.off()


#### make the Deseq analysis on all feces samples and calprotectinelevel on Genus level all feces samples, no filter for foldchange ####
variable_of_interest = "calprotectinelevel"
name = "feces"
 level = "Genus"

res = dffilteredfeces %>%
  phyloseq::subset_samples( calprotectinelevel != "") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run","pays","stunted", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrcalprotectinelevelgenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Genus_bycalprotectinelevel_corr_gender_age_feces_LRTnofilterfoldchangeallsamples.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrcalprotectinelevelgenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by alphaantitrypsin level (high calprotectinelevel vs. low calprotectinelevel), controlled gender, age, LRT model")
dev.off()

# make graph plot
corrcalprotectinelevelgenderage_reduced_prevabdeseq$taxonomy<- paste0(corrcalprotectinelevelgenderage_reduced_prevabdeseq$Family, "/", corrcalprotectinelevelgenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_genera by calpro status-merged dataset.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 12, # define plot width and height. completely up to user.
  height = 8)
ggplot(data=corrcalprotectinelevelgenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
 geom_bar(position="dodge",stat="identity", color="black") +
 coord_flip() +
 theme(axis.text.x = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 xlab("")+
 theme(axis.title.y = element_text(size=14))+
 theme(title = element_text(size=16, face="bold"))+
 ylab("log2 fold change")+
 ggtitle("Significantly different genera by fecal calprotectine levels, merged dataset")
dev.off()

#- 5011
#### make the Deseq analysis on feces from Bangui samples and calprotectinelevel on Genus level feces from Bangui samples, no filter for foldchange ####
variable_of_interest = "calprotectinelevel"
name = "feces_Bangui"
 level = "Genus"

res = dffilteredfecesBangui %>%
  phyloseq::subset_samples( calprotectinelevel != "") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrcalprotectinelevelgenderage_reduced_prevabdeseq <- payscorr
#make the actual plot

pdf("differentially_present_taxa_Genus_bycalprotectinelevel_corr_gender_age_feces_LRTnofilterfoldchangeallsamplesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrcalprotectinelevelgenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by alphaantitrypsin level (high calprotectinelevel vs. low calprotectinelevel), controlled gender, age, LRT model")
dev.off()

# make graph plot
corrcalprotectinelevelgenderage_reduced_prevabdeseq$taxonomy<- paste0(corrcalprotectinelevelgenderage_reduced_prevabdeseq$Family, "/", corrcalprotectinelevelgenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_genera by calpro status-CAR.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 12, # define plot width and height. completely up to user.
  height = 8)
ggplot(data=corrcalprotectinelevelgenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
 geom_bar(position="dodge",stat="identity", color="black") +
 coord_flip() +
 theme(axis.text.x = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 xlab("")+
 theme(axis.title.y = element_text(size=14))+
 theme(title = element_text(size=16, face="bold"))+
 ylab("log2 fold change")+
 ggtitle("Significantly different genera by fecal calprotectine levels, CAR")
dev.off()


#### make the Deseq analysis on feces from Mada samples and calprotectinelevel on Genus level feces from Mada samples, no filter for foldchange ####
#- 3920

variable_of_interest = "calprotectinelevel"
name = "feces_Mada"
 level = "Genus"

res = dffilteredfecesMada %>%
  phyloseq::subset_samples( calprotectinelevel != "") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrcalprotectinelevelgenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Genus_bycalprotectinelevel_corr_gender_age_feces_LRTnofilterfoldchangeallsamplesMada.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 10, # define plot width and height. completely up to user.
  height = 8)
ggplot(corrcalprotectinelevelgenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
 scale_fill_manual(values = c("red", "blue", "green")) +
 theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by alphaantitrypsin level (high calprotectinelevel vs. low calprotectinelevel), controlled gender, age, LRT model")
dev.off()

# make graph plot
corrcalprotectinelevelgenderage_reduced_prevabdeseq$taxonomy<- paste0(corrcalprotectinelevelgenderage_reduced_prevabdeseq$Family, "/", corrcalprotectinelevelgenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_genera by calpro status-Mada.pdf", #name of file to print. can also include relative or absolute path before filename.
  width = 12, # define plot width and height. completely up to user.
  height = 8)
ggplot(data=corrcalprotectinelevelgenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
 geom_bar(position="dodge",stat="identity", color="black") +
 coord_flip() +
 theme(axis.text.x = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 theme(axis.text.y = element_text(size=12))+
 xlab("")+
 theme(axis.title.y = element_text(size=14))+
 theme(title = element_text(size=16, face="bold"))+
 ylab("log2 fold change")+
 ggtitle("Significantly different genera by fecal calprotectine levels, Madagascar")
dev.off()

#- 3972 (5116)
#### make the Deseq analysis on feces from Mada samples and calprotectinelevel on Genus level feces from Mada samples including current breastfeeding, no filter for foldchange ####
#- 3920

variable_of_interest = "calprotectinelevel"
name = "feces_Mada"
level = "Genus"

res = dffilteredfecesMada %>%
  phyloseq::subset_samples( calprotectinelevel != "") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
               , full_model = c("age","sexe","run", "current_breastfeeding", variable_of_interest)
               , variable_of_interest = variable_of_interest
               , sig_level = 0.05
  )
merged   = prevalence_abundance( dataset = res$glom
                                 , level = level , name = name
                                 , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrcalprotectinelevelgenderage_reduced_prevabdeseq <- payscorr

#make the actual plot

pdf("differentially_present_taxa_Genus_bycalprotectinelevel_corr_gender_age_feces_LRTnofilterfoldchangeallsamplesMadaincbf.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 10, # define plot width and height. completely up to user.
    height = 8)
ggplot(corrcalprotectinelevelgenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
  scale_fill_manual(values = c("red", "blue", "green")) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by alphaantitrypsin level (high calprotectinelevel vs. low calprotectinelevel), controlled gender, age, LRT model")
dev.off()

# make graph plot
corrcalprotectinelevelgenderage_reduced_prevabdeseq$taxonomy<- paste0(corrcalprotectinelevelgenderage_reduced_prevabdeseq$Family, "/", corrcalprotectinelevelgenderage_reduced_prevabdeseq$Genus)
pdf("differentially_present_genera by calpro status-Madaincbf.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 12, # define plot width and height. completely up to user.
    height = 8)
ggplot(data=corrcalprotectinelevelgenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
  geom_bar(position="dodge",stat="identity", color="black") +
  coord_flip() +
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  xlab("")+
  theme(axis.title.y = element_text(size=14))+
  theme(title = element_text(size=16, face="bold"))+
  ylab("log2 fold change")+
  ggtitle("Significantly different genera by fecal calprotectine levels including current breastfeeding as covariable, Madagascar")
dev.off()

#- 3972 (5116)
#### make the Deseq analysis on all feces samples and sibo on Genus level all feces samples, no filter for foldchange ####

variable_of_interest = "sibo"
name = "feces"
 level = "Genus"

res = dffilteredfeces %>%
  phyloseq::subset_samples( sibo != "") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrsibogenderage_reduced_prevabdeseq <- payscorr

      #make the actual plot

      pdf("differentially_present_taxa_Genus_bysibo_corr_gender_age_feces_LRTnofilterfoldchangeallsamples.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
      ggplot(corrsibogenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
       scale_fill_manual(values = c("red", "blue", "green")) +
       theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by sibo status (yes vs. no) in feces samples, controlled gender, age, stunted LRT model")
      dev.off()

      corrsibogenderage_reduced_prevabdeseq$taxonomy<- paste0(corrsibogenderage_reduced_prevabdeseq$Family, "/", corrsibogenderage_reduced_prevabdeseq$Genus)


      pdf("differentially_present_genera in feces by sibo status-merged.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 12, # define plot width and height. completely up to user.
        height = 8)
      ggplot(data=corrsibogenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
       geom_bar(position="dodge",stat="identity", color="black") +
       coord_flip() +
       theme(axis.text.x = element_text(size=12))+
       theme(axis.text.y = element_text(size=12))+
       theme(axis.text.y = element_text(size=12))+
       xlab("")+
       theme(axis.title.y = element_text(size=14))+
       theme(title = element_text(size=16, face="bold"))+
       ylab("log2 fold change")+
       ggtitle("Significantly different genera in feces by SIBO status, merged dataset")
      dev.off()



#### make the Deseq analysis on feces from Bangui samples and sibo on Genus level feces from Bangui samples, no filter for foldchange ####
variable_of_interest = "sibo"
name = "feces_Bangui"
 level = "Genus"

res = dffilteredfeces %>%
  phyloseq::subset_samples( sibo != "") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrsibogenderage_reduced_prevabdeseq <- payscorr

      #make the actual plot

      pdf("differentially_present_taxa_Genus_bysibo_corr_gender_age_feces_LRTnofilterfoldchangeallsamplesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
      ggplot(corrsibogenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
       scale_fill_manual(values = c("red", "blue", "green")) +
       theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by sibo status (yes vs. no) in feces samples, controlled gender, age, stunted LRT model")
      dev.off()

      # make graph plot
      corrsibogenderage_reduced_prevabdeseq$taxonomy<- paste0(corrsibogenderage_reduced_prevabdeseq$Family, "/", corrsibogenderage_reduced_prevabdeseq$Genus)
      pdf("differentially_present_genera in feces by sibo status-CAR.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 12, # define plot width and height. completely up to user.
        height = 8)
      ggplot(data=corrsibogenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
       geom_bar(position="dodge",stat="identity", color="black") +
       coord_flip() +
       theme(axis.text.x = element_text(size=12))+
       theme(axis.text.y = element_text(size=12))+
       theme(axis.text.y = element_text(size=12))+
       xlab("")+
       theme(axis.title.y = element_text(size=14))+
       theme(title = element_text(size=16, face="bold"))+
       ylab("log2 fold change")+
       ggtitle("Significantly different genera in feces by SIBO status, CAR")
      dev.off()

#### make the Deseq analysis on feces from Mada samples and sibo on Genus level feces from Mada samples, no filter for foldchange ####

variable_of_interest = "sibo"
name = "feces_Mada"
 level = "Genus"

res = dffilteredfecesMada %>%
  phyloseq::subset_samples( sibo != "") %>%
  filter_prune( ) %>%
  deseq_cooks( level = level
             , full_model = c("age","sexe","run", variable_of_interest)
             , variable_of_interest = variable_of_interest
             , sig_level = 0.05
             )
merged   = prevalence_abundance( dataset = res$glom
                               , level = level , name = name
                               , group = variable_of_interest, write_csv = TRUE )
payscorr = deseq_names( res, merged, name = name
                        , fold_change = 0
                        , sig_level = 0.05
                        , write_csv = TRUE ) %>% max_reorder_results()

corrsibogenderage_reduced_prevabdeseq <- payscorr

      #make the actual plot

      pdf("differentially_present_taxa_Genus_bysibo_corr_gender_age_feces_LRTnofilterfoldchangeallsamplesMada.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
      ggplot(corrsibogenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=2) +
       scale_fill_manual(values = c("red", "blue", "green")) +
       theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8, margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Differences by sibo status (yes vs. no) in feces samples, controlled gender, age, stunted LRT model")
      dev.off()

      corrsibogenderage_reduced_prevabdeseq$taxonomy<- paste0(corrsibogenderage_reduced_prevabdeseq$Family, "/", corrsibogenderage_reduced_prevabdeseq$Genus)


      pdf("differentially_present_genera in feces by sibo status-Mada.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 12, # define plot width and height. completely up to user.
        height = 8)
      ggplot(data=corrsibogenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
       geom_bar(position="dodge",stat="identity", color="black") +
       coord_flip() +
       theme(axis.text.x = element_text(size=12))+
       theme(axis.text.y = element_text(size=12))+
       theme(axis.text.y = element_text(size=12))+
       xlab("")+
       theme(axis.title.y = element_text(size=14))+
       theme(title = element_text(size=16, face="bold"))+
       ylab("log2 fold change")+
       ggtitle("Significantly different genera in feces by SIBO status, Mada")
      dev.off()

#- 5276
#### calprotectine and AAT by stunting and sibo in small intestinal aspirates ####
      sampledata<-data.frame(phyloseq::sample_data(dffilteredduodenal))
      sampledata2 = dplyr::filter(sampledata, sibo!="")
      ggplot(sampledata2, aes(x=pays, y=calprotectinesurliquideduodn, fill=sibo)) +
       geom_boxplot()

      siboyes = dplyr::filter(sampledata2, sampledata2$sibo=="yes")
      sibono = dplyr::filter(sampledata2, sampledata2$sibo=="no")

      wilcox.test(siboyes$calprotectinesurliquideduodn, sibono$calprotectinesurliquideduodn) # p-value = 0.07436

      ggplot(sampledata2, aes(x=pays, y=aatsurldmgml, fill=sibo)) +
       geom_boxplot()

      siboyes2 = dplyr::filter(sampledata2, sampledata2$sibo=="yes")
      sibono2 = dplyr::filter(sampledata2, sampledata2$sibo=="no")

      wilcox.test(siboyes2$aatsurldmgml, sibono2$aatsurldmgml) # p-value = ns

#- 5296
#### Is there correlation of duodenal and fecal levels of calprotectine and AAT as well as cfu/ml and inflammation? ####
      dffilteredduodenal_inflammation = phyloseq::subset_samples(dffilteredduodenal, !is.na(calprotectinesurliquideduodn))
      dffilteredduodenal_inflammation = phyloseq::subset_samples(dffilteredduodenal_inflammation, !is.na(aatsurldmgml))

      sampledata_inflammation<-as.data.frame(phyloseq::sample_data(dffilteredduodenal_inflammation))
      sampledata_inflammation$logAATfeces<-log10(sampledata_inflammation$aatmggdeps)
      sampledata_inflammation$logAATduodenal<-log10(sampledata_inflammation$aatsurldmgml)
      
      sampledata_inflammation<-filter(sampledata_inflammation, current_breastfeeding!="yes")

      cor.test(sampledata_inflammation$logAATfeces, sampledata_inflammation$logAATduodenal, method = c("spearman"), use = "complete.obs")
      cor.test(sampledata_inflammation$aatsurldmgml, sampledata_inflammation$aatmggdeps, method = c("spearman"), use = "complete.obs")

      ggpubr::ggscatter(sampledata_inflammation, x = "aatsurldmgml", y = "aatmggdeps",
           add = "reg.line", conf.int = TRUE,
           cor.coef = TRUE, cor.method = "spearman",
           xlab = "AAT level in duodenum (mg/ml)", ylab = "AAT levels in feces (mg/g de poids sec)")

      sampledata_inflammation$logcalprofeces<-log10(sampledata_inflammation$calprotectineggdeps)
      sampledata_inflammation$logcalproduodenal<-log10(sampledata_inflammation$calprotectinesurliquideduodn)

      cor.test(sampledata_inflammation$logcalprofeces, sampledata_inflammation$logcalproduodenal, method = c("spearman"), use = "complete.obs")
      cor.test(sampledata_inflammation$calprotectineggdeps, sampledata_inflammation$calprotectinesurliquideduodn, method = c("spearman"), use = "complete.obs")


      ggpubr::ggscatter(sampledata_inflammation, x = "calprotectinesurliquideduodn", y = "calprotectineggdeps",
           add = "reg.line", conf.int = TRUE,
           cor.coef = TRUE, cor.method = "spearman",
           xlab = "calprotectine levels in duodenum (mg/ml)", ylab = "Calprotecine levels in feces (mg/g de poids sec)")

      cor.test(sampledata_inflammation$calprotectineggdeps, sampledata_inflammation$ufcml, method = c("spearman"), use = "complete.obs")


      cor.test(sampledata_inflammation$aatmggdeps, sampledata_inflammation$ufcml, method = c("spearman"), use = "complete.obs")


      cor.test(sampledata_inflammation$calprotectinesurliquideduodn, sampledata_inflammation$ufcml, method = c("spearman"), use = "complete.obs")


      ggpubr::ggscatter(sampledata_inflammation, x = "calprotectinesurliquideduodn", y = "ufcml",
           add = "reg.line", conf.int = TRUE,
           cor.coef = TRUE, cor.method = "spearman",
           xlab = "calprotectine levels in duodenum (mg/ml)", ylab = "cfu per ml in duodenum")

      cor.test(sampledata_inflammation$aatsurldmgml, sampledata_inflammation$ufcml, method = c("spearman"), use = "complete.obs")


      ggpubr::ggscatter(sampledata_inflammation, x = "aatsurldmgml", y = "ufcml",
           add = "reg.line", conf.int = TRUE,
           cor.coef = TRUE, cor.method = "spearman",
           xlab = "AAT levels in duodenum (mg/ml)", ylab = "cfu per ml in duodenum")


##### Analysis of cytokines in duodenal samples #####
### Make heatmaps of co-occuring cytokines in duodenum  ####
      # first import the full cytokine dataset
      datacytokine <- read.csv("Cytokinedatapaper.csv", header=TRUE, sep=";")

      # Define data sets to cross-correlate with each other
      Dataduodenal_wilcox = (datacytokine[, c(10:39)])
      corrz = microbiome::associate( Dataduodenal_wilcox , method = "spearman", p.adj.threshold = 0.05, n.signif = 1)

      pdf("heatmapcytokinesduodenum.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10)
      p = microbiome::heat(corrz, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)
      print(p)
      dev.off()


### Make heatmaps of co-occuring cytokines in duodenum CAR ####

      # Define data sets to cross-correlate with each other
      dffilteredduodenal_inflammation_RCA = dplyr::filter(datacytokine, pays=="RCA")
      Dataduodenal_wilcox <- (phyloseq::sample_data(dffilteredduodenal_inflammation_RCA)[, c(10:39)]) ## choose the columns you want to include
      Dataduodenal_wilcox<-as.matrix(Dataduodenal_wilcox)

      x<- Dataduodenal_wilcox
      y<- Dataduodenal_wilcox

      correlations = microbiome::associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
      correlation.table = microbiome::associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
      pdf("heatmapcytokinesduodenumCAR.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10)
      p = microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)
      print(p)
      dev.off()



### Make heatmaps of co-occuring cytokines in duodenum Mada ####

      # Define data sets to cross-correlate with each other
      dffilteredduodenal_inflammation_Mada = dplyr::filter(datacytokine, pays=="Madagascar")
      Dataduodenal_wilcox <- (phyloseq::sample_data(dffilteredduodenal_inflammation_Mada)[, c(10:39)]) ## choose the columns you want to include
      Dataduodenal_wilcox<-as.matrix(Dataduodenal_wilcox)

      x<- Dataduodenal_wilcox
      y<- Dataduodenal_wilcox

      correlations = microbiome::associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
      correlation.table = microbiome::associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
      pdf("heatmapcytokinesduodenumMada.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10)
      p = microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)
      print(p)
      dev.off()


### Make heatmaps of co-occuring cytokines in duodenum by SIBO status ####

      # Define data sets to cross-correlate with each other
      dffilteredduodenal_inflammation_SIBO = dplyr::filter(datacytokine, sibo=="yes")
      Dataduodenal_wilcox <- (phyloseq::sample_data(dffilteredduodenal_inflammation_SIBO)[, c(10:39)]) ## choose the columns you want to include
      Dataduodenal_wilcox<-as.matrix(Dataduodenal_wilcox)

      x<- Dataduodenal_wilcox
      y<- Dataduodenal_wilcox

      correlations = microbiome::associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
      correlation.table = microbiome::associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
      pdf("heatmapcytokinesduodenumSIBO.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10)
      p = microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)
      print(p)
      dev.off()

      dffilteredduodenal_inflammation_noSIBO = dplyr::filter( datacytokine , sibo=="no")
      Dataduodenal_wilcox <- (phyloseq::sample_data(dffilteredduodenal_inflammation_noSIBO)[, c(9:38)]) ## choose the columns you want to include
      Dataduodenal_wilcox<-as.matrix(Dataduodenal_wilcox)

      x<- Dataduodenal_wilcox
      y<- Dataduodenal_wilcox

      correlations = microbiome::associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
      correlation.table = microbiome::associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
      pdf("heatmapcytokinesduodenumnoSIBO.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10)
      p = microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)
      print(p)
      dev.off()

### Make heatmaps of co-occuring cytokines in duodenum with calprotectin and alphaantitrypsin in feces ####

      # Define data sets to cross-correlate with each other
      Dataduodenal_wilcox <- datacytokine[, c(5:6, 10:39)] ## choose the columns you want to include
      Dataduodenal_wilcox<-as.matrix(Dataduodenal_wilcox)

      x<- Dataduodenal_wilcox
      y<- Dataduodenal_wilcox

      correlations = microbiome::associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
      correlation.table = microbiome::associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)

      pdf("heatmapcytokinesduodenumcomparedtoAATCalproinfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10)
      p = microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)
      print(p)
      dev.off()



### Make heatmaps of co-occuring cytokines in duodenum with calprotectin and alphaantitrypsin in feces RCA ####

      # Define data sets to cross-correlate with each other
      Dataduodenal_wilcox <- dffilteredduodenal_inflammation_RCA[, c(5:6, 10:39)] ## choose the columns you want to include
      Dataduodenal_wilcox<-as.matrix(Dataduodenal_wilcox)

      x<- Dataduodenal_wilcox
      y<- Dataduodenal_wilcox

      correlations = microbiome::associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
      correlation.table = microbiome::associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)

      pdf("heatmapcytokinesduodenumcomparedtoAATCalproinfecesRCA.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10)
      p = microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)
      print(p)
      dev.off()



### Make heatmaps of co-occuring cytokines in duodenum with calprotectin and alphaantitrypsin in feces Mada ####

      # Define data sets to cross-correlate with each other
      Dataduodenal_wilcox <- dffilteredduodenal_inflammation_Mada[, c(5:6, 10:39)] ## choose the columns you want to include
      Dataduodenal_wilcox<-as.matrix(Dataduodenal_wilcox)

      x<- Dataduodenal_wilcox
      y<- Dataduodenal_wilcox

      correlations = microbiome::associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
      correlation.table = microbiome::associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
      pdf("heatmapcytokinesduodenumcomparedtoAATCalproinfecesMada.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10)
      p = microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)
      print(p)
      dev.off()


#- 5504
### Make heatmaps of co-occuring cytokines in duodenum with calprotectin and alphaantitrypsin in duodenum ####

      # Define data sets to cross-correlate with each other
      Dataduodenal_wilcox <- datacytokine[, c(7:8, 10:39)] ## choose the columns you want to include
      dim(Dataduodenal_wilcox)
      Dataduodenal_wilcox2<-Dataduodenal_wilcox[complete.cases(Dataduodenal_wilcox)==TRUE, ]
      Dataduodenal_wilcox2[, 2]<-Dataduodenal_wilcox2[, 2]*1000  # as it cannot use values below 10
      Dataduodenal_wilcox2<-Dataduodenal_wilcox2[complete.cases(Dataduodenal_wilcox2)==TRUE, ]
      Dataduodenal_wilcox2<-as.matrix(Dataduodenal_wilcox2)

      x<- Dataduodenal_wilcox2
      y<- Dataduodenal_wilcox2

      correlations = microbiome::associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
      correlation.table = microbiome::associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)

      pdf("heatmapcytokinesduodenumcomparedtoAATCalproinduodenum.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10)
      p = microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)
      print(p)
      dev.off()

      colnames(Dataduodenal_wilcox2)
      dim(Dataduodenal_wilcox2)

### Make heatmaps of co-occuring cytokines in duodenum with calprotectin and alphaantitrypsin in duodenum RCA ####

      # Define data sets to cross-correlate with each other
      Dataduodenal_wilcox <- dffilteredduodenal_inflammation_RCA[, c(7:8, 10:39)] ## choose the columns you want to include
      dim(Dataduodenal_wilcox)
      Dataduodenal_wilcox2<-Dataduodenal_wilcox[complete.cases(Dataduodenal_wilcox)==TRUE, ]
      Dataduodenal_wilcox2[, 2]<-Dataduodenal_wilcox2[, 2]*1000  # as it cannot use values below 10
      Dataduodenal_wilcox2<-Dataduodenal_wilcox2[complete.cases(Dataduodenal_wilcox2)==TRUE, ]
      Dataduodenal_wilcox2<-as.matrix(Dataduodenal_wilcox2)

      x<- Dataduodenal_wilcox2
      y<- Dataduodenal_wilcox2

      correlations = microbiome::associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
      correlation.table = microbiome::associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)

      pdf("heatmapcytokinesduodenumcomparedtoAATCalproinduodenumRCA.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10)
      p = microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)
      print(p)
      dev.off()


colnames(tax)
### Make heatmaps of co-occuring cytokines in duodenum with calprotectin and alphaantitrypsin in duodenum Mada ####

      # Define data sets to cross-correlate with each other
      Dataduodenal_wilcox <- dffilteredduodenal_inflammation_RCA[, c(7:8, 10:39)] ## choose the columns you want to include
      dim(Dataduodenal_wilcox)
      Dataduodenal_wilcox2<-Dataduodenal_wilcox[complete.cases(Dataduodenal_wilcox)==TRUE, ]
      Dataduodenal_wilcox2[, 2]<-Dataduodenal_wilcox2[, 2]*1000  # as it cannot use values below 10
      Dataduodenal_wilcox2<-Dataduodenal_wilcox2[complete.cases(Dataduodenal_wilcox2)==TRUE, ]
      Dataduodenal_wilcox2<-as.matrix(Dataduodenal_wilcox2)

      x<- Dataduodenal_wilcox2
      y<- Dataduodenal_wilcox2

      correlations = microbiome::associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
      correlation.table = microbiome::associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)

      pdf("heatmapcytokinesduodenumcomparedtoAATCalproinduodenumMada.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10)
      p = microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)
      print(p)
      dev.off()


### PCA analysis on cytokines ####
      dim(datacytokine)
      data<-datacytokine[, c(10:39, 9, 2)]
      data$plate<-as.factor(data$plate)
      data<-data[complete.cases(data)==TRUE, ]
      data_log<-log2(data[, c(1:30)])

      res.pca <- prcomp(data_log, scale = FALSE, center = FALSE)
      get_eig(res.pca)

      pdf("PlotindividusPcoACytokinesellipse.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      factoextra::fviz_pca_ind(res.pca,
                   geom.ind = "point",
                   col.ind = data[, 31], # colorer by groups
                   palette = c("#00AFBB", "#E7B800", "#FC4E07", "blue", "grey"),
                   addEllipses = TRUE, # Ellipses de concentration
                   legend.title = "Groups"
      )
      dev.off()

      # This shows us that there is a problem with plate 4 -> take it out!

      pdf("PlotindividusPcoACytokinespaysellipse.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      factoextra::fviz_pca_ind(res.pca,
                   geom.ind = "point", # Montre les points seulement (mais pas le "text")
                   col.ind = as.factor(data[, 32]), # colorer by groups
                   palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                   addEllipses = TRUE, # Ellipses de concentration
                   legend.title = "Groups"
      )
      dev.off()


#### make models for SIBO status ####
#### Step 1: Pick Cytokines Based on Individual Relative Signifiance with no correction for cofactors using Wilcox Test for SIBO without plate 4. ####

      data_no4 = dplyr::filter(datacytokine, plate!=4)
      df_wilcox <- as.data.frame(data_no4[, c(10:39)])
      meta_wilcox <- as.data.frame(data_no4) #take metadata

      dim(df_wilcox)
      df_wilcox=t(df_wilcox)
      dim(meta_wilcox)

      meta_wilcox$sibo<-as.factor(as.character(meta_wilcox$sibo))

      MW.p = apply(df_wilcox,1,
                   function(x) wilcox.test(c(x)~meta_wilcox$sibo)$p.value)
      p.res = data.frame(row.names(df_wilcox),MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni")
      #export results
      write.csv(p.res,"Duodenal.sibo.noplate4.csv")
      #View(p.res) ## multiple are surviving, even Bonferoni!

#### Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 1,
                     function(x) cor.test(c(x), meta_wilcox$ufcml, method="pearson")$p.value)
      estcorr <- apply(df_wilcox,1,
                       function(x) cor.test(c(x), meta_wilcox$ufcml, method="pearson")$estimate)
      corres = data.frame(row.names(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      #export results
      write.csv(corres,"duodenal.sibo.noplate4.csv")
      #View(corres) ## none is surviving multiple correction

#### Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
#### prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox

      #add other categorical factors, etc.
      df <- data.frame(t(data1))
      dim(df)
      df$Country <- meta_wilcox$pays
      df$batch <- meta_wilcox$plate
      df$sibo <- meta_wilcox$sibo

#### make a loop for logistic regressions ####
      library(broom)
      long = reshape2::melt(df, id.vars = c("sibo",  "batch", "Country"))
      long = dplyr::filter(long, sibo!="") ## keep only the ones with valid data for sibo

      logresults<- long %>%
        dplyr::group_by(variable) %>%
        dplyr::do(tidy(glm(as.factor(sibo) ~ value + batch + Country, .,  family=binomial))) %>%
        dplyr::filter(term == "value") %>%
        dplyr::mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>%
        dplyr::ungroup()%>%
        dplyr::select(variable, Beta, SE, "p.value") %>%
        as.data.frame()

      #View(logresults)
      logresults <- logresults[order(logresults$p.value), ]
      logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
      write.csv(logresults,"Duodenal.sibo.multi.noplate4.csv")

#### make box plots of the ones showing a trend ####
      df = dplyr::filter(df, sibo!="") ## keep only the ones with valid data for stunted

      par(mfrow=c(1,1))
      boxplot(df$fgfbasicconc~df$sibo, data= df, main="Relative abundance of FGF basic by SIBO status",
              xlab="SIBO no/yes", ylab="Concentration in duodenum")

      boxplot(df$mip1alphaconc~df$sibo, data= df, main="Relative abundance of MIP1 alpha by SIBO status",
              xlab="SIBO no/yes", ylab="Concentration in duodenum")

      boxplot(df$il1bconc~df$sibo, data= df, main="Relative abundance of IL1beta by SIBO status",
              xlab="SIBO no/yes", ylab="Concentration in duodenum")

      boxplot(df$vegfconc~df$sibo, data= df, main="Relative abundance of VEGF by SIBO status",
              xlab="SIBO no/yes", ylab="Concentration in duodenum")

      boxplot(df$rantesconc~df$sibo, data= df, main="Relative abundance of Rantes by SIBO status",
              xlab="SIBO no/yes", ylab="Concentration in duodenum")

      boxplot(df$il12conc~df$sibo, data= df, main="Relative abundance of IL2 by SIBO status",
              xlab="SIBO no/yes", ylab="Concentration in duodenum")

      boxplot(df$il6conc~df$sibo, data= df, main="Relative abundance of IL6 by SIBO status",
              xlab="SIBO no/yes", ylab="Concentration in duodenum")

      boxplot(df$egfconc~df$sibo, data= df, main="Relative abundance of EGF alpha by SIBO status",
              xlab="SIBO no/yes", ylab="Concentration in duodenum")

      boxplot(df$mcp1conc~df$sibo, data= df, main="Relative abundance of MCP1 by SIBO status",
              xlab="SIBO no/yes", ylab="Concentration in duodenum")

      boxplot(df$il1raconc~df$sibo, data= df, main="Relative abundance of IL1 receptor by SIBO status",
              xlab="SIBO no/yes", ylab="Concentration in duodenum")


#### calculate the fold change  ####
      # first calculate the mean for each bile acid
      data_no4$sibo<-as.factor(data_no4$sibo)
      sibono<- dplyr::filter(data_no4, sibo=="no")
      siboyes<- dplyr::filter(data_no4)

      sibono<-sibono[, c(10:39)]
      siboyes<-siboyes[, c(10:39)]

      mean_no<- colMeans(sibono)
      mean_yes<- colMeans(siboyes)

      # now calculate the fold change
      foldchange<- as.matrix( gtools::foldchange(mean_yes, mean_no))
      no  <- as.matrix(mean_no)
      yes <- as.matrix(mean_yes)

      results<- cbind(mean_no, mean_yes, foldchange)
      colnames(results)<-c("mean no", "mean yes", "foldchange")
      results<-as.data.frame(results)
      results$variable <-rownames(results)

#### make a plot with the significant cytokines and their fold change ####
      # now take only the significant results and merge the fold change with the results
      logresults_sig = dplyr::filter(logresults, rel.fdr <= 0.05)
      results_sig    = dplyr::filter(results, variable %in% logresults_sig$variable)
#- 4600

      dim(results_sig) # check if dimensions are ok
      dim(logresults_sig) # check if dimensions are ok

      merged<- merge(results_sig, logresults_sig)
      merged<- merged[order(merged$foldchange ), ]

      #make the graph
      pdf("differentially_present_cytokines_duodenum_sibo.noplate4.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 4, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=merged, aes(x = reorder(variable, foldchange), y=foldchange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() +
        theme(axis.text.x = element_text(size=12))+
        theme(axis.text.y = element_text(size=12))+
        theme(axis.text.y = element_text(size=12))+
        xlab("")+
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("fold change")+
        ggtitle("Significantly different cytokines by SIBO status in duodenum, plate4 excluded")
      dev.off()




#- 4628
#### Is there correlation between microbiota composition in the SIBO and pro-inflammatory cytokines? ####

      samp = correlation_setup( dataset = dffilteredduodenal
                                , variables = c("calprotectineggdeps", "aatmggdeps" )
                                , pattern = "conc"
                                , level = "Species")

      #-- correlation at the species level
      associate_plot( x = samp$cyto, y = samp$bact, mode = "matrix" )

      ## is there correlation between the bacterial taxa?
      bact.table = associate_plot( x = samp$ bact, mode = "table", pdfname = "heatmap_bacs_duodenum.pdf" )

      ## is there correlation between the cytokines?
      cyto.table = associate_plot( x = samp$ cyto, mode = "table", pdfname = "heatmap_cyto_duodenum.pdf" )


#- 4646
#### duodenal samples: restrict the analysis on Bangui only and Species, Genus, Class, Order ####
  for (level in c("Species","Genus","Class","Order")){

       samp = correlation_setup( dataset = dffilteredduodenalBangui
                                , variables = c("calprotectineggdeps", "aatmggdeps" )
                                , pattern = "conc"
                                , level = level )

      #-- correlation at the species level
      associate_plot( x = samp$cyto, y = samp$bact, mode = "matrix" )

      ## is there correlation between the bacterial taxa?
      bact.table = associate_plot( x = samp$ bact, mode = "table"
                                   , pdfname = paste0("heatmap_bacs", level,"_duodenum_Bangui.pdf" ))
      bact.table = associate_plot( x = samp$ bact_unrestrict, mode = "table"
                                   , pdfname = paste0("heatmap_bacs", level,"_duodenum_Bangui_all.pdf" ))

      ## is there correlation between the cytokines?
      cyto.table = associate_plot( x = samp$ cyto, mode = "table"
                                   , pdfname = paste0("heatmap_cyto", level,"_duodenum_Bangui.pdf" ))

  }# end loop

#- 4670
#### duodenal samples: restrict the analysis on Mada only and Species####
   for (level in c("Species","Genus","Class" )){

       samp = correlation_setup( dataset = dffilteredduodenalMada
                                , variables = c("calprotectineggdeps", "aatmggdeps" )
                                , pattern = "conc"
                                , level = level )

      #-- correlation at the species level
      associate_plot( x = samp$cyto, y = samp$bact, mode = "matrix" )

      ## is there correlation between the bacterial taxa?
      bact.table = associate_plot( x = samp$ bact, mode = "table"
                                   , pdfname = paste0("heatmap_bacs", level,"_duodenum_Mada.pdf" ))
      bact.table = associate_plot( x = samp$ bact_unrestrict, mode = "table"
                                   , pdfname = paste0("heatmap_bacs", level,"_duodenum_Mada_all.pdf" ))

      ## is there correlation between the cytokines?
      cyto.table = associate_plot( x = samp$ cyto, mode = "table"
                                   , pdfname = paste0("heatmap_cyto", level,"_duodenum_Mada.pdf" ))

  }# end loop

#- 4694
#### duodenal samples: Bacterial co-occurrence by SIBO status ####
  for (sibo in c( "yes", "no" )){ #sibo = "yes"

       samp = correlation_setup( dataset =
                        dffilteredduodenal %>% phyloseq::subset_samples( sibo == sibo )
                                , variables = c("calprotectineggdeps", "aatmggdeps" )
                                , pattern = "conc"
                                , level = "Species" )

      #-- correlation at the species level
      associate_plot( x = samp$cyto, y = samp$bact, mode = "matrix" )

      ## is there correlation between the bacterial taxa?
      bact.table = associate_plot( x = samp$ bact, mode = "table"
                                   , pdfname = paste0("heatmap_bacs", level,"_duodenum_sibo_", sibo,".pdf" ))
      # bact.table = associate_plot( x = samp$ bact_unrestrict, mode = "table"
      #                              , pdfname = paste0("heatmap_bacs", level,"_duodenum_sibo_", sibo,"_all.pdf" ))
      # ## is there correlation between the cytokines?
      # cyto.table = associate_plot( x = samp$ cyto, mode = "table"
      #                              , pdfname = paste0("heatmap_cyto", level,"_duodenum_sibo_", sibo,".pdf" ))

  }# end loop



#- 4720
#### duodenal samples: Bacterial co-occurence by haz status on Species level, and country of origin ####
       for (haz in c( "severely stunted", "moderately stunted" )){ # haz = "moderately stunted"

       samp = correlation_setup( dataset =
                        dffilteredduodenal %>%
                          phyloseq::subset_samples( haz == haz )
                                , variables = c("calprotectineggdeps", "aatmggdeps" )
                                , pattern = "conc"
                                , level = "Species" )

      #-- correlation at the species level
      associate_plot( x = samp$cyto, y = samp$bact, mode = "matrix" )

      ## is there correlation between the bacterial taxa?
      bact.table = associate_plot( x = samp$ bact, mode = "table"
              , pdfname = paste0("heatmap_bacs", level,"_duodenum_haz_", gsub(pattern = " ", "", haz), ".pdf" ))
      # bact.table = associate_plot( x = samp$ bact_unrestrict, mode = "table"
      #                              , pdfname = paste0("heatmap_bacs", level,"_duodenum_sibo_", sibo,"_all.pdf" ))
      # ## is there correlation between the cytokines?
      # cyto.table = associate_plot( x = samp$ cyto, mode = "table"
      #                              , pdfname = paste0("heatmap_cyto", level,"_duodenum_sibo_", sibo,".pdf" ))

  }# end loop


#- 4746
#### duodenal samples: Bacterial co-occurence by haz status on Species level by country of origin ####
   for( pays in c("Madagascar", "RCA")){ #pays = "RCA"
       for (haz in c( "severely stunted", "moderately stunted" )){ # haz = "moderately stunted"

       samp = correlation_setup( dataset =
                        dffilteredduodenal %>%
                          phyloseq::subset_samples( haz == haz ) %>%
                          phyloseq::subset_samples( pays == pays )
                                , variables = c("calprotectineggdeps", "aatmggdeps" )
                                , pattern = "conc"
                                , level = "Species" )

      #-- correlation at the species level
      associate_plot( x = samp$cyto, y = samp$bact, mode = "matrix" )

      ## is there correlation between the bacterial taxa?
      bact.table = associate_plot( x = samp$ bact, mode = "table"
              , pdfname = paste0("heatmap_bacs", level,"_duodenum_haz_", gsub(pattern = " ", "", haz), pays, ".pdf" ))
      # bact.table = associate_plot( x = samp$ bact_unrestrict, mode = "table"
      #                              , pdfname = paste0("heatmap_bacs", level,"_duodenum_sibo_", sibo,"_all.pdf" ))
      # ## is there correlation between the cytokines?
      # cyto.table = associate_plot( x = samp$ cyto, mode = "table"
      #                              , pdfname = paste0("heatmap_cyto", level,"_duodenum_sibo_", sibo,".pdf" ))

  }# end loop
      }# end country loop



#- 4776
##### rel. abundance of specific taxa, nutritional status, SIBO and inflammation markers on merged dataset in duodenum ####
      dffilteredduodenalallsibo = phyloseq::subset_samples(dffilteredduodenal, sibo!="")
      dffilteredduodenalallsiboGenus = phyloseq::tax_glom(dffilteredduodenalallsibo, "Genus")

      dffilteredduodenalallsiboGenusrel<- microbiome::transform(dffilteredduodenalallsiboGenus, "compositional")

      sampledata<-data.frame(t(phyloseq::otu_table(dffilteredduodenalallsiboGenusrel)))
      metadata<-data.frame(phyloseq::sample_data(dffilteredduodenalallsiboGenus))
      sampletaxa<-data.frame(phyloseq::tax_table(dffilteredduodenalallsiboGenus))
      which(duplicated(sampletaxa$Genus))
      Genusnames<-as.vector(sampletaxa$Genus)
      length(Genusnames)

      Genusnames[111]<-"uncultured Elysipelotrichaceae"
      Genusnames[131]<-"uncultured Lachnospiraceae"
      Genusnames[164]<-"uncultured Coriobacteriaceae"
      Genusnames[167]<-"uncultured Porphyromonadaceae"
      Genusnames[168]<-"uncultured Rhodospirillaceae"
      Genusnames[207]<-"uncultured Ruminococcaceae"
      Genusnames[211]<-"uncultured Family_XIII"

      row.names(sampletaxa)<-Genusnames # Make a check to see if newly assigned annotation is right

      colnames(sampledata)<-Genusnames
#- 4801
      # Streptococcus
      pdf("SIBOStreptoBoxplot.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 5, # define plot width and height. completely up to user.
        height = 10)
      ggplot(sampledata, aes(x=metadata$pays, y=sampledata$Streptococcus, fill=metadata$sibo)) +
       geom_boxplot()
      dev.off()

      for( bacteria in c("Fusobacterium", 'Streptococcus'
                         , "Campylobacter","Rothia","Veillonella"
                         , "Neisseria", "Haemophilus", "Actinobacillus"
                         , "Porphyromonas", "Granulicatella"
                         , "Prevotella")){

        # boxplots
       ggplot(sampledata, aes(x=metadata$pays, y=sampledata[, bacteria], fill=metadata$sibo)) +
       geom_boxplot() + labs( y = bacteria , x = "pays")

       ggplot(sampledata, aes(x=metadata$pays, y=sampledata[, bacteria], fill=metadata$haz)) +
       geom_boxplot() + labs( y = bacteria , x = "pays")

       # wilcox test on group variable across strata
        out = wilcox_strata_group(metadata = metadata, sampledata = sampledata
                          , bacteria = bacteria
                          , group = 'sibo', strata = 'pays')
      print(out)
      # correlation test
      out = cor.test(sampledata[, bacteria], as.numeric(metadata$haz_cont), method = "spearman", use = "complete.obs")
      print(out)
      # ggscatter plots for different controls
      out = make_ggscatter_plots( sampledata = sampledata, metadata = metadata
                            , bacteria = bacteria , controls = c("haz_cont","ufcml") )
      print( out )
      ## ggscatter plots for different controls
      temp = metadata %>%
        dplyr::rename(  calprotectine_duodenum = calprotectinesurliquideduodn, AAT_duodenum = aatsurldmgml )

      out = make_ggscatter_plots( sampledata = sampledata
                  , metadata = temp
                  , bacteria = bacteria
                  , controls = c("calprotectine_duodenum","AAT_duodenum") )

      print(out )
      }# end loop



      ## Proteobacteria

      dffilteredduodenalallsiboPhylum = phyloseq::tax_glom(dffilteredduodenalallsibo, "Phylum")
      dffilteredduodenalallsiboPhylumrel<- microbiome::transform(dffilteredduodenalallsiboPhylum, "compositional")
      sampledata<-data.frame(t(phyloseq::otu_table(dffilteredduodenalallsiboPhylumrel)))
      metadata<-data.frame(phyloseq::sample_data(dffilteredduodenalallsiboPhylumrel))
      sampletaxa<-data.frame(phyloseq::tax_table(dffilteredduodenalallsiboPhylumrel))
      which(duplicated(sampletaxa$Phylum))
      colnames(sampledata)<-sampletaxa$Phylum


#- 4860
    for( bacteria in c("Proteobacteria" )){

        # boxplots
       ggplot(sampledata, aes(x=metadata$pays, y=sampledata[, bacteria], fill=metadata$sibo)) +
       geom_boxplot() + labs( y = bacteria , x = "pays")

       ggplot(sampledata, aes(x=metadata$pays, y=sampledata[, bacteria], fill=metadata$haz)) +
       geom_boxplot()+ labs( y = bacteria , x = "pays")

       # wilcox test on group variable across strata
       out = wilcox_strata_group(metadata = metadata, sampledata = sampledata
                          , bacteria = bacteria
                          , group = 'sibo', strata = 'pays')
      print(out)
      # correlation test
      out = cor.test(sampledata[, bacteria], as.numeric(metadata$haz_cont), method = "spearman", use = "complete.obs")
      print(out)

      # ggscatter plots for different controls
      out = make_ggscatter_plots( sampledata = sampledata, metadata = metadata
                            , bacteria = bacteria , controls = c("haz_cont","ufcml") )

      print(out)
      ## ggscatter plots for different controls
      temp = metadata %>%
        dplyr::rename(  calprotectine_duodenum = calprotectinesurliquideduodn, AAT_duodenum = aatsurldmgml )

      out = make_ggscatter_plots( sampledata = sampledata
                  , metadata = temp
                  , bacteria = bacteria
                  , controls = c("calprotectine_duodenum","AAT_duodenum") )
      print(out)

      }# end loop
      
#### Association of SIBO with Helicobacter and gastric pH ####
      dffilteredgastricallsibo = phyloseq::subset_samples(dffilteredgastric, sibo!="")
      dffilteredgastricallsiboGenus = phyloseq::tax_glom(dffilteredgastricallsibo, "Genus")
      
      dffilteredgastricallsiboGenusrel<- microbiome::transform(dffilteredgastricallsiboGenus, "compositional")
      
      
      sampledata<-data.frame(t(phyloseq::otu_table(dffilteredgastricallsiboGenusrel)))
      metadata<-data.frame(phyloseq::sample_data(dffilteredgastricallsiboGenus))
      sampletaxa<-data.frame(phyloseq::tax_table(dffilteredgastricallsiboGenus))
      View(sampletaxa)
      which(duplicated(sampletaxa$Genus))
      Correspondance<-as.data.frame(sampletaxa$Genus, row.names(sampletaxa))
      Correspondance<-t(Correspondance)
      sampledata2<-rbind(sampledata, Correspondance)
      
      dim(sampledata2)
      colnames(sampledata2)<-sampledata2[112, ]
      sampledata2<-sampledata2[-112, ]
      sampledata3<-data.frame(sampledata2)
      sampledata3$Helicobacter<-as.numeric(sampledata3$Helicobacter)
      
      
      
      
      # Helicobacter
      pdf("SIBOHelicobacterBoxplot.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 5, # define plot width and height. completely up to user.
          height = 10)
      ggplot(sampledata3, aes(x=metadata$sibo, y=sampledata3[, "Helicobacter"], fill=metadata$sibo)) +
        geom_boxplot()
      dev.off()
      
      siboyes = dplyr::filter(sampledata3, metadata$sibo=="yes")
      sibono = dplyr::filter(sampledata3, metadata$sibo=="no")
      
      dim(siboyes)
      dim(sibono)
      
      wilcox.test(siboyes$Helicobacter, sibono$Helicobacter) # p-value = 0.752
      
      result = cor(sampledata3[, "Helicobacter"], metadata$ufcml, method = "spearman")
      result  #0.004181051
      
      result = cor(sampledata3[, "Helicobacter"], metadata$ph_estomac, method = "spearman")
      result  #-0.1473176
      
      result = cor(metadata$ufcml, metadata$ph_estomac, method = "spearman")
      result  #0.1661987
      
      
      test<-as.data.frame(as.matrix(cbind(sampledata3[, "Helicobacter"], metadata$ufcml, metadata$ph_estomac)))
      colnames(test)<-c("Helicobacter", "ufcml", "ph_estomac")
      
      
      ggpubr::ggscatter(test, x = "Helicobacter", y = "ufcml",
                        add = "reg.line", conf.int = TRUE,
                        cor.coef = TRUE, cor.method = "spearman",
                        xlab = "Relative abundance of Helicobacter spp.", ylab = "CFU/ml in the small intestine")
      
      ggpubr::ggscatter(test, x = "Helicobacter", y = "ph_estomac",
                        add = "reg.line", conf.int = TRUE,
                        cor.coef = TRUE, cor.method = "spearman",
                        xlab = "Relative abundance of Helicobacter spp.", ylab = "pH in stomach")
      
      ggpubr::ggscatter(test, x = "ufcml", y = "ph_estomac",
                        add = "reg.line", conf.int = TRUE,
                        cor.coef = TRUE, cor.method = "spearman",
                        xlab = "CFU/ml in the small intestine", ylab = "pH in stomach")
      
      
      
      # Now assess also for presence/absence of Helicobacter
      
      prevalence = function(x){ #this only returns prevalence counts per genus
        x[x >= 1] <- 1
        return(x)
      }
      
      dffilteredgastricallsiboGenuspres<- dffilteredgastricallsiboGenus %>% #this produces prevalence "counts" for each Genus
        transform_sample_counts(fun = prevalence)
      
      Helicobacter<-t(otu_table(dffilteredgastricallsiboGenuspres))
      colnames(Helicobacter)
      Helicobacter<-Helicobacter[, "ASV21"]
      length(Helicobacter)
      Helicobacter<-as.character(Helicobacter)
      Helicobacter
      
      
      sibo<-sample_data(dffilteredgastricallsiboGenuspres)$sibo
      length(sibo)
      prop.table(table(sibo,Helicobacter), 1)
      chisq.test(table(sibo,Helicobacter))
      
      Helicobacter<-as.data.frame(Helicobacter)
      row.names(Helicobacter)<-colnames((otu_table(dffilteredgastricallsiboGenuspres)))
      
      helicobacteryes = dplyr::filter(data.frame(sample_data(dffilteredgastricallsiboGenuspres)), Helicobacter$Helicobacter=="1")
      helicobacterno = dplyr::filter(data.frame(sample_data(dffilteredgastricallsiboGenuspres)), Helicobacter$Helicobacter=="0")
      
      dim(helicobacteryes)
      dim(helicobacterno)
      
      test2<-data.frame(sample_data(dffilteredgastricallsiboGenuspres))
      test2
      
      wilcox.test(helicobacteryes$ufcml, helicobacterno$ufcml) # p-value = 0.83
      wilcox.test(helicobacteryes$ph_estomac, helicobacterno$ph_estomac) # p-value = 0.115
      
      
      
      
      

#### generate output for Piecrust 2 DOES NOT WORK YET ####
      library(biomformat)
      library(dada2)
      View(head(otu_table(dffiltered4)))
      otu<-as(otu_table(dffiltered4), "matrix")
      otu_biom<-make_biom(data=otu)
      write_biom(otu_biom, "otu_biom")
      
      # include the sequences into the sequence table
      otutable<-otu_table(dffiltered4)
      otutable$ASV<-row.names(otutable)
      colnames<-colnames(otu_table(dffiltered4))
      otutable<-as.data.frame(otutable)
      colnames(otutable)<-c(colnames, "ASV")
      
      taxtable<-tax_table(dffiltered4)
      taxtable$ASV<-row.names(taxtable)
      taxtable<-as.data.frame(taxtable)
      colnames(taxtable)<-c(tax_table(dffiltered4), "ASV")
      
      forpiecrust<-cbind(otutable, taxtable, by="ASV")
      dim(forpiecrust)
      dim(taxtable)
      View(head(taxtable))
      
      View(head(forpiecrust))
      colnames(forpiecrust)
      row.names(forpiecrust)<-forpiecrust[, 933]
      forpiecrust<-forpiecrust[, 1:924]
      forpiecrust$ASV<-row.names(forpiecrust)

      forpiecrust<-t(forpiecrust)
      View(head(forpiecrust))
      
      seqs <- colnames(forpiecrust)
      row.names(forpiecrust)
      
      write.table(seqs, "dada_seqs.txt",quote=FALSE)
      write.table(otab, "dada_table.txt",quote=FALSE,sep="\t")


      # next line to be performed on Terminal and replace seq by ASV
      grep -v '^x' dada_seqs.txt | awk '{print ">seq"$1"\n"$2}' > dada_seqs.fa; rm dada_seqs.txt
      
      