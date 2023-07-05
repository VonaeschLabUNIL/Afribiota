###### Analysis of dataset 18S Afribiota using a human blocking primer  #######
#### Vonaesch et al., February 2023 #######
##### install all necessary packages ####
library(tidyr)
library(DESeq2)
library(qualpalr)
library(ggplot2)
library(biomformat)
library(plyr)
library(dplyr)
library(microbiome)
library(data.table)
library(reshape2)
library(readstata13)
require(devtools)
library(vegan)
library(phyloseq)
library(binom)
library(dplyr)
library(tidyverse)
packageDescription("phyloseq")$Version 
packageDescription("vegan")$Version 
version # to get the version from R

##### set the work directory and set working environment ####
  setwd("~/Users/admin/Desktop/Papers and manuscripts/under revision/Eukaryome paper/Revised version October 22 including last datasets/Older things/Test")
  rm(list=ls(all=TRUE)) # go get a clean workspace 

##### read in the files  #####

seqs<-read.table("Table S16.txt", header=TRUE )
row.names(seqs)<-seqs$row_names
seqs<-seqs[, -1]
min(rowSums(seqs))
max(rowSums(seqs))
mean(rowSums(seqs))

asv<-otu_table(seqs, taxa_are_rows = FALSE)

# import the tax table
taxtable<- read.table("Table S15.csv", sep = ",", quote = "", header=TRUE) 
View(head(taxtable))
taxtable2<-filter(taxtable, taxtable$dataset=="Dalhousie")
taxtable2=taxtable2[, -2]
dim(taxtable2)
View(head(taxtable2))
taxonomy=tax_table(taxtable2)
taxonomy<-gsub("DAL", "", taxonomy, fixed = TRUE)
row.names(taxonomy)<-taxonomy[, 1]
View(head(taxonomy))
colnames(taxonomy)<-colnames(taxtable2)
taxonomy=taxonomy[,-1]

physeq= phyloseq(asv, taxonomy)
sample_names(physeq)

metadata<-read.csv("Table S17.csv", sep=",")
View(metadata)
row.names(metadata)<-metadata$SampleID
metadata<-sample_data(metadata)

data<-merge_phyloseq(physeq, metadata)
data
View(head(tax_table(data)))

# take out reads that are not assigned beyond a very high rank
data2<-subset_taxa(data, rank1=="Eukaryota")
data2<-subset_taxa(data2, rank7!="Unassigned_Eukaryota")
data2<-subset_taxa(data2, rank7!="Unassigned_Metazoa")

#### Prune singlets ####      
prune_singlets= function(x) {x[x <= 1] <- 0
return(x)}

df <- data2 %>% #this produces a maximum possible prevalence count per species per ASV
  transform_sample_counts(fun = prune_singlets)

#### clean metadata and tax-table, filter your samples according to the blanks and analyze what leads to differences in total sample read count #####
sample_data(df)$SampleType<-as.factor(sample_data(df)$SampleType)
levels(sample_data(df)$SampleType)
sum(sample_sums(df))
sample_data(df)$read_count<-as.numeric(sample_sums(df))

##clean up your metadata for full afribiota dataset
#create a variable age in years
sample_data(df)$ageyears<-cut(as.numeric(as.character(sample_data(df)$age)), c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
which(is.na(sample_data(df)$ageyears)) # the controls have no assocaited age year
levels(sample_data(df)$ageyears) <- c("2-3 years", "3-4 years", "4-5 years")
levels(sample_data(df)$ageyears)

sample_data(df)$haz<-as.factor(sample_data(df)$haz)
levels(sample_data(df)$haz) 

sample_data(df)$pays<-as.factor(sample_data(df)$pays)
levels(sample_data(df)$pays) 

sample_data(df)$stunted<-as.factor(sample_data(df)$stunted)
levels(sample_data(df)$stunted)  # we need to transform this to factor and give this labels

levels(sample_data(df)$stunted) <- c("non-stunted", "stunted")

sample_data(df)$whz<-as.factor(sample_data(df)$whz)
levels(sample_data(df)$whz)

sample_data(df)$whz_cont<-as.numeric(sample_data(df)$whz_cont)

sampledata=as.data.frame(sample_data(df))
sum=sum(sampledata$read_count)
sum ##2365060

mean=mean(sampledata$read_count) ##7604.695 this is the mean sample sums
mean

median=median(sampledata$read_count) ##6057 this is the median sample sums
median

# kick out samples that do not match the criteria
sample_data(df)$ph_estomac<-as.numeric(sample_data(df)$ph_estomac)
sample_data(df)$ph_intestin<-as.numeric(sample_data(df)$ph_intestin)

df <- df %>%
  subset_samples(raison_hospi=="Recrutement communautaire") # to filter out samples which are not recruited in the community

df <- df %>%
  subset_samples(whz_cont <=2)

df <- df %>%
  subset_samples(whz_cont >=-2)  # to keep only samples with no acute undernutrition or obesity

dfa <- df %>%
  subset_samples(SampleType=="feces")

dfc <- df %>%
  subset_samples(SampleType == "duodenal")
dfc <- dfc %>%
  subset_samples(ph_intestin>=5)

df<-merge_phyloseq(dfa, dfc)

df<-subset_samples(df, pays!="")

df # we are left with 287 samples

length(which(sample_data(df)$SampleType=="duodenal")) # we start with 26 samples
length(which(sample_data(df)$SampleType=="feces")) # we start with 261 samples 


# now filter out low sample sums samples
dff = prune_samples((sample_data(df)$read_count>=500)==TRUE, df) # only keep samples with more than 500 reads
dim(sample_data(dff))
dim(tax_table(dff)) #940 ASV's

length(which(sample_data(dff)$SampleType=="duodenal")) # we start with 23 samples
length(which(sample_data(dff)$SampleType=="feces")) # we start with 241 samples 

# filter out low abundance taxa, as they are likely contaminants. Here, we decided that an ASV needs to have at least 50 sequences in at least 1% of the samples
dff=filter_taxa(dff, function(x) sum(x > 50) > (0.01*length(x)), TRUE)
dim(tax_table(dff)) #173 ASV's. Note: we lost a lot of taxa, mainly due to the fact that they are spuriously distributed in between the samples!

#track the evolution of your sample numbers
table(sample_data(dff)$SampleType) # 23 duodenal and 241 feces

# now make subsets of samples according to SampleType
df_feces=subset_samples(dff, sample_data(dff)$SampleType=="feces") 
df_duodenal=subset_samples(dff, sample_data(dff)$SampleType=="duodenal")

sum_feces=sum(sample_sums(df_feces))
sum_feces ## 1671100

sum_duodenal=sum(sample_sums(df_duodenal))
sum_duodenal ## 271837


#### Look what is influencing the overall sequence count for feces: run, Sample Type, age, stunting, country ####
    sampledata= as.data.frame(sample_data(dff))
    sampledata$SampleType<-as.factor(sampledata$SampleType)
    sampledata$ageyears<-as.factor(sampledata$ageyears)
    
    boxplot(sampledata$read_count~sampledata$SampleType, main= "Sample Sums", xlab="Sample Type", ylab= "")
    kruskal.test(sampledata$read_count~sampledata$SampleType) #  p-value = 2.15e-05
    
    sampledata$row.names=row.names(sampledata)
    sampledata=as.data.frame(sampledata)
    sampledata=unclass(sampledata)
    
    # now make a multivariate model to see who is contributing independently 
    
    samplesums_total <- lm(read_count ~ calprotectinelevel + haz  + ageyears + pays + SampleType, data= sampledata, na.action=na.omit)
    summary(samplesums_total) # to see the results
    
    # redo the analysis just for the feces
    sampledata= sample_data(df_feces)
    
    boxplot(sampledata$read_count~sampledata$pays, main= "Sample Sums", xlab="pays of origin", ylab= "")
    kruskal.test(sampledata$read_count~sampledata$pays) # p-value = 1.701e-13
    
    boxplot(sampledata$read_count~sampledata$stunted + sampledata$pays, main= "Sample Sums", xlab="Stunting", ylab= "")
    kruskal.test(sampledata$read_count~sampledata$stunted) # p-value = 0.2081
    
    boxplot(sampledata$read_count~sampledata$ageyears, main= "Sample Sums", xlab="Age in years", ylab= "")
    kruskal.test(sampledata$read_count~sampledata$ageyears) # p-value = 0.1645
    
    boxplot(sampledata$read_count~sampledata$calprotectinelevel, main= "Sample Sums", xlab="Calprotectinelevel", ylab= "")
    kruskal.test(sampledata$read_count~sampledata$calprotectinelevel) # p-value = 0.03721
    
    sampledata$row.names=row.names(sampledata)
    sampledata=as.data.frame(sampledata)
    sampledata=unclass(sampledata)
    
    #now make a multivariate model to see who is contributing independently
    
    samplesums_total <- lm(read_count ~  calprotectinelevel + haz  + ageyears + pays, data=sampledata, na.action=na.omit)
    summary(samplesums_total) # to see the results
    
    
    #### describe your dataset in terms of sequences contributed by plants, mammals etc. ####
    
    df_clean<-dff 
    df_clean## 264 samples, 173 taxa 
    get_taxa_unique(df_clean, "rank1") # to see how many kingdoms we have, we have one kingdom: Eukaryota
    get_taxa_unique(df_clean, "rank4") # 
    
    get_taxa_unique(df_clean, "rank3") ## 10 unique taxa
    tax_table(df_clean)[, 3][which(tax_table(df_clean)[, 3] == "")] <- NA
    tax_table(df_clean)[, 3][which(tax_table(df_clean)[, 3] == "uncultured eukaryote")] <- NA
    
     #how many reads are contributed by mammals?
    read_count<-sample_data(df_clean)$read_count
    length(read_count)
    
    mammal_reads <- df_clean %>%
      subset_taxa(rank4 == "Mammalia")
    mammal_reads
    
    mammal_sums <- sample_sums(mammal_reads)
    sample_data(df_clean)$mammal_sums <- sample_sums(mammal_reads)
    mammal_sums
    plot(sort(mammal_sums))
    length(mammal_sums)
    
    sampledata_c=sample_data(df_clean)
    sampledata_c$row.names=row.names(sampledata_c)
    sampledata_c=as.data.frame(sampledata_c)
    sampledata_c=unclass(sampledata_c)
    
    boxplot(sample_data(df_clean)$mammal_sums~sample_data(df_clean)$SampleType, main= "Mammal Sums", xlab="SampleType", ylab= "", data=sampledata_c)
    kruskal.test(data=sampledata_c, sample_data(df_clean)$mammal_sums~sample_data(df_clean)$SampleType)  #p-value = 1.97e-13
    
    #how many reads are contributed by vertebrates?
    vertebrate_reads = df_clean %>%
      subset_taxa((rank4== "Vertebrata") | (rank5== "Vertebrata") |(rank6== "Vertebrata"))
    
    vertebrate_sums <- sample_sums(vertebrate_reads)
    sample_data(df_clean)$vertebrate_sums <- sample_sums(vertebrate_reads)
    length(vertebrate_sums)
    
    #how many reads are contributed by plants?
    plant_reads= df_clean %>%
      subset_taxa(rank2=="Archaeplastida")
    
    plant_sums <- sample_sums(plant_reads)
    sample_data(df_clean)$plant_sums <- sample_sums(plant_reads)
    length(plant_sums)
    
    #how many reads are contributed by microeukaryotes?
    microeuk_reads= df_clean %>%
      subset_taxa(rank4!="Vertebrata" & rank5!="Vertebrata" & rank6!="Vertebrata" & rank2!="Archaeplastida" & rank4!="Mammalia" & rank2!="Unassigned" ) 
    
    microeuk_reads_feces= df_feces %>%
      subset_taxa(rank4!="Vertebrata" & rank5!="Vertebrata" & rank6!="Vertebrata" & rank2!="Archaeplastida" & rank4!="Mammalia" & rank2!="Unassigned" ) 
    
    microeuk_reads_duodenal= df_duodenal %>%
      subset_taxa(rank4!="Vertebrata" & rank5!="Vertebrata" & rank6!="Vertebrata" & rank2!="Archaeplastida" & rank4!="Mammalia" & rank2!="Unassigned" ) 
    
    microeuk_sums <- sample_sums(microeuk_reads)
    sample_data(df_clean)$microeuk_sums <- sample_sums(microeuk_reads)
    mean(sample_data(df_clean)$microeuk_sums)
    plot(sort(microeuk_sums))
    length(microeuk_sums)
    
    microeuk_sums_feces <- sample_sums(microeuk_reads_feces)
    sample_data(df_feces)$microeuk_reads_feces <- sample_sums(microeuk_reads_feces)
    mean(sample_data(df_feces)$microeuk_reads_feces) ## 2879.029
    plot(sort(microeuk_sums_feces))
    
    microeuk_sums_duodenal <- sample_sums(microeuk_reads_duodenal)
    sample_data(df_duodenal)$microeuk_reads_duodenal <- sample_sums(microeuk_reads_duodenal)
    mean(sample_data(df_duodenal)$microeuk_reads_duodenal) # 49.21739
    plot(sort(microeuk_sums_duodenal)) ## note, there are much less sequences in duodenal!!
    
    
    # how much of which kind of sequences in the different sample types?
    o<-ordered(sample_data(df_clean)$SampleType, levels=c( "duodenal", "feces"))
    
    pdf("SampleTypevsmicroeukreads.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    boxplot(sample_data(df_clean)$microeuk_sums~o, main= "Microeukaryotes Sums", xlab="Sample Type", ylab= "", data=sampledata_c)
    dev.off()
    
    boxplot(sample_data(df_clean)$plant_sums~o, main= "Plant Sums", xlab="Sample Type", ylab= "", data=sampledata_c)
    boxplot(sample_data(df_clean)$mammal_sums~o, main= "Mammal Sums", xlab="Sample Type", ylab= "", data=sampledata_c)
    
    
    # now look a bit closer in the overall dataset and then in fecal samples
    microeuk_reads_feces= df_feces %>%
      subset_taxa(rank4!="Vertebrata" & rank5!="Vertebrata" & rank6!="Vertebrata" & rank2!="Archaeplastida" & rank4!="Mammalia" & rank2!="Unassigned" ) 
    
    sample_data(df_feces)$microeuk_reads <- sample_sums(microeuk_reads_feces)
    
    
    # export table of results for contribution by different reads
    length(read_count)
    length(mammal_sums)
    length(plant_sums)
    length(vertebrate_sums)
    length(microeuk_sums)
    
    ssums <- cbind(mammal_sums, plant_sums, vertebrate_sums, microeuk_sums,  read_count)
    ssums
    View(ssums)
    
    write.csv(ssums, "MammalPlantVertebratemicroeukunassignedreadsDAL.csv")
    
    mammal_percentage<-mammal_sums/read_count*100
    length(mammal_percentage)
    View(mammal_percentage)
    
    plant_percentage<- plant_sums/read_count*100
    length(plant_percentage)
    
    microeuk_percentage<-microeuk_sums/read_count*100
    length(microeuk_percentage)
    
    percentages<-cbind(mammal_percentage, plant_percentage, microeuk_percentage, read_count)
    View(percentages)
    write.csv(percentages, "MammalPlantmicroeukunassignedreadpercentagesDAL.csv")

    ####  make a graphic representation of percentage of plants, mammals and microeuks ####
    percentages<-as.data.frame(percentages)
    View(percentages)
    
    
    pdf("percentageplantsDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(percentages, aes(plant_percentage)) + 
      geom_histogram(color = "black", fill = "blue", binwidth = 1) +
      ggtitle("Distribution of plant percentage") + 
      xlab("percentage of plant sequences/total reads") + theme(axis.text = element_text(size = 16), axis.title = element_text(size=18, face="bold"), plot.title=element_text(size=22, face="bold"))
    dev.off() 
    
    pdf("percentagemammalsDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(percentages, aes(mammal_percentage)) + 
      geom_histogram(color = "black", fill = "red", binwidth = 1) +
      ggtitle("Distribution of mammal percentage") + 
      xlab("percentage of mammal sequences/total reads") + theme(axis.text = element_text(size = 16), axis.title = element_text(size=18, face="bold"), plot.title=element_text(size=22, face="bold")) +
      theme(axis.title.y = element_blank())
    dev.off() 
    
    pdf("percentagemicroeukDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(percentages, aes(x = microeuk_percentage)) + 
      geom_histogram(color = "black", fill = "purple", binwidth = 1) +
      ggtitle("Distribution of microeukaroytes percentage") + 
      xlab("percentage of microeukaryotic sequences/total reads") + theme(axis.text = element_text(size = 16), axis.title = element_text(size=18, face="bold"), plot.title=element_text(size=22, face="bold")) +
      theme(axis.title.y = element_blank())
    dev.off() 
    
 
    #### filter your dataset to kick out Vertebrates and Plants and make subset for feces and duodenal samples #####                 
    # Filter your data to get get rid of plants and mammals
    
    #remove mammal and plant and other vertebrate sequences as well as unassigned Eukaryotes
    dfclean <- df_clean %>%
      subset_taxa(rank4!= "Mammalia") %>%
      subset_taxa(rank5!= "Mammalia") %>%
      subset_taxa(rank6!= "Mammalia") %>%
      subset_taxa(rank4!= "Vertebrata") %>%
      subset_taxa(rank5!= "Vertebrata") %>%
      subset_taxa(rank6!= "Vertebrata") %>%
      subset_taxa(rank2!= "Archaeplastida") %>% #plants and algae
      subset_taxa(rank4!= "Porifera") %>% #sponges
      subset_taxa(rank4!= "Arthropoda") %>% #arthropods
      subset_taxa(rank5!= "Arthropoda") %>% #arthropods
      subset_taxa(rank6!= "Arthropoda") %>% #arthropods
      subset_taxa(rank4!= "Mollusca") %>% #mussels
      subset_taxa(rank2!="Unassigned")
    
      dfcleaninclunassigned<- df_clean %>%
      subset_taxa(rank4!= "Mammalia") %>%
      subset_taxa(rank5!= "Mammalia") %>%
      subset_taxa(rank6!= "Mammalia") %>%
      subset_taxa(rank4!= "Vertebrata") %>%
      subset_taxa(rank5!= "Vertebrata") %>%
      subset_taxa(rank6!= "Vertebrata") %>%
      subset_taxa(rank2!= "Archaeplastida") %>% #plants and algae
      subset_taxa(rank4!= "Porifera") %>% #sponges
      subset_taxa(rank4!= "Arthropoda") %>% #arthropods
      subset_taxa(rank5!= "Arthropoda") %>% #arthropods
      subset_taxa(rank6!= "Arthropoda") %>% #arthropods
      subset_taxa(rank4!= "Mollusca") #mussels
    
    dfclean # 127 taxa, 264 samples
    sample_data(dfclean)$read_count=sample_sums(dfclean)
    
    sum(sample_sums(dfclean)) 
    mean(sample_sums(dfclean)) 
    median(sample_sums(dfclean)) 
    max(sample_sums(dfclean)) 
    min(sample_sums(dfclean)) 
  
    length(which(sample_data(dfclean)$SampleType=="duodenal")) # we start with 23 samples
    length(which(sample_data(dfclean)$SampleType=="feces")) # we start with 241 samples 
    
    get_taxa_unique(dfclean, "rank2") ## 6 unique taxa 
    get_taxa_unique(dfclean, "rank3") ## 9 unique taxa 
    get_taxa_unique(dfclean, "rank4") ## 10 unique taxa 
    get_taxa_unique(dfclean, "rank5") ## 12 unique taxa 
    get_taxa_unique(dfclean, "rank6") ## 15 unique taxa 
    get_taxa_unique(dfclean, "rank7") ## 19 unique taxa 
    get_taxa_unique(dfclean, "rank8") ## 25 unique taxa 
    
    dfclean #127

    dfcleanfeces<-subset_samples(dfclean, SampleType=="feces")
    dfcleanfeces # 241 fecal samples
    table(sample_data(dfcleanfeces)$pays)
    
    dfcleanduodenal<-subset_samples(dfclean, SampleType=="duodenal")
    dfcleanduodenal # 23 duodenal samples
    table(sample_data(dfcleanduodenal)$pays)
    
    length(which(sample_data(dfclean)$SampleType=="duodenal")) # we start with 23 samples
    length(which(sample_data(dfclean)$SampleType=="feces")) # we start with 241 samples 
    
#### get information on the overall characteristics of the study group####
    table(sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$stunted, sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$haz, sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$sexe, sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$ageyears, sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$calprotectinelevel, sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$alphaantitrypsinlevel, sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$anemie2, sample_data(dfcleanfeces)$pays)
    
    
### compare data from microscopy to 18S data ####
    table(sample_data(dfcleanfeces)$Parasitology_performed, sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$Ascaris, sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$Trichuris, sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$Giardia, sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$Entamoeba , sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$Chilomastix, sample_data(dfcleanfeces)$pays)
    table(sample_data(dfcleanfeces)$Enterobius, sample_data(dfcleanfeces)$pays)
    
    # now calculate the prevalence by 18S for the same parasites
    prevalence = function(x){ #this only returns prevalence counts 
      x[x >= 5] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dfcleanfeces2<-subset_samples(dfcleanfeces, sample_data(dfcleanfeces)$Parasitology_performed==1)
    dfcleanfeces_Genus= tax_glom(dfcleanfeces2, "rank7")
    

    prev_counts.dffiltered_feces_s <- dfcleanfeces_Genus %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_s) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_feces_possible <- dfcleanfeces_Genus %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_s) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_feces_s/prev_counts.dffiltered_feces_possible)*100
    

    TAX1 = as(tax_table(dfcleanfeces_Genus), "matrix")
    tax_table<-as.data.frame(TAX1)
    tax_table$Genus<-row.names((tax_table(dfcleanfeces_Genus)))
    
    merge.prev_counts.dffiltered_feces= merge(tax_table, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_feces, "Prevalence_parasitosubset_Genus.csv")
    
    
    
    
#### calculate prevalence for ASV in fecal dataset filtered  ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_feces= dfcleanfeces
    
    sample_data(dffiltered_feces)$collapse<-"SampleType"
    sample_data(dffiltered_feces)$collapse<-as.factor(sample_data(dffiltered_feces)$collapse)
    
    prev_counts.dffiltered_feces_s <- dffiltered_feces %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("SampleType") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_s) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_feces_possible <- dffiltered_feces %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("collapse") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_s) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_feces_s/prev_counts.dffiltered_feces_possible)*100
    
    TAX1 = as(tax_table(dfcleanfeces), "matrix")
    tax_table<-as.data.frame(TAX1)
    tax_table$ASV<-row.names((tax_table(dfcleanfeces)))
    
    merge.prev_counts.dffiltered_feces_s_pays= merge(tax_table, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_feces_s_pays, "merge.prev_counts.dffiltered_feces_ASV_withlowabundance.csv")
    
    #### export abundance tables for ASV overall and fuse them with prevalence tables to get a single one  ####

    Afribiota_feces_abundance_ASV <- microbiome::transform(dffiltered_feces, "compositional")
    
    Afribiota_feces_abundance_ASV <- merge_samples(Afribiota_feces_abundance_ASV, "SampleType")
    
    Afribiota_feces_abundance_ASV= Afribiota_feces_abundance_ASV %>%
      otu_table() %>%
      as.data.frame()
    
    Afribiota_feces_abundance_ASV<-t(Afribiota_feces_abundance_ASV)
    
    colnames(Afribiota_feces_abundance_ASV) <- paste("abundance", colnames(Afribiota_feces_abundance_ASV), sep = ".") #add something to distinguish between relative abundance and prevalence

    row.names(merge.prev_counts.dffiltered_feces_s_pays) <- merge.prev_counts.dffiltered_feces_s_pays$Row.names
    
    merge.abunanceprev.feces.ASV.pays= merge(merge.prev_counts.dffiltered_feces_s_pays, Afribiota_feces_abundance_ASV, by="row.names")
    
    test=merge.abunanceprev.feces.ASV.pays[,1]==merge.abunanceprev.feces.ASV.pays[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.feces.ASV.pays <- merge.abunanceprev.feces.ASV.pays[,-1:-2]
    
    write.csv(merge.abunanceprev.feces.ASV.pays, "merge.abunanceprev.feces.ASV_withlowabundance.csv")
    
    
    
    #### calculate prevalence for species for feces dataset for whole dataset merged  ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_feces= tax_glom(dfcleanfeces, "rank8")
    
    sample_data(dffiltered_feces)$collapse<-"collapse"
    
    prev_counts.dffiltered_feces_s <- dffiltered_feces %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("collapse") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_s) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_feces_possible <- dffiltered_feces %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("collapse") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_s) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_feces_s/prev_counts.dffiltered_feces_possible)*100
    
    TAX1 = as(tax_table(dfcleanfeces), "matrix")
    tax_table<-as.data.frame(TAX1)
    
    merge.prev_counts.dffiltered_feces_s_pays= merge(tax_table, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_feces_s_pays, "merge.prev_counts.dffiltered_feces_specieswithlowabundance.csv")
    
    #### export abundance tables for species overall and fuse them with prevalence tables to get a single one with low abundance ####
    
    Afribiota_feces_abundance_species <- microbiome::transform(dffiltered_feces, "compositional")
    
    Afribiota_feces_abundance_species <- merge_samples(Afribiota_feces_abundance_species, "collapse")
    
    Afribiota_feces_abundance_species=  transform_sample_counts(Afribiota_feces_abundance_species, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_species= Afribiota_feces_abundance_species %>%
      otu_table() %>%
      as.data.frame()
    
    Afribiota_feces_abundance_species <- t(Afribiota_feces_abundance_species)
    
    
    colnames(Afribiota_feces_abundance_species) <- paste("abundance", colnames(Afribiota_feces_abundance_species), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    row.names(merge.prev_counts.dffiltered_feces_s_pays) <- merge.prev_counts.dffiltered_feces_s_pays$Row.names
    
    merge.abunanceprev.feces.species.pays= merge(merge.prev_counts.dffiltered_feces_s_pays, Afribiota_feces_abundance_species, by="row.names")
    
    test=merge.abunanceprev.feces.species.pays[,1]==merge.abunanceprev.feces.species.pays[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.feces.species.pays <- merge.abunanceprev.feces.species.pays[,-1:-2]
    
    write.csv(merge.abunanceprev.feces.species.pays, "merge.abunanceprev.feces.specieswithlowabundance.csv")
    
    
    #### calculate prevalence for ASV in fecal dataset by country of origin ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_feces= dfcleanfeces
    
    
    prev_counts.dffiltered_feces_s <- dffiltered_feces %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_s) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_feces_possible <- dffiltered_feces %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_s) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_feces_s/prev_counts.dffiltered_feces_possible)*100
    
    TAX1 = as(tax_table(dfcleanfeces), "matrix")
    tax_table<-as.data.frame(TAX1)
    tax_table.fecespays$ASV<-row.names((tax_table(dfcleanfeces)))
    
    merge.prev_counts.dffiltered_feces_s_pays= merge(tax_table.fecespays, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_feces_s_pays, "merge.prev_counts.dffiltered_feces_ASV_bypays_withlowabundance.csv")
    
    #### export abundance tables for ASV by country of origin  ####
    
    Afribiota_feces_abundance_ASV <- microbiome::transform(dffiltered_feces, "compositional")
    
    Afribiota_feces_abundance_ASV <- merge_samples(Afribiota_feces_abundance_ASV, "pays")
    
    Afribiota_feces_abundance_ASV=  transform_sample_counts(Afribiota_feces_abundance_ASV, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_ASV= Afribiota_feces_abundance_ASV %>%
      otu_table() %>%
      as.data.frame()
    
    Afribiota_feces_abundance_ASV <- t(Afribiota_feces_abundance_ASV)
    
    colnames(Afribiota_feces_abundance_ASV) <- paste("abundance", colnames(Afribiota_feces_abundance_ASV), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    row.names(merge.prev_counts.dffiltered_feces_s_pays) <- merge.prev_counts.dffiltered_feces_s_pays$Row.names
    
    merge.abunanceprev.feces.ASV.pays= merge(merge.prev_counts.dffiltered_feces_s_pays, Afribiota_feces_abundance_ASV, by="row.names")
    
    test=merge.abunanceprev.feces.ASV.pays[,1]==merge.abunanceprev.feces.ASV.pays[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.feces.ASV.pays <- merge.abunanceprev.feces.ASV.pays[,-1:-2]
    
    write.csv(merge.abunanceprev.feces.ASV.pays, "merge.abunanceprev.feces.ASV_byPays_withlowabundance.csv")
    
   #### calculate prevalence for species for filtered dataset in feces for pays of origin  ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_feces= tax_glom(dfcleanfeces, "rank8")
    dffiltered_feces<-subset_samples(dffiltered_feces, pays!="")
    sample_data(dffiltered_feces)$pays<-as.factor(sample_data(dffiltered_feces)$pays)
    levels(sample_data(dffiltered_feces)$pays)
    
    prev_counts.dffiltered_feces_s <- dffiltered_feces %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_s) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_feces_possible <- dffiltered_feces %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_s) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_feces_s/prev_counts.dffiltered_feces_possible)*100
    
    TAX1 = as(tax_table(dffiltered_feces), "matrix")
    tax_table<-as.data.frame(TAX1)
    
    merge.prev_counts.dffiltered_feces_s_pays= cbind(as.matrix(tax_table), as.matrix(test.prev))
    
    
    write.csv(merge.prev_counts.dffiltered_feces_s_pays, "merge.prev_counts.dffiltered_feces_s_payswithlowabundance.csv")
    
    #### export abundance tables for species overall and fuse them with prevalence tables to get a single one ####
    Afribiota_feces_abundance_species <- dffiltered_feces %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("pays") %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      t() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_species) <- paste("abundance", colnames(Afribiota_feces_abundance_species), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    dim(merge.prev_counts.dffiltered_feces_s_pays)
    dim(Afribiota_feces_abundance_species)
    
    merge.abunanceprev.feces.species.pays= cbind(as.matrix(merge.prev_counts.dffiltered_feces_s_pays), as.matrix(Afribiota_feces_abundance_species))

    test=merge.abunanceprev.feces.species.pays[,1]==merge.abunanceprev.feces.species.pays[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.feces.species.pays <- merge.abunanceprev.feces.species.pays[,-1:-2]
    
    write.csv(merge.abunanceprev.feces.species.pays, "merge.abunanceprev.feces.species.payswithlowabundance.csv")
   
    #### calculate prevalence for genus for filtered dataset in feces and pays  ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_feces= tax_glom(dfcleanfeces, "rank7")
    
    
    prev_counts.dffiltered_feces_g <- dffiltered_feces %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_g) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_g), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_feces_possible <- dffiltered_feces %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_g) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_feces_g/prev_counts.dffiltered_feces_possible)*100
    row.names(test.prev)
    dim(test.prev)
    
    TAX1 = as(tax_table(dffiltered_feces), "matrix")
    tax_table<-as.data.frame(TAX1)
    dim(tax_table)
    
    merge.prev_counts.dffiltered_feces_g_pays= cbind(as.matrix(tax_table), as.matrix(test.prev))
    
    write.csv(merge.prev_counts.dffiltered_feces_g_pays, "merge.prev_counts.dffiltered_feces_gwithlowabundance.csv")
    
    #### export abundance tables for genus and fuse them with prevalence tables to get a single one  ####
    
    Afribiota_feces_abundance_genus <- dffiltered_feces %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("pays") %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      t() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_genus) <- paste("abundance", colnames(Afribiota_feces_abundance_genus), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    dim(merge.prev_counts.dffiltered_feces_g_pays)
    dim(Afribiota_feces_abundance_genus)
    
    merge.abunanceprev.feces.genus.pays= cbind(as.matrix(merge.prev_counts.dffiltered_feces_g_pays), as.matrix(Afribiota_feces_abundance_genus))
    
    test=merge.abunanceprev.feces.genus.pays[,1]==merge.abunanceprev.feces.genus.pays[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.feces.genus.pays <- merge.abunanceprev.feces.genus.pays[,-1:-2]
    
    write.csv(merge.abunanceprev.feces.genus.pays, "merge.abunanceprev.feces.genuswithlowabundance.csv")
    View(merge.abunanceprev.feces.genus.pays)
    
    #### calculate prevalence for species for filtered dataset in  for SampleType only for samples from stunted children ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    dfcleanstunted<-subset_samples(dfclean, stunted=="stunted")
    dffiltered= tax_glom(dfcleanstunted, "rank8")
    
    prev_counts.dffiltered_s <- dffiltered %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("SampleType") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_s) <- paste("prevalence", colnames(prev_counts.dffiltered_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_possible <- dffiltered %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("SampleType") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_s) <- paste("prevalence", colnames(prev_counts.dffiltered_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_s/prev_counts.dffiltered_possible)*100
    
    TAX1 = as(tax_table(dffiltered), "matrix")
    tax_table<-as.data.frame(TAX1)
    
    merge.prev_counts.dffiltered_s_SampleType= cbind(as.matrix(tax_table), as.matrix(test.prev))
    
    write.csv(merge.prev_counts.dffiltered_s_SampleType, "merge.prev_counts.dffilteredspecies_SampleType.csv")
    
    
    #### export abundance tables for species and SampleType of origin and fuse them with prevalence tables to get a single one  ####
    Afribiota__abundance_species <- dffiltered %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("SampleType") %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      t() %>%
      as.data.frame()
    
    
    colnames(Afribiota__abundance_species) <- paste("abundance", colnames(Afribiota__abundance_species), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    merge.abunanceprevspecies.SampleType= cbind(as.matrix(merge.prev_counts.dffiltered_s_SampleType), as.matrix(Afribiota__abundance_species))
    
    test=merge.abunanceprevspecies.SampleType[,1]==merge.abunanceprevspecies.SampleType[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprevspecies.SampleType <- merge.abunanceprevspecies.SampleType[,-1:-2]
    
    write.csv(merge.abunanceprevspecies.SampleType, "merge.abunanceprevspecies.SampleType.csv")
    
  
        
#### now make a plot of the prevalence of genus by country of origin ####
    ## first add genus level data of taxa to metadata ##
    dffiltered_feces= tax_glom(dfcleanfeces, "rank7")
    dffiltered_prev<-transform_sample_counts(dffiltered_feces, fun = prevalence) 
    otu_table_g=otu_table(dffiltered_prev) %>%
      unclass() %>%
      as.data.frame()
    
    View(otu_table_g)
    otu_table_g<-t(otu_table_g)
    dim(otu_table_g)
    otu_table_g<-as.data.frame(otu_table_g)
    
    dataframe<-as.data.frame(tax_table(dffiltered_prev)[, 7])
    colnames(dataframe)<-"taxonomy"
    View(dataframe)
    dim(dataframe)
    
    otu_table_g$taxonomy<-dataframe$taxonomy
    otu_table_g<-t(otu_table_g)
    dim(otu_table_g)
    View(otu_table_g)
    colnames(otu_table_g)<-otu_table_g[242,]
    
    sample_data_g=as.data.frame(sample_data(dffiltered_prev))
    Genusdata18S=merge(otu_table_g, sample_data_g, by="row.names", all=TRUE)
    View(Genusdata18S)
    dim(Genusdata18S)
    Genusdata18S = Genusdata18S[-242,] # to take away the line with taxonomy text
    Genusdata18S$Row.names=as.character(Genusdata18S$Row.names)
    write.csv(Genusdata18S, "Genusdata18SfecesDAL.csv")
    View(Genusdata18S)
    
    ## now make the actual graph
    library(binom)
    
    percentage= function(x) {per <- (sum(as.numeric(as.character(x)), na.rm=T)/length(x))*100 
    return (per)}
    
    
    Genusdata18SCAR <-as.data.frame(Genusdata18S[which(Genusdata18S$pays=="RCA"), ]) # to generate a list of all the ones that are from CAR
    Genusdata18SMada <-as.data.frame(Genusdata18S[which(Genusdata18S$pays=="Madagascar"), ]) # to generate a list of all the ones that are from Madagascar
    
    nrow(Genusdata18SCAR) # 91
    nrow(Genusdata18SMada) # 150
    
    Genusdata18SCAR2<-Genusdata18SCAR[, 2:20] 
    Genusdata18SCAR2<- as.matrix(sapply(Genusdata18SCAR2, as.numeric))
    measured<-as.data.frame(colSums(Genusdata18SCAR2))
    colnames(measured)<-"Positive"
    dim(measured)
    
    theoretical<-numeric( length = 19)
    theoretical[1:19] <- 91
    
    test<-as.data.frame(cbind(measured, theoretical))
    
    CICAR <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Genusdata18SMada2<-Genusdata18SMada[, 2:20] 
    Genusdata18SMada2<- as.matrix(sapply(Genusdata18SMada2, as.numeric))
    measured<-colSums(Genusdata18SMada2)
    length(measured)
    theoretical<-numeric( length = 19)
    theoretical[1:19] <- 150
    
    test<-as.data.frame(cbind(measured, theoretical))
    
    CIMada <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Lower=rbind(CICAR$lower*100, CIMada$lower*100)
    row.names(Lower)=c("CAR", "Madagascar")
    colnames(Lower)<-colnames(Genusdata18SCAR2)
    
    Upper=rbind(CICAR$upper*100, CIMada$upper*100)
    row.names(Upper)=c("CAR", "Madagascar")
    colnames(Upper)<-colnames(Genusdata18SMada2)
    
    Genusdata18SCARcollapsed=colwise(percentage)(Genusdata18SCAR)
    Genusdata18SMadacollapsed=colwise(percentage)(Genusdata18SMada)
    
    PercentageGenus=rbind(Genusdata18SCARcollapsed, Genusdata18SMadacollapsed)
    
    row.names(PercentageGenus)=c("CAR", "Madagascar")
    View(PercentageGenus)
    ncol(PercentageGenus)
    
    PercentageGenus <- PercentageGenus[,-(21:850), drop=FALSE] # kick-out the metadata
    PercentageGenus <- PercentageGenus[,-(1),drop=FALSE] # kick-out the metadata
    
    write.csv(PercentageGenus, "PercentageGenuspaysfeces.csv")
    
    PercentageGenus2=t(PercentageGenus)
    Upper2=t(Upper)
    Lower2=t(Lower)
    
    PercentageGenus2toplot=melt(PercentageGenus2, value.name="PercentageGenus2", varnames=c("Taxon", "pays"))
    Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "pays"))
    Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "pays"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Uppertoplot, by=c("Taxon", "pays"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Lowertoplot, by=c("Taxon", "pays"))
    
    PercentageGenus2toplot<-PercentageGenus2toplot[order(-PercentageGenus2),]
    
    pdf("generaaccordingtopaysfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 5)
    ggplot(data=PercentageGenus2toplot, aes(x=reorder(Taxon, -PercentageGenus2), y=PercentageGenus2,  color=pays, fill=pays)) +
      geom_errorbar(data=PercentageGenus2toplot, mapping=aes(x=reorder(Taxon, -PercentageGenus2), ymin=Lower, ymax=Upper),
                    size=0.3,
                    width=0.9,  
                    color="black", 
                    position=position_dodge(.9)) +
      geom_bar(stat="identity", position="dodge") +
      theme_bw() +
      theme(strip.background=element_rect(fill="white")) +
      theme(strip.text=element_text(color="black")) +
      scale_color_manual(values=c("blue", "red", "green")) +
      scale_fill_manual(values=c("blue", "red", "green")) +
      theme(axis.text.x=element_text(angle=-45, vjust=1, hjust=0)) +
      labs(title="Percentage of 18S genera in fecal samples according to country of origin")
    dev.off()
    
    
    # now filter out the low abundance taxa
    PercentageGenus2 <- PercentageGenus[, (PercentageGenus[1, ] > 10) | (PercentageGenus[2, ] > 10)] # now kick-out taxa with less than 10% prevalence in any of the two countries
    
    PercentageGenus2=t(PercentageGenus2)
    
    PercentageGenus2toplot=melt(PercentageGenus2, value.name="PercentageGenus2", varnames=c("Taxon", "pays"))
    Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "pays"))
    Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "pays"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Uppertoplot, by=c("Taxon", "pays"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Lowertoplot, by=c("Taxon", "pays"))
    
    PercentageGenus2toplot<-PercentageGenus2toplot[order(-PercentageGenus2),]
    
    pdf("generaaccordingtopaysfecesfiltered.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 5)
    ggplot(data=PercentageGenus2toplot, aes(x=reorder(Taxon, -PercentageGenus2), y=PercentageGenus2,  color=pays, fill=pays)) +
      geom_errorbar(data=PercentageGenus2toplot, mapping=aes(x=reorder(Taxon, -PercentageGenus2), ymin=Lower, ymax=Upper),
                    size=0.3,
                    width=0.9,  
                    color="black", 
                    position=position_dodge(.9)) +
      geom_bar(stat="identity", position="dodge") +
      theme_bw() +
      theme(strip.background=element_rect(fill="white")) +
      theme(strip.text=element_text(color="black")) +
      scale_color_manual(values=c("blue", "red", "green")) +
      scale_fill_manual(values=c("blue", "red", "green")) +
      theme(axis.text.x=element_text(angle=-45, vjust=1, hjust=0)) +
      labs(title="Percentage of 18S genera in fecal samples according to country of origin")
    dev.off()  
    
    
    
    #### now make a plot of the prevalence of genus by stunting ####
    ## first add genus level data of taxa to metadata ##
    dffiltered_feces= tax_glom(dfcleanfeces, "rank7")
    dffiltered_prev<-transform_sample_counts(dffiltered_feces, fun = prevalence) 
    otu_table_g=otu_table(dffiltered_prev) %>%
      unclass() %>%
      as.data.frame()
    
    View(otu_table_g)
    dim(otu_table_g)
    otu_table_g<-as.data.frame(t(otu_table_g))
    
    dataframe<-as.data.frame(tax_table(dffiltered_prev)[, 7])
    colnames(dataframe)<-"taxonomy"
    View(dataframe)
    
    dataframe$taxonomy <- sub("uncultured", "Unassigned_Saccharomycetales_2", dataframe$taxonomy) # replace by a better name
    colnames(dataframe)<-"taxonomy"
    
    otu_table_g$taxonomy<-dataframe$taxonomy
    otu_table_g<-t(otu_table_g)
    dim(otu_table_g)
    colnames(otu_table_g)<-otu_table_g[242,]
    
    sample_data_g=as.data.frame(sample_data(dffiltered_prev))
    Genusdata18S=merge(otu_table_g, sample_data_g, by="row.names", all=TRUE)
    View(Genusdata18S)
    dim(Genusdata18S)
    Genusdata18S = Genusdata18S[-242,] # to take away the line with taxonomy text
    
    Genusdata18S$Row.names=as.character(Genusdata18S$Row.names)
    write.csv(Genusdata18S, "Genusdata18Sfeces.csv")
    View(Genusdata18S)
    
    ## now make the actual graph
    library(binom)
    
    percentage= function(x) {per <- (sum(as.numeric(as.character(x)), na.rm=T)/length(x))*100 
    return (per)}
    
    Genusdata18Sstunted <-as.data.frame(Genusdata18S[which(Genusdata18S$stunted=="stunted"), ]) # to generate a list of all the ones that are from CAR
    Genusdata18Snonstunted <-as.data.frame(Genusdata18S[which(Genusdata18S$stunted=="non-stunted"), ]) # to generate a list of all the ones that are from Madagascar
    
    nrow(Genusdata18Sstunted) # 126
    nrow(Genusdata18Snonstunted) # 115
    colnames(Genusdata18Sstunted)
    
    Genusdata18Sstunted2<-Genusdata18Sstunted[, 2:20] 
    Genusdata18Sstunted2<- as.matrix(sapply(Genusdata18Sstunted2, as.numeric))
    measured<-as.data.frame(colSums(Genusdata18Sstunted2))
    colnames(measured)<-"Positive"
    dim(measured)
    
    theoretical<-numeric( length = 19)
    theoretical[1:19] <- 126
    
    test<-as.data.frame(cbind(measured, theoretical))
    
    CIstunted <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Genusdata18Snonstunted2<-Genusdata18Snonstunted[, 2:20] 
    Genusdata18Snonstunted2<- as.matrix(sapply(Genusdata18Snonstunted2, as.numeric))
    measured<-colSums(Genusdata18Snonstunted2)
    length(measured)
    theoretical<-numeric( length = 19)
    theoretical[1:19] <- 115
    
    test<-as.data.frame(cbind(measured, theoretical))
    
    CInonstunted <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Lower=rbind(CIstunted$lower*100, CInonstunted$lower*100)
    row.names(Lower)=c("stunted", "non-stunted")
    colnames(Lower)<-colnames(Genusdata18Sstunted2)
    
    Upper=rbind(CIstunted$upper*100, CInonstunted$upper*100)
    row.names(Upper)=c("stunted", "non-stunted")
    colnames(Upper)<-colnames(Genusdata18Snonstunted2)
    
    Genusdata18Sstuntedcollapsed=colwise(percentage)(Genusdata18Sstunted)
    Genusdata18Snonstuntedcollapsed=colwise(percentage)(Genusdata18Snonstunted)
    
    PercentageGenus=rbind(Genusdata18Sstuntedcollapsed, Genusdata18Snonstuntedcollapsed)
    
    row.names(PercentageGenus)=c("stunted", "non-stunted")
    View(PercentageGenus)
    ncol(PercentageGenus)
    
    PercentageGenus <- PercentageGenus[,-(21:855), drop=FALSE] # kick-out the metadata
    PercentageGenus <- PercentageGenus[,-(1),drop=FALSE] # kick-out the metadata
    
    write.csv(PercentageGenus, "PercentageGenusstuntedfeces.csv")
    
    PercentageGenus2=t(PercentageGenus)
    Upper2=t(Upper)
    Lower2=t(Lower)
    
    PercentageGenus2toplot=melt(PercentageGenus2, value.name="PercentageGenus2", varnames=c("Taxon", "stunted"))
    Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "stunted"))
    Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "stunted"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Uppertoplot, by=c("Taxon", "stunted"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Lowertoplot, by=c("Taxon", "stunted"))
    
    PercentageGenus2toplot<-PercentageGenus2toplot[order(-PercentageGenus2),]
    
    pdf("generaaccordingtostuntedfecesDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 5)
    ggplot(data=PercentageGenus2toplot, aes(x=reorder(Taxon, -PercentageGenus2), y=PercentageGenus2,  color=stunted, fill=stunted)) +
      geom_errorbar(data=PercentageGenus2toplot, mapping=aes(x=reorder(Taxon, -PercentageGenus2), ymin=Lower, ymax=Upper),
                    size=0.3,
                    width=0.9,  
                    color="black", 
                    position=position_dodge(.9)) +
      geom_bar(stat="identity", position="dodge") +
      theme_bw() +
      theme(strip.background=element_rect(fill="white")) +
      theme(strip.text=element_text(color="black")) +
      scale_color_manual(values=c("blue", "red", "green")) +
      scale_fill_manual(values=c("blue", "red", "green")) +
      theme(axis.text.x=element_text(angle=-45, vjust=1, hjust=0)) +
      labs(title="Percentage of 18S genera in fecal samples according to stunting status")
    dev.off()
    
    
    # now filter out the low abundance taxa
    PercentageGenus2 <- PercentageGenus[, (PercentageGenus[1, ] > 10) | (PercentageGenus[2, ] > 10)] # now kick-out taxa with less than 10% prevalence in any of the two countries
    
    PercentageGenus2=t(PercentageGenus2)
    
    PercentageGenus2toplot=melt(PercentageGenus2, value.name="PercentageGenus2", varnames=c("Taxon", "stunted"))
    Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "stunted"))
    Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "stunted"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Uppertoplot, by=c("Taxon", "stunted"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Lowertoplot, by=c("Taxon", "stunted"))
    
    PercentageGenus2toplot<-PercentageGenus2toplot[order(-PercentageGenus2),]
    
    pdf("generaaccordingtostuntedfecesfilteredDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 5)
    ggplot(data=PercentageGenus2toplot, aes(x=reorder(Taxon, -PercentageGenus2), y=PercentageGenus2,  color=stunted, fill=stunted)) +
      geom_errorbar(data=PercentageGenus2toplot, mapping=aes(x=reorder(Taxon, -PercentageGenus2), ymin=Lower, ymax=Upper),
                    size=0.3,
                    width=0.9,  
                    color="black", 
                    position=position_dodge(.9)) +
      geom_bar(stat="identity", position="dodge") +
      theme_bw() +
      theme(strip.background=element_rect(fill="white")) +
      theme(strip.text=element_text(color="black")) +
      scale_color_manual(values=c("blue", "red", "green")) +
      scale_fill_manual(values=c("blue", "red", "green")) +
      theme(axis.text.x=element_text(angle=-45, vjust=1, hjust=0)) +
      labs(title="Percentage of 18S genera in fecal samples according to stunting status")
    dev.off()  
    
    
    #### which clinical variables are associated with rel abundance of genus? ####
    
    dfcleanfecesgenus<-tax_glom(dfcleanfeces, "rank7")
    
    eukaryome.rel.genus <- microbiome::transform(dfcleanfecesgenus, "compositional")
    eukaryome.rel.genus<-subset_samples(eukaryome.rel.genus, pays!="")
    
    eukaryome.rel.genus=subset_samples(eukaryome.rel.genus, sample_sums(eukaryome.rel.genus)!="0")
    eukaryome.rel.genus = filter_taxa(eukaryome.rel.genus, function(x) sum(x) > 0, TRUE)
    
    df_wilcox <- as.matrix((otu_table(eukaryome.rel.genus))) #take rel abund and Wilcoxon rank-sum
    meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.genus)) #take metadata
    
    dim(df_wilcox)
    dim(meta_wilcox)
    
    ASVnames<-row.names(tax_table(eukaryome.rel.genus))
    
    # assess for association with country of origin
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$pays)$p.value)
    Genusnames<-colnames(df_wilcox)
    p.res = data.frame(ASVnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.genus)
    p.res = cbind(as.matrix(p.res),as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.pays.GenusWilcoxDAL.csv")
    View(p.res) 
    
    
    # assess for association with stunting
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$stunted)$p.value)
    Genusnames<-colnames(df_wilcox)
    p.res = data.frame(ASVnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.genus)
    p.res = cbind(as.matrix(p.res),as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.stunted.GenusWilcoxDAL.csv")
    View(p.res) 
    
    # assess for age (use younger and older than three years)
    meta_wilcox$ageyears2<-meta_wilcox$age
    meta_wilcox$ageyears2[meta_wilcox$age>=36]<-5
    meta_wilcox$ageyears2[meta_wilcox$age<36]<-2
    meta_wilcox$ageyears2
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$ageyears2)$p.value)
    Genusnames<-colnames(df_wilcox)
    p.res = data.frame(ASVnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.genus)
    p.res = cbind(as.matrix(p.res),as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.ageyears.GenusWilcoxDAL.csv")
    View(p.res) 
    
    # assess for calprotectin level
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$calprotectinelevel)$p.value)
    Genusnames<-colnames(df_wilcox)
    p.res = data.frame(ASVnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.genus)
    p.res = cbind(as.matrix(p.res),as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.calpro.GenusWilcoxDAL.csv")
    View(p.res) 
    
    # assess for AAT level
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$alphaantitrypsinlevel)$p.value)
    Genusnames<-colnames(df_wilcox)
    p.res = data.frame(ASVnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.genus)
    p.res = cbind(as.matrix(p.res),as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.AAT.GenusWilcoxDAL.csv")
    View(p.res) 
    
    # assess for anemia level
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$anemie2)$p.value)
    Genusnames<-colnames(df_wilcox)
    p.res = data.frame(ASVnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.genus)
    p.res = cbind(as.matrix(p.res),as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.anemie.GenusWilcoxDAL.csv")
    View(p.res) 
    
    #### which clinical variables are associated with presence/ absence of genera? ####
    dfcleanfecesgenus<-tax_glom(dfcleanfeces, "rank7")
    
    fecesprevalence <- dfcleanfecesgenus%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence))) #take rel abund and Wilcoxon rank-sum
    meta_chi2 <- as.data.frame(sample_data(fecesprevalence)) #take metadata
    taxa_chi2<- tax_table(fecesprevalence) #take taxtable
    
    meta_chi2$stunted<-as.factor(meta_chi2$stunted)
    meta_chi2$pays<-as.factor(meta_chi2$pays)
    meta_chi2$haz<-as.factor(meta_chi2$haz)
    meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
    
    # on country of origin
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value)
    rank7names<-taxa_chi2[, 7]
    chi2results = data.frame(rank7names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    
    write.csv(chi2results_sig,"Chi2resultsFeces.GenuspaysDAL.csv")
    
    # on stunting status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
    rank7names<-taxa_chi2[, 7]
    chi2results = data.frame(rank7names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    
    write.csv(chi2results_sig,"Chi2resultsFeces.GenusstuntedDAL.csv")
    
    # on age 
    meta_chi2$ageyears2<-meta_chi2$age
    meta_chi2$ageyears2[meta_chi2$age>=36]<-5
    meta_chi2$ageyears2[meta_chi2$age<36]<-2
    meta_chi2$ageyears2
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears2), simulate.p.value = TRUE)$p.value)
    rank7names<-taxa_chi2[, 7]
    chi2results = data.frame(rank7names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    
       
    # on calprotectine status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevels), simulate.p.value = TRUE)$p.value)
    rank7names<-colnames(df_chi2)
    chi2results = data.frame(rank7names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    
    write.csv(chi2results_sig,"Chi2resultsFeces.GenuscalproDAL.csv")
                            
                            
    # on AAT status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value)
    rank7names<-taxa_chi2[, 7]
    chi2results = data.frame(rank7names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
                            
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
                            
    write.csv(chi2results_sig,"Chi2resultsFeces.GenusaatDAL.csv")
                            
    # on anemia status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
    chi2results = data.frame(rank7names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
                            
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.GenusanemiaDAL.csv")
                            
                            
                            
                            
    #### which clinical variables are associated with presence/ absence of species? ####
    dfcleanfecesspecies<-tax_glom(dfcleanfeces, "rank8")
    
    fecesprevalence <- dfcleanfecesspecies%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence))) #take rel abund and Wilcoxon rank-sum
    meta_chi2 <- as.data.frame(sample_data(fecesprevalence)) #take metadata
    taxa_chi2<- tax_table(fecesprevalence) #take taxtable
    
    meta_chi2$stunted<-as.factor(meta_chi2$stunted)
    meta_chi2$pays<-as.factor(meta_chi2$pays)
    meta_chi2$haz<-as.factor(meta_chi2$haz)
    meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
    meta_chi2$alphaantitrypsinlevel<-as.factor(meta_chi2$alphaantitrypsinlevel)
    meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
    
    # on country of origin
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value)
    rank8names<-taxa_chi2[, 8]
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.speciespaysDAL.csv")
    
    # on stunting status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.speciesstuntedDAL.csv")
    
    # on age 
    meta_chi2$ageyears2<-meta_chi2$age
    meta_chi2$ageyears2[meta_chi2$age>=36]<-5
    meta_chi2$ageyears2[meta_chi2$age<36]<-2
    meta_chi2$ageyears2
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears2), simulate.p.value = TRUE)$p.value)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.speciesageyears2DAL.csv")
    
    # on calprotectine status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.speciescalproDAL.csv")
    
    # on AAT status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.speciesaatDAL.csv")
    
    # on anemia status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.speciesanemiaDAL.csv")
    
   ####  now make a plot of the prevalence of different Blasto subtypes by country of origin ####
    ## first add Species level data of taxa to metadata ##
    dffiltered_feces= tax_glom(dfcleanfeces, "rank8")
    dffiltered_prev<-transform_sample_counts(dffiltered_feces, fun = prevalence) 
    otu_table_g=otu_table(dffiltered_prev) %>%
      unclass() %>%
      as.data.frame()
    
    View(otu_table_g)
    otu_table_g<-t(otu_table_g)
    dim(otu_table_g)
    otu_table_g<-as.data.frame(otu_table_g)
    
    dataframe<-as.data.frame(tax_table(dffiltered_prev)[, 8])
    colnames(dataframe)<-"taxonomy"
    View(dataframe)
    
    otu_table_g$taxonomy<-dataframe$taxonomy
    otu_table_g<-t(otu_table_g)
    dim(otu_table_g)
    colnames(otu_table_g)<-otu_table_g[242,]
    
    sample_data_g=as.data.frame(sample_data(dffiltered_prev))
    Speciesdata18S=merge(otu_table_g, sample_data_g, by="row.names", all=TRUE)
    View(Speciesdata18S)
    dim(Speciesdata18S)
    Speciesdata18S = Speciesdata18S[-242,] # to take away the line with taxonomy text
    
    Speciesdata18S$Row.names=as.character(Speciesdata18S$Row.names)
    write.csv(Speciesdata18S, "Speciesdata18Sfeces.csv")
    colnames(Speciesdata18S)
    
    # filter out the Blasto columns
    Speciesdata18S_Blasto<-Speciesdata18S[, -4]
    
    
    ## now make the actual graph
    library(binom)
    
    percentage= function(x) {per <- (sum(as.numeric(as.character(x)), na.rm=T)/length(x))*100 
    return (per)}
    
    
    Speciesdata18S_BlastoCAR <-as.data.frame(Speciesdata18S_Blasto[which(Speciesdata18S_Blasto$pays=="RCA"), ]) # to generate a list of all the ones that are from CAR
    Speciesdata18S_BlastoMada <-as.data.frame(Speciesdata18S_Blasto[which(Speciesdata18S_Blasto$pays=="Madagascar"), ]) # to generate a list of all the ones that are from Madagascar
    
    nrow(Speciesdata18S_BlastoCAR) # 91
    nrow(Speciesdata18S_BlastoMada) # 150
    
    Speciesdata18S_BlastoCAR2<-Speciesdata18S_BlastoCAR[, 2:4] 
    Speciesdata18S_BlastoCAR2<- as.matrix(sapply(Speciesdata18S_BlastoCAR2, as.numeric))
       measured<-as.data.frame(colSums(Speciesdata18S_BlastoCAR2))
    colnames(measured)<-"Positive"
    dim(measured)
    
    theoretical<-numeric( length = 3)
    theoretical[1:3] <- 91
    
    test<-as.data.frame(cbind(measured, theoretical))
    
    CICAR <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Speciesdata18S_BlastoMada2<-Speciesdata18S_BlastoMada[, 2:4] 
    Speciesdata18S_BlastoMada2<- as.matrix(sapply(Speciesdata18S_BlastoMada2, as.numeric))
    measured<-colSums(Speciesdata18S_BlastoMada2)
    length(measured)
    theoretical<-numeric( length = 3)
    theoretical[1:3] <- 150
    
    test<-as.data.frame(cbind(measured, theoretical))
    
    CIMada <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Lower=rbind(CICAR$lower*100, CIMada$lower*100)
    row.names(Lower)=c("CAR", "Madagascar")
    colnames(Lower)<-colnames(Speciesdata18S_BlastoCAR2)
    
    Upper=rbind(CICAR$upper*100, CIMada$upper*100)
    row.names(Upper)=c("CAR", "Madagascar")
    colnames(Upper)<-colnames(Speciesdata18S_BlastoMada2)
    
    Speciesdata18S_BlastoCARcollapsed=colwise(percentage)(Speciesdata18S_BlastoCAR)
    Speciesdata18S_BlastoMadacollapsed=colwise(percentage)(Speciesdata18S_BlastoMada)
    
    PercentageBlasto=rbind(Speciesdata18S_BlastoCARcollapsed, Speciesdata18S_BlastoMadacollapsed)
    
    row.names(PercentageBlasto)=c("CAR", "Madagascar")
    View(PercentageBlasto)
    dim(PercentageBlasto)
    
    PercentageBlasto <- PercentageBlasto[,-(5:56), drop=FALSE] # kick-out the metadata
    PercentageBlasto <- PercentageBlasto[,-(1),drop=FALSE] # kick-out the metadata
    
    write.csv(PercentageBlasto, "PercentageBlastosubtypespaysfeces.csv")
    
    PercentageBlasto2=t(PercentageBlasto)
    Upper2=t(Upper)
    Lower2=t(Lower)
    
    PercentageBlasto2toplot=melt(PercentageBlasto2, value.name="PercentageBlasto2", varnames=c("Taxon", "pays"))
    Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "pays"))
    Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "pays"))
    PercentageBlasto2toplot<-merge(PercentageBlasto2toplot, Uppertoplot, by=c("Taxon", "pays"))
    PercentageBlasto2toplot<-merge(PercentageBlasto2toplot, Lowertoplot, by=c("Taxon", "pays"))
    
    PercentageBlasto2toplot<-PercentageBlasto2toplot[order(-PercentageBlasto2),]
    
    pdf("Blastosubtypesaccordingtopaysfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 5)
    ggplot(data=PercentageBlasto2toplot, aes(x=reorder(Taxon, -PercentageBlasto2), y=PercentageBlasto2,  color=pays, fill=pays)) +
      geom_errorbar(data=PercentageBlasto2toplot, mapping=aes(x=reorder(Taxon, -PercentageBlasto2), ymin=Lower, ymax=Upper),
                    size=0.3,
                    width=0.9,  
                    color="black", 
                    position=position_dodge(.9)) +
      geom_bar(stat="identity", position="dodge") +
      theme_bw() +
      theme(strip.background=element_rect(fill="white")) +
      theme(strip.text=element_text(color="black")) +
      scale_color_manual(values=c("blue", "red", "green")) +
      scale_fill_manual(values=c("blue", "red", "green")) +
      theme(axis.text.x=element_text(angle=-45, vjust=1, hjust=0)) +
      labs(title="Percentage of different Blastocytis subtypes in fecal samples according to country of origin")
    dev.off()
    
    
    #### find out which samples have Blastocystis and if they are having one or several subtypes cutoff 5 seqs ####
    
    # first subset samples to have only info on Blasto
    dfcleanfeces
    dfcleanfeces_blasto=subset_taxa(dfcleanfeces, rank7=="Blastocystis")
    dfcleanfeces_blasto_s=tax_glom(dfcleanfeces_blasto, "rank8") # this gives us a otu_table with how many subtypes per sample
    
    prevalenceyesno = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function, but need to be expressed as a percentage only for samples belonging to a particular species
      x[x >= 5] <- 1
      x[x >= 1] <- 1
      return(x)}
    
    
    dfcleanfeces_blasto_s_pres <- dfcleanfeces_blasto_s %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalenceyesno)
    
    
    # transform real counts into table with header
    otu_table=as.data.frame(t(otu_table(dfcleanfeces_blasto_s)))
    dim(otu_table)
    View(otu_table)
    
    tax_table=as.data.frame(t(tax_table(dfcleanfeces_blasto_s)))
    tax_table=tax_table[-(1:7), ]
    tax_table=tax_table[-(2:8), ]
    dim(tax_table)
    View(tax_table)
    
    tax_table<-t(tax_table)
    
    row.names(otu_table)<-tax_table[,1]
    View(otu_table)
    
    otu_table$sumrow <- rowSums(otu_table, na.rm = TRUE)
    
    write.csv(otu_table, "BlastorealcountsAfribiota.csv")
    
    # transform real counts into table with header at least five sequences
    
    otu_table=as.data.frame((otu_table(dfcleanfeces_blasto_s_pres)))
    dim(otu_table)
    View(otu_table)
    colnames(otu_table)<-tax_table[, 1]
    View(otu_table)
    
    otu_table$sumrow <- rowSums(otu_table, na.rm = TRUE)
    
    write.csv(otu_table, "BlastoatleastfiveAfribiota.csv")
    
    # summarize your Blasto positive samples
    dim(otu_table) #241
    
    sum(otu_table$Blastocystis_ST1 == 1) # 111 samples have Blasto ST1
    sum(otu_table$Blastocystis_ST2 == 1) # 79 samples have Blasto ST2
    sum(otu_table$Blastocystis_ST3 == 1) # 128 samples have Blasto ST3
    sum(otu_table$sumrow >=1) #189 samples have Blasto count of at least 1
    100/241*189  ## hence 78% are positive
    100/241*111 # 46%
    100/241*79 # 33%
    100/241*128 # 53%
    sum(otu_table$sumrow ==0) # 52 samples have no Blastocystis
    
    # this would also work for all numeric columns in a single command
    res <- unlist(lapply(otu_table, function(x) if(is.numeric(x)) sum(x, na.rm=T)))
    res
    
    # summarize your Blasto positive samples which have only one Blasto
    sum(otu_table$Blastocystis_ST1 == 0 & otu_table$Blastocystis_ST2 == 0 & otu_table$Blastocystis_ST3 == 0) # 0 samples have only Blasto hominis
    
    sum(otu_table$Blastocystis_ST1 == 1 & otu_table$Blastocystis_ST2 == 0 & otu_table$Blastocystis_ST3 == 0) 
    listST1=rownames(otu_table)[otu_table$Blastocystis_ST1 == 1 & otu_table$Blastocystis_ST2 == 0 & otu_table$Blastocystis_ST3 == 0] 
    listST1 # 28 have only  Blastocystis ST1
    
    sum(otu_table$Blastocystis_ST1 == 0 & otu_table$Blastocystis_ST2 == 1 & otu_table$Blastocystis_ST3 == 0) 
    listST2=rownames(otu_table)[otu_table$Blastocystis_ST1 == 0 & otu_table$Blastocystis_ST2 == 1 & otu_table$Blastocystis_ST3 == 0] 
    listST2 # 23 sample have only Blasto ST2
    
    sum(otu_table$Blastocystis_ST1 == 0 & otu_table$Blastocystis_ST2 == 0 & otu_table$Blastocystis_ST3 == 1) 
    listST3=rownames(otu_table)[otu_table$Blastocystis_ST1 == 0 & otu_table$Blastocystis_ST2 == 0 & otu_table$Blastocystis_ST3 == 1] 
    listST3 # 42 sample have only Blasto ST3
    
    # summarize your Blasto positive samples which have several Blasto, first make a column saying how many different Blasto species
    sum(otu_table$sumrow == 1) # 93/241 samples have only one Blasto 
    100/241*93 # 38%
    100/189*93 #49%
    
    sum(otu_table$sumrow == 2) # 63/241 samples have two Blasto
    100/241*63 # 26%
    
    sum(otu_table$sumrow == 3) # 33/241 samples have tree Blasto
    100/241*33 # 14%
    
    sum(otu_table$sumrow >1) # 96/241 samples have several Blasto
    100/241*96 # 40%
    100/189*96 #51%
    
    
    ####  now make a plot of the prevalence of Blastocystis Cluster by country of origin ####
    ## first add Cluster level data of taxa to metadata #
    dffiltered_feces_cluster<-dfcleanfeces
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -12] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    percentage= function(x) {per <- (sum(as.numeric(as.character(x)), na.rm=T)/length(x))*100 
    return (per)}
    
    dffiltered_prev<-transform_sample_counts(dffiltered_feces_cluster, fun = prevalence) 
    otu_table_g=otu_table(dffiltered_prev) %>%
      unclass() %>%
      as.data.frame()
    
    View(tax_table(dffiltered_feces_cluster))
    
    View(otu_table_g)
    otu_table_g<-t(otu_table_g)
    dim(otu_table_g)
    otu_table_g<-as.data.frame(otu_table_g)
    
    dataframe<-as.data.frame(tax_table(dffiltered_prev)[, 12])
    colnames(dataframe)<-"taxonomy"
    View(dataframe)
    dim(dataframe)
    
    otu_table_g$taxonomy<-dataframe$taxonomy
    otu_table_g<-t(otu_table_g)
    dim(otu_table_g)
    View(otu_table_g)
    colnames(otu_table_g)<-otu_table_g[242,]
    
    sample_data_g=as.data.frame(sample_data(dffiltered_prev))
    Clusterdata18S=merge(otu_table_g, sample_data_g, by="row.names", all=TRUE)
    View(Clusterdata18S)
    dim(Clusterdata18S)
    Clusterdata18S = Clusterdata18S[-242,] # to take away the line with taxonomy text
    
    Clusterdata18S$Row.names=as.character(Clusterdata18S$Row.names)
    write.csv(Clusterdata18S, "Clusterdata18SfecesDAL.csv")
    View(Clusterdata18S)
    
    ## now make the actual graph
    library(binom)
    
    percentage= function(x) {per <- (sum(as.numeric(as.character(x)), na.rm=T)/length(x))*100 
    return (per)}
    
    
    Clusterdata18SCAR <-as.data.frame(Clusterdata18S[which(Clusterdata18S$pays=="RCA"), ]) # to generate a list of all the ones that are from CAR
    Clusterdata18SMada <-as.data.frame(Clusterdata18S[which(Clusterdata18S$pays=="Madagascar"), ]) # to generate a list of all the ones that are from Madagascar
    
    nrow(Clusterdata18SCAR) # 91
    nrow(Clusterdata18SMada) # 150
    
    Clusterdata18SCAR2<-Clusterdata18SCAR[, 2:17] 
    Clusterdata18SCAR2<- as.matrix(sapply(Clusterdata18SCAR2, as.numeric))
    measured<-as.data.frame(colSums(Clusterdata18SCAR2))
    colnames(measured)<-"Positive"
    dim(measured)
    
    theoretical<-numeric( length = 16)
    theoretical[1:16] <- 91
    
    test<-as.data.frame(cbind(measured, theoretical))
    
    CICAR <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Clusterdata18SMada2<-Clusterdata18SMada[, 2:17] 
    Clusterdata18SMada2<- as.matrix(sapply(Clusterdata18SMada2, as.numeric))
    measured<-colSums(Clusterdata18SMada2)
    length(measured)
    theoretical<-numeric( length = 16)
    theoretical[1:16] <- 150
    
    test<-as.data.frame(cbind(measured, theoretical))
    
    CIMada <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Lower=rbind(CICAR$lower*100, CIMada$lower*100)
    row.names(Lower)=c("CAR", "Madagascar")
    colnames(Lower)<-colnames(Clusterdata18SCAR2)
    
    Upper=rbind(CICAR$upper*100, CIMada$upper*100)
    row.names(Upper)=c("CAR", "Madagascar")
    colnames(Upper)<-colnames(Clusterdata18SMada2)
    
    Clusterdata18SCARcollapsed=colwise(percentage)(Clusterdata18SCAR)
    Clusterdata18SMadacollapsed=colwise(percentage)(Clusterdata18SMada)
    
    PercentageCluster=rbind(Clusterdata18SCARcollapsed, Clusterdata18SMadacollapsed)
    
    row.names(PercentageCluster)=c("CAR", "Madagascar")
    View(PercentageCluster)
    ncol(PercentageCluster)
    
    PercentageCluster <- PercentageCluster[,-(18:48), drop=FALSE] # kick-out the metadata
    PercentageCluster <- PercentageCluster[,-(1),drop=FALSE] # kick-out the metadata
    
    write.csv(PercentageCluster, "PercentageClusterpaysfeces.csv")
    
    PercentageCluster2=t(PercentageCluster)
    Upper2=t(Upper)
    Lower2=t(Lower)
    
    PercentageCluster2toplot=melt(PercentageCluster2, value.name="PercentageCluster2", varnames=c("Taxon", "pays"))
    Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "pays"))
    Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "pays"))
    PercentageCluster2toplot<-merge(PercentageCluster2toplot, Uppertoplot, by=c("Taxon", "pays"))
    PercentageCluster2toplot<-merge(PercentageCluster2toplot, Lowertoplot, by=c("Taxon", "pays"))
    
    PercentageCluster2toplot<-PercentageCluster2toplot[order(-PercentageCluster2),]
    
    pdf("Blastoclusteraccordingtopaysfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 5)
    ggplot(data=PercentageCluster2toplot, aes(x=reorder(Taxon, -PercentageCluster2), y=PercentageCluster2,  color=pays, fill=pays)) +
      geom_errorbar(data=PercentageCluster2toplot, mapping=aes(x=reorder(Taxon, -PercentageCluster2), ymin=Lower, ymax=Upper),
                    size=0.3,
                    width=0.9,  
                    color="black", 
                    position=position_dodge(.9)) +
      geom_bar(stat="identity", position="dodge") +
      theme_bw() +
      theme(strip.background=element_rect(fill="white")) +
      theme(strip.text=element_text(color="black")) +
      scale_color_manual(values=c("blue", "red", "green")) +
      scale_fill_manual(values=c("blue", "red", "green")) +
      theme(axis.text.x=element_text(angle=-45, vjust=1, hjust=0)) +
      labs(title="Percentage of Blastocystis cluster in fecal samples according to country of origin")
    dev.off()
    pro

#### which clinical variables are associated with rel. abundance of Blasto at cluster level? ####
    dffiltered_feces_cluster<-dfcleanfeces
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -12] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    eukaryome.rel.cluster <- microbiome::transform(dffiltered_feces_cluster, "compositional")
    eukaryome.rel.cluster<-subset_samples(eukaryome.rel.cluster, pays!="")
    
    eukaryome.rel.cluster=subset_samples(eukaryome.rel.cluster, sample_sums(eukaryome.rel.cluster)!="0")
    eukaryome.rel.cluster = filter_taxa(eukaryome.rel.cluster, function(x) sum(x) > 0, TRUE)
    
    df_wilcox <- as.matrix((otu_table(eukaryome.rel.cluster))) #take rel abund and Wilcoxon rank-sum
    meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.cluster)) #take metadata
    
    dim(df_wilcox)
    dim(meta_wilcox)
    
    # assess for association with country of origin
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$pays)$p.value)
    clusternames<-colnames(df_wilcox)
    p.res = data.frame(clusternames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.cluster)
    dim(p.res)
    dim(Tax_corr)
    
    p.res = cbind(as.matrix(p.res), as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.pays.clusterWilcoxDAL.csv")
    View(p.res) 
    
    
    # assess for association with stunting
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$stunted)$p.value)
    clusternames<-colnames(df_wilcox)
    p.res = data.frame(clusternames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.cluster)
    p.res = cbind(as.matrix(p.res), as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.stunted.clusterWilcoxDAL.csv")
    View(p.res) 
    
    # assess for age (use younger and older than three years)
    meta_wilcox$ageyears2<-meta_wilcox$age
    meta_wilcox$ageyears2[meta_wilcox$age>=36]<-5
    meta_wilcox$ageyears2[meta_wilcox$age<36]<-2
    meta_wilcox$ageyears2
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$ageyears2)$p.value)
    clusternames<-colnames(df_wilcox)
    p.res = data.frame(clusternames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.cluster)
    p.res = cbind(as.matrix(p.res), as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.ageyears.clusterWilcoxDAL.csv")
    View(p.res) 
    
    # assess for calprotectin level
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$calprotectinelevel)$p.value)
    clusternames<-colnames(df_wilcox)
    p.res = data.frame(clusternames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.cluster)
    p.res = cbind(as.matrix(p.res), as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.calpro.clusterWilcoxDAL.csv")
    View(p.res) 
    
    # assess for AAT level
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$alphaantitrypsinlevel)$p.value)
    clusternames<-colnames(df_wilcox)
    p.res = data.frame(clusternames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.cluster)
    p.res = cbind(as.matrix(p.res),as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.AAT.clusterWilcoxDAL.csv")
    View(p.res) 
    
    # assess for anemia level
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$anemie2)$p.value)
    clusternames<-colnames(df_wilcox)
    p.res = data.frame(clusternames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.cluster)
    p.res = cbind(as.matrix(p.res), as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.anemie.clusterWilcoxDAL.csv")
    View(p.res) 
    
    #### which clinical variables are associated with presence/ absence of cluster? ####
    fecesprevalence <- dffiltered_feces_cluster%>% #this produces prevalence "counts" for each cluster, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence))) #take rel abund and Wilcoxon rank-sum
    meta_chi2 <- as.data.frame(sample_data(fecesprevalence)) #take metadata
    taxa_chi2<- tax_table(fecesprevalence) #take taxtable
    
    meta_chi2$stunted<-as.factor(meta_chi2$stunted)
    meta_chi2$pays<-as.factor(meta_chi2$pays)
    meta_chi2$haz<-as.factor(meta_chi2$haz)
    meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
    meta_chi2$alphaantitrypsinlevel<-as.factor(meta_chi2$alphaantitrypsinlevel)
    meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
    
    clusternames<-taxa_chi2[, 12]
    
    # on country of origin
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value)
    chi2results = data.frame(clusternames,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.clusterpaysDAL.csv")
    
    # on stunting status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
    chi2results = data.frame(clusternames,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.clusterstuntedDAL.csv")
    
    # on age 
    meta_chi2$ageyears2<-meta_chi2$age
    meta_chi2$ageyears2[meta_chi2$age>=36]<-5
    meta_chi2$ageyears2[meta_chi2$age<36]<-2
    meta_chi2$ageyears2
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears2), simulate.p.value = TRUE)$p.value)
    chi2results = data.frame(clusternames,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.clusterageyears2DAL.csv")
    
    # on calprotectine status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears2), simulate.p.value = TRUE)$p.value)
    chi2results = data.frame(clusternames,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.clustercalproDAL.csv")
    
    
    # on AAT status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value)
    chi2results = data.frame(clusternames,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.clusteraatDAL.csv")
    
    # on anemia status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
    chi2results = data.frame(clusternames,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.clusteranemiaDAL.csv")
    
    
#### now make a plot of the prevalence of different Entamoeba subtypes by country of origin ####
    ## first add Species level data of taxa to metadata ##
    dffiltered_feces= tax_glom(dfcleanfeces, "rank8")
    dffiltered_prev<-transform_sample_counts(dffiltered_feces, fun = prevalence) 
    otu_table_g=otu_table(dffiltered_prev) %>%
      unclass() %>%
      as.data.frame()
    
    View(otu_table_g)
    dim(otu_table_g)
    otu_table_g<-as.data.frame(t(otu_table_g))
    
    dataframe<-as.data.frame(tax_table(dffiltered_prev)[, 8])
    colnames(dataframe)<-"taxonomy"
    View(dataframe)
    
    otu_table_g$taxonomy<-dataframe$taxonomy
    otu_table_g<-t(otu_table_g)
    dim(otu_table_g)
    colnames(otu_table_g)<-otu_table_g[242,]
    
    sample_data_g=as.data.frame(sample_data(dffiltered_prev))
    Speciesdata18S=merge(otu_table_g, sample_data_g, by="row.names", all=TRUE)
    View(Speciesdata18S)
    dim(Speciesdata18S)
    Speciesdata18S = Speciesdata18S[-242,] # to take away the line with taxonomy text
    
    Speciesdata18S$Row.names=as.character(Speciesdata18S$Row.names)
    write.csv(Speciesdata18S, "Speciesdata18Sfeces.csv")
    View(Speciesdata18S)
    
    # filter out the Entamoeba columns
    Speciesdata18S_Entamoeba<-Speciesdata18S[, -(2:9)]
    Speciesdata18S_Entamoeba<-Speciesdata18S_Entamoeba[, -(3:4)]
    Speciesdata18S_Entamoeba<-Speciesdata18S_Entamoeba[, -(5:9)]
    Speciesdata18S_Entamoeba<-Speciesdata18S_Entamoeba[, -(6)]
    
    ## now make the actual graph
    library(binom)
    
    percentage= function(x) {per <- (sum(as.numeric(as.character(x)), na.rm=T)/length(x))*100 
    return (per)}
    
    
    Speciesdata18S_EntamoebaCAR <-as.data.frame(Speciesdata18S_Entamoeba[which(Speciesdata18S_Entamoeba$pays=="RCA"), ]) # to generate a list of all the ones that are from CAR
    Speciesdata18S_EntamoebaMada <-as.data.frame(Speciesdata18S_Entamoeba[which(Speciesdata18S_Entamoeba$pays=="Madagascar"), ]) # to generate a list of all the ones that are from Madagascar
    
    nrow(Speciesdata18S_EntamoebaCAR) # 91
    nrow(Speciesdata18S_EntamoebaMada) # 150
    
    Speciesdata18S_EntamoebaCAR2<-Speciesdata18S_EntamoebaCAR[, 2:6] 
    Speciesdata18S_EntamoebaCAR2<- as.matrix(sapply(Speciesdata18S_EntamoebaCAR2, as.numeric))
    measured<-as.data.frame(colSums(Speciesdata18S_EntamoebaCAR2))
    colnames(measured)<-"Positive"
    dim(measured)
    
    theoretical<-numeric( length = 5)
    theoretical[1:5] <- 91
    
    test<-as.data.frame(cbind(measured, theoretical))
    
    CICAR <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Speciesdata18S_EntamoebaMada2<-Speciesdata18S_EntamoebaMada[, 2:6] 
    Speciesdata18S_EntamoebaMada2<- as.matrix(sapply(Speciesdata18S_EntamoebaMada2, as.numeric))
    measured<-colSums(Speciesdata18S_EntamoebaMada2)
    length(measured)
    theoretical<-numeric( length = 5)
    theoretical[1:5] <- 150
    
    test<-as.data.frame(cbind(measured, theoretical))
    
    CIMada <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Lower=rbind(CICAR$lower*100, CIMada$lower*100)
    row.names(Lower)=c("CAR", "Madagascar")
    colnames(Lower)<-colnames(Speciesdata18S_EntamoebaCAR2)
    
    Upper=rbind(CICAR$upper*100, CIMada$upper*100)
    row.names(Upper)=c("CAR", "Madagascar")
    colnames(Upper)<-colnames(Speciesdata18S_EntamoebaMada2)
    
    Speciesdata18S_EntamoebaCARcollapsed=colwise(percentage)(Speciesdata18S_EntamoebaCAR)
    Speciesdata18S_EntamoebaMadacollapsed=colwise(percentage)(Speciesdata18S_EntamoebaMada)
    
    PercentageEntamoeba=rbind(Speciesdata18S_EntamoebaCARcollapsed, Speciesdata18S_EntamoebaMadacollapsed)
    
    row.names(PercentageEntamoeba)=c("CAR", "Madagascar")
    View(PercentageEntamoeba)
    ncol(PercentageEntamoeba)
    
    PercentageEntamoeba <- PercentageEntamoeba[,-(7:42), drop=FALSE] # kick-out the metadata
    PercentageEntamoeba <- PercentageEntamoeba[,-(1),drop=FALSE] # kick-out the metadata
    
    write.csv(PercentageEntamoeba, "PercentageEntamoebasubtypespaysfeces.csv")
    
    PercentageEntamoeba2=t(PercentageEntamoeba)
    Upper2=t(Upper)
    Lower2=t(Lower)
    
    PercentageEntamoeba2toplot=melt(PercentageEntamoeba2, value.name="PercentageEntamoeba2", varnames=c("Taxon", "pays"))
    Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "pays"))
    Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "pays"))
    PercentageEntamoeba2toplot<-merge(PercentageEntamoeba2toplot, Uppertoplot, by=c("Taxon", "pays"))
    PercentageEntamoeba2toplot<-merge(PercentageEntamoeba2toplot, Lowertoplot, by=c("Taxon", "pays"))
    
    PercentageEntamoeba2toplot<-PercentageEntamoeba2toplot[order(-PercentageEntamoeba2),]
    
    pdf("Entamoebasubtypesaccordingtopaysfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 6, # define plot width and height. completely up to user.
        height = 5)
    ggplot(data=PercentageEntamoeba2toplot, aes(x=reorder(Taxon, -PercentageEntamoeba2), y=PercentageEntamoeba2,  color=pays, fill=pays)) +
      geom_errorbar(data=PercentageEntamoeba2toplot, mapping=aes(x=reorder(Taxon, -PercentageEntamoeba2), ymin=Lower, ymax=Upper),
                    size=0.3,
                    width=0.9,  
                    color="black", 
                    position=position_dodge(.9)) +
      geom_bar(stat="identity", position="dodge") +
      theme_bw() +
      theme(strip.background=element_rect(fill="white")) +
      theme(strip.text=element_text(color="black")) +
      scale_color_manual(values=c("blue", "red", "green")) +
      scale_fill_manual(values=c("blue", "red", "green")) +
      theme(axis.text.x=element_text(angle=-45, vjust=1, hjust=0)) +
      labs(title="Percentage of different Entamoeba subtypes in fecal samples according to country of origin")
    dev.off()
    
    
    
    #### now make a plot of the prevalence of different Entamoeba subtypes by stunting status ####
    ## first add Species level data of taxa to metadata ##
    dffiltered_feces= tax_glom(dfcleanfeces, "rank8")
    dffiltered_prev<-transform_sample_counts(dffiltered_feces, fun = prevalence) 
    otu_table_g=otu_table(dffiltered_prev) %>%
      unclass() %>%
      as.data.frame()
    
    View(otu_table_g)
    dim(otu_table_g)
    otu_table_g<-as.data.frame(t(otu_table_g))
    
    dataframe<-as.data.frame(tax_table(dffiltered_prev)[, 8])
    colnames(dataframe)<-"taxonomy"
    View(dataframe)
    
    otu_table_g$taxonomy<-dataframe$taxonomy
    otu_table_g<-t(otu_table_g)
    dim(otu_table_g)
    colnames(otu_table_g)<-otu_table_g[242,]
    
    sample_data_g=as.data.frame(sample_data(dffiltered_prev))
    Speciesdata18S=merge(otu_table_g, sample_data_g, by="row.names", all=TRUE)
    View(Speciesdata18S)
    dim(Speciesdata18S)
    Speciesdata18S = Speciesdata18S[-242,] # to take away the line with taxonomy text
    
    Speciesdata18S$Row.names=as.character(Speciesdata18S$Row.names)
    write.csv(Speciesdata18S, "Speciesdata18Sfeces.csv")
    View(Speciesdata18S)
    
    Speciesdata18S_Entamoeba<-Speciesdata18S[, -(2:9)]
    Speciesdata18S_Entamoeba<-Speciesdata18S_Entamoeba[, -(3:4)]
    Speciesdata18S_Entamoeba<-Speciesdata18S_Entamoeba[, -(5:9)]
    Speciesdata18S_Entamoeba<-Speciesdata18S_Entamoeba[, -(6)]
    
    
        ## now make the actual graph
    library(binom)
    
    percentage= function(x) {per <- (sum(as.numeric(as.character(x)), na.rm=T)/length(x))*100 
    return (per)}
    
    
    Speciesdata18S_Entamoebastunted <-as.data.frame(Speciesdata18S_Entamoeba[which(Speciesdata18S_Entamoeba$stunted=="stunted"), ]) # to generate a list of all the ones that are from stunted
    Speciesdata18S_Entamoebanonstunted <-as.data.frame(Speciesdata18S_Entamoeba[which(Speciesdata18S_Entamoeba$stunted=="non-stunted"), ]) # to generate a list of all the ones that are from nonstunted
    
    nrow(Speciesdata18S_Entamoebastunted) # 126
    nrow(Speciesdata18S_Entamoebanonstunted) # 115
    
    Speciesdata18S_Entamoebastunted2<-Speciesdata18S_Entamoebastunted[, 2:6] 
    Speciesdata18S_Entamoebastunted2<- as.matrix(sapply(Speciesdata18S_Entamoebastunted2, as.numeric))
    measured<-as.data.frame(colSums(Speciesdata18S_Entamoebastunted2))
    colnames(measured)<-"Positive"
    dim(measured)
    
    theoretical<-numeric( length = 5)
    theoretical[1:5] <- 126
    
    test<-as.data.frame(cbind(measured, theoretical))
    
    CIstunted <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Speciesdata18S_Entamoebanonstunted2<-Speciesdata18S_Entamoebanonstunted[, 2:6] 
    Speciesdata18S_Entamoebanonstunted2<- as.matrix(sapply(Speciesdata18S_Entamoebanonstunted2, as.numeric))
    measured<-colSums(Speciesdata18S_Entamoebanonstunted2)
    length(measured)
    theoretical<-numeric( length = 5)
    theoretical[1:5] <- 115
    
    test<-as.data.frame(cbind(measured, theoretical))
    
    CInonstunted <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Lower=rbind(CIstunted$lower*100, CInonstunted$lower*100)
    row.names(Lower)=c("stunted", "nonstunted")
    colnames(Lower)<-colnames(Speciesdata18S_Entamoebastunted2)
    
    Upper=rbind(CIstunted$upper*100, CInonstunted$upper*100)
    row.names(Upper)=c("stunted", "nonstunted")
    colnames(Upper)<-colnames(Speciesdata18S_Entamoebanonstunted2)
    
    Speciesdata18S_Entamoebastuntedcollapsed=colwise(percentage)(Speciesdata18S_Entamoebastunted)
    Speciesdata18S_Entamoebanonstuntedcollapsed=colwise(percentage)(Speciesdata18S_Entamoebanonstunted)
    
    PercentageEntamoeba=rbind(Speciesdata18S_Entamoebastuntedcollapsed, Speciesdata18S_Entamoebanonstuntedcollapsed)
    
    row.names(PercentageEntamoeba)=c("stunted", "nonstunted")
    View(PercentageEntamoeba)
    ncol(PercentageEntamoeba)
    dim(PercentageEntamoeba)
    
    write.csv(PercentageEntamoeba, "PercentageEntamoebasubtypesstuntedfeces.csv")
    
    PercentageEntamoeba2=t(PercentageEntamoeba)
    Upper2=t(Upper)
    Lower2=t(Lower)
    
    PercentageEntamoeba2toplot=melt(PercentageEntamoeba2, value.name="PercentageEntamoeba2", varnames=c("Taxon", "stunted"))
    Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "stunted"))
    Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "stunted"))
    PercentageEntamoeba2toplot<-merge(PercentageEntamoeba2toplot, Uppertoplot, by=c("Taxon", "stunted"))
    PercentageEntamoeba2toplot<-merge(PercentageEntamoeba2toplot, Lowertoplot, by=c("Taxon", "stunted"))
    
    PercentageEntamoeba2toplot<-PercentageEntamoeba2toplot[order(-PercentageEntamoeba2),]
    
    pdf("Entamoebasubtypesaccordingtostuntedfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 6, # define plot width and height. completely up to user.
        height = 5)
    ggplot(data=PercentageEntamoeba2toplot, aes(x=reorder(Taxon, -PercentageEntamoeba2), y=PercentageEntamoeba2,  color=stunted, fill=stunted)) +
      geom_errorbar(data=PercentageEntamoeba2toplot, mapping=aes(x=reorder(Taxon, -PercentageEntamoeba2), ymin=Lower, ymax=Upper),
                    size=0.3,
                    width=0.9,  
                    color="black", 
                    position=position_dodge(.9)) +
      geom_bar(stat="identity", position="dodge") +
      theme_bw() +
      theme(strip.background=element_rect(fill="white")) +
      theme(strip.text=element_text(color="black")) +
      scale_color_manual(values=c("green", "red")) +
      scale_fill_manual(values=c("green", "red")) +
      theme(axis.text.x=element_text(angle=-45, vjust=1, hjust=0)) +
      labs(title="Percentage of different Entamoeba subtypes in fecal samples according to stunting status")
    dev.off()
    
    
 #### find out which samples have Entamoeba and if they are having one or several subtypes cutoff 2 or 5 seqs ####
    
    # first subset samples to have only info on Entamoeba
    dfcleanfeces
    dfcleanfeces_Entamoeba=subset_taxa(dfcleanfeces, rank7=="Entamoeba")
    dfcleanfeces_Entamoeba_s=tax_glom(dfcleanfeces_Entamoeba, "rank8") # this gives us a otu_table with how many subtypes per sample
    
    prevalenceyesno2 = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function, but need to be expressed as a percentage only for samples belonging to a particular species
      x[x >= 2] <- 1
      return(x)}
    
    prevalenceyesno = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function, but need to be expressed as a percentage only for samples belonging to a particular species
      x[x >= 5] <- 1
      x[x >= 1] <- 1
      return(x)}
    
    dfcleanfeces_Entamoeba_s_pres <- dfcleanfeces_Entamoeba_s %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalenceyesno2)
    
    dfcleanfeces_Entamoeba_s_pres5 <- dfcleanfeces_Entamoeba_s %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalenceyesno)
    
    
    # transform real counts into table with header
    otu_table=as.data.frame(t(otu_table(dfcleanfeces_Entamoeba_s)))
    dim(otu_table)
    View(otu_table)
    
    tax_table=as.data.frame(t(tax_table(dfcleanfeces_Entamoeba_s)))
    tax_table=tax_table[-(1:7), ]
    tax_table=tax_table[-(2:4), ]
    dim(tax_table)
    View(tax_table)
    
    tax_table<-t(tax_table)
    
    row.names(otu_table)<-tax_table[,1]
    View(otu_table)
    
    otu_table$sumrow <- rowSums(otu_table, na.rm = TRUE)
    
    write.csv(otu_table, "EntamoebarealcountsAfribiota.csv")
    
    # transform real counts into table with header at least two or five sequences
    
    otu_table=as.data.frame(t(otu_table(dfcleanfeces_Entamoeba_s_pres5)))
    dim(otu_table)
    View(otu_table)
    row.names(otu_table)<-tax_table[, 1]
    View(otu_table)
    
    otu_table$sumrow <- rowSums(otu_table, na.rm = TRUE)
    
    write.csv(otu_table, "EntamoebaatleastfiveAfribiota.csv")
    otu_table<-as.data.frame(t(otu_table))
    
    # summarize your Entamoeba positive samples
    
    sum(otu_table$Entamoeba_coli == 1) # 54 samples have Entamoeba coli
    sum(otu_table$Entamoeba_dispar == 1) # 50 samples have Entamoeba dispar
    sum(otu_table$Entamoeba_hartmanni == 1) # 53 samples have Entamoeba hartmanni
    sum(otu_table$Entamoeba_polecki == 1) # 20 samples have Entamoeba polecki
    sum(otu_table$Entamoeba_bovis == 1) # 25 samples have Entamoeba bovis
    
    otu_table$sumrow<-rowSums(otu_table)
    sum(otu_table$sumrow >=1) #123 samples have Entamoeba count of at least 1
    sum(otu_table$sumrow >=2) #61 samples have Entamoeba count of at least 2¨
    dim(otu_table)
    
    100/242*123  ## hence 51% are positive
    100/242*54  ## hence 22% are positive
    100/242*50  ## hence 21% are positive
    100/242*53  ## hence 22% are positive
    100/242*20  ## hence 8% are positive
    100/242*25  ## hence 10% are positive
    
    100/242*61  ## hence 25% have more than one species
    
    sum(otu_table$sumrow ==0) # 119 samples have no Entamoeba
    
    # this would also work for all numeric columns in a single command
    res <- unlist(lapply(otu_table, function(x) if(is.numeric(x)) sum(x, na.rm=T)))
    res
    
    #  summarize your Entamoeba positive samples which have only one Entamoeba
  
    sum(otu_table$Entamoeba_coli == 1 & otu_table$Entamoeba_dispar == 0 & otu_table$Entamoeba_hartmanni == 0) 
    listcoli=rownames(otu_table)[otu_table$Entamoeba_coli == 1 & otu_table$Entamoeba_dispar == 0 & otu_table$Entamoeba_hartmanni == 0] 
    listcoli # 28 have only one Entamoeba coli
    
    sum(otu_table$Entamoeba_coli == 0 & otu_table$Entamoeba_dispar == 1 & otu_table$Entamoeba_hartmanni == 0) 
    listdispar=rownames(otu_table)[otu_table$Entamoeba_coli == 0 & otu_table$Entamoeba_dispar == 1 & otu_table$Entamoeba_hartmanni == 0] 
    listdispar # 24 sample have only Entamoeba dispar
    
    sum(otu_table$Entamoeba_coli == 0 & otu_table$Entamoeba_dispar == 0 & otu_table$Entamoeba_hartmanni == 1) 
    listhartmanni=rownames(otu_table)[otu_table$Entamoeba_coli == 0 & otu_table$Entamoeba_dispar == 0 & otu_table$Entamoeba_hartmanni == 1] 
    listhartmanni # 19 sample have only Entamoeba hartmanni
    
    # summarize your Entamoeba positive samples which have several Entamoeba, first make a column saying how many different Entamoeba species
    sum(otu_table$sumrow == 1) # 62/242 samples have only one Entamoeba 
    100/242*62 # 26%
    
    sum(otu_table$sumrow == 2) # 43/242 samples have two Entamoeba
    100/242*43 # 18%
    
    sum(otu_table$sumrow == 3) # 15/242 samples have tree Entamoeba
    100/242*15 # 6%
    
    sum(otu_table$sumrow >1) # 61/242 samples have several Entamoeba
    100/242*61 # 25%
    100/123*61 # corresponding to 50% of these which have at all Entamoeba
    
    sum(otu_table$sumrow==0) # 119/242 have no Entamoeba at all
    
 
 #### which clinical variables are associated with presence/ absence  at entamoeba subtype level? ####
   
    dffiltered_feces= tax_glom(dfcleanfeces, "rank8")
    dffiltered_prev<-transform_sample_counts(dffiltered_feces, fun = prevalence) 
    otu_table_g=otu_table(dffiltered_prev) %>%
      unclass() %>%
      as.data.frame()
    
    dffiltered_prev_entamoeba<-subset_taxa(dffiltered_prev, rank7=="Entamoeba")
    
    Speciesotu<-otu_table(dffiltered_prev_entamoeba)
    Metadataotu<-sample_data(dffiltered_prev_entamoeba)
    
    rank8names<-tax_table(dffiltered_prev_entamoeba)[, 8]
    
    df_chi2 <- as.matrix(Speciesotu) #take Entamoeba pres abs data
    meta_chi2 <- as.data.frame(Metadataotu) #take metadata
    
    # on country of origin
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value)
    
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebapaysDAL.csv")
    
    # on stunting status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebastuntedDAL.csv")
    
    # on age 
    meta_chi2$ageyears2<-meta_chi2$age
    meta_chi2$ageyears2[meta_chi2$age>=36]<-5
    meta_chi2$ageyears2[meta_chi2$age<36]<-2
    meta_chi2$ageyears2
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears2), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.Entamoebaageyears2DAL.csv")
    
    # on calprotectine status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebacalproDAL.csv")
    
    
    # on AAT status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebaaatDAL.csv")
    
    # on anemia status
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebaanemiaDAL.csv")
    

 #### Make rarefaction curves  (see if the number still increases when samples added) and rarefy your dataset in an iterative manner to avoid bias ####
    library(vegan)
    raremax <- min(rowSums(t(data.frame(otu_table(dfclean)))))
    raremax
    dfclean
    
    col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
    lty <- c("solid", "dashed", "longdash", "dotdash")
    pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
    
    out <- with(pars[1:57, ],
                rarecurve(t(data.frame(otu_table(dfclean))), step = 20, sample = 10000, col = "red",
                          lty = "solid", label = FALSE))
    
    # we can also make rarefaction curves colored by SampleType
    rarecurve(t(data.frame(otu_table(dfclean))), step=100, col=sample_data(dfclean)$SampleType, lwd=3, ylab="ASVs", label=F)
    abline(v=(min(rowSums(t(data.frame(otu_table(dfclean)))))))
    
    
    
    #### rarefy data to 1000 sequences ####
    # Now make a loop and rarefy to even depth in repetitive manner to avoid inducing bias. It is always important to set a seed when you subsample so your result is replicable 
    
    set.seed(3)
    for (i in 1:100) {
      # Subsample
      dfrar1000 <- rarefy_even_depth(dfcleanfeces, sample.size = 1000, verbose = FALSE, replace = TRUE)}
    
    length(which(sample_data(dfrar1000)$pays=="RCA")) #54 are from RCA
    length(which(sample_data(dfrar1000)$pays=="Madagascar")) #98 are from Mada
    
    # kick out samples that do not match the criteria for dataset with singlets
    sample_data(data2)$ph_estomac<-as.numeric(sample_data(data2)$ph_estomac)
    sample_data(data2)$ph_intestin<-as.numeric(sample_data(data2)$ph_intestin)
    
    data2 <- data2 %>%
      subset_samples(raison_hospi=="Recrutement communautaire") # to filter out samples which are not recruited in the community
    
    
    data2 <- data2 %>%
      subset_samples((SampleType == "feces") | ph_intestin>=5 | ph_estomac<=4) # to keep only samples with good pH range (<=4 for gastric, >=5 for duodenal))
    
    data2 <- data2 %>%
      subset_samples(whz_cont <=2)
    max(sample_data(df)$whz_cont)
    
    data2 <- data2 %>%
      subset_samples(whz_cont >=-2)
    
    min(sample_data(data2)$whz_cont) # to keep only samples with no acute undernutrition or obesity
    
    data2 # we are left with 289 samples
    
    # now filter out low sample sums samples. We have chosen 500, as this is the average in the blanks
    sample_data(data2)$read_count<-as.numeric(sample_sums(data2))
    data2 = prune_samples((sample_data(data2)$read_count>=500)==TRUE, data2) # only keep samples with more than 500 reads
    dim(sample_data(data2)) ## we are left with 266 samples
    
    # filter out low abundance taxa, as they are likely contaminants. Here, we decided that an ASV needs to have at least 50 sequences in at least 1% of the samples
    data2=filter_taxa(data2, function(x) sum(x > 50) > (0.01*length(x)), TRUE)
    dim(tax_table(data2)) #174 ASV's. Note: we lost a lot of taxa, mainly due to the fact that they are spuriously distributed in between the samples!
    
    # now also kick out the non microeuk reads
    data2 <- data2 %>%
      subset_taxa(rank4!= "Mammalia") %>%
      subset_taxa(rank5!= "Mammalia") %>%
      subset_taxa(rank6!= "Mammalia") %>%
      subset_taxa(rank4!= "Vertebrata") %>%
      subset_taxa(rank5!= "Vertebrata") %>%
      subset_taxa(rank6!= "Vertebrata") %>%
      subset_taxa(rank2!= "Archaeplastida") %>% #plants and algae
      subset_taxa(rank4!= "Porifera") %>% #sponges
      subset_taxa(rank4!= "Arthropoda") %>% #arthropods
      subset_taxa(rank5!= "Arthropoda") %>% #arthropods
      subset_taxa(rank6!= "Arthropoda") %>% #arthropods
      subset_taxa(rank4!= "Mollusca") %>% #mussels
      subset_taxa(rank2!="Unassigned")
    
    data2 # 127 taxa are left
    
    set.seed(3)
    for (i in 1:100) {
      # Subsample
      data2_rar <- rarefy_even_depth(data2, sample.size = 1000, verbose = FALSE, replace = TRUE)}
    
    
#### Alpha diversity measures on untrimmed, rarefied datasets ####
    
    df_SampleType<-subset_samples(data2_rar, (SampleType=="gastric" | SampleType=="duodenal" | SampleType=="feces"))
    sample_data(df_SampleType)$SampleType<- factor(sample_data(df_SampleType)$SampleType,
                                                   levels = c("gastric", "duodenal", "feces"))
    
    pdf("alphadiversityASVlevelSampleTypefeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_SampleType, x= "SampleType", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to Sample Type")
    p + geom_boxplot(data = p$data, aes(x = SampleType, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Sample Type")
    dev.off()
    
    
    data2_feces<-subset_samples(data2_rar, SampleType=="feces") # do it on untrimmed dataset and filter out the unwanted samples. Do it on rarefied dataset
    data2_feces<-subset_samples(data2_feces, pays!="")
    data2_feces #152 samples
    
    data2_duodenal<-subset_samples(data2_rar, SampleType=="duodenal") # do it on untrimmed dataset and filter out the unwanted samples. Do it on rarefied dataset
    data2_duodenal<-subset_samples(data2_duodenal, pays!="")
    data2_duodenal #only one sample is left
    
    
    # on country of origin
    pdf("alphadiversityASVlevelpaysfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_feces, x= "pays", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to Country of origin")
    p + geom_boxplot(data = p$data, aes(x = pays, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Country of origin")
    dev.off()
     
    # calculate if there is a significant difference in Indeces
    results = estimate_richness(data2_feces, measures="Shannon")
    d = sample_data(data2_feces)
    NN = results[d[,'haz'] == 'normonutris',]
    MCM = results[d[,'haz'] == 'malnutris chronique modere',]
    MCS = results[d[,'haz'] == 'malnutris chronique severe',]
    wilcox.test(NN, MCM)
    wilcox.test(NN, MCS)
    wilcox.test(MCM, MCS)
    
    Mada=results[d[,'pays'] == 'Madagascar',]
    CAR=results[d[,'pays'] == 'RCA',]
    wilcox.test(Mada, CAR) # p-value = 4.665e-06
    
    
    results = estimate_richness(data2_feces, measures="Chao1")
    d = sample_data(data2_feces)
    Mada=results[d[,'pays'] == 'Madagascar',]
    CAR=results[d[,'pays'] == 'RCA',]
    wilcox.test(Mada$Chao1, CAR$Chao1) # p-value = 0.02726
    
    results = estimate_richness(data2_feces, measures="Observed")
    d = sample_data(data2_feces)
    Mada=results[d[,'pays'] == 'Madagascar',]
    CAR=results[d[,'pays'] == 'RCA',]
    wilcox.test(Mada, CAR) # p-value = 0.02516
    
    results = estimate_richness(data2_feces, measures="Simpson")
    d = sample_data(data2_feces)
    Mada=results[d[,'pays'] == 'Madagascar',]
    CAR=results[d[,'pays'] == 'RCA',]
    wilcox.test(Mada, CAR) # p-value = 1.628e-06
    

    # on stunting status
    pdf("alphadiversityASVlevelhazfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_feces, x= "haz", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to height-for-age z-score")
    p + geom_boxplot(data = p$data, aes(x = haz, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Height-for-age z-score")
    dev.off()
    
    # calculate if there is a significant difference in Indeces
    results = estimate_richness(data2_feces, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_feces)
    NN = results[d[,'haz'] == 'normonutris',]
    MCM = results[d[,'haz'] == 'malnutris chronique modere',]
    MCS = results[d[,'haz'] == 'malnutris chronique severe',]
    wilcox.test(NN$Shannon, MCM$Shannon)
    wilcox.test(NN$Shannon, MCS$Shannon)
    wilcox.test(MCM$Shannon, MCS$Shannon)
    
    wilcox.test(NN$Chao1, MCM$Chao1)
    wilcox.test(NN$Chao1, MCS$Chao1)
    wilcox.test(MCM$Chao1, MCS$Chao1)
    
    wilcox.test(NN$Simpson, MCM$Simpson)
    wilcox.test(NN$Simpson, MCS$Simpson)
    wilcox.test(MCM$Simpson, MCS$Simpson)
    
    wilcox.test(NN$Observed, MCM$Observed)
    wilcox.test(NN$Observed, MCS$Observed)
    wilcox.test(MCM$Observed, MCS$Observed)
    
    statutnut<-d[,'haz']
    my_data<-cbind(statutnut, results)
    res.aov <- aov(Shannon ~ haz, data = my_data)
    summary(res.aov)
    

    # on stunting status Bangui
    data2_fecesB<-subset_samples(data2_feces, pays=="RCA")
    data2_fecesB # 54 samples
    
    pdf("alphadiversityASVlevelhazfecesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_fecesB, x= "haz", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to height-for-age z-score")
    p + geom_boxplot(data = p$data, aes(x = haz, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Height-for-age z-score")
    dev.off()
    
    # calculate if there is a significant difference Indeces
    results = estimate_richness(data2_fecesB, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_fecesB)
    NN = results[d[,'haz'] == 'normonutris',]
    MCM = results[d[,'haz'] == 'malnutris chronique modere',]
    MCS = results[d[,'haz'] == 'malnutris chronique severe',]
    wilcox.test(NN$Shannon, MCM$Shannon)
    wilcox.test(NN$Shannon, MCS$Shannon)
    wilcox.test(MCM$Shannon, MCS$Shannon)
    
    wilcox.test(NN$Chao1, MCM$Chao1)
    wilcox.test(NN$Chao1, MCS$Chao1)
    wilcox.test(MCM$Chao1, MCS$Chao1)
    
    wilcox.test(NN$Simpson, MCM$Simpson)
    wilcox.test(NN$Simpson, MCS$Simpson)
    wilcox.test(MCM$Simpson, MCS$Simpson)
    
    wilcox.test(NN$Observed, MCM$Observed)
    wilcox.test(NN$Observed, MCS$Observed)
    wilcox.test(MCM$Observed, MCS$Observed)
    
    # on stunting status Tana
    data2_fecesA<-subset_samples(data2_feces, pays=="Madagascar")
    data2_fecesA #98 samples
    
    pdf("alphadiversityASVlevelhazfecesTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_fecesA, x= "haz", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to height-for-age z-score")
    p + geom_boxplot(data = p$data, aes(x = haz, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Height-for-age z-score")
    dev.off()
    
    # calculate if there is a significant difference in Indeces
    results = estimate_richness(data2_fecesA, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_fecesA)
    NN = results[d[,'haz'] == 'normonutris',]
    MCM = results[d[,'haz'] == 'malnutris chronique modere',]
    MCS = results[d[,'haz'] == 'malnutris chronique severe',]
    wilcox.test(NN$Shannon, MCM$Shannon)
    wilcox.test(NN$Shannon, MCS$Shannon)
    wilcox.test(MCM$Shannon, MCS$Shannon)
    
    wilcox.test(NN$Chao1, MCM$Chao1)
    wilcox.test(NN$Chao1, MCS$Chao1)
    wilcox.test(MCM$Chao1, MCS$Chao1)
    
    wilcox.test(NN$Simpson, MCM$Simpson)
    wilcox.test(NN$Simpson, MCS$Simpson)
    wilcox.test(MCM$Simpson, MCS$Simpson)
    
    wilcox.test(NN$Observed, MCM$Observed)
    wilcox.test(NN$Observed, MCS$Observed)
    wilcox.test(MCM$Observed, MCS$Observed)
    
    # on calprotectine levels
    sample_data(data2_feces)$calprotectinelevel<-as.factor(sample_data(data2_feces)$calprotectinelevel)
    data2_fecescalpro<-subset_samples(data2_feces, !is.na(calprotectinelevel))
    data2_fecescalpro #141 sampples
    
    pdf("alphadiversityASVlevelcalprotectinelevelfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_fecescalpro, x= "calprotectinelevel", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to Calprotectine levels")
    p + geom_boxplot(data = p$data, aes(x = calprotectinelevel, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Calprotectine levels")
    dev.off()
    
    # calculate if there is a significant difference in Indeces
    results = estimate_richness(data2_fecescalpro, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_fecescalpro)
    normal = results[d[,'calprotectinelevel'] == '0',]
    high = results[d[,'calprotectinelevel'] == '1',]
    
    wilcox.test(normal$Shannon, high$Shannon)
    wilcox.test(normal$Chao1, high$Chao1)
    wilcox.test(normal$Simpson, high$Simpson)
    wilcox.test(normal$Observed, high$Observed)
    
    
    # on calprotectine levels Bangui
    data2_fecescalproB<-subset_samples(data2_fecescalpro, pays=="RCA")
    data2_fecescalproB #53 samples
    
    pdf("alphadiversityASVlevelcalprotectinelevelfecesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_fecescalproB, x= "calprotectinelevel", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to Calprotectine levels")
    p + geom_boxplot(data = p$data, aes(x = calprotectinelevel, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Calprotectine levels")
    dev.off()
    
    # calculate if there is a significant difference Indeces
    results = estimate_richness(data2_fecescalpro, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_fecescalpro)
    normal = results[d[,'calprotectinelevel'] == '0',]
    high = results[d[,'calprotectinelevel'] == '1',]
    
    wilcox.test(normal$Shannon, high$Shannon)
    wilcox.test(normal$Chao1, high$Chao1)
    wilcox.test(normal$Simpson, high$Simpson)
    wilcox.test(normal$Observed, high$Observed)
    
    # on calprotectine levels Tana
    data2_fecescalproA<-subset_samples(data2_fecescalpro, pays=="Madagascar")
    data2_fecescalproA #88 samples
    
    pdf("alphadiversityASVlevelcalprotectinelevelfecesTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_fecescalproA, x= "calprotectinelevel", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to Calprotectine levels")
    p + geom_boxplot(data = p$data, aes(x = calprotectinelevel, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Calprotectine levels")
    dev.off()
    
    # calculate if there is a significant difference in Indeces
    results = estimate_richness(data2_fecescalpro, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_fecescalpro)
    normal = results[d[,'calprotectinelevel'] == '0',]
    high = results[d[,'calprotectinelevel'] == '1',]
    
    wilcox.test(normal$Shannon, high$Shannon)
    wilcox.test(normal$Chao1, high$Chao1)
    wilcox.test(normal$Simpson, high$Simpson)
    wilcox.test(normal$Observed, high$Observed)
    
    # on alphaantitrypsin levels
    sample_data(data2_feces)$alphaantitrypsinlevel<-as.factor(sample_data(data2_feces)$alphaantitrypsinlevel)
    data2_fecesaat<-subset_samples(data2_feces, !is.na(alphaantitrypsinlevel))
    sample_data(data2_feces)$alphaantitrypsinlevel<-as.factor(sample_data(data2_feces)$alphaantitrypsinlevel)
    levels(sample_data(data2_feces)$alphaantitrypsinlevel)
    data2_fecesaat #141 samples
    
    pdf("alphadiversityASVlevelalphaantitrypsinlevelfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_fecesaat, x= "alphaantitrypsinlevel", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to alphaantitrypsin levels")
    p + geom_boxplot(data = p$data, aes(x = alphaantitrypsinlevel, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("alphaantitrypsin levels")
    dev.off()
    
    # calculate if there is a significant difference in Indeces
    results = estimate_richness(data2_fecesaat, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_fecesaat)
    normal = results[d[,'alphaantitrypsinlevel'] == '0',]
    high = results[d[,'alphaantitrypsinlevel'] == '2',]
    
    wilcox.test(normal$Shannon, high$Shannon)
    wilcox.test(normal$Chao1, high$Chao1)
    wilcox.test(normal$Simpson, high$Simpson)
    wilcox.test(normal$Observed, high$Observed)
    
    
    # on alphaantitrypsin levels Bangui
    data2_fecesaatB<-subset_samples(data2_fecesaat, pays=="RCA")
    data2_fecesaatB #53 samples
    
    pdf("alphadiversityASVlevelalphaantitrypsinlevelfecesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_fecesaatB, x= "alphaantitrypsinlevel", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to alphaantitrypsin levels")
    p + geom_boxplot(data = p$data, aes(x = alphaantitrypsinlevel, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("alphaantitrypsin levels")
    dev.off()
    
    # calculate if there is a significant difference Indeces
    results = estimate_richness(data2_fecesaat, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_fecesaat)
    normal = results[d[,'alphaantitrypsinlevel'] == '0',]
    high = results[d[,'alphaantitrypsinlevel'] == '2',]
    
    wilcox.test(normal$Shannon, high$Shannon)
    wilcox.test(normal$Chao1, high$Chao1)
    wilcox.test(normal$Simpson, high$Simpson)
    wilcox.test(normal$Observed, high$Observed)
    
    # on alphaantitrypsin levels Tana
    data2_fecesaatA<-subset_samples(data2_fecesaat, pays=="Madagascar")
    data2_fecesaatA #88 samples
    
    pdf("alphadiversityASVlevelalphaantitrypsinlevelfecesTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_fecesaatA, x= "alphaantitrypsinlevel", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to alphaantitrypsin levels")
    p + geom_boxplot(data = p$data, aes(x = alphaantitrypsinlevel, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("alphaantitrypsin levels")
    dev.off()
    
    # calculate if there is a significant difference in Indeces
    results = estimate_richness(data2_fecesaat, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_fecesaat)
    normal = results[d[,'alphaantitrypsinlevel'] == '0',]
    high = results[d[,'alphaantitrypsinlevel'] == '2',]
    
    wilcox.test(normal$Shannon, high$Shannon)
    wilcox.test(normal$Chao1, high$Chao1)
    wilcox.test(normal$Simpson, high$Simpson)
    wilcox.test(normal$Observed, high$Observed)
    
    # on anemie
    sample_data(data2_feces)$anemie2<-as.factor(sample_data(data2_feces)$anemie2)
    data2_fecesanemie<-subset_samples(data2_feces, !is.na(anemie2))
    sample_data(data2_feces)$anemie2<-as.factor(sample_data(data2_feces)$anemie2)
    levels(sample_data(data2_feces)$anemie2)
    
    data2_fecesanemie #142 samples
    
    pdf("alphadiversityASVlevelanemie2feces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_fecesanemie, x= "anemie2", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to anemie")
    p + geom_boxplot(data = p$data, aes(x = anemie2, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("anemie")
    dev.off()
    
    # calculate if there is a significant difference in Indeces
    results = estimate_richness(data2_fecesanemie, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_fecesanemie)
    non = results[d[,'anemie2'] == 'Non',]
    oui = results[d[,'anemie2'] == 'Oui',]
    
    wilcox.test(non$Shannon, oui$Shannon)
    wilcox.test(non$Chao1, oui$Chao1)
    wilcox.test(non$Simpson, oui$Simpson)
    wilcox.test(non$Observed, oui$Observed)
    
    
    # on anemie Bangui
    data2_fecesanemieB<-subset_samples(data2_fecesanemie, pays=="RCA")
    
    pdf("alphadiversityASVlevelanemie2fecesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_fecesanemieB, x= "anemie2", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to anemie")
    p + geom_boxplot(data = p$data, aes(x = anemie2, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("anemie")
    dev.off()
    
    # calculate if there is a significant difference Indeces
    results = estimate_richness(data2_fecesanemie, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_fecesanemie)
    non = results[d[,'anemie2'] == 'Non',]
    oui = results[d[,'anemie2'] == 'Oui',]
    
    wilcox.test(non$Shannon, oui$Shannon)
    wilcox.test(non$Chao1, oui$Chao1)
    wilcox.test(non$Simpson, oui$Simpson)
    wilcox.test(non$Observed, oui$Observed)
    
    # on anemie Tana
    data2_fecesanemieA<-subset_samples(data2_fecesanemie, pays=="Madagascar")
    
    pdf("alphadiversityASVlevelanemie2fecesTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_fecesanemieA, x= "anemie2", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to anemie")
    p + geom_boxplot(data = p$data, aes(x = anemie2, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("anemie")
    dev.off()
    
    # calculate if there is a significant difference in Indeces
    results = estimate_richness(data2_fecesanemie, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_fecesanemie)
    non = results[d[,'anemie2'] == 'Non',]
    oui = results[d[,'anemie2'] == 'Oui',]
    
    wilcox.test(non$Shannon, oui$Shannon)
    wilcox.test(non$Chao1, oui$Chao1)
    wilcox.test(non$Simpson, oui$Simpson)
    wilcox.test(non$Observed, oui$Observed)
    
    # on age in years
    sample_data(data2_feces)$ageyears<-cut(as.numeric(as.character(sample_data(data2_feces)$age)), c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
    which(is.na(sample_data(data2_feces)$ageyears)) # the controls have no assocaited age year
    levels(sample_data(data2_feces)$ageyears) <- c("2-3 years", "3-4 years", "4-5 years")
    levels(sample_data(data2_feces)$ageyears)
    
    pdf("alphadiversityASVlevelageyearsfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(data2_feces, x= "ageyears", measures=c("Observed", "Chao1", "Simpson",  "Shannon"), title = "Alpha Diversity according to height-for-age z-score")
    p + geom_boxplot(data = p$data, aes(x = ageyears, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Age in years")
    dev.off()
    
    
    # calculate if there is a significant difference in Indeces
    results = estimate_richness(data2_feces, measures= c("Observed", "Shannon", "Chao1", "Simpson"))
    d = sample_data(data2_feces)
    young = results[d[,'ageyears'] == '2-3 years',]
    middle = results[d[,'ageyears'] == '3-4 years',]
    old = results[d[,'ageyears'] == '4-5 years',]
    wilcox.test(young$Shannon, middle$Shannon)
    wilcox.test(young$Shannon, old$Shannon) #ns
    wilcox.test(middle$Shannon, old$Shannon) #ns
    
    wilcox.test(young$Chao1, middle$Chao1)
    wilcox.test(young$Chao1, old$Chao1)
    wilcox.test(middle$Chao1, old$Chao1)
    
    wilcox.test(young$Simpson, middle$Simpson)
    wilcox.test(young$Simpson, old$Simpson) #ns
    wilcox.test(middle$Simpson, old$Simpson) #ns
    
    wilcox.test(young$Observed, middle$Observed)
    wilcox.test(young$Observed, old$Observed)
    wilcox.test(middle$Observed, old$Observed)
    
    ageyears<-d[,'ageyears']
    my_data<-cbind(ageyears, results)
    res.aov <- aov(Shannon ~ ageyears, data = my_data)
    summary(res.aov)

    #### Make a PERMANOVA to look for factors influencing  dispersion of fecal samples ####
    dffiltered2<-subset_samples(dfrar1000, calprotectinelevel!="")
    dffiltered2<-subset_samples(dffiltered2, anemie2!="")
    dffiltered2<-subset_samples(dffiltered2, (rowSums(otu_table(dffiltered2))!=0))
    dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
    dffiltered2transf<-transform(dffiltered2, "log10")
    
    project_bray <- phyloseq::distance(dffiltered2transf, method = "bray")
    sample_df <- data.frame(sample_data(dffiltered2transf))
    
    
    #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
    res.adonis <- adonis(project_bray ~  haz + anemie2 + calprotectinelevel + ageyears + read_count + pays , data=sample_df, method="bray")
    res.adonis 
    results<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results)
    write.csv(results, "Disperionstestallfecalsamples.csv")
    
    
    #change order to account for confounding
    res.adonis <- adonis(project_bray ~  read_count + pays + haz + anemie2 + calprotectinelevel + ageyears  , data=sample_df, method="bray")
    res.adonis 
    results2<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results2)
    write.csv(results2, "Disperionstestallfecalsamples2.csv")
    
    #change order to account for confounding
    res.adonis <- adonis(project_bray ~  read_count + pays + ageyears+ haz + anemie2 + calprotectinelevel   , data=sample_df, method="bray")
    res.adonis 
    results3<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results3)
    write.csv(results3, "Disperionstestallfecalsamples3.csv")
    
    results3$variable<-row.names(results3)
    results3$Contribution<-(results$R2*100)
    results_filt<-filter(results, results$Pr..F.<0.05)
    
    # make  graph plot
    pdf("Contributionfulldatasetfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 6, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=results_filt, aes(x = reorder(variable, Contribution), y=Contribution)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("% Contribution")+
      ggtitle("Sig. contribution to beta-dispersion")
    dev.off()
    
    # STOPPED HERE now in Bangui dataset only
    dfrar1000B<-subset_samples(dffiltered2transf, pays=="RCA")
    dffiltered2<-subset_samples(dfrar1000B, (rowSums(otu_table(dffiltered2))!=0))
    dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
    
    project_bray <- phyloseq::distance(dffiltered2, method = "bray")
    sample_df <- data.frame(sample_data(dffiltered2))
    
    
    #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
    res.adonis <- adonis(project_bray ~  haz + anemie2 + calprotectinelevel + ageyears + read_count , data=sample_df, method="bray")
    res.adonis 
    results<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results)
    write.csv(results, "DisperionstestallfecalsamplesBangui.csv")
    
    #change order to account for confounding
    res.adonis <- adonis(project_bray ~  read_count + haz + anemie2 + calprotectinelevel + ageyears , data=sample_df, method="bray")
    res.adonis 
    results<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results)
    write.csv(results, "DisperionstestallfecalsamplesBangui2.csv")
    
    results$variable<-row.names(results)
    results$Contribution<-(results$R2*100)
    results_filt<-filter(results, results$Pr..F.<0.05)
    
    # make  graph plot
    pdf("ContributionfulldatasetfecesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 6, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=results_filt, aes(x = reorder(variable, Contribution), y=Contribution)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("% Contribution")+
      ggtitle("Sig. contribution to beta-dispersion")
    dev.off()
    
    # now in Tana dataset only
    dfrar1000A<-subset_samples(dffiltered2transf, pays=="Madagascar")
    dffiltered2<-subset_samples(dfrar1000A, (rowSums(otu_table(dffiltered2))!=0))
    dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
    
    
    project_bray <- phyloseq::distance(dffiltered2, method = "bray")
    sample_df <- data.frame(sample_data(dffiltered2))
    
    
    #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
    res.adonis <- adonis(project_bray ~  haz + anemie2 + calprotectinelevel + ageyears + read_count , data=sample_df, method="bray")
    res.adonis 
    results<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results)
    write.csv(results, "DisperionstestallfecalsamplesTana.csv")
    
    #change order to account for confounding
    res.adonis <- adonis(project_bray ~  read_count + haz + anemie2 + calprotectinelevel + ageyears , data=sample_df, method="bray")
    res.adonis 
    results<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results)
    write.csv(results, "DisperionstestallfecalsamplesTana2.csv")
    
    results$variable<-row.names(results)
    results$Contribution<-(results$R2*100)
    results_filt<-filter(results, results$Pr..F.<0.05)
    
    # make  graph plot
    pdf("ContributionfulldatasetfecesTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 6, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=results_filt, aes(x = reorder(variable, Contribution), y=Contribution)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("% Contribution")+
      ggtitle("Sig. contribution to beta-dispersion")
    dev.off()
    
    
    
    #### Make a PERMANOVA to look for factors influencing  dispersion of fecal samples presence/absence and Jaccard ####
    dfrar1000presabs<-transform_sample_counts(dfrar1000, fun = prevalence)
    dffiltered2<-subset_samples(dfrar1000presabs, calprotectinelevel!="")
    dffiltered2<-subset_samples(dffiltered2, anemie2!="")
    dffiltered2<-subset_samples(dffiltered2, SampleType=="feces")
    dffiltered2<-subset_samples(dffiltered2, (rowSums(otu_table(dffiltered2))!=0))
    dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
    
    project_jaccard <- phyloseq::distance(dffiltered2, method = "jaccard")
    sample_df <- data.frame(sample_data(dffiltered2))
    
    
    #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
    res.adonis <- adonis(project_jaccard ~  haz + anemie2 + calprotectinelevel + ageyears + read_count + pays , data=sample_df, method="jaccard")
    res.adonis 
    results<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results)
    write.csv(results, "Disperionstestallfecalsamplespresabs.csv")
    
    #change order to correct for confounding
    res.adonis <- adonis(project_jaccard ~ read_count + pays +  haz + anemie2 + calprotectinelevel + ageyears , data=sample_df, method="jaccard")
    res.adonis 
    results2<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results2)
    write.csv(results, "Disperionstestallfecalsamplespresabs2.csv")
    
    results$variable<-row.names(results)
    results$Contribution<-(results$R2*100)
    results_filt<-filter(results, results$Pr..F.<0.05)
    
    # make  graph plot
    pdf("Contributionfulldatasetfecespresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 6, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=results_filt, aes(x = reorder(variable, Contribution), y=Contribution)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("% Contribution")+
      ggtitle("Sig. contribution to beta-dispersion")
    dev.off()
    
    # now in Bangui dataset only
    dfrar1000B<-subset_samples(dffiltered2, pays=="RCA")
    dffiltered2<-subset_samples(dfrar1000B, (rowSums(otu_table(dffiltered2))!=0))
    dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
    
    project_jaccard <- phyloseq::distance(dffiltered2, method = "jaccard")
    sample_df <- data.frame(sample_data(dffiltered2))
    
    
    #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
    res.adonis <- adonis(project_jaccard ~  haz + anemie2 + calprotectinelevel + ageyears + read_count , data=sample_df, method="jaccard")
    res.adonis 
    results<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results)
    write.csv(results, "DisperionstestallfecalsamplespresabsBangui.csv")
    
    #now change to make sure we account for confounding
    res.adonis <- adonis(project_jaccard ~  read_count +  haz + anemie2 + calprotectinelevel + ageyears , data=sample_df, method="jaccard")
    res.adonis 
    results2<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results2)
    write.csv(results, "DisperionstestallfecalsamplespresabsBangui2.csv")
    
    
    results$variable<-row.names(results)
    results$Contribution<-(results$R2*100)
    results_filt<-filter(results, results$Pr..F.<0.05)
    
    # make  graph plot
    pdf("ContributionfulldatasetfecesBanguifecespresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 6, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=results_filt, aes(x = reorder(variable, Contribution), y=Contribution)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("% Contribution")+
      ggtitle("Sig. contribution to beta-dispersion")
    dev.off()
    
    # now in Tana dataset only
    dfrar1000A<-subset_samples(dffiltered2transf, pays=="Madagascar")
    dffiltered2<-subset_samples(dfrar1000A, (rowSums(otu_table(dffiltered2))!=0))
    dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
    
    
    project_jaccard <- phyloseq::distance(dffiltered2, method = "jaccard")
    sample_df <- data.frame(sample_data(dffiltered2))
    
    
    #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
    res.adonis <- adonis(project_jaccard ~  haz + anemie2 + calprotectinelevel + ageyears + read_count , data=sample_df, method="jaccard")
    res.adonis 
    results<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results)
    write.csv(results, "DisperionstestallfecalsamplespresabsTana.csv")
    
    #now control for confounding by changing position
    res.adonis <- adonis(project_jaccard ~  read_count + haz + anemie2 + calprotectinelevel + ageyears , data=sample_df, method="jaccard")
    res.adonis 
    results2<-data.frame(as.data.frame(res.adonis$aov.tab))
    View(results2)
    write.csv(results2, "DisperionstestallfecalsamplespresabsTana.csv")
    
    results2$variable<-row.names(results)
    results$Contribution<-(results2$R2*100)
    results_filt<-filter(results2, results$Pr..F.<0.05)
    
    # make  graph plot
    pdf("ContributionfulldatasetfecesTanapresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 6, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=results_filt, aes(x = reorder(variable, Contribution), y=Contribution)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("% Contribution")+
      ggtitle("Sig. contribution to beta-dispersion")
    dev.off()
    
    
    
    
    #### PCoA illustrations for rarefied dataset rel. abundance and presence absence on ASV level ####
    dfrar1000rel<-transform_sample_counts(dfrar1000, function(x) x / sum(x) *100)
    dfrar1000reltransf<- transform(dfrar1000rel, 'log10')
    
    ## pays of origin
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000reltransf, "PCoA", "bray")
    
    pdf("DalhousiefecesbraycurtisPCoATaxapays.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000reltransf, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoApays.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000reltransf, GP.ord, type="samples", color="pays") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000prel<-dfrar1000 %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000prel, "PCoA", "jaccard")
    
    pdf("DalhousiefecesbraycurtisPCoApayspresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000prel, GP.ord, type="samples", color="pays") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoATaxapayspresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000prel, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    ## calprotectinelevel
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000reltransf, "PCoA", "bray")
    sample_data(dfrar1000reltransf)$calprotectinelevel<-as.factor(sample_data(dfrar1000reltransf)$calprotectinelevel)
    sample_data(dfrar1000)$calprotectinelevel<-as.factor(sample_data(dfrar1000)$calprotectinelevel)
    
    pdf("DalhousiefecesbraycurtisPCoATaxacalprotectinelevel.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000reltransf, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoAcalprotectinelevel.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000reltransf, GP.ord, type="samples", color="calprotectinelevel") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000prel<-dfrar1000 %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000prel, "PCoA", "jaccard")
    
    pdf("DalhousiefecesbraycurtisPCoAcalprotectinelevelpresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000prel, GP.ord, type="samples", color="calprotectinelevel") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoATaxacalprotectinelevelpresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000prel, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    ## anemie2
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000reltransf, "PCoA", "bray")
    sample_data(dfrar1000reltransf)$anemie2<-as.factor(sample_data(dfrar1000reltransf)$anemie2)
    sample_data(dfrar1000)$anemie2<-as.factor(sample_data(dfrar1000)$anemie2)
    
    pdf("DalhousiefecesbraycurtisPCoATaxaanemie2.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000reltransf, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoAanemie2.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000reltransf, GP.ord, type="samples", color="anemie2") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000prel<-dfrar1000 %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000prel, "PCoA", "jaccard")
    
    pdf("DalhousiefecesbraycurtisPCoAanemie2presabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000prel, GP.ord, type="samples", color="anemie2") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoATaxaanemie2presabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000prel, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    ## calprotectinelevel Bangui
    dfrar1000reltransfB<-subset_samples(dfrar1000reltransf, pays="RCA")
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000reltransfB, "PCoA", "bray")
    sample_data(dfrar1000reltransfB)$calprotectinelevel<-as.factor(sample_data(dfrar1000reltransfB)$calprotectinelevel)
    sample_data(dfrar1000)$calprotectinelevel<-as.factor(sample_data(dfrar1000)$calprotectinelevel)
    
    pdf("DalhousiefecesbraycurtisPCoATaxacalprotectinelevelBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000reltransfB, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoAcalprotectinelevelBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000reltransfB, GP.ord, type="samples", color="calprotectinelevel") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000B<-subset_samples(dfrar1000, pays="RCA")
    dfrar1000prelB<-dfrar1000B %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000prelB, "PCoA", "jaccard")
    
    pdf("DalhousiefecesbraycurtisPCoAcalprotectinelevelBanguipresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000prelB, GP.ord, type="samples", color="calprotectinelevel") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoATaxacalprotectinelevelBanguipresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000prelB, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    ## anemie2 Bangui
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000reltransfB, "PCoA", "bray")
    sample_data(dfrar1000reltransfB)$anemie2<-as.factor(sample_data(dfrar1000reltransfB)$anemie2)
    sample_data(dfrar1000)$anemie2<-as.factor(sample_data(dfrar1000)$anemie2)
    
    pdf("DalhousiefecesbraycurtisPCoATaxaanemie2Bangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000reltransfB, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoAanemie2Bangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000reltransfB, GP.ord, type="samples", color="anemie2") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000prelB<-dfrar1000 %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000prelB, "PCoA", "jaccard")
    
    pdf("DalhousiefecesbraycurtisPCoAanemie2presabsBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000prelB, GP.ord, type="samples", color="anemie2") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoATaxaanemie2Banguipresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000prelB, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    ## calprotectinelevel Tana
    dfrar1000reltransfA<-subset_samples(dfrar1000reltransf, pays="Madagascar")
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000reltransfA, "PCoA", "bray")
    sample_data(dfrar1000reltransfA)$calprotectinelevel<-as.factor(sample_data(dfrar1000reltransfA)$calprotectinelevel)
    sample_data(dfrar1000)$calprotectinelevel<-as.factor(sample_data(dfrar1000)$calprotectinelevel)
    
    pdf("DalhousiefecesbraycurtisPCoATaxacalprotectinelevelTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000reltransfA, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoAcalprotectinelevelTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000reltransfA, GP.ord, type="samples", color="calprotectinelevel") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000A<-subset_samples(dfrar1000, pays="Madagascar")
    dfrar1000prelA<-dfrar1000A %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000prelA, "PCoA", "jaccard")
    
    pdf("DalhousiefecesbraycurtisPCoAcalprotectinelevelTanapresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000prelA, GP.ord, type="samples", color="calprotectinelevel") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoATaxacalprotectinelevelTanapresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000prelA, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    ## anemie2 Tana
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000reltransfA, "PCoA", "bray")
    sample_data(dfrar1000reltransfA)$anemie2<-as.factor(sample_data(dfrar1000reltransfA)$anemie2)
    sample_data(dfrar1000)$anemie2<-as.factor(sample_data(dfrar1000)$anemie2)
    
    pdf("DalhousiefecesbraycurtisPCoATaxaanemie2Tana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000reltransfA, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoAanemie2Tana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000reltransfA, GP.ord, type="samples", color="anemie2") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000prelA<-dfrar1000 %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000prelA, "PCoA", "jaccard")
    
    pdf("DalhousiefecesbraycurtisPCoAanemie2presabsTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000prelA, GP.ord, type="samples", color="anemie2") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousiefecesbraycurtisPCoATaxaanemie2Tanapresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000prelA, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
  
    
    #### Plot Ordination with significant taxa as explanatory arrows on Genus level ####
    dfrar2genus<-tax_glom(dfrar1000rel, "rank7")
    OTU1 = as(otu_table(dfrar2genus), "matrix")
    ft.df.c<-as.data.frame(OTU1)
    
    TAX1 = as(tax_table(dfrar2genus), "matrix")
    ft.df.taxa<-as.data.frame(TAX1)
    
    META1= as(sample_data(dfrar2genus), "matrix")
    ft.df.meta<-as.data.frame(META1)
    
    nrow(ft.df.taxa)
    nrow(ft.df.c)
    
    test<-cbind(as.matrix(ft.df.taxa), as.matrix(ft.df.c))
    colnames(test)
    row.names(test)
    
    ##### Make genus-level OTU table with Sample IDs as row names 
    test2<-test[, -(1:6)]
    test2<-test2[, -(2:9)]
    View(test2)
    
    test2[, 1]<-as.character(test2[, 1])
    
    test2[, 1]<-make.unique(test2[, 1])
    
    test3<-t(test2)
    colnames(test3)<-test3[1, ]
    test3<-test3[-1, ]
    class(test3) <- "numeric"
    test3<-as.data.frame(test3)
    test3[is.na(test3)] <- 0 #convert NAs to zeros
    which(is.na(test3))==TRUE
    
    test3$sampleid<-row.names(test3)
    dim(test3)
    
    test4<-filter(test3, (rowSums(test3[, 1:19])!=0))
    genus.dat <- colnames(test4[, 1:19])
    
    dim(test3)
    dim(test4)
    dim(ft.df.meta)
   
    # PCoA is not included in vegan. 
    # We will use the ape package instead
    library(ape)
    PCOA <- pcoa(dist)
    
    # plot the eigenvalues and interpret
    barplot(PCOA$values$Relative_eig[1:10])
    biplot.pcoa(PCOA)
    biplot.pcoa(PCOA, test4[, 1:19], scale.=F, center=T)
    
    vectorPCOA<-cbind(PCOA$vectors[, 1:2])
    
    #### Include genus level Abundances in Ordination 
    env.fit <- envfit(vectorPCOA, test4[, 1:19], perm = 999) 
    
    # look at p-values of different taxa
    env.fit 
    env.fit.sig  <- as.data.frame(scores(env.fit, display = "vectors"))
    
    #### Extracting significant pvalues from envfit taxa
    #shortcutting ef$vectors
    A <- as.list(env.fit$vectors)
    #creating the dataframe
    pvals<-as.data.frame(A$pvals)
    arrows<-as.data.frame(A$arrows*sqrt(A$r))
    C<-cbind(arrows, pvals)
    #subset
    Cred<-subset(C,pvals<=0.001)
    Cred <- cbind(Cred, genus = rownames(Cred))
    
    #### Format taxa scores for plotting
    df_envfit<-scores(env.fit,display=c("vectors"))
    df_envfit<-df_envfit*vegan:::ordiArrowMul(df_envfit)
    df_envfit<-as.data.frame(df_envfit)
    df_envfit$genus <- rownames(df_envfit)
    df_envfit <- df_envfit %>% filter(genus %in% Cred$genus)
    
    #### Get PCOA scores for axis 1 and 2
    nrow(ft.df.meta)
    nrow(vectorPCOA) # we kicked-out 0 samples in the process. need to filter metadata file
    
    ft.scores <- as.data.frame(vectorPCOA)
    ft.scores$sampleid <- rownames(test4)
    
    # Paste metadata to PCoA scores for plotting
    ft.scores$pays <- ft.df.meta$pays
    ft.scores$stunted <- ft.df.meta$stunted
    
    # get distinct rows
    ft.scores.d <- distinct(ft.scores)
    
    
    #### Make ordination plot with Families as explanatory variables
    ft.shapes <- c(0,1,2,5,16)
    
    pdf(file="~/Desktop/GenuspaysDAL_explantory_arrowsrelabun.pdf",
        width = 8,
        height = 5)
    p<-ggplot() 
    p+geom_point(aes(ft.scores.d$Axis.1,ft.scores.d$Axis.2, colour=ft.scores.d$pays)) + 
      scale_color_manual(values = c("cyan3","mediumblue","darkorange3","goldenrod1")) +
      scale_shape_manual(values = ft.shapes) +
      geom_segment(aes(x = 0, y = 0, xend = df_envfit$Axis.1*0.5, yend = df_envfit$Axis.2*0.5), arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)+
      geom_text(aes(df_envfit$Axis.1*0.5, df_envfit$Axis.2*0.5, label = df_envfit$genus),color="#808080",alpha=0.5, vjust=2.5, hjust=.5) +
      theme_classic()
    dev.off()
    
    
    #### Make PCoA illustrations for rarefied dataset rel. abundance and presence absence on Species level ####
    dfrar1000Species<-tax_glom(dfrar1000, "rank8")
    dfrar1000Speciesrel<-transform_sample_counts(dfrar1000Species, function(x) x / sum(x) *100)
    dfrar1000Speciesreltransf<- transform(dfrar1000Speciesrel, 'log10')
    
    ## pays of origin
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000Speciesreltransf, "PCoA", "bray")
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxapays.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000Speciesreltransf, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoApays.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000Speciesreltransf, GP.ord, type="samples", color="pays") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000Speciesprel<-dfrar1000Species %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000Speciesprel, "PCoA", "jaccard")
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoApayspresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000Speciesprel, GP.ord, type="samples", color="pays") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxapayspresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000Speciesprel, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    ## calprotectinelevel
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000Speciesreltransf, "PCoA", "bray")
    sample_data(dfrar1000Speciesreltransf)$calprotectinelevel<-as.factor(sample_data(dfrar1000Speciesreltransf)$calprotectinelevel)
    sample_data(dfrar1000Species)$calprotectinelevel<-as.factor(sample_data(dfrar1000Species)$calprotectinelevel)
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxacalprotectinelevel.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000Speciesreltransf, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoAcalprotectinelevel.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000Speciesreltransf, GP.ord, type="samples", color="calprotectinelevel") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000Speciesprel<-dfrar1000Species %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000Speciesprel, "PCoA", "jaccard")
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoAcalprotectinelevelpresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000Speciesprel, GP.ord, type="samples", color="calprotectinelevel") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxacalprotectinelevelpresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000Speciesprel, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    ## anemie2
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000Speciesreltransf, "PCoA", "bray")
    sample_data(dfrar1000Speciesreltransf)$anemie2<-as.factor(sample_data(dfrar1000Speciesreltransf)$anemie2)
    sample_data(dfrar1000Species)$anemie2<-as.factor(sample_data(dfrar1000Species)$anemie2)
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxaanemie2.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000Speciesreltransf, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoAanemie2.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000Speciesreltransf, GP.ord, type="samples", color="anemie2") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000Speciesprel<-dfrar1000Species %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000Speciesprel, "PCoA", "jaccard")
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoAanemie2presabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000Speciesprel, GP.ord, type="samples", color="anemie2") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxaanemie2presabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000Speciesprel, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    ## calprotectinelevel Bangui
    dfrar1000SpeciesreltransfB<-subset_samples(dfrar1000Speciesreltransf, pays="RCA")
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000SpeciesreltransfB, "PCoA", "bray")
    sample_data(dfrar1000SpeciesreltransfB)$calprotectinelevel<-as.factor(sample_data(dfrar1000SpeciesreltransfB)$calprotectinelevel)
    sample_data(dfrar1000Species)$calprotectinelevel<-as.factor(sample_data(dfrar1000Species)$calprotectinelevel)
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxacalprotectinelevelBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000SpeciesreltransfB, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoAcalprotectinelevelBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000SpeciesreltransfB, GP.ord, type="samples", color="calprotectinelevel") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000SpeciesB<-subset_samples(dfrar1000Species, pays="RCA")
    dfrar1000SpeciesprelB<-dfrar1000SpeciesB %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000SpeciesprelB, "PCoA", "jaccard")
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoAcalprotectinelevelBanguipresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000SpeciesprelB, GP.ord, type="samples", color="calprotectinelevel") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxacalprotectinelevelBanguipresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000SpeciesprelB, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    ## anemie2 Bangui
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000SpeciesreltransfB, "PCoA", "bray")
    sample_data(dfrar1000SpeciesreltransfB)$anemie2<-as.factor(sample_data(dfrar1000SpeciesreltransfB)$anemie2)
    sample_data(dfrar1000Species)$anemie2<-as.factor(sample_data(dfrar1000Species)$anemie2)
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxaanemie2Bangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000SpeciesreltransfB, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoAanemie2Bangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000SpeciesreltransfB, GP.ord, type="samples", color="anemie2") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000SpeciesprelB<-dfrar1000Species %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000SpeciesprelB, "PCoA", "jaccard")
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoAanemie2presabsBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000SpeciesprelB, GP.ord, type="samples", color="anemie2") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxaanemie2Banguipresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000SpeciesprelB, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    ## calprotectinelevel Tana
    dfrar1000SpeciesreltransfA<-subset_samples(dfrar1000Speciesreltransf, pays="Madagascar")
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000SpeciesreltransfA, "PCoA", "bray")
    sample_data(dfrar1000SpeciesreltransfA)$calprotectinelevel<-as.factor(sample_data(dfrar1000SpeciesreltransfA)$calprotectinelevel)
    sample_data(dfrar1000Species)$calprotectinelevel<-as.factor(sample_data(dfrar1000Species)$calprotectinelevel)
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxacalprotectinelevelTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000SpeciesreltransfA, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoAcalprotectinelevelTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000SpeciesreltransfA, GP.ord, type="samples", color="calprotectinelevel") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000SpeciesA<-subset_samples(dfrar1000Species, pays="Madagascar")
    dfrar1000SpeciesprelA<-dfrar1000SpeciesA %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000SpeciesprelA, "PCoA", "jaccard")
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoAcalprotectinelevelTanapresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000SpeciesprelA, GP.ord, type="samples", color="calprotectinelevel") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxacalprotectinelevelTanapresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000SpeciesprelA, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    ## anemie2 Tana
    
    # make ordination using PCoA log 10 transformed
    GP.ord <- ordinate(dfrar1000SpeciesreltransfA, "PCoA", "bray")
    sample_data(dfrar1000SpeciesreltransfA)$anemie2<-as.factor(sample_data(dfrar1000SpeciesreltransfA)$anemie2)
    sample_data(dfrar1000Species)$anemie2<-as.factor(sample_data(dfrar1000Species)$anemie2)
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxaanemie2Tana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000SpeciesreltransfA, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoAanemie2Tana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000SpeciesreltransfA, GP.ord, type="samples", color="anemie2") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    # now also make one with prensence/absence data
    dfrar1000SpeciesprelA<-dfrar1000Species %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    GP.ord <- ordinate(dfrar1000SpeciesprelA, "PCoA", "jaccard")
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoAanemie2presabsTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p2 = plot_ordination(dfrar1000SpeciesprelA, GP.ord, type="samples", color="anemie2") 
    p2 + geom_point(size=1) + ggtitle("samples")
    dev.off()
    
    pdf("DalhousieSpeciesfecesbraycurtisPCoATaxaanemie2Tanapresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    p1 = plot_ordination(dfrar1000SpeciesprelA, GP.ord, type="taxa", color="rank3", title="taxa")
    print(p1)
    dev.off()
    
    
    
    #### Plot Ordination with significant taxa as explanatory arrows on Genus level pres abs ####
    dfrar2genus<-tax_glom(dfrar1000rel, "rank7")
    dfrar2genuspresabs<-transform_sample_counts(dfrar2genus, fun = prevalence)
    OTU1 = as(otu_table(dfrar2genuspresabs), "matrix")
    ft.df.c<-as.data.frame(t(OTU1))
    
    TAX1 = as(tax_table(dfrar2genuspresabs), "matrix")
    ft.df.taxa<-as.data.frame(TAX1)
    
    META1= as(sample_data(dfrar2genuspresabs), "matrix")
    ft.df.meta<-as.data.frame(META1)
    
    dim(ft.df.taxa)
    dim(ft.df.c)
    
    test<-cbind(as.matrix(ft.df.taxa), as.matrix(ft.df.c))
    colnames(test)
    row.names(test)
    
    ##### Make genus-level OTU table with Sample IDs as row names 
    test2<-test[, -(1:6)]
    test2<-test2[, -(2:9)]
    View(test2)
    
    test2[, 1]<-as.character(test2[, 1])
    
    test2[, 1]<-make.unique(test2[, 1])
    
    test3<-t(test2)
    colnames(test3)<-test3[1, ]
    test3<-test3[-1, ]
    class(test3) <- "numeric"
    test3<-as.data.frame(test3)
    test3[is.na(test3)] <- 0 #convert NAs to zeros
    which(is.na(test3))==TRUE
    
    test3$sampleid<-row.names(test3)
    dim(test3)
    
    test4<-filter(test3, (rowSums(test3[, 1:19])!=0))
    genus.dat <- colnames(test4[, 1:19])
    
    dim(test3)
    dim(test4)
    dim(ft.df.meta)
    
    #### use PCoA 
    # First step is to calculate a distance matrix. 
    # Here we use Jaccard distance metric
    
    dist <- vegdist(test4[, 1:19],  method = "jaccard")
    
    # PCoA is not included in vegan. 
    # We will use the ape package instead
    library(ape)
    PCOA <- pcoa(dist)
    
    # plot the eigenvalues and interpret
    barplot(PCOA$values$Relative_eig[1:10])
    biplot.pcoa(PCOA)
    biplot.pcoa(PCOA, test4[, 1:19], scale.=F, center=T)
    
    vectorPCOA<-cbind(PCOA$vectors[, 1:2])
    
    #### Include genus level Abundances in Ordination 
    env.fit <- envfit(vectorPCOA, test4[, 1:19], perm = 999) 
    
    # look at p-values of different taxa
    env.fit 
    env.fit.sig  <- as.data.frame(scores(env.fit, display = "vectors"))
    
    #### Extracting significant pvalues from envfit taxa
    #shortcutting ef$vectors
    A <- as.list(env.fit$vectors)
    #creating the dataframe
    pvals<-as.data.frame(A$pvals)
    arrows<-as.data.frame(A$arrows*sqrt(A$r))
    C<-cbind(arrows, pvals)
    #subset
    Cred<-subset(C,pvals<=0.001)
    Cred <- cbind(Cred, genus = rownames(Cred))
    
    #### Format taxa scores for plotting
    df_envfit<-scores(env.fit,display=c("vectors"))
    df_envfit<-df_envfit*vegan:::ordiArrowMul(df_envfit)
    df_envfit<-as.data.frame(df_envfit)
    df_envfit$genus <- rownames(df_envfit)
    df_envfit <- df_envfit %>% filter(genus %in% Cred$genus)
    
    #### Get PCOA scores for axis 1 and 2
    nrow(ft.df.meta)
    nrow(vectorPCOA) # we kicked-out 0 samples in the process. need to filter metadata file
    
    ft.scores <- as.data.frame(vectorPCOA)
    ft.scores$sampleid <- rownames(test4)
    
    # Paste metadata to PCoA scores for plotting
    ft.scores$pays <- ft.df.meta$pays
    ft.scores$stunted <- ft.df.meta$stunted
    
    # get distinct rows
    ft.scores.d <- distinct(ft.scores)
    
    
    #### Make ordination plot with genera as explanatory variables
    ft.shapes <- c(0,1,2,5,16)
    
    pdf(file="~/Desktop/GenuspaysDAL_explantory_arrowspresabs.pdf",
        width = 8,
        height = 5)
    p<-ggplot() 
    p+geom_point(aes(ft.scores.d$Axis.1,ft.scores.d$Axis.2, colour=ft.scores.d$pays)) + 
      scale_color_manual(values = c("cyan3","mediumblue","darkorange3","goldenrod1")) +
      scale_shape_manual(values = ft.shapes) +
      geom_segment(aes(x = 0, y = 0, xend = df_envfit$Axis.1*0.2, yend = df_envfit$Axis.2*0.2), arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)+
      geom_text(aes(df_envfit$Axis.1*0.2, df_envfit$Axis.2*0.2, label = df_envfit$genus),color="#808080",alpha=0.5, vjust=2.5, hjust=.5) +
      theme_classic()
    dev.off()
    
    
    
    #### Plot Ordination with significant taxa as explanatory arrows on Species level ####
    dfrar2species<-tax_glom(dfrar1000rel, "rank8")
    OTU1 = as(otu_table(dfrar2species), "matrix")
    ft.df.c<-as.data.frame(t(OTU1))
    
    TAX1 = as(tax_table(dfrar2species), "matrix")
    ft.df.taxa<-as.data.frame(TAX1)
    
    META1= as(sample_data(dfrar2species), "matrix")
    ft.df.meta<-as.data.frame(META1)
    
    nrow(ft.df.taxa)
    nrow(ft.df.c)
    
    test<-cbind(as.matrix(ft.df.taxa), as.matrix(ft.df.c))
    colnames(test)
    row.names(test)
    
    ##### Make species-level OTU table with Sample IDs as row names 
    test2<-test[, -(1:7)]
    test2<-test2[, -(2:8)]
    View(test2)
    
    test2[, 1]<-as.character(test2[, 1])
    
    test2[, 1]<-make.unique(test2[, 1])
    
    test3<-t(test2)
    colnames(test3)<-test3[1, ]
    test3<-test3[-1, ]
    class(test3) <- "numeric"
    test3<-as.data.frame(test3)
    test3[is.na(test3)] <- 0 #convert NAs to zeros
    which(is.na(test3))==TRUE
    
    test3$sampleid<-row.names(test3)
    dim(test3)
    
    test4<-filter(test3, (rowSums(test3[, 1:25])!=0))
    species.dat <- colnames(test4[, 1:25])
    
    dim(test3)
    dim(test4)
    dim(ft.df.meta)
    
    #### NMDS ordination 
    ft.bc <- metaMDS(test4[, 1:25], distance="bray")
    
    #### does not find convergence!
    
    #### use PCoA instead
    # First step is to calculate a distance matrix. 
    # Here we use Bray-Curtis distance metric
    test4_log<-transform(test4[, 1:25], 'log10')
    dist <- vegdist(test4_log[, 1:25],  method = "bray")
    
    # PCoA is not included in vegan. 
    # We will use the ape package instead
    library(ape)
    PCOA <- pcoa(dist)
    
    # plot the eigenvalues and interpret
    barplot(PCOA$values$Relative_eig[1:10])
    biplot.pcoa(PCOA)
    biplot.pcoa(PCOA, test4[, 1:25], scale.=F, center=T)
    
    vectorPCOA<-cbind(PCOA$vectors[, 1:2])
    
    #### Include species level Abundances in Ordination 
    env.fit <- envfit(vectorPCOA, test4[, 1:25], perm = 999) 
    
    # look at p-values of different taxa
    env.fit 
    env.fit.sig  <- as.data.frame(scores(env.fit, display = "vectors"))
    
    #### Extracting significant pvalues from envfit taxa
    #shortcutting ef$vectors
    A <- as.list(env.fit$vectors)
    #creating the dataframe
    pvals<-as.data.frame(A$pvals)
    arrows<-as.data.frame(A$arrows*sqrt(A$r))
    C<-cbind(arrows, pvals)
    #subset
    Cred<-subset(C,pvals<=0.001)
    Cred <- cbind(Cred, species = rownames(Cred))
    
    #### Format taxa scores for plotting
    df_envfit<-scores(env.fit,display=c("vectors"))
    df_envfit<-df_envfit*vegan:::ordiArrowMul(df_envfit)
    df_envfit<-as.data.frame(df_envfit)
    df_envfit$species <- rownames(df_envfit)
    df_envfit <- df_envfit %>% filter(species %in% Cred$species)
    
    #### Get PCOA scores for axis 1 and 2
    nrow(ft.df.meta)
    nrow(vectorPCOA) # we kicked-out 0 samples in the process. need to filter metadata file
    
    ft.scores <- as.data.frame(vectorPCOA)
    ft.scores$sampleid <- rownames(test4)
    
    # Paste metadata to PCoA scores for plotting
    ft.scores$pays <- ft.df.meta$pays
    ft.scores$stunted <- ft.df.meta$stunted
    
    # get distinct rows
    ft.scores.d <- distinct(ft.scores)
    
    
    #### Make ordination plot with Families as explanatory variables
    ft.shapes <- c(0,1,2,5,16)
    
    pdf(file="~/Desktop/speciespaysDAL_explantory_arrowsrelabun.pdf",
        width = 8,
        height = 5)
    p<-ggplot() 
    p+geom_point(aes(ft.scores.d$Axis.1,ft.scores.d$Axis.2, colour=ft.scores.d$pays)) + 
      scale_color_manual(values = c("red","blue")) +
      scale_shape_manual(values = ft.shapes) +
      geom_segment(aes(x = 0, y = 0, xend = df_envfit$Axis.1*0.003, yend = df_envfit$Axis.2*0.003), arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)+
      geom_text(aes(df_envfit$Axis.1*0.003, df_envfit$Axis.2*0.003, label = df_envfit$species),color="#808080",alpha=0.5, vjust=2.5, hjust=.5) +
      theme_classic()
    dev.off()
    
    
    
    
    #### Plot Ordination with significant taxa as explanatory arrows on species level pres abs ####
    dfrar2species<-tax_glom(dfrar1000rel, "rank8")
    dfrar2speciespresabs<-transform_sample_counts(dfrar2species, fun = prevalence)
    OTU1 = as(otu_table(dfrar2speciespresabs), "matrix")
    ft.df.c<-as.data.frame(t(OTU1))
    
    TAX1 = as(tax_table(dfrar2speciespresabs), "matrix")
    ft.df.taxa<-as.data.frame(TAX1)
    
    META1= as(sample_data(dfrar2speciespresabs), "matrix")
    ft.df.meta<-as.data.frame(META1)
    
    nrow(ft.df.taxa)
    nrow(ft.df.c)
    
    test<-cbind(as.matrix(ft.df.taxa), as.matrix(ft.df.c))
    colnames(test)
    row.names(test)
    
    ##### Make species-level OTU table with Sample IDs as row names 
    test2<-test[, -(1:7)]
    test2<-test2[, -(2:8)]
    View(test2)
    
    test2[, 1]<-as.character(test2[, 1])
    
    test2[, 1]<-make.unique(test2[, 1])
    
    test3<-t(test2)
    colnames(test3)<-test3[1, ]
    test3<-test3[-1, ]
    class(test3) <- "numeric"
    test3<-as.data.frame(test3)
    test3[is.na(test3)] <- 0 #convert NAs to zeros
    which(is.na(test3))==TRUE
    
    test3$sampleid<-row.names(test3)
    dim(test3)
    
    test4<-filter(test3, (rowSums(test3[, 1:25])!=0))
    species.dat <- colnames(test4[, 1:25])
    
    dim(test3)
    dim(test4)
    dim(ft.df.meta)
    
    #### NMDS ordination 
    ft.bc <- metaMDS(test4[, 1:25], distance="bray")
    
    #### does not find convergence!
    
    #### use PCoA instead
    # First step is to calculate a distance matrix. 
    # Here we use Bray-Curtis distance metric
    test4_log<-transform(test4[, 1:25], 'log10')
    dist <- vegdist(test4_log[, 1:25],  method = "bray")
    
    # PCoA is not included in vegan. 
    # We will use the ape package instead
    library(ape)
    PCOA <- pcoa(dist)
    
    # plot the eigenvalues and interpret
    barplot(PCOA$values$Relative_eig[1:10])
    biplot.pcoa(PCOA)
    biplot.pcoa(PCOA, test4[, 1:25], scale.=F, center=T)
    
    vectorPCOA<-cbind(PCOA$vectors[, 1:2])
    
    #### Include species level Abundances in Ordination 
    env.fit <- envfit(vectorPCOA, test4[, 1:25], perm = 999) 
    
    # look at p-values of different taxa
    env.fit 
    env.fit.sig  <- as.data.frame(scores(env.fit, display = "vectors"))
    
    #### Extracting significant pvalues from envfit taxa
    #shortcutting ef$vectors
    A <- as.list(env.fit$vectors)
    #creating the dataframe
    pvals<-as.data.frame(A$pvals)
    arrows<-as.data.frame(A$arrows*sqrt(A$r))
    C<-cbind(arrows, pvals)
    #subset
    Cred<-subset(C,pvals<=0.001)
    Cred <- cbind(Cred, species = rownames(Cred))
    
    #### Format taxa scores for plotting
    df_envfit<-scores(env.fit,display=c("vectors"))
    df_envfit<-df_envfit*vegan:::ordiArrowMul(df_envfit)
    df_envfit<-as.data.frame(df_envfit)
    df_envfit$species <- rownames(df_envfit)
    df_envfit <- df_envfit %>% filter(species %in% Cred$species)
    
    #### Get PCOA scores for axis 1 and 2
    nrow(ft.df.meta)
    nrow(vectorPCOA) # we kicked-out 0 samples in the process. need to filter metadata file
    
    ft.scores <- as.data.frame(vectorPCOA)
    ft.scores$sampleid <- rownames(test4)
    
    # Paste metadata to PCoA scores for plotting
    ft.scores$pays <- ft.df.meta$pays
    ft.scores$stunted <- ft.df.meta$stunted
    
    # get distinct rows
    ft.scores.d <- distinct(ft.scores)
    
    
    #### Make ordination plot with Families as explanatory variables
    ft.shapes <- c(0,1,2,5,16)
    
    pdf(file="~/Desktop/speciespaysDAL_explantory_arrowspresabs.pdf",
        width = 8,
        height = 5)
    p<-ggplot() 
    p+geom_point(aes(ft.scores.d$Axis.1,ft.scores.d$Axis.2, colour=ft.scores.d$pays)) + 
      scale_color_manual(values = c("red","blue")) +
      scale_shape_manual(values = ft.shapes) +
      geom_segment(aes(x = 0, y = 0, xend = df_envfit$Axis.1*0.35, yend = df_envfit$Axis.2*0.35), arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)+
      geom_text(aes(df_envfit$Axis.1*0.35, df_envfit$Axis.2*0.35, label = df_envfit$species),color="#808080",alpha=0.5, vjust=2.5, hjust=.5) +
      theme_classic()
    dev.off()
    
    
    
    
#### Generate core eukaryomes for given taxon level in feces on filtered dataset ####    
    ## Anything conserved on ASV level?
    eukaryome.rel <- microbiome::transform(dfcleanfeces, "compositional")
    eukaryome.rel
    View(head(otu_table(eukaryome.rel)))
    
    core.taxa <- core(eukaryome.rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
    # no taxon is present in 90% of all samples
    
    core.taxa <- core(eukaryome.rel, detection = 0.0001, prevalence = 50/100, include.lowest = TRUE) 
    # no taxon is present in 50% of all samples
    
    core.taxa <- core(eukaryome.rel, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE) 
    core.taxa # 1 Taxon
    View(tax_table(core.taxa)) # Blastocystis ST3
    
    core.taxa <- core(eukaryome.rel, detection = 0.00000001, prevalence = 10/100, include.lowest = TRUE)
    core.taxa # 16 taxa 
    View(tax_table(core.taxa)) # different Blastocystis, Ascaris, Saccharomyces, Entamoeba, Giardia
    
    ## ASV level stratified by country?
    eukaryome.rel.B=subset_samples(eukaryome.rel, pays=="RCA")
    
    core.taxa <- core(eukaryome.rel.B, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
    # no taxon is present in 90% of all samples
    
    core.taxa <- core(eukaryome.rel.B, detection = 0.00001, prevalence = 50/100, include.lowest = TRUE) 
    # 0 are conserved
   
    core.taxa <- core(eukaryome.rel.B, detection = 0.00000001, prevalence = 50/100, include.lowest = TRUE) 
    ##0 taxa  
    
    core.taxa <- core(eukaryome.rel.B, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE)
    core.taxa # 1 taxon
    View(tax_table(core.taxa)) ##Blasto ST3 
    
    core.taxa <- core(eukaryome.rel.B, detection = 0.00000001, prevalence = 10/100, include.lowest = TRUE)
    core.taxa # 14 taxa 
    View(tax_table(core.taxa))
    
    eukaryome.rel.A=subset_samples(eukaryome.rel, pays=="Madagascar")
    
    core.taxa <- core(eukaryome.rel.A, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
    # no taxon is present in 90% of all samples
    
    core.taxa <- core(eukaryome.rel.A, detection = 0.00001, prevalence = 50/100, include.lowest = TRUE) 
    # no taxon is present in 50% of all samples
    
    core.taxa <- core(eukaryome.rel.A, detection = 0.00000001, prevalence = 50/100, include.lowest = TRUE) 
    # no taxon is present in 50% of all samples
    
    core.taxa <- core(eukaryome.rel.A, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE)
    core.taxa # 2 taxa   
    View(tax_table(core.taxa)) # Blasto ST3
    
    core.taxa <- core(eukaryome.rel.A, detection = 0.00000001, prevalence = 10/100, include.lowest = TRUE)
    core.taxa # 25 taxa 
    View(tax_table(core.taxa))
    
    ## Anything conserved on rank 4?
    dfcleanfeces_rank4<-tax_glom(dfcleanfeces, "rank4")
    eukaryomerank4.rel <- microbiome::transform(dfcleanfeces_rank4, "compositional")
    core.taxa <- core(eukaryomerank4.rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank4.rel, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank4.rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
    core.taxa
    View(tax_table(core.taxa)) #Incerta sedis
    
    eukaryomerank4.rel.A=subset_samples(eukaryomerank4.rel, pays=="Madagascar")
    eukaryomerank4.rel.B=subset_samples(eukaryomerank4.rel, pays=="RCA")
    
    core.taxa <- core(eukaryomerank4.rel.A, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank4.rel.A, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank4.rel.A, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank4.rel.A, detection = 0.0001, prevalence = 65/100, include.lowest = TRUE)
    core.taxa
    View(tax_table(core.taxa)) #"Incerta sedis
    
    core.taxa <- core(eukaryomerank4.rel.B, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank4.rel.B, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank4.rel.B, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE) 
    core.taxa
    View(tax_table(core.taxa)) #Incerta sedis
    
    ## Anything conserved on rank 8?
    dfcleanfeces_rank8<-tax_glom(dfcleanfeces, "rank8")
    eukaryomerank8.rel <- microbiome::transform(dfcleanfeces_rank8, "compositional")
    core.taxa <- core(eukaryomerank8.rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank8.rel, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank8.rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
    core.taxa <- core(eukaryomerank8.rel, detection = 0.0001, prevalence = 50/100, include.lowest = TRUE)
    core.taxa
    View(tax_table(core.taxa)) #1 species are conserved in 50 percent of samples; Blasto ST3
    
    core.taxa <- core(eukaryomerank8.rel, detection = 0.0001, prevalence = 25/100, include.lowest = TRUE)
    core.taxa
    View(tax_table(core.taxa)) #4 species are conserved in 25 percent of samples
    
    core.taxa <- core(eukaryomerank8.rel, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE)
    core.taxa
    tax_table(core.taxa) #4 species are conserved in 25 percent of samples
    
    core.taxa <- core(eukaryomerank8.rel, detection = 0.0001, prevalence = 10/100, include.lowest = TRUE) 
    core.taxa
    View(tax_table(core.taxa)) #12 taxa
    
       
    
    eukaryomerank8.rel.A=subset_samples(eukaryomerank8.rel, pays=="Madagascar")
    eukaryomerank8.rel.B=subset_samples(eukaryomerank8.rel, pays=="RCA")
    
    core.taxa <- core(eukaryomerank8.rel.A, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank8.rel.A, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank8.rel.A, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank8.rel.A, detection = 0.0001, prevalence = 65/100, include.lowest = TRUE)
    core.taxa <- core(eukaryomerank8.rel.A, detection = 0.0001, prevalence = 50/100, include.lowest = TRUE)
    core.taxa
    tax_table(core.taxa) #Blasto ST3
    
    core.taxa <- core(eukaryomerank8.rel.A, detection = 0.0001, prevalence = 10/100, include.lowest = TRUE) 
    core.taxa
    View(tax_table(core.taxa)) #12 taxa
    
    core.taxa <- core(eukaryomerank8.rel.B, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank8.rel.B, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank8.rel.B, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE) #3 fungal species
    core.taxa <- core(eukaryomerank8.rel.B, detection = 0.0001, prevalence = 50/100, include.lowest = TRUE)  
    core.taxa <- core(eukaryomerank8.rel.B, detection = 0.0001, prevalence = 30/100, include.lowest = TRUE)  
    core.taxa # 3 taxa, Blasto ST3 and Blasto ST2 and Blasto ST1
    View(tax_table(core.taxa) )
    
    core.taxa <- core(eukaryomerank8.rel.B, detection = 0.0001, prevalence = 10/100, include.lowest = TRUE) 
    core.taxa
    View(tax_table(core.taxa)) #12 taxa
    
    ## Anything conserved on rank 7?
    dfcleanfeces_rank7<-tax_glom(dfcleanfeces, "rank7")
    eukaryomerank7.rel <- microbiome::transform(dfcleanfeces_rank7, "compositional")
    core.taxa <- core(eukaryomerank7.rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank7.rel, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank7.rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
    core.taxa <- core(eukaryomerank7.rel, detection = 0.00000001, prevalence = 75/100, include.lowest = TRUE)
    core.taxa # 1 taxon are shared
    View(tax_table(core.taxa)) # Blastocystis
    core.taxa <- core(eukaryomerank7.rel, detection = 0.00000001, prevalence = 50/100, include.lowest = TRUE)
    core.taxa # 2 taxa are shared
    View(tax_table(core.taxa)) # Blastocystis and Entamoeba
    
    core.taxa <- core(eukaryomerank7.rel, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE)
    core.taxa # 3 taxa are shared
    View(tax_table(core.taxa)) # Entamoeba, Blastocystis, Saccharomyces
    
    eukaryomerank7.rel.A=subset_samples(eukaryomerank7.rel, pays=="Madagascar")
    eukaryomerank7.rel.B=subset_samples(eukaryomerank7.rel, pays=="RCA")
    
    core.taxa <- core(eukaryomerank7.rel.A, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank7.rel.A, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank7.rel.A, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank7.rel.A, detection = 0.0001, prevalence = 65/100, include.lowest = TRUE)
    core.taxa #1 taxon
    core.taxa <- core(eukaryomerank7.rel.A, detection = 0.00000001, prevalence = 50/100, include.lowest = TRUE)
    core.taxa # 2 taxa
    View(tax_table(core.taxa))
    
    core.taxa <- core(eukaryomerank7.rel.A, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE)
    core.taxa # 3 taxa, Blasto, Entamoeba, Saccharomyces
    View(tax_table(core.taxa))
    
    core.taxa <- core(eukaryomerank7.rel.B, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank7.rel.B, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank7.rel.B, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE) 
    core.taxa <- core(eukaryomerank7.rel.B, detection = 0.0001, prevalence = 50/100, include.lowest = TRUE) 
    core.taxa
    View(tax_table(core.taxa)) #Blasto
    
    core.taxa <- core(eukaryomerank7.rel.B, detection = 0.00000001, prevalence = 50/100, include.lowest = TRUE)
    core.taxa # 1 taxon
    tax_table(core.taxa) # Blasto
    
    core.taxa <- core(eukaryomerank7.rel.B, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE)
    core.taxa # 3 taxa
    tax_table(core.taxa) # Blasto, Saccharomyces,  Entamoeba
    
    
    #### Generate tables of rel. abundance by different categories ####    
    eukaryome.rel <- microbiome::transform(dfcleanfeces, "compositional")
    eukaryome.rel<-subset_samples(eukaryome.rel, pays!="")
    
    # Plot Phylum content
    
    eukaryome.rel.p <- tax_glom(eukaryome.rel, "rank2")
    
    eukaryome.rel.p = transform_sample_counts(eukaryome.rel.p, function(x) 100 * x/sum(x))
    eukaryome.rel.p = prune_samples(sample_sums(eukaryome.rel.p)>0, eukaryome.rel.p)
    
       
    eukaryome.rel.p= eukaryome.rel.p  %>%
      psmelt() %>%                                         # Melt to long format
      arrange(rank2) 
    
    pdf("AbundancetablePhylumpersampleDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 40, # define plot width and height. completely up to user.
        height = 10) 
    ggplot(eukaryome.rel.p, aes(x = Sample, y = Abundance, fill = rank2)) + 
      facet_wrap(pays~ageyears, scales = "free_x") +
      theme(strip.text = element_text(colour = 'black', size=18), axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.title.y = element_text(size=18), title =element_text(size=25, face='bold')) +
      geom_bar(stat = "identity") +
      #scale_fill_manual(values = phylum_colors) +
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
    
    pdf("AbundancetablePhylumpersamplepaysstuntedDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 40, # define plot width and height. completely up to user.
        height = 10) 
    ggplot(eukaryome.rel.p, aes(x = Sample, y = Abundance, fill = rank2)) + 
      facet_wrap(pays~stunted, scales = "free_x") +
      theme(strip.text = element_text(colour = 'black', size=18), axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.title.y = element_text(size=18), title =element_text(size=25, face='bold')) +
      geom_bar(stat = "identity") +
      #scale_fill_manual(values = phylum_colors) +
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
    
    
    # Plot rank4 content
    
    eukaryome.rel.4 <- tax_glom(eukaryome.rel, "rank4")
    
    eukaryome.rel.4 = transform_sample_counts(eukaryome.rel.4, function(x) 100 * x/sum(x))
    
    
    eukaryome.rel.4= eukaryome.rel.4  %>%
      psmelt() %>%                                         # Melt to long format
      arrange(rank4) 
    
    pdf("Abundancetablerank4persampleDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 40, # define plot width and height. completely up to user.
        height = 10) 
    ggplot(eukaryome.rel.4, aes(x = Sample, y = Abundance, fill = rank4)) + 
      facet_wrap(pays~ageyears, scales = "free_x") +
      theme(strip.text = element_text(colour = 'black', size=18), axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.title.y = element_text(size=18), title =element_text(size=25, face='bold')) +
      geom_bar(stat = "identity") +
      #scale_fill_manual(values = phylum_colors) +
      scale_x_discrete(
        drop = FALSE
      ) +
      # Remove x axis title
      theme(axis.title.x = element_blank()) + 
      theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
      ylab("Relative Abundance") +
      ggtitle("rank4 Composition of Afribiota samples") 
    dev.off()
    
    pdf("Abundancetablerank4persamplepaysstuntedDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 40, # define plot width and height. completely up to user.
        height = 10) 
    ggplot(eukaryome.rel.4, aes(x = Sample, y = Abundance, fill = rank4)) + 
      facet_wrap(pays~stunted, scales = "free_x") +
      theme(strip.text = element_text(colour = 'black', size=18), axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.title.y = element_text(size=18), title =element_text(size=25, face='bold')) +
      geom_bar(stat = "identity") +
      #scale_fill_manual(values = phylum_colors) +
      scale_x_discrete(
        drop = FALSE
      ) +
      # Remove x axis title
      theme(axis.title.x = element_blank()) + 
      theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
      ylab("Relative Abundance") +
      ggtitle("rank4 Composition of Afribiota samples") 
    dev.off()
    
    # plot rank4 content merged by country
    
    eukaryome.rel.p = merge_samples(eukaryome.rel, "pays")
    
    eukaryome.rel.p_p <- tax_glom(eukaryome.rel.p, "rank4")
    
    eukaryome.rel.p_p = transform_sample_counts(eukaryome.rel.p_p, function(x) 100 * x/sum(x))
    
    eukaryome.rel.p_p_s= eukaryome.rel.p_p  %>%
      psmelt() %>%                                         # Melt to long format
      arrange(rank4) 
    
    class_colors <- c(
      "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
      "#AD6F3B", "#673770","black", "#652926", "#C84248", 
      "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
    )
    
    pdf("Abundancetablerank4paysnopruningDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 11, # define plot width and height. completely up to user.
        height = 8) 
    ggplot(eukaryome.rel.p_p_s, aes(x =Sample, y = Abundance, fill = rank4)) + 
      theme_bw() +
      theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
      theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
      geom_bar(stat = "identity", width = 1.0) +
      scale_y_continuous(expand = c(0.01,0.01)) +
      theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
      ylab("Relative Abundance") +
      xlab("Pays") + ggtitle ("Relative Abundance of phyla according to country of origin")
    dev.off()
    
    # plot rank2 content merged by country
    
    eukaryome.rel.p = merge_samples(eukaryome.rel, "pays")
    
    eukaryome.rel.p_p <- tax_glom(eukaryome.rel.p, "rank2")
    
    eukaryome.rel.p_p = transform_sample_counts(eukaryome.rel.p_p, function(x) 100 * x/sum(x))
    
    eukaryome.rel.p_p_s= eukaryome.rel.p_p  %>%
      psmelt() %>%                                         # Melt to long format
      arrange(rank2) 
    
    class_colors <- c(
      "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
      "#AD6F3B", "#673770","black", "#652926", "#C84248", 
      "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
    )
    
    pdf("AbundancetablePhylumpaysnopruningDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 11, # define plot width and height. completely up to user.
        height = 8) 
    ggplot(eukaryome.rel.p_p_s, aes(x =Sample, y = Abundance, fill = rank2)) + 
      theme_bw() +
      theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
      theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
      geom_bar(stat = "identity", width = 1.0) +
      scale_y_continuous(expand = c(0.01,0.01)) +
      theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
      ylab("Relative Abundance") +
      xlab("Pays") + ggtitle ("Relative Abundance of phyla according to country of origin")
    dev.off()
    
    # plot rank3 content merged by country
    
    eukaryome.rel.p = merge_samples(eukaryome.rel, "pays")
    
    eukaryome.rel.p_p <- tax_glom(eukaryome.rel.p, "rank3")
    
    eukaryome.rel.p_p = transform_sample_counts(eukaryome.rel.p_p, function(x) 100 * x/sum(x))
    
    eukaryome.rel.p_p_s= eukaryome.rel.p_p  %>%
      psmelt() %>%                                         # Melt to long format
      arrange(rank3) 
    
    class_colors <- c(
      "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
      "#AD6F3B", "#673770","black", "#652926", "#C84248", 
      "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
    )
    
    pdf("Abundancetablerank3paysnopruningDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 11, # define plot width and height. completely up to user.
        height = 8) 
    ggplot(eukaryome.rel.p_p_s, aes(x =Sample, y = Abundance, fill = rank3)) + 
      theme_bw() +
      theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
      theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
      geom_bar(stat = "identity", width = 1.0) +
      scale_y_continuous(expand = c(0.01,0.01)) +
      theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
      ylab("Relative Abundance") +
      xlab("Pays") + ggtitle ("Relative Abundance of rank3 according to country of origin")
    dev.off()
    
    
    
    #### Generate tables of rel. abundance for origin of taxa  ####   
    
    # Plot origin of taxa in feces
    eukaryome.rel <- microbiome::transform(dfcleanfeces, "compositional")
    eukaryome.rel.get <- tax_glom(eukaryome.rel, "gut_environmental_diet")
    
    eukaryome.rel.get.red = transform_sample_counts(eukaryome.rel.get, function(x) 100 * x/sum(x))
    eukaryome.rel.get.red = prune_samples(sample_sums(eukaryome.rel.get.red)>0, eukaryome.rel.get.red)
    
    eukaryome.rel.get.red<-subset_samples(eukaryome.rel.get.red, pays!="")
    eukaryome.rel.get.red<-subset_samples(eukaryome.rel.get.red, stunted!="")
    
    eukaryome.rel.get.red= eukaryome.rel.get.red  %>%
      psmelt() %>%                                         # Melt to long format
      arrange(gut_environmental_diet) 
    
    pdf("AbundancetableOriginpersamplefecesDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 40, # define plot width and height. completely up to user.
        height = 10) 
    ggplot(eukaryome.rel.get.red, aes(x = Sample, y = Abundance, fill = gut_environmental_diet)) + 
      facet_wrap(pays~stunted, scales = "free_x") +
      theme(strip.text = element_text(colour = 'black', size=18), axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.title.y = element_text(size=18), title =element_text(size=25, face='bold')) +
      geom_bar(stat = "identity") +
      #scale_fill_manual(values = gut_environmental_diet_colors) +
      scale_x_discrete(
        drop = FALSE
      ) +
      # Remove x axis title
      theme(axis.title.x = element_blank()) + 
      theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
      ylab("Relative Abundance") +
      ggtitle("Origin of taxa in Afribiota fecal samples") 
    dev.off()
    
    
    # Plot origin of taxa in all sample types
    
    eukaryome.rel.get <- tax_glom(eukaryome.rel, "gut_environmental_diet")
    eukaryome.rel.get<-subset_samples(eukaryome.rel.get, (SampleType=="gastric" | SampleType=="duodenal" | SampleType=="feces"))
    
    eukaryome.rel.get = transform_sample_counts(eukaryome.rel.get, function(x) 100 * x/sum(x))
    
    
    eukaryome.rel.get.red<-subset_samples(eukaryome.rel.get, pays!="")
    eukaryome.rel.get.red<-subset_samples(eukaryome.rel.get.red, stunted!="")
    
    eukaryome.rel.get.red= eukaryome.rel.get.red  %>%
      psmelt() %>%                                         # Melt to long format
      arrange(gut_environmental_diet) 
    
    pdf("AbundancetableOriginpersampletypeDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 40, # define plot width and height. completely up to user.
        height = 10) 
    ggplot(eukaryome.rel.get.red, aes(x = Sample, y = Abundance, fill = gut_environmental_diet)) + 
      facet_wrap(pays~SampleType, scales = "free_x") +
      theme(strip.text = element_text(colour = 'black', size=18), axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.title.y = element_text(size=18), title =element_text(size=25, face='bold')) +
      geom_bar(stat = "identity") +
      #scale_fill_manual(values = gut_environmental_diet_colors) +
      scale_x_discrete(
        drop = FALSE
      ) +
      # Remove x axis title
      theme(axis.title.x = element_blank()) + 
      theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
      ylab("Relative Abundance") +
      ggtitle("Origin of taxa in Afribiota samples") 
    dev.off()
    
    # plot origin content merged by country for feces
    sample_data(eukaryome.rel)$pays<-as.factor(sample_data(eukaryome.rel)$pays)
    eukaryome.rel<-subset_samples(eukaryome.rel, pays!="")
    which(is.na(sample_data(eukaryome.rel)$pays))
    
    eukaryome.rel.pays = merge_samples(eukaryome.rel, "pays")
    
    eukaryome.rel.pays_p <- tax_glom(eukaryome.rel.pays, "gut_environmental_diet")
    
    eukaryome.rel.pays_p = transform_sample_counts(eukaryome.rel.pays_p, function(x) 100 * x/sum(x))
    
    eukaryome.rel.pays_p_s= eukaryome.rel.pays_p  %>%
      psmelt() %>%                                         # Melt to long format
      arrange(gut_environmental_diet) 
    
    class_colors <- c(
      "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
      "#AD6F3B", "#673770","black", "#652926", "#C84248", 
      "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
    )
    
    pdf("Abundancetablegut_environmental_dietpaysnopruningDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 11, # define plot width and height. completely up to user.
        height = 8) 
    ggplot(eukaryome.rel.pays_p_s, aes(x =Sample, y = Abundance, fill = gut_environmental_diet)) + 
      theme_bw() +
      theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
      theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
      geom_bar(stat = "identity", width = 1.0) +
      scale_y_continuous(expand = c(0.01,0.01)) +
      theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
      ylab("Relative Abundance") +
      xlab("Pays") + ggtitle ("Relative Abundance according to origin of taxa in feces by country of origin")
    dev.off()
    
    
    
    # plot origin content merged by stunting for feces
    sample_data(eukaryome.rel)$stunted<-as.factor(sample_data(eukaryome.rel)$stunted)
    eukaryome.rel<-subset_samples(eukaryome.rel, stunted!="")
    which(is.na(sample_data(eukaryome.rel)$stunted))
    
    eukaryome.rel.stunted = merge_samples(eukaryome.rel, "stunted")
    
    eukaryome.rel.stunted_p <- tax_glom(eukaryome.rel.stunted, "gut_environmental_diet")
    
    eukaryome.rel.stunted_p = transform_sample_counts(eukaryome.rel.stunted_p, function(x) 100 * x/sum(x))
    
    eukaryome.rel.stunted_p_s= eukaryome.rel.stunted_p  %>%
      psmelt() %>%                                         # Melt to long format
      arrange(gut_environmental_diet) 
    
    class_colors <- c(
      "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
      "#AD6F3B", "#673770","black", "#652926", "#C84248", 
      "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
    )
    
    pdf("Abundancetablegut_environmental_dietstuntednopruningDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 13, # define plot width and height. completely up to user.
        height = 8) 
    ggplot(eukaryome.rel.stunted_p_s, aes(x =Sample, y = Abundance, fill = gut_environmental_diet)) + 
      theme_bw() +
      theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
      theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
      geom_bar(stat = "identity", width = 1.0) +
      scale_y_continuous(expand = c(0.01,0.01)) +
      theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
      ylab("Relative Abundance") +
      xlab("stunted") + ggtitle ("Relative Abundance according to stunting status in feces by country of origin")
    dev.off()
    
    
    
    #### Generate tables of rel. abundance by different categories on custom ranks ####    
    dffiltered_feces_custom= tax_glom(dfcleanfeces, "Rank_taxa_plots")
    tax_table(dffiltered_feces_custom)<-tax_table(dffiltered_feces_custom)[, 1:10]
    dffiltered_feces_custom= tax_glom(dffiltered_feces_custom, "Rank_taxa_plots")
    
    eukaryome.rel <- microbiome::transform(dffiltered_feces_custom, "compositional")
    eukaryome.rel<-subset_samples(eukaryome.rel, pays!="")
    
    # Plot customtaxa content
    eukaryome.rel.p <- tax_glom(eukaryome.rel, "Rank_taxa_plots")
    
    eukaryome.rel.p = transform_sample_counts(eukaryome.rel.p, function(x) 100 * x/sum(x))
    eukaryome.rel.p = prune_samples(sample_sums(eukaryome.rel.p)>0, eukaryome.rel.p)
    
    eukaryome.rel.p= eukaryome.rel.p  %>%
      psmelt() %>%                                         # Melt to long format
      arrange(Rank_taxa_plots) 
    
    
    
    pdf("AbundancetablecustomtaxapersamplepaysstuntedDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 40, # define plot width and height. completely up to user.
        height = 10) 
    ggplot(eukaryome.rel.p, aes(x = Sample, y = Abundance, fill = Rank_taxa_plots)) + 
      facet_wrap(pays~stunted, scales = "free_x") +
      theme(strip.text = element_text(colour = 'black', size=18), axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.title.y = element_text(size=18), title =element_text(size=25, face='bold')) +
      geom_bar(stat = "identity") +
      #scale_fill_manual(values = customtaxa_colors) +
      scale_x_discrete(
        drop = FALSE
      ) +
      # Remove x axis title
      theme(axis.title.x = element_blank()) + 
      theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
      ylab("Relative Abundance") +
      ggtitle("Custom taxa composition of Afribiota samples") 
    dev.off()
    
    
    
    # plot Rank_taxa_plots content merged by country
    
    eukaryome.rel.p = merge_samples(eukaryome.rel, "pays")
    
    eukaryome.rel.p_p <- tax_glom(eukaryome.rel.p, "Rank_taxa_plots")
    
    eukaryome.rel.p_p = transform_sample_counts(eukaryome.rel.p_p, function(x) 100 * x/sum(x))
    
    eukaryome.rel.p_p_s= eukaryome.rel.p_p  %>%
      psmelt() %>%                                         # Melt to long format
      arrange(Rank_taxa_plots) 
    
    class_colors <- c(
      "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
      "#AD6F3B", "#673770","black", "#652926", "#C84248", 
      "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
    )
    
    pdf("Abundancetablerank4paysnopruningDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 11, # define plot width and height. completely up to user.
        height = 8) 
    ggplot(eukaryome.rel.p_p_s, aes(x =Sample, y = Abundance, fill = Rank_taxa_plots)) + 
      theme_bw() +
      theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
      theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
      geom_bar(stat = "identity", width = 1.0) +
      scale_y_continuous(expand = c(0.01,0.01)) +
      theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
      ylab("Relative Abundance") +
      xlab("Pays") + ggtitle ("Relative Abundance of phyla according to country of origin")
    dev.off()
    
    
    
    
    #### Correlations at rank 2  ####
    
    dffiltered_feces_rank2<-dfcleanfeces
    dffiltered_feces_rank2<-tax_glom(dffiltered_feces_rank2, "rank2")
    eukaryome.rel.s<-microbiome::transform(dffiltered_feces_rank2, "compositional")
    eukaryomerelforcorrelation <- otu_table(eukaryome.rel.s)
    dim(eukaryomerelforcorrelation) ## 6 taxa and 241 samples are still in the run. Comment: on species and genus level, almost nothing is conserved!
    
    i <- (colSums(eukaryomerelforcorrelation, na.rm=T) != 0) # T if colSum is not 0, F otherwise
    eukaryomerelforcorrelation_nozero <- eukaryomerelforcorrelation[, i] # all the non-zero columns
    eukaryomerelforcorrelation_nozero<-as.matrix(eukaryomerelforcorrelation_nozero)
    eukaryomerelforcorrelation_nozero[is.na(eukaryomerelforcorrelation_nozero)] <- 0
    dim(eukaryomerelforcorrelation_nozero) ## 241
    View(eukaryomerelforcorrelation_nozero)
    
    View(tax_table(eukaryome.rel.s))
    
    toplot<-as.data.frame(cbind(as.numeric(eukaryomerelforcorrelation_nozero[, 3]) , as.numeric(eukaryomerelforcorrelation_nozero[, 1])))
    
    scatter_plot <- ggplot(toplot, aes(V1, V2))
    scatter_plot + geom_point() + labs(x = "Opisthokonta", y = "Stramenopiles") + geom_smooth(method="lm")
    
    p<- cor.test(toplot$V1, toplot$V2, method = c("spearman"))
    p # p-value = 1.933e-13 //-0.4503847
    
    # now Ophistokonta and Amoebazoa
    
    toplot<-as.data.frame(cbind(as.numeric(eukaryomerelforcorrelation_nozero[, 3]) , as.numeric(eukaryomerelforcorrelation_nozero[, 4])))
    
    scatter_plot <- ggplot(toplot, aes(V1, V2))
    scatter_plot + geom_point() + labs(x = "Opisthokonta", y = "Amoebozoa") + geom_smooth(method="lm")
    
    p<- cor.test(toplot$V1, toplot$V2, method = c("spearman"))
    p # p-value = 0.05772 //rho -0.1224253
    
    # now Ophistokonta and Excavata
     toplot<-as.data.frame(cbind(as.numeric(eukaryomerelforcorrelation_nozero[, 3]) , as.numeric(eukaryomerelforcorrelation_nozero[, 2])))
    
    scatter_plot <- ggplot(toplot, aes(V1, V2))
    scatter_plot + geom_point() + labs(x = "Opisthokonta", y = "Excavata") + geom_smooth(method="lm")
    
    p<- cor.test(toplot$V1, toplot$V2, method = c("spearman"))
    p # p-value = 0.4426 //rho -0.04969142
  
    
    #### Correlations at rank 3  ####
    
    dffiltered_feces_rank3<-dfcleanfeces
    dffiltered_feces_rank3<-tax_glom(dffiltered_feces_rank3, "rank3")
    eukaryome.rel.s<-microbiome::transform(dffiltered_feces_rank3, "compositional")
    eukaryomerelforcorrelation <- otu_table(eukaryome.rel.s)
    dim(eukaryomerelforcorrelation) ## 9 taxa and 241 samples are still in the run. Comment: on species and genus level, almost nothing is conserved!
    
    i <- (colSums(eukaryomerelforcorrelation, na.rm=T) != 0) # T if colSum is not 0, F otherwise
    eukaryomerelforcorrelation_nozero <- eukaryomerelforcorrelation[, i] # all the non-zero columns
    eukaryomerelforcorrelation_nozero<-as.matrix(eukaryomerelforcorrelation_nozero)
    eukaryomerelforcorrelation_nozero[is.na(eukaryomerelforcorrelation_nozero)] <- 0
    dim(eukaryomerelforcorrelation_nozero) ## 240
    View(eukaryomerelforcorrelation_nozero)
    
    View(tax_table(eukaryome.rel.s))
    
    # now Fungi and Incerta sedis
    toplot<-as.data.frame(cbind(as.numeric(eukaryomerelforcorrelation_nozero[, 3]) , as.numeric(eukaryomerelforcorrelation_nozero[, 1])))
    
    scatter_plot <- ggplot(toplot, aes(V1, V2))
    scatter_plot + geom_point() + labs(x = "Fungi", y = "Incertae_Sedis/Stramenopiles") + geom_smooth(method="lm")
    
    p<- cor.test(toplot$V1, toplot$V2, method = c("spearman"))
    p # p-value = 7.789e-10 //rho -0.3828951   
    
    # now Fungi and Conosa
     toplot<-as.data.frame(cbind(as.numeric(eukaryomerelforcorrelation_nozero[, 3]) , as.numeric(eukaryomerelforcorrelation_nozero[, 5])))
    
    scatter_plot <- ggplot(toplot, aes(V1, V2))
    scatter_plot + geom_point() + labs(x = "Fungi", y = "Conosa") + geom_smooth(method="lm")
    
    p<- cor.test(toplot$V1, toplot$V2, method = c("spearman"))
    p # p-value = 0.05909 //rho -0.1217648
    
    # now Fungi and Metazoa
    
    toplot<-as.data.frame(cbind(as.numeric(eukaryomerelforcorrelation_nozero[, 3]) , as.numeric(eukaryomerelforcorrelation_nozero[, 4])))
    
    scatter_plot <- ggplot(toplot, aes(V1, V2))
    scatter_plot + geom_point() + labs(x = "Fungi", y = "Metazoa") + geom_smooth(method="lm")
    
    p<- cor.test(toplot$V1, toplot$V2, method = c("spearman"))
    p # p-value = 0.9145 //rho 0.006950346 
    
    ####  Are there co-occuring eukaryotes in feces based on rank3 ? ####
    # Define data sets to cross-correlate with each other and only keep taxa with given threshold of 0.1%
    dffiltered_feces_species<-tax_glom(dfcleanfeces, "rank3")
    eukaryome.rel.s<-microbiome::transform(dffiltered_feces_species, "compositional")
    core.taxa <- core(eukaryome.rel.s, detection = 0.01, prevalence = 10/100, include.lowest = TRUE) 
    eukaryomerelforcorrelation <- otu_table(core.taxa)
    dim(eukaryomerelforcorrelation) ## 5 taxa and 241 samples are still in the run. Comment: on species and genus level, almost nothing is conserved!
    
    i <- (colSums(eukaryomerelforcorrelation, na.rm=T) != 0) # T if colSum is not 0, F otherwise
    eukaryomerelforcorrelation_nozero <- eukaryomerelforcorrelation[, i] # all the non-zero columns
    eukaryomerelforcorrelation_nozero<-as.matrix(eukaryomerelforcorrelation_nozero)
    eukaryomerelforcorrelation_nozero[is.na(eukaryomerelforcorrelation_nozero)] <- 0
    dim(eukaryomerelforcorrelation_nozero) ## 241
    
    correspondance<-tax_table(core.taxa)
    
    colnames(eukaryomerelforcorrelation_nozero)<-correspondance[, 3]
    
    
    dim(eukaryomerelforcorrelation_nozero)
    x<- eukaryomerelforcorrelation_nozero[1:241, ]
    y<- eukaryomerelforcorrelation_nozero[1:241, ]
    
    correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
    correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
    
    pdf("heatmapfecesDalhousierank3.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 10)
    p <- heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05) 
    print(p)
    dev.off()
    View(tax_table(core.taxa))
    
    
    ####  Are there co-occuring eukaryotes in feces based on custom rank? ####
    # Define data sets to cross-correlate with each other and only keep taxa with given threshold of 0.1%
    dffiltered_feces_custom<-dfcleanfeces
    tax_table(dffiltered_feces_custom)<-tax_table(dffiltered_feces_custom)[, -(1:9)]
    dffiltered_feces_custom<-tax_glom(dffiltered_feces_custom, "Rank_taxa_plots")
    eukaryome.rel.s<-microbiome::transform(dffiltered_feces_custom, "compositional")
    core.taxa <- core(eukaryome.rel.s, detection = 0.01, prevalence = 10/100, include.lowest = TRUE) 
    eukaryomerelforcorrelation <- otu_table(core.taxa)
    dim(eukaryomerelforcorrelation) ## 6 taxa and 241 samples are still in the run. Comment: on species and genus level, almost nothing is conserved!
    
    i <- (colSums(eukaryomerelforcorrelation, na.rm=T) != 0) # T if colSum is not 0, F otherwise
    eukaryomerelforcorrelation_nozero <- eukaryomerelforcorrelation[, i] # all the non-zero columns
    eukaryomerelforcorrelation_nozero<-as.matrix(eukaryomerelforcorrelation_nozero)
    eukaryomerelforcorrelation_nozero[is.na(eukaryomerelforcorrelation_nozero)] <- 0
    dim(eukaryomerelforcorrelation_nozero) ## 241
    
    dim(eukaryomerelforcorrelation_nozero)
    x<- eukaryomerelforcorrelation_nozero[1:241, ]
    y<- eukaryomerelforcorrelation_nozero[1:241, ]
    
    correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
    correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
    pdf("heatmapfecesDalhousiecustom.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 10)
    p <- heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05) 
    print(p)
    dev.off()
    
    
    View(tax_table(dffiltered_feces_custom))
    
   
    
    ####  Are there co-occuring eukaryotes in feces based on rank8 ? ####
    # Define data sets to cross-correlate with each other and only keep taxa with given threshold of 0.1%
    dffiltered_feces_species<-tax_glom(dfcleanfeces, "rank8")
    eukaryome.rel.s<-microbiome::transform(dffiltered_feces_species, "compositional")
    core.taxa <- core(eukaryome.rel.s, detection = 0.01, prevalence = 10/100, include.lowest = TRUE) 
    eukaryomerelforcorrelation <- otu_table(core.taxa)
    dim(eukaryomerelforcorrelation) ## 12 taxa and 241 samples are still in the run. Comment: on species and genus level, almost nothing is conserved!
    
    i <- (colSums(eukaryomerelforcorrelation, na.rm=T) != 0) # T if colSum is not 0, F otherwise
    eukaryomerelforcorrelation_nozero <- eukaryomerelforcorrelation[, i] # all the non-zero columns
    eukaryomerelforcorrelation_nozero<-as.matrix(eukaryomerelforcorrelation_nozero)
    eukaryomerelforcorrelation_nozero[is.na(eukaryomerelforcorrelation_nozero)] <- 0
    dim(eukaryomerelforcorrelation_nozero) ## 241
    
    dim(eukaryomerelforcorrelation_nozero)
    x<- eukaryomerelforcorrelation_nozero[1:241, ]
    y<- eukaryomerelforcorrelation_nozero[1:241, ]
    
    correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
    correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
    pdf("heatmapfecesDalhousiespecies.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 10)
    p <- heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05) 
    print(p)
    dev.off()
    View(tax_table(core.taxa))
    
    
    
    
    ####  Are there co-occuring eukaryotes in feces based on rank7 ? ####
    # Define data sets to cross-correlate with each other and only keep taxa with given threshold of 0.1%
    dffiltered_feces_species<-tax_glom(dfcleanfeces, "rank7")
    eukaryome.rel.s<-microbiome::transform(dffiltered_feces_species, "compositional")
    core.taxa <- core(eukaryome.rel.s, detection = 0.01, prevalence = 10/100, include.lowest = TRUE) 
    eukaryomerelforcorrelation <- otu_table(core.taxa)
    dim(eukaryomerelforcorrelation) ## 6 taxa and 241 samples are still in the run. Comment: on species and genus level, almost nothing is conserved!
    
    i <- (colSums(eukaryomerelforcorrelation, na.rm=T) != 0) # T if colSum is not 0, F otherwise
    eukaryomerelforcorrelation_nozero <- eukaryomerelforcorrelation[, i] # all the non-zero columns
    eukaryomerelforcorrelation_nozero<-as.matrix(eukaryomerelforcorrelation_nozero)
    eukaryomerelforcorrelation_nozero[is.na(eukaryomerelforcorrelation_nozero)] <- 0
    dim(eukaryomerelforcorrelation_nozero) ## 241
    
    dim(eukaryomerelforcorrelation_nozero)
    x<- eukaryomerelforcorrelation_nozero[1:241, ]
    y<- eukaryomerelforcorrelation_nozero[1:241, ]
    
    correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
    correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
    pdf("heatmapfecesDalhousiesgenus.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 10)
    p <- heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05) 
    print(p)
    dev.off()
    
   
    #### now make the actual correlation plots from the ones that seem to be more or less correlated ####
    # comment from heatmaps:custom level: ASV15 (Saccharomycetales) and ASV5 (Blastocystis) seem co-excluding in the custom rank clustered data
    # species level: ASV10 (Blastocystis ST2) and ASV15 (Saccharomyces pastorianus) are weakly co-excluding
    # genus level: ASV5 (Blastocystis) and ASV15 (Saccharomyces) are co-excluding
    # genus level: ASV5 (Blastocystis) and ASV34 (Ascaris) are weakly co-excluding
    
    # plot custom level correlation
    dffiltered_feces_custom<-dfcleanfeces
    tax_table(dffiltered_feces_custom)<-tax_table(dffiltered_feces_custom)[, -(1:9)]
    dffiltered_feces_custom<-tax_glom(dffiltered_feces_custom, "Rank_taxa_plots")
    eukaryome.rel.s<-microbiome::transform(dffiltered_feces_custom, "compositional")
    
    plot(as.numeric(otu_table(eukaryome.rel.s)[, 1]),as.numeric(otu_table(eukaryome.rel.s)[, 3]))
    
    toplot<-as.data.frame(cbind(as.numeric(otu_table(eukaryome.rel.s)[, 1]),as.numeric(otu_table(eukaryome.rel.s)[, 3])))
    
    scatter_plot <- ggplot(toplot, aes(V1, V2))
    scatter_plot + geom_point() + labs(x = "Blastocystis", y = "Saccharomycetes") + geom_smooth(method="lm")
    
    p<- cor.test(otu_table(eukaryome.rel.s)[, 1], otu_table(eukaryome.rel.s)[, 3], method = c("spearman"))
    p # p-value = 1.115e-09 //rho -0.3803482
    
    # now for Species level dataset
    dffiltered_feces_species<-dfcleanfeces
    dffiltered_feces_species<-tax_glom(dffiltered_feces_species, "rank8")
    eukaryome.rel.s<-microbiome::transform(dffiltered_feces_species, "compositional")
    
    View(otu_table(eukaryome.rel.s))
    
   
    # now asssess correlation ASV10 and ASV15
    
    plot(as.numeric(otu_table(eukaryome.rel.s)[, 2]),as.numeric(otu_table(eukaryome.rel.s)[,5]))
    
    toplot<-as.data.frame(cbind(as.numeric(otu_table(eukaryome.rel.s)[, 2]),as.numeric(otu_table(eukaryome.rel.s)[, 5])))
    
    scatter_plot <- ggplot(toplot, aes(V1, V2))
    scatter_plot + geom_point() + labs(x = "Blastocystis ST2", y = "Saccharomyces pastorianus") + geom_smooth(method="lm")
    
    p<- cor.test(otu_table(eukaryome.rel.s)[, 2], otu_table(eukaryome.rel.s)[, 5], method = c("spearman"))
    
    # now for Genus level dataset
    dffiltered_feces_species<-dfcleanfeces
    dffiltered_feces_species<-tax_glom(dffiltered_feces_species, "rank7")
    eukaryome.rel.s<-microbiome::transform(dffiltered_feces_species, "compositional")
    
    View(otu_table(eukaryome.rel.s))
    
    # now asssess correlation ASV5 and ASV15
    
    plot(as.numeric(otu_table(eukaryome.rel.s)[, 1]),as.numeric(otu_table(eukaryome.rel.s)[, 3]))
    
    toplot<-as.data.frame(cbind(as.numeric(otu_table(eukaryome.rel.s)[, 1]),as.numeric(otu_table(eukaryome.rel.s)[, 3])))
    
    scatter_plot <- ggplot(toplot, aes(V1, V2))
    scatter_plot + geom_point() + labs(x = "Blastocystis", y = "Saccharomyces") + geom_smooth(method="lm")
    
    p<- cor.test(otu_table(eukaryome.rel.s)[, 1], otu_table(eukaryome.rel.s)[, 3], method = c("spearman"))
    p # p-value = 2.758e-06 //rho -0.2967019  
    
    # now asssess correlation ASV5 and ASV34
    
    plot(as.numeric(otu_table(eukaryome.rel.s)[, 1]),as.numeric(otu_table(eukaryome.rel.s)[,4]))
    
    toplot<-as.data.frame(cbind(as.numeric(otu_table(eukaryome.rel.s)[, 1]),as.numeric(otu_table(eukaryome.rel.s)[, 4])))
    
    scatter_plot <- ggplot(toplot, aes(V1, V2))
    scatter_plot + geom_point() + labs(x = "Blastocystis", y = "Ascaris") + geom_smooth(method="lm")
    
    p<- cor.test(otu_table(eukaryome.rel.s)[, 1], otu_table(eukaryome.rel.s)[, 4], method = c("spearman"))
    p # p-value = 0.008488 //rho -0.169199 
    
   
 #### Are given taxa associated with the nutritional status/stunting? logistic regression in loop####
 #### cluster level: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test.-stunting ####
    eukaryome.rel<-microbiome::transform(dfcleanfeces, "compositional")
    eukaryome.rel.cluster<-eukaryome.rel
    tax_table(eukaryome.rel.cluster)= tax_table(eukaryome.rel.cluster)[, -12] # need to take away accession number as this does not let collapse correctly
    eukaryome.rel.cluster= tax_glom(eukaryome.rel.cluster, "cluster")
    eukaryome.rel.cluster<-subset_taxa(eukaryome.rel.cluster, cluster!="")
    eukaryome.rel.cluster<-microbiome::transform(eukaryome.rel.cluster, "compositional")
    eukaryome.rel.cluster=subset_samples(eukaryome.rel.cluster, sample_sums(eukaryome.rel.cluster)!="0")
    eukaryome.rel.cluster = filter_taxa(eukaryome.rel.cluster, function(x) sum(x) > 0, TRUE)
    
    df_wilcox <- as.matrix((otu_table(eukaryome.rel.cluster))) #take rel abund and Wilcoxon rank-sum
    meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.cluster)) #take metadata
    tax_wilcox <- as.matrix(tax_table(eukaryome.rel.cluster))
    
    dim(df_wilcox)
    dim(meta_wilcox)
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$stunted)$p.value)
    clusternames<-tax_wilcox[, 12]
    p.res = data.frame(clusternames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.cluster)
    p.res = cbind(as.matrix(p.res),as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.stunted.clusterWilcoxDAL.csv")
    View(p.res) ## signficiant: none after multiple testing!
    
    #### cluster level: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables -stunting ####
    mcorr <- apply(df_wilcox, 2,
                   function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$p.value)
    estcorr <- apply(df_wilcox,2,
                     function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$estimate)
    corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
    
    # Perform multiple comparison correction using a given method of choice
    corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
    corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
    # Merge with tax info
    corres = cbind(as.matrix(corres), as.matrix(Tax_corr))
    
    #export results
    write.csv(corres,"Feces.stuntedCorrelationDAL.csv")
    View(corres) ## several are associated, but they do not survive multiple testing!
    
    #### cluster level: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test.-country of origin ####
    eukaryome.rel.cluster<-subset_taxa(eukaryome.rel.cluster, cluster!="")
    eukaryome.rel.cluster<-microbiome::transform(eukaryome.rel.cluster, "compositional")
    eukaryome.rel.cluster=subset_samples(eukaryome.rel.cluster, sample_sums(eukaryome.rel.cluster)!="0")
    eukaryome.rel.cluster = filter_taxa(eukaryome.rel.cluster, function(x) sum(x) > 0, TRUE)
    
    df_wilcox <- as.matrix((otu_table(eukaryome.rel.cluster))) #take rel abund and Wilcoxon rank-sum
    meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.cluster)) #take metadata
    
    dim(df_wilcox)
    dim(meta_wilcox)
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$pays)$p.value)
    clusternames<-as.matrix(tax_table(eukaryome.rel.cluster))[, 12]
    p.res = data.frame(clusternames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-as.data.frame(tax_table(eukaryome.rel.cluster))
    p.res = cbind(as.matrix(p.res), as.matrix(Tax_corr))
    #export results
    write.csv(p.res,"Feces.pays.clusterWilcoxDAL.csv")
    View(p.res) 
    
    #### cluster level:  continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables -country of origin ####
    mcorr <- apply(df_wilcox, 2,
                   function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$p.value)
    estcorr <- apply(df_wilcox,2,
                     function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$estimate)
    corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
    # Perform multiple comparison correction using a given method of choice
    corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
    corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
    # Merge with tax info
    corres = cbind(as.matrix(corres), as.matrix(Tax_corr))
    
    #export results
    write.csv(corres,"Feces.paysCorrelationDAL.csv")
    View(corres)  
    
    
    #### cluster level:  LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
    #### cluster level: prepare your data ####
    dim(df_wilcox)
    dim(meta_wilcox)
    data1 <- df_wilcox
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_wilcox$age
    df$sexe <- meta_wilcox$sexe
    df$Country <- meta_wilcox$pays
    df$stunted <- meta_wilcox$stunted
    df$calpro <- meta_wilcox$calprotectinelevel
    df$aat <- meta_wilcox$alphaantitrypsinlevel
    df$anemie<-meta_wilcox$anemie2
    df$totalreads<-meta_wilcox$read_count
    
    
    View(df) ## numbers look ok now!
    
    #### cluster level: make a loop for logistic regressions for stunting with inflammation ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "calpro", "aat", "Country", "totalreads",  "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value+ totalreads+ calpro + Country + anemie , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info 
    Tax_corr<-tax_table(eukaryome.rel.cluster)
    Tax_corr_log<-cbind (logresults, Tax_corr[, 12])
    write.csv(Tax_corr_log,"LogresultsFeces.stuntedclusterwithinflaDAL.csv")
    
    #### cluster level: make a loop for logistic regressions for country with inflammation ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "calpro", "aat", "Country", "totalreads",  "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, Country!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(Country ~ value + anemie + calpro + stunted  + totalreads , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    Tax_corr<-tax_table(eukaryome.rel.cluster)
    Tax_corr_log<-cbind (logresults, Tax_corr[, 12])
    write.csv(Tax_corr_log,"LogresultsFeces.paysclusterwithinflaDAL.csv")
    
    
#### ASV level: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
    eukaryome.rel=subset_samples(eukaryome.rel, sample_sums(eukaryome.rel)!="0")
    eukaryome.rel = filter_taxa(eukaryome.rel, function(x) sum(x) > 0, TRUE)
    
    df_wilcox <- as.matrix((otu_table(eukaryome.rel))) #take rel abund and Wilcoxon rank-sum
    meta_wilcox <- as.data.frame(sample_data(eukaryome.rel)) #take metadata
    
    dim(df_wilcox)
    dim(meta_wilcox)
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$stunted)$p.value)
    ASVnames<-colnames(df_wilcox)
    p.res = data.frame(ASVnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-as.data.frame(tax_table(eukaryome.rel))
    p.res = merge(p.res,Tax_corr, by="row.names")
    #export results
    write.csv(p.res,"Feces.stunted.ASVWilcoxDAL.csv")
    View(p.res) ## signficiant: none after multiple testing!
    
    #### ASV level: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
    mcorr <- apply(df_wilcox, 2,
                   function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$p.value)
    estcorr <- apply(df_wilcox,2,
                     function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$estimate)
    corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
    # Perform multiple comparison correction using a given method of choice
    corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
    corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
    # Merge with tax info
    corres = merge(corres,Tax_corr, by="row.names")
    
    #export results
    write.csv(corres,"Feces.stuntedCorrelationDAL.csv")
    View(corres) ## several are associated, but they do not survive multiple testing!
    
    #### ASV level: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
    #### ASV level: prepare your data ####
    dim(df_wilcox)
    dim(meta_wilcox)
    data1 <- df_wilcox
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_wilcox$age
    df$sexe <- meta_wilcox$sexe
    df$Country <- meta_wilcox$pays
    df$stunted <- meta_wilcox$stunted
    df$calpro <- meta_wilcox$calprotectinelevel
    df$aat <- meta_wilcox$alphaantitrypsinlevel
    df$anemie<-meta_wilcox$anemie2
    df$totalreads<-meta_wilcox$read_count
    
    microeuk_reads= dfcleanfeces %>%
      subset_taxa(rank4!="Vertebrata" & rank5!="Vertebrata" & rank6!="Vertebrata" & rank2!="Archaeplastida" & rank4!="Mammalia" & rank2!="Unassigned" ) 
    microeuk_reads_sub<-as.data.frame(otu_table(microeuk_reads))
    microeuk_reads_sub$rownames<-row.names(microeuk_reads_sub)
    microeuk_reads_sub<-filter(microeuk_reads_sub, microeuk_reads_sub$rownames %in% row.names(df))
    dim(microeuk_reads_sub)
    microeuk_reads_sub<-microeuk_reads_sub[, -128]
    df$microeukreads<-rowSums(microeuk_reads_sub)
    df$totalreads<-sample_data(eukaryome.rel)$read_count
    
    View(df) ## numbers look ok now!
    
    #### ASV level: make a loop for logistic regressions for stunting with inflammation ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "calpro", "aat", "Country", "totalreads", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value+ totalreads+ calpro + Country + anemie , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.stuntedASVwithinflaDAL.csv")
    
    #### ASV level: make a loop for logistic regressions for country with inflammation ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "calpro", "aat", "Country", "totalreads", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, Country!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(Country ~ value + anemie + calpro + stunted , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.paysASVwithinflaDAL.csv")
    
 #### Species level: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
    eukaryome.rel.s <- tax_glom(eukaryome.rel, "rank8")
    eukaryome.rel.s=subset_samples(eukaryome.rel.s, sample_sums(eukaryome.rel.s)!="0")
    eukaryome.rel.s = filter_taxa(eukaryome.rel.s, function(x) sum(x) > 0, TRUE)
    
    df_wilcox <- as.matrix((otu_table(eukaryome.rel.s))) #take rel abund and Wilcoxon rank-sum
    meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.s)) #take metadata
    
    dim(df_wilcox)
    dim(meta_wilcox)
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$stunted)$p.value)
    Speciesnames<-colnames(df_wilcox)
    p.res = data.frame(Speciesnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-as.data.frame(tax_table(eukaryome.rel))
    p.res = merge(p.res,Tax_corr, by="row.names")
    #export results
    write.csv(p.res,"Feces.stunted.SpeciesWilcoxDAL.csv")
    View(p.res) ## signficiant: none after multiple testing!
    
    #### Species level: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
    mcorr <- apply(df_wilcox, 2,
                   function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$p.value)
    estcorr <- apply(df_wilcox,2,
                     function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$estimate)
    corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
    # Perform multiple comparison correction using a given method of choice
    corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
    corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
    # Merge with tax info
    corres = merge(corres,Tax_corr, by="row.names")
    
    #export results
    write.csv(corres,"Feces.stuntedSpeciesCorrelationDAL.csv")
    View(corres) ## several are associated, but they do not survive multiple testing!
    
    #### Species level: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
    #### Species level: prepare your data ####
    dim(df_wilcox)
    dim(meta_wilcox)
    data1 <- df_wilcox
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_wilcox$age
    df$sexe <- meta_wilcox$sexe
    df$Country <- meta_wilcox$pays
    df$stunted <- meta_wilcox$stunted
    df$calpro <- meta_wilcox$calprotectinelevel
    df$aat <- meta_wilcox$alphaantitrypsinlevel
    df$anemie<-meta_wilcox$anemie2
    df$totalreads<-meta_wilcox$read_count
    
    #### Species level: make a loop for logistic regressions for stunting with inflammation ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "Country", "totalreads", "aat", "calpro", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value+ calpro + anemie + Country + totalreads , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.stuntedSpecieswithinflaDAL.csv")
   
    #### Species level: make a loop for logistic regressions for country with inflammation ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "calpro", "aat", "Country", "totalreads", "microeukreads", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, Country!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(Country ~ value+ calpro + stunted + totalreads + microeukreads + anemie, .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.paysSpecieswithinflaDAL.csv")
    
    #### Species level: prepare data for without inflammation ####
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_wilcox$age
    df$sexe <- meta_wilcox$sexe
    df$Country <- meta_wilcox$pays
    df$stunted <- meta_wilcox$stunted
    df$anemie<-meta_wilcox$anemie2
    df$totalreads<-meta_wilcox$read_count
    
    #### Species level: make a loop for logistic regressions for stunting without inflammation ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "Country", "totalreads", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value + anemie + Country , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.stuntedSpecieswithoutinflaDAL.csv")

      
    #### Species level: make a loop for logistic regressions for country without inflammation ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "Country", "totalreads", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, Country!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(Country ~ value+ stunted + totalreads + anemie, .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.paysSpecieswithoutinflaDAL.csv")
    
    
 #### Species level: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. For Bangui only ####
    eukaryome.rel.s <- tax_glom(eukaryome.rel, "rank8")
    eukaryome.rel.sb <- subset_samples(eukaryome.rel.s, pays=="RCA")
    eukaryome.rel.sb=subset_samples(eukaryome.rel.sb, sample_sums(eukaryome.rel.s)!="0")
    eukaryome.rel.sb = filter_taxa(eukaryome.rel.sb, function(x) sum(x) > 0, TRUE)
    
    df_wilcox <- as.matrix((otu_table(eukaryome.rel.sb))) #take rel abund and Wilcoxon rank-sum
    meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.sb)) #take metadata
    
    dim(df_wilcox)
    dim(meta_wilcox)
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$stunted)$p.value)
    Speciesnames<-colnames(df_wilcox)
    p.res = data.frame(Speciesnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-as.data.frame(tax_table(eukaryome.rel))
    p.res = merge(p.res,Tax_corr, by="row.names")
    #export results
    write.csv(p.res,"Feces.stunted.SpeciesWilcoxDALBangui.csv")
    View(p.res) ## 
    
    #### Species level: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
    mcorr <- apply(df_wilcox, 2,
                   function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$p.value)
    estcorr <- apply(df_wilcox,2,
                     function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$estimate)
    corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
    # Perform multiple comparison correction using a given method of choice
    corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
    corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
    # Merge with tax info
    corres = merge(corres,Tax_corr, by="row.names")
    
    #export results
    write.csv(corres,"Feces.stuntedSpeciesCorrelationDALBangui.csv")
    View(corres) ## several are associated, but they do not survive multiple testing!
    
    #### Species level: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
    #### Species level: prepare your data ####
    dim(df_wilcox)
    dim(meta_wilcox)
    data1 <- df_wilcox
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_wilcox$age
    df$sexe <- meta_wilcox$sexe
    df$stunted <- meta_wilcox$stunted
    df$calpro <- meta_wilcox$calprotectinelevel
    df$aat <- meta_wilcox$alphaantitrypsinlevel
    df$anemie<-meta_wilcox$anemie2
    
    df$totalreads<-meta_wilcox$read_count
    
    #### Species level: make a loop for logistic regressions for stunting with inflammation ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "totalreads", "aat", "calpro", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value+ totalreads + calpro + anemie , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.stuntedSpecieswithinflaDALBangui.csv")
    
    
    #### Species level: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. For Tana only ####
    eukaryome.rel.s <- tax_glom(eukaryome.rel, "rank8")
    eukaryome.rel.st <- subset_samples(eukaryome.rel.s, pays=="Madagascar")
    eukaryome.rel.st=subset_samples(eukaryome.rel.st, sample_sums(eukaryome.rel.s)!="0")
    eukaryome.rel.st = filter_taxa(eukaryome.rel.st, function(x) sum(x) > 0, TRUE)
    
    df_wilcox <- as.matrix((otu_table(eukaryome.rel.st))) #take rel abund and Wilcoxon rank-sum
    meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.st)) #take metadata
    
    dim(df_wilcox)
    dim(meta_wilcox)
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$stunted)$p.value)
    Speciesnames<-colnames(df_wilcox)
    p.res = data.frame(Speciesnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-as.data.frame(tax_table(eukaryome.rel))
    p.res = merge(p.res,Tax_corr, by="row.names")
    #export results
    write.csv(p.res,"Feces.stunted.SpeciesWilcoxDALTana.csv")
    View(p.res) ## signficiant: none
    
    #### Species level: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
    mcorr <- apply(df_wilcox, 2,
                   function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$p.value)
    estcorr <- apply(df_wilcox,2,
                     function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$estimate)
    corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
    # Perform multiple comparison correction using a given method of choice
    corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
    corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
    # Merge with tax info
    corres = merge(corres,Tax_corr, by="row.names")
    
    #export results
    write.csv(corres,"Feces.stuntedSpeciesCorrelationDALTana.csv")
    View(corres) ## several are associated, but they do not survive multiple testing!
    
    #### Species level: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
    #### Species level: prepare your data ####
    dim(df_wilcox)
    dim(meta_wilcox)
    data1 <- df_wilcox
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_wilcox$age
    df$sexe <- meta_wilcox$sexe
    df$stunted <- meta_wilcox$stunted
    df$calpro <- meta_wilcox$calprotectinelevel
    df$aat <- meta_wilcox$alphaantitrypsinlevel
    df$anemie <- meta_wilcox$anemie2
    
    df$totalreads<-meta_wilcox$read_count
    
    #### Species level: make a loop for logistic regressions for stunting with inflammation ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "totalreads", "aat", "calpro", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value + anemie + totalreads + calpro , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.stuntedSpecieswithinflaDALTana.csv")
    
    
 #### Genus level: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
    eukaryome.rel.g <- tax_glom(eukaryome.rel, "rank7")
    eukaryome.rel.g=subset_samples(eukaryome.rel.g, sample_sums(eukaryome.rel.g)!="0")
    eukaryome.rel.g = filter_taxa(eukaryome.rel.g, function(x) sum(x) > 0, TRUE)
    
    df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.g))) #take rel abund and Wilcoxon rank-sum
    meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.g)) #take metadata
    
    dim(df_wilcox)
    dim(meta_wilcox)
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$stunted)$p.value)
    Genusnames<-colnames(df_wilcox)
    p.res = data.frame(Genusnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-as.data.frame(tax_table(eukaryome.rel))
    p.res = merge(p.res,Tax_corr, by="row.names")
    #export results
    write.csv(p.res,"Feces.stunted.GenusWilcoxDAL.csv")
    View(p.res) ## signficiant: none after multiple testing!
    
    #### Genus level: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
    mcorr <- apply(df_wilcox, 2,
                   function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$p.value)
    estcorr <- apply(df_wilcox,2,
                     function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$estimate)
    corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
    # Perform multiple comparison correction using a given method of choice
    corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
    corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
    # Merge with tax info
    corres = merge(corres,Tax_corr, by="row.names")
    
    plot(as.numeric(meta_wilcox$haz_cont),as.numeric(df_wilcox[, 1]))
    p<- cor.test(meta_wilcox$haz_cont, df_wilcox[, 1], method = c("spearman"))
    p # p-value = 0.0002664, rho 0.2907506 
    
    df_wilcox2<-as.matrix(df_wilcox)
    df_wilcox2<-subset(df_wilcox2, df_wilcox2[, 1] > 0)
    
    meta_wilcox2<-filter(meta_wilcox, row.names(meta_wilcox) %in% row.names(df_wilcox2) )
    plot(as.numeric(meta_wilcox2$haz_cont),as.numeric(df_wilcox2[, 1]))
    toplot<-as.data.frame(cbind(meta_wilcox2$haz_cont, df_wilcox2[, 1]))
    
    scatter_plot <- ggplot(toplot, aes(V1, ASV5))
    scatter_plot + geom_point() + labs(x = "Rel. abundance of Blastocystis", y = "Height-for-age z-score") + geom_smooth(method="lm")
    
    p<- cor.test(meta_wilcox2$haz_cont, df_wilcox2[, 1], method = c("spearson"))
    p # p-value = 0.0002664
    
    #export results
    write.csv(corres,"Feces.stuntedGenusCorrelationDAL.csv")
    View(corres) ## several are associated, but they do not survive multiple testing!
    
    #### Genus level: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
    #### Genus level: prepare your data ####
    dim(df_wilcox)
    dim(meta_wilcox)
    data1 <- df_wilcox
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_wilcox$age
    df$sexe <- meta_wilcox$sexe
    df$Country <- meta_wilcox$pays
    df$stunted <- meta_wilcox$stunted
    df$calpro <- meta_wilcox$calprotectinelevel
    df$aat <- meta_wilcox$alphaantitrypsinlevel
    df$anemie <- meta_wilcox$anemie2
    df$totalreads<-meta_wilcox$read_count
    
    #### Genus level: make a loop for logistic regressions ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "calpro", "aat", "Country", "totalreads", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value + anemie + Country , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.stuntedGenuswithoutinflaDAL.csv")
    
    #### Genus level: make a loop for logistic regressions for country with inflammation ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "calpro", "aat", "Country", "totalreads", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, Country!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(Country ~ value+ anemie + stunted + calpro , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.paysGenuswithinflaDAL.csv")
    
    #### Genus level: prepare data for without inflammation ####
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_wilcox$age
    df$sexe <- meta_wilcox$sexe
    df$Country <- meta_wilcox$pays
    df$stunted <- meta_wilcox$stunted
    df$totalreads<-meta_wilcox$read_count
    df$run<-meta_wilcox$run
    df$anemie<-meta_wilcox$anemie2
    
    #### Genus level: make a loop for logistic regressions stunting without inflammation ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "Country", "totalreads", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value + anemie + Country , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.stuntedGenuswithoutinflaDAL.csv")
    
    
    
 #### Genus level: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. Bangui only ####
    eukaryome.rel.g <- tax_glom(eukaryome.rel, "rank7")
    eukaryome.rel.gb <- subset_samples(eukaryome.rel.g, pays=="RCA")
    eukaryome.rel.gb=subset_samples(eukaryome.rel.gb, sample_sums(eukaryome.rel.gb)!="0")
    eukaryome.rel.gb = filter_taxa(eukaryome.rel.gb, function(x) sum(x) > 0, TRUE)
    
    df_wilcox <- as.matrix((otu_table(eukaryome.rel.gb))) #take rel abund and Wilcoxon rank-sum
    meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.gb)) #take metadata
    
    dim(df_wilcox)
    dim(meta_wilcox)
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$stunted)$p.value)
    Genusnames<-colnames(df_wilcox)
    p.res = data.frame(Genusnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-as.data.frame(tax_table(eukaryome.rel))
    p.res = merge(p.res,Tax_corr, by="row.names")
    #export results
    write.csv(p.res,"Feces.stunted.GenusWilcoxDALBangui.csv")
    View(p.res) ## signficiant: none after multiple testing!
    
    #### Genus level: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
    mcorr <- apply(df_wilcox, 2,
                   function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$p.value)
    estcorr <- apply(df_wilcox,2,
                     function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$estimate)
    corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
    # Perform multiple comparison correction using a given method of choice
    corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
    corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
    # Merge with tax info
    corres = merge(corres,Tax_corr, by="row.names")
    
    df_wilcox2<-as.matrix(df_wilcox)
    df_wilcox2<-subset(df_wilcox2, df_wilcox2[, 1] > 0)
    
    meta_wilcox2<-filter(meta_wilcox, as.vector(row.names(meta_wilcox)) %in% as.vector(row.names(df_wilcox2) ))
    
    plot(as.numeric(meta_wilcox$haz_cont),as.numeric(df_wilcox[, 1]))
    plot(as.numeric(meta_wilcox2$haz_cont),as.numeric(df_wilcox2[, 1]))
    
    toplot<-as.data.frame(cbind(meta_wilcox$haz_cont, df_wilcox[, 1]))
    
    scatter_plot <- ggplot(toplot, aes(V1, ASV5))
    scatter_plot + geom_point() + labs(x = "Rel. abundance of Blastocystis", y = "Height-for-age z-score") + geom_smooth(method="lm")
    
    p<- cor.test(meta_wilcox$haz_cont, df_wilcox[, 1], method = c("spearman"))
    p # p-value = 0.8924 //0.0150799
    
    df_wilcox2<-as.matrix(df_wilcox)
    df_wilcox2<-subset(df_wilcox2, df_wilcox2[, 1] > 0)
    
    meta_wilcox2<-filter(meta_wilcox, meta_wilcox$row_names %in% row.names(df_wilcox2) )
    
    toplot<-as.data.frame(cbind(meta_wilcox2$haz_cont, df_wilcox2[, 1]))
    
    scatter_plot <- ggplot(toplot, aes(V1, ASV5))
    scatter_plot + geom_point() + labs(x = "Rel. abundance of Blastocystis", y = "Height-for-age z-score") + geom_smooth(method="lm")
    
    plot(as.numeric(meta_wilcox2$haz_cont),as.numeric(df_wilcox2[, 1]))
    p<- cor.test(meta_wilcox2$haz_cont, df_wilcox2[, 1], method = c("spearman"))
    p # p-value = 0.8924 //0.0150799
    
    #export results
    write.csv(corres,"Feces.stuntedGenusCorrelationDALBangui.csv")
    View(corres) ## several are associated, but they do not survive multiple testing!
    
    #### Genus level: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
    #### Genus level: prepare your data ####
    dim(df_wilcox)
    dim(meta_wilcox)
    data1 <- df_wilcox
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_wilcox$age
    df$sexe <- meta_wilcox$sexe
    df$stunted <- meta_wilcox$stunted
    df$calpro <- meta_wilcox$calprotectinelevel
    df$aat <- meta_wilcox$alphaantitrypsinlevel
    df$anemie <- meta_wilcox$anemie2
    
    df$totalreads<-meta_wilcox$read_count
    
    #### Genus level: make a loop for logistic regressions for stunting ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "calpro", "aat", "totalreads", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value+ anemie  , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.stuntedGenuswithinflaDALBangui.csv")
    
    
    
 #### Genus level: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. Tana only ####
    eukaryome.rel.g <- tax_glom(eukaryome.rel, "rank7")
    eukaryome.rel.gt <- subset_samples(eukaryome.rel.g, pays=="Madagascar")
    eukaryome.rel.gt=subset_samples(eukaryome.rel.gt, sample_sums(eukaryome.rel.gt)!="0")
    eukaryome.rel.gt = filter_taxa(eukaryome.rel.gt, function(x) sum(x) > 0, TRUE)
    
    df_wilcox <- as.matrix((otu_table(eukaryome.rel.gt))) #take rel abund and Wilcoxon rank-sum
    meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.gt)) #take metadata
    
    dim(df_wilcox)
    dim(meta_wilcox)
    
    MW.p = apply(df_wilcox,2,
                 function(x) wilcox.test(c(x)~meta_wilcox$stunted)$p.value)
    Genusnames<-colnames(df_wilcox)
    p.res = data.frame(Genusnames,MW.p)
    # Perform multiple comparison correction using a given method of choice
    p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
    p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
    # Merge with tax info
    Tax_corr<-as.data.frame(tax_table(eukaryome.rel))
    p.res = merge(p.res,Tax_corr, by="row.names")
    #export results
    write.csv(p.res,"Feces.stunted.GenusWilcoxDALTana.csv")
    View(p.res) ## signficiant: none after multiple testing!
    
    #### Genus level: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
    mcorr <- apply(df_wilcox, 2,
                   function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$p.value)
    estcorr <- apply(df_wilcox,2,
                     function(x) cor.test(c(x), meta_wilcox$haz_cont, method="spearman")$estimate)
    corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
    # Perform multiple comparison correction using a given method of choice
    corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
    corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
    # Merge with tax info
    corres = merge(corres,Tax_corr, by="row.names")
    
    df_wilcox2<-as.matrix(df_wilcox)
    df_wilcox2<-subset(df_wilcox2, df_wilcox2[, 1] > 0)
    
    meta_wilcox2<-filter(meta_wilcox, row.names(meta_wilcox) %in% row.names(df_wilcox2) )
    
    plot(as.numeric(meta_wilcox$haz_cont),as.numeric(df_wilcox[, 1]))
    plot(as.numeric(meta_wilcox2$haz_cont),as.numeric(df_wilcox2[, 1]))
    
    toplot<-as.data.frame(cbind(meta_wilcox$haz_cont, df_wilcox[, 1]))
    
    scatter_plot <- ggplot(toplot, aes(V1, ASV5))
    scatter_plot + geom_point() + labs(x = "Rel. abundance of Blastocystis", y = "Height-for-age z-score") + geom_smooth(method="lm")
    
    p<- cor.test(meta_wilcox$haz_cont, df_wilcox[, 1], method = c("spearman"))
    p # p-value = 0.001828 //0.2602695
    
    df_wilcox2<-as.matrix(df_wilcox)
    df_wilcox2<-subset(df_wilcox2, df_wilcox2[, 1] > 0)
    
    meta_wilcox2<-filter(meta_wilcox, meta_wilcox$row_names %in% row.names(df_wilcox2) )
    
    toplot<-as.data.frame(cbind(meta_wilcox2$haz_cont, df_wilcox2[, 1]))
    
    scatter_plot <- ggplot(toplot, aes(V1, ASV5))
    scatter_plot + geom_point() + labs(x = "Rel. abundance of Blastocystis", y = "Height-for-age z-score") + geom_smooth(method="lm")
    
    plot(as.numeric(meta_wilcox2$haz_cont),as.numeric(df_wilcox2[, 1]))
    p<- cor.test(meta_wilcox2$haz_cont, df_wilcox2[, 1], method = c("spearman"))
    p # p-value = 0.001828 //0.2602695 
    
    #export results
    write.csv(corres,"Feces.stuntedGenusCorrelationDALTana.csv")
    View(corres) 
    #### Genus level: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
    #### Genus level: prepare your data ####
    dim(df_wilcox)
    dim(meta_wilcox)
    data1 <- df_wilcox
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_wilcox$age
    df$sexe <- meta_wilcox$sexe
    df$stunted <- meta_wilcox$stunted
    df$calpro <- meta_wilcox$calprotectinelevel
    df$aat <- meta_wilcox$alphaantitrypsinlevel
    df$anemie <- meta_wilcox$anemie2

    
    df$totalreads<-meta_wilcox$read_count
    
    #### Genus level: make a loop for logistic regressions ####
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "calpro", "aat", "totalreads", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value+ anemie + calpro , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    # Merge with tax info
    row.names(logresults)<-logresults$variable
    Tax_corr$variable<-row.names(Tax_corr)
    Tax_corr_log<-filter(Tax_corr, row.names(Tax_corr) %in% logresults$variable)
    logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFeces.stuntedGenuswithinflaDALTana.csv")
    
    
    
    
    
#### DeSeq2 ####  
    #### on ASV level and stunting  ####
    dfcleanfeces<-subset_samples(dfcleanfeces, stunted!="")
    dfcleanfeces<-subset_samples(dfcleanfeces, pays!="")
    dfcleanfeces<-subset_samples(dfcleanfeces, ageyears!="")
    dfcleanfeces<-subset_samples(dfcleanfeces, sexe!="")
    dfcleanfeces<-subset_samples(dfcleanfeces, anemie2!="")
    dfcleanfeces<-subset_samples(dfcleanfeces, calprotectinelevel!="")
    dfcleanfeces = filter_taxa(dfcleanfeces, function(x) sum(x) > 0, TRUE)
    dfcleanfeces<-subset_samples(dfcleanfeces, (rowSums(otu_table(dfcleanfeces))!=0))
    dfcleanfeces<-subset_samples(dfcleanfeces, (colSums(otu_table(dfcleanfeces))!=0))
    
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces, ~  read_count + pays + calprotectinelevel + anemie2+ stunted)
    
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count+ pays + calprotectinelevel + anemie2)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    res_p_ordered_filt_2
    
    # Make an abundance and prevalence table at ASV level 
    #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
    prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
      x[x >= 2] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
      x[x >= 0] <- 1
      return(x)
    }
    
    prev_counts.stunted <- dfcleanfeces %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.stunted <- dfcleanfeces %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.stunted/prev_possible.stunted)*100
    
    tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces))
    
    merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "stunted")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.stunted.prev) <-merge.stunted.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.stunted.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    stuntedcorrstuntedgenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    View(stuntedcorrstuntedgenderage_reduced_prevabdeseq)
    write.csv(stuntedcorrstuntedgenderage_reduced_prevabdeseq, "stuntedcorrstuntedgenderage_reduced_prevabdeseqDAL.csv")
    
    #plot
    
    # Order order
    x = tapply(stuntedcorrstuntedgenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(stuntedcorrstuntedgenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank5, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank5 = factor(as.character(stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank5), levels=names(x))
    # Family order
    x = tapply(stuntedcorrstuntedgenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank6, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank6 = factor(as.character(stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank6), levels=names(x))
    # Genus order
    x = tapply(stuntedcorrstuntedgenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank7, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank7 = factor(as.character(stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank7), levels=names(x))
    # Species order
    x = tapply(stuntedcorrstuntedgenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank8, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank8 = factor(as.character(stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank8), levels=names(x))
    # Accession order
    x = tapply(stuntedcorrstuntedgenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrstuntedgenderage_reduced_prevabdeseq$accession, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrstuntedgenderage_reduced_prevabdeseq$Accession = factor(as.character(stuntedcorrstuntedgenderage_reduced_prevabdeseq$accession), levels=names(x))
    
    #make the actual plot
    
    pdf("differentially_present_taxa_ASV_stunted_corr_stunted_sexe_age_feces_LRTDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(stuntedcorrstuntedgenderage_reduced_prevabdeseq, aes(x=accession, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
    dev.off()
    
    # make  graph plot
    stuntedcorrstuntedgenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank6, "/", stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank7, "/", stuntedcorrstuntedgenderage_reduced_prevabdeseq$rank8)
    pdf("differentially_present_ASVbystuntedDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 12, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=stuntedcorrstuntedgenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different ASV by stunted ")
    dev.off()
    
    
    
    #### on ASV level and country correcting for inflammation ####
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces, ~  read_count+ stunted + calprotectinelevel + anemie2 + pays)
    
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count+ stunted + calprotectinelevel + anemie2)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") #log2( Madagascar/RCA ) #so when we say things are "upregulated" we mean more prevalent in children from RCA )
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) # 
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    res_p_ordered_filt_2
    
    # Make an abundance and prevalence table at ASV level 
    #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
    prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
      x[x >= 2] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
      x[x >= 0] <- 1
      return(x)
    }
    
    prev_counts.pays <- dfcleanfeces %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.pays <- dfcleanfeces %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.pays/prev_possible.pays)*100
    
    tax_table.pays =  as.data.frame(tax_table(dfcleanfeces))
    
    merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "pays")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.pays.prev) <-merge.pays.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.pays.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    payscorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    View(payscorrcountrygenderage_reduced_prevabdeseq)
    write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrcountrygenderage_reduced_prevabdeseqDAL.csv")
    
    #plot
    
    # Order order
    x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$rank5, function(x) max(x))
    x = sort(x, TRUE)
    payscorrcountrygenderage_reduced_prevabdeseq$rank5 = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$rank5), levels=names(x))
    # Family order
    x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$rank6, function(x) max(x))
    x = sort(x, TRUE)
    payscorrcountrygenderage_reduced_prevabdeseq$rank6 = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$rank6), levels=names(x))
    # Genus order
    x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$rank7, function(x) max(x))
    x = sort(x, TRUE)
    payscorrcountrygenderage_reduced_prevabdeseq$rank7 = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$rank7), levels=names(x))
    # Species order
    x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$rank8, function(x) max(x))
    x = sort(x, TRUE)
    payscorrcountrygenderage_reduced_prevabdeseq$rank8 = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$rank8), levels=names(x))
    # Accession order
    x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$accession, function(x) max(x))
    x = sort(x, TRUE)
    payscorrcountrygenderage_reduced_prevabdeseq$Accession = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$accession), levels=names(x))
    
    #make the actual plot
    
    pdf("differentially_present_taxa_ASV_pays_corr_country_sexe_age_feces_LRTDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(payscorrcountrygenderage_reduced_prevabdeseq, aes(x=accession, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
    dev.off()
    
    # make  graph plot
    payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$rank6, "/", payscorrcountrygenderage_reduced_prevabdeseq$rank7, "/", payscorrcountrygenderage_reduced_prevabdeseq$rank8)
    pdf("differentially_present_ASVbypaysDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 12, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different ASV by country of origin")
    dev.off()
    
    
    
    #### on Species/rank8 level and stunting  ####
    dfcleanfeces_s<-tax_glom(dfcleanfeces, "rank8")
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, stunted!="")
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, pays!="")
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, ageyears!="")
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, sexe!="")
    dfcleanfeces_s = filter_taxa(dfcleanfeces_s, function(x) sum(x) > 0, TRUE)
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, (rowSums(otu_table(dfcleanfeces))!=0))
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s, ~ read_count + pays + anemie2 + age + stunted)
    
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + pays + anemie2 + age)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    View(res_p_ordered_filt_2)
    
    
    
    
    #### on Species/rank8 level and stunting, correcting for alphaantitrypsin  ####
    dfcleanfeces_s<-tax_glom(dfcleanfeces, "rank8")
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, stunted!="")
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, pays!="")
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, ageyears!="")
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, sexe!="")
    dfcleanfeces_s = filter_taxa(dfcleanfeces_s, function(x) sum(x) > 0, TRUE)
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, alphaantitrypsinlevel!="")
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s2, (rowSums(otu_table(dfcleanfeces))!=0))
    
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + pays + anemie2 + alphaantitrypsinlevel + stunted)
    
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + pays + anemie2 + alphaantitrypsinlevel)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    View(res_p_ordered_filt_2)
    
    
    #### on Species/rank8 level and stunting, correcting for calprotectine  ####
    dfcleanfeces_s<-tax_glom(dfcleanfeces, "rank8")
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, stunted!="")
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, pays!="")
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, ageyears!="")
    dfcleanfeces_s<-subset_samples(dfcleanfeces_s, sexe!="")
    dfcleanfeces_s = filter_taxa(dfcleanfeces_s, function(x) sum(x) > 0, TRUE)
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, calprotectinelevel!="")
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s2, (rowSums(otu_table(dfcleanfeces))!=0))
    
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + pays + anemie2 + calprotectinelevel + stunted)
    
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + pays  + anemie2 + calprotectinelevel)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) # 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    View(res_p_ordered_filt_2)
    
    
    #### on Species/rank8 level and stunting, not corrected for inflammation, Mada only  ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, pays!="RCA")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + anemie2 + stunted)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + anemie2 )
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) # 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    res_p_ordered_filt_2
    
    
    
    #### on Species/rank8 level and stunting, not corrected for inflammation, CAR only  ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, pays!="Madagascar")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + anemie2 + stunted)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + anemie2)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    res_p_ordered_filt_2
    
    prev_counts.stunted <- dfcleanfeces %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.stunted <- dfcleanfeces %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.stunted/prev_possible.stunted)*100
    
    tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces))
    
    merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "stunted")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.stunted.prev) <-merge.stunted.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.stunted.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    stuntedcorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    View(stuntedcorrcountrygenderage_reduced_prevabdeseq)
    write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "stuntdcorrcountrygenderage_reduced_prevabdeseqDALBAngui.csv")
    
    #plot
    
    # Order order
    x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$rank5, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrcountrygenderage_reduced_prevabdeseq$rank5 = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$rank5), levels=names(x))
    # Family order
    x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$rank6, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrcountrygenderage_reduced_prevabdeseq$rank6 = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$rank6), levels=names(x))
    # Genus order
    x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$rank7, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrcountrygenderage_reduced_prevabdeseq$rank7 = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$rank7), levels=names(x))
    # Species order
    x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$rank8, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrcountrygenderage_reduced_prevabdeseq$rank8 = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$rank8), levels=names(x))
    # Accession order
    x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$accession, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrcountrygenderage_reduced_prevabdeseq$Accession = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$accession), levels=names(x))
    
    #make the actual plot
    
    pdf("differentially_present_taxa_Species_stunted_corr_country_sexe_age_feces_LRTDALBAngui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=accession, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
    dev.off()
    
    # make  graph plot
    stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$rank6, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$rank7, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$rank8)
    pdf("differentially_present_SpeciesbystuntedDALBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 12, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different Species by stunting status, Bangui dataset only")
    dev.off()
    
    #### on Species/rank8 level and stunting, correcting for calprotectine, Mada only   ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, calprotectinelevel!="")
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s2, pays!="RCA")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + anemie2 + calprotectinelevel + stunted)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + anemie2 + calprotectinelevel)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) # 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 
    
    prev_counts.stunted <- dfcleanfeces %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.stunted <- dfcleanfeces %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.stunted/prev_possible.stunted)*100
    
    tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces))
    
    merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "stunted")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.stunted.prev) <-merge.stunted.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.stunted.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    stuntedcorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    View(stuntedcorrcountrygenderage_reduced_prevabdeseq)
    write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "stuntdcorrcountrygenderage_reduced_prevabdeseqDALBAngui.csv")
    
    #plot
    
    # Order order
    x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$rank5, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrcountrygenderage_reduced_prevabdeseq$rank5 = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$rank5), levels=names(x))
    # Family order
    x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$rank6, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrcountrygenderage_reduced_prevabdeseq$rank6 = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$rank6), levels=names(x))
    # Genus order
    x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$rank7, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrcountrygenderage_reduced_prevabdeseq$rank7 = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$rank7), levels=names(x))
    # Species order
    x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$rank8, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrcountrygenderage_reduced_prevabdeseq$rank8 = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$rank8), levels=names(x))
    # Accession order
    x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$accession, function(x) max(x))
    x = sort(x, TRUE)
    stuntedcorrcountrygenderage_reduced_prevabdeseq$Accession = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$accession), levels=names(x))
    
    #make the actual plot
    
    pdf("differentially_present_taxa_Species_stunted_corr_country_sexe_age_feces_LRTDALTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=accession, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
    dev.off()
    
    # make  graph plot
    stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$rank6, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$rank7, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$rank8)
    pdf("differentially_present_SpeciesbystuntedDALTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 12, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different Species by stunting status, Antananarivo dataset only")
    dev.off()
    
   
    #### on Species/rank8 level and stunting, correcting for calprotectine, CAR only ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, calprotectinelevel!="")
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s2, pays!="Madagascar")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + anemie2 + calprotectinelevel + stunted)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + anemie2 + calprotectinelevel)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    # show results preserving taxa table
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_s2)[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    # Make an abundance and prevalence table at rank 8 /Species level
    #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
    prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
      x[x >= 2] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
      x[x >= 0] <- 1
      return(x)
    }
    
    prev_counts.stunted <- dfcleanfeces_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.stunted <- dfcleanfeces_s2 %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.stunted/prev_possible.stunted)*100
    
    tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces_s2))
    
    merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces_s2 %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "stunted")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.stunted.prev) <-merge.stunted.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.stunted.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    stuntedcorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "CARonlystuntedcorrcountrygenderagecalprotectine_reduced_prevabdeseqDAL.csv")
    
    #plot
    
    # Order order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5), levels=names(x))
    # Family order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6), levels=names(x))
    # Genus order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7), levels=names(x))
    # Species order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8), levels=names(x))
    
    
    #make the actual plot
    
    pdf("CARonlydifferentially_present_taxa_Species_stunted_corr_country_gender_age_calpro_feces_LRTDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=rank8, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunting status in feces from CAR, controlled for  gender, age & calpro,  LRT model")
    dev.off()
    
    # make a boxplot of rel. abundance
    dfcleanfeces_s2.rel <- microbiome::transform(dfcleanfeces_s2, "compositional")
    data<-data.frame(otu_table(dfcleanfeces_s2.rel))
    data$stunted<-sample_data(dfcleanfeces_s2.rel)$stunted
    
    # make  graph plot
    stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$rank6, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$rank7, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$rank8)
    pdf("CARonlydifferentially_present_SpeciesbystuntingcorrcalproDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different Species by stunting status in CAR, corrected for age, country of origin, calprotectine level")
    dev.off()
    
    
    #### on Species/rank8 level and stunting, correcting for alphaantitrypsin, Mada only  ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, alphaantitrypsinlevel!="")
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s2, pays!="RCA")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + anemie2 + alphaantitrypsinlevel + stunted)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + anemie2 + alphaantitrypsinlevel)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 
    
    #### on Species/rank8 level and stunting, correcting for alphaantitrypsin, CAR only  ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, alphaantitrypsinlevel!="")
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s2, pays!="Madagascar")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + anemie2 + alphaantitrypsinlevel + stunted)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + anemie2 + alphaantitrypsinlevel)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    # show results preserving taxa table
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_s2)[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    # Make an abundance and prevalence table at rank 8 /Species level
    #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
    prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
      x[x >= 2] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
      x[x >= 0] <- 1
      return(x)
    }
    
    prev_counts.stunted <- dfcleanfeces_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.stunted <- dfcleanfeces_s2 %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.stunted/prev_possible.stunted)*100
    
    tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces_s2))
    
    merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces_s2 %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "stunted")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.stunted.prev) <-merge.stunted.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.stunted.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    stuntedcorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "CARonlystuntedcorrcountrygenderagealphaantitrypsine_reduced_prevabdeseqDAL.csv")
    
    #plot
    
    # Order order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5), levels=names(x))
    # Family order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6), levels=names(x))
    # Genus order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7), levels=names(x))
    # Species order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8), levels=names(x))
    
    
    #make the actual plot
    
    pdf("CARonlydifferentially_present_taxa_Species_stunted_corr_country_gender_age_calpro_feces_LRTDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=rank8, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunting status in feces from CAR, controlled for  gender, age & calpro,  LRT model")
    dev.off()
    
    # make a boxplot of rel. abundance
    dfcleanfeces_s2.rel <- microbiome::transform(dfcleanfeces_s2, "compositional")
    data<-data.frame(otu_table(dfcleanfeces_s2.rel))
    data$stunted<-sample_data(dfcleanfeces_s2.rel)$stunted
    
    # make  graph plot
    stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$rank6, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$rank7, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$rank8)
    pdf("CARonlydifferentially_present_SpeciesbystuntingcorrcalproDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different Species by stunting status in CAR, corrected for age, country of origin, alphaantitrypsine level")
    dev.off()
    
    
    
    
    
    #### on Species/rank8 level and country, not correcting for inflammation ####
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s, ~  read_count + stunted + anemie2 + pays)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + stunted + anemie2)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) # 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    # show results preserving taxa table
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_s)[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    # Make an abundance and prevalence table at rank 8 /Species level
    #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
    prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
      x[x >= 2] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
      x[x >= 0] <- 1
      return(x)
    }
    
    prev_counts.pays <- dfcleanfeces_s %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.pays <- dfcleanfeces_s %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.pays/prev_possible.pays)*100
    
    tax_table.pays =  as.data.frame(tax_table(dfcleanfeces_s))
    
    merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces_s %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "pays")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.pays.prev) <-merge.pays.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.pays.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    payscorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrcountrygenderage_reduced_prevabdeseqDAL.csv")
    
    #plot
    
    # Order order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5), levels=names(x))
    # Family order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6), levels=names(x))
    # Genus order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7), levels=names(x))
    # Species order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8), levels=names(x))
    
    
    #make the actual plot
    
    pdf("differentially_present_taxa_Species_pays_corr_country_gender_age_feces_LRTDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=rank8, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
    dev.off()
    
    # make a boxplot of rel. abundance
    dfcleanfeces_s.rel <- microbiome::transform(dfcleanfeces_s, "compositional")
    data<-data.frame(otu_table(dfcleanfeces_s.rel))
    data$pays<-sample_data(dfcleanfeces_s.rel)$pays
    
    # make  graph plot
    payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$rank6, "/", payscorrcountrygenderage_reduced_prevabdeseq$rank7, "/", payscorrcountrygenderage_reduced_prevabdeseq$rank8)
    pdf("differentially_present_SpeciesbypaysDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different Species by country of origin")
    dev.off()
    
    
    #### on Species/rank8 level and country, correcting for calprotectine ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, calprotectinelevel!="")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~ read_count + stunted + anemie2 + calprotectinelevel + pays )
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + stunted + anemie2 + calprotectinelevel)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    # show results preserving taxa table
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_s2)[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    # Make an abundance and prevalence table at rank 8 /Species level
    #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
    prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
      x[x >= 2] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
      x[x >= 0] <- 1
      return(x)
    }
    
    prev_counts.pays <- dfcleanfeces_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.pays <- dfcleanfeces_s2 %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.pays/prev_possible.pays)*100
    
    tax_table.pays =  as.data.frame(tax_table(dfcleanfeces))
    
    merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces_s2 %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "pays")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.pays.prev) <-merge.pays.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.pays.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    payscorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrcountrygenderagecalprotectine_reduced_prevabdeseqDAL.csv")
    
    #plot
    
    # Order order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5), levels=names(x))
    # Family order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6), levels=names(x))
    # Genus order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7), levels=names(x))
    # Species order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8), levels=names(x))
    
    
    #make the actual plot
    
    pdf("differentially_present_taxa_Species_pays_corr_country_gender_age_calpro_feces_LRTDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=rank8, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & calpro,  LRT model")
    dev.off()
    
    # make a boxplot of rel. abundance
    dfcleanfeces_s2.rel <- microbiome::transform(dfcleanfeces_s2, "compositional")
    data<-data.frame(otu_table(dfcleanfeces_s2.rel))
    data$pays<-sample_data(dfcleanfeces_s2.rel)$pays
    
    # make  graph plot
    payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$rank6, "/", payscorrcountrygenderage_reduced_prevabdeseq$rank7, "/", payscorrcountrygenderage_reduced_prevabdeseq$rank8)
    pdf("differentially_present_SpeciesbypayscorrcalproDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 8, # define plot width and height. completely up to user.
        height = 4)
    ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different Species by country of origin, corrected for stunting status, calprotectine level")
    dev.off()
    
    
    #### on Species/rank8 level and country, correcting for alphaantitrypsin ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, alphaantitrypsinlevel!="")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + stunted + anemie2 + alphaantitrypsinlevel + pays )
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + stunted + anemie2 + alphaantitrypsinlevel)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    # show results preserving taxa table
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_s2)[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    # Make an abundance and prevalence table at rank 8 /Species level
    #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
    prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
      x[x >= 2] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
      x[x >= 0] <- 1
      return(x)
    }
    
    prev_counts.pays <- dfcleanfeces_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.pays <- dfcleanfeces_s2 %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.pays/prev_possible.pays)*100
    
    tax_table.pays =  as.data.frame(tax_table(dfcleanfeces))
    
    merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces_s2 %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "pays")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.pays.prev) <-merge.pays.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.pays.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    payscorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrcountrygenderagealphaantitrypsin_reduced_prevabdeseqDAL.csv")
    
    #plot
    
    # Order order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5), levels=names(x))
    # Family order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6), levels=names(x))
    # Genus order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7), levels=names(x))
    # Species order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8), levels=names(x))
    
    
    #make the actual plot
    
    pdf("differentially_present_taxa_Species_pays_corr_country_gender_age_aat_feces_LRTDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=rank8, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & calpro,  LRT model")
    dev.off()
    
    # make a boxplot of rel. abundance
    dfcleanfeces_s2.rel <- microbiome::transform(dfcleanfeces_s2, "compositional")
    data<-data.frame(otu_table(dfcleanfeces_s2.rel))
    data$pays<-sample_data(dfcleanfeces_s2.rel)$pays
    
    # make  graph plot
    payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$rank6, "/", payscorrcountrygenderage_reduced_prevabdeseq$rank7, "/", payscorrcountrygenderage_reduced_prevabdeseq$rank8)
    pdf("differentially_present_SpeciesbypayscorraatDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different Species by country of origin, corrected for stunting status, alphaantitrypsin level")
    dev.off()
    
    
    
    #### on Species/rank8 level and calprotectinelevel  ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, calprotectinelevel!="")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + stunted + anemie2 + pays + calprotectinelevel)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + stunted + anemie2 + pays)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "calprotectinelevel") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    res_p_ordered_filt_2
    
    # Make an abundance and prevalence table at rank 8 /Species level
    #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
    prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
      x[x >= 2] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
      x[x >= 0] <- 1
      return(x)
    }
    
    prev_counts.calprotectinelevel <- dfcleanfeces_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("calprotectinelevel") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.calprotectinelevel) <- paste("prevalence", colnames(prev_counts.calprotectinelevel), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.calprotectinelevel <- dfcleanfeces_s2 %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("calprotectinelevel") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.calprotectinelevel) <- paste("prevalence", colnames(prev_possible.calprotectinelevel), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.calprotectinelevel/prev_possible.calprotectinelevel)*100
    
    tax_table.calprotectinelevel =  as.data.frame(tax_table(dfcleanfeces))
    
    merge.calprotectinelevel.prev= merge(tax_table.calprotectinelevel, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces_s2 %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "calprotectinelevel")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.calprotectinelevel.prev) <-merge.calprotectinelevel.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.calprotectinelevel.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    write.csv(calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq, "calprotectinelevelcorrcountrygenderagealphaantitrypsin_reduced_prevabdeseqDAL.csv")
    
    #plot
    
    # Order order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5), levels=names(x))
    # Family order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6), levels=names(x))
    # Genus order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7), levels=names(x))
    # Species order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8), levels=names(x))
    
    
    #make the actual plot
    
    pdf("differentially_present_taxa_Species_calprotectinelevel_corr_country_gender_age_aat_feces_LRTDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=rank8, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & calpro,  LRT model")
    dev.off()
    
    # make a boxplot of rel. abundance
    dfcleanfeces_s2.rel <- microbiome::transform(dfcleanfeces_s2, "compositional")
    data<-data.frame(otu_table(dfcleanfeces_s2.rel))
    data$calprotectinelevel<-sample_data(dfcleanfeces_s2.rel)$calprotectinelevel
    
    # make  graph plot
    calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq$rank6, "/", calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq$rank7, "/", calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq$rank8)
    pdf("differentially_present_SpeciesbycalprotectinelevelcorraatDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different Species by calprotectine level, corrected for anemie, stunting status, country of origin")
    dev.off()
    
    
    #### on Species/rank8 level and calprotectinelevel, Madagascar only  ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, calprotectinelevel!="")
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s2, pays!="RCA")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + stunted + anemie2 + calprotectinelevel)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + stunted + anemie2)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "calprotectinelevel") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    res_p_ordered_filt_2
    
    # Make an abundance and prevalence table at rank 8 /Species level
    #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
    prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
      x[x >= 2] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
      x[x >= 0] <- 1
      return(x)
    }
    
    prev_counts.calprotectinelevel <- dfcleanfeces_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("calprotectinelevel") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.calprotectinelevel) <- paste("prevalence", colnames(prev_counts.calprotectinelevel), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.calprotectinelevel <- dfcleanfeces_s2 %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("calprotectinelevel") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.calprotectinelevel) <- paste("prevalence", colnames(prev_possible.calprotectinelevel), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.calprotectinelevel/prev_possible.calprotectinelevel)*100
    
    tax_table.calprotectinelevel =  as.data.frame(tax_table(dfcleanfeces))
    
    merge.calprotectinelevel.prev= merge(tax_table.calprotectinelevel, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces_s2 %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "calprotectinelevel")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.calprotectinelevel.prev) <-merge.calprotectinelevel.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.calprotectinelevel.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    write.csv(calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq, "calprotectinelevelcorrcountrygenderagealphaantitrypsin_reduced_prevabdeseqDALMada.csv")
    
    #plot
    
    # Order order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5), levels=names(x))
    # Family order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6), levels=names(x))
    # Genus order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7), levels=names(x))
    # Species order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8), levels=names(x))
    
    
    #make the actual plot
    
    pdf("differentially_present_taxa_Species_calprotectinelevel_corr_country_gender_age_aat_feces_LRTDALMada.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=rank8, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & calpro,  LRT model")
    dev.off()
    
    # make a boxplot of rel. abundance
    dfcleanfeces_s2.rel <- microbiome::transform(dfcleanfeces_s2, "compositional")
    data<-data.frame(otu_table(dfcleanfeces_s2.rel))
    data$calprotectinelevel<-sample_data(dfcleanfeces_s2.rel)$calprotectinelevel
    
    # make  graph plot
    calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq$rank6, "/", calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq$rank7, "/", calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq$rank8)
    pdf("differentially_present_SpeciesbycalprotectinelevelcorraatDALMada.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(rank8, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different Species by calprotectine level, corrected for age, stunting status, Madagascar dataset")
    dev.off()
    
    
    
    
    
    #### on Species/rank8 level and calprotectinelevel, CAR only  ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, calprotectinelevel!="")
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s2, pays!="Madagascar")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + stunted + anemie2 + calprotectinelevel)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + stunted + anemie2)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "calprotectinelevel") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    res_p_ordered_filt_2
    
    # Make an abundance and prevalence table at rank 8 /Species level
    #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
    prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
      x[x >= 2] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
      x[x >= 0] <- 1
      return(x)
    }
    
    prev_counts.calprotectinelevel <- dfcleanfeces_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("calprotectinelevel") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.calprotectinelevel) <- paste("prevalence", colnames(prev_counts.calprotectinelevel), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.calprotectinelevel <- dfcleanfeces_s2 %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("calprotectinelevel") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.calprotectinelevel) <- paste("prevalence", colnames(prev_possible.calprotectinelevel), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.calprotectinelevel/prev_possible.calprotectinelevel)*100
    
    tax_table.calprotectinelevel =  as.data.frame(tax_table(dfcleanfeces))
    
    merge.calprotectinelevel.prev= merge(tax_table.calprotectinelevel, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces_s2 %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "calprotectinelevel")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.calprotectinelevel.prev) <-merge.calprotectinelevel.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.calprotectinelevel.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    write.csv(calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq, "calprotectinelevelcorrcountrygenderagealphaantitrypsin_reduced_prevabdeseqDALCAR.csv")
    
    #plot
    
    # Order order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5), levels=names(x))
    # Family order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6), levels=names(x))
    # Genus order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7), levels=names(x))
    # Species order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8), levels=names(x))
    
    
    #make the actual plot
    
    pdf("differentially_present_taxa_Species_calprotectinelevel_corr_country_gender_age_aat_feces_LRTDALCAR.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=rank8, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & calpro,  LRT model")
    dev.off()
    
    # make a boxplot of rel. abundance
    dfcleanfeces_s2.rel <- microbiome::transform(dfcleanfeces_s2, "compositional")
    data<-data.frame(otu_table(dfcleanfeces_s2.rel))
    data$calprotectinelevel<-sample_data(dfcleanfeces_s2.rel)$calprotectinelevel
    
    # make  graph plot
    calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq$rank6, "/", calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq$rank7, "/", calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq$rank8)
    pdf("differentially_present_SpeciesbycalprotectinelevelcorraatDALCAR.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=calprotectinelevelcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(rank8, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different Species by calprotectine level, corrected for age, stunting status, CAR dataset")
    dev.off()
    
    
    #### on Species/rank8 level and AAT level  ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, alphaantitrypsinlevel!="")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + stunted + anemie2 + pays + alphaantitrypsinlevel)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + stunted + anemie2 + pays)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "alphaantitrypsinlevel") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) # 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    res_p_ordered_filt_2
       
    
    #### on Species/rank8 level and AAT level, Madagascar only  ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, alphaantitrypsinlevel!="")
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s2, pays!="RCA")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + stunted + anemie2 + alphaantitrypsinlevel)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + stunted + anemie2)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "alphaantitrypsinlevel") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    res_p_ordered_filt_2
   
    #### on Species/rank8 level and AAT level, RCA only ####
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s, alphaantitrypsinlevel!="")
    dfcleanfeces_s2<-subset_samples(dfcleanfeces_s2, pays!="Madagascar")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_s2, ~  read_count + stunted + anemie2 + alphaantitrypsinlevel)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + stunted + anemie2)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "alphaantitrypsinlevel") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01)
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    res_p_ordered_filt_2
    # show results preserving taxa table
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_s)[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    # Make an abundance and prevalence table at rank 8 /Species level
    #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
    prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
      x[x >= 2] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
      x[x >= 0] <- 1
      return(x)
    }
    
    prev_counts.pays <- dfcleanfeces_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("alphaantitrypsinlevel") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.pays <- dfcleanfeces_s2 %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("alphaantitrypsinlevel") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.pays/prev_possible.pays)*100
    
    tax_table.pays =  as.data.frame(tax_table(dfcleanfeces))
    
    merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces_s2 %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "alphaantitrypsinlevel")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.pays.prev) <-merge.pays.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.pays.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    payscorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "CARAATlevelcorrcountrygenderage_reduced_prevabdeseqDAL.csv")
    
    #plot
    
    # Order order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5), levels=names(x))
    # Family order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6), levels=names(x))
    # Genus order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7), levels=names(x))
    # Species order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8), levels=names(x))
    
    
    #make the actual plot
    
    pdf("CARonlydifferentially_present_taxa_Species_AATlevel_corr_country_gender_age_feces_LRTDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=rank8, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("High vs. low AAT levels in feces from CAR, controlled for stunting, gender, age & country of origin LRT model")
    dev.off()
    
   
        
    # make  graph plot
    payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$rank6, "/", payscorrcountrygenderage_reduced_prevabdeseq$rank7, "/", payscorrcountrygenderage_reduced_prevabdeseq$rank8)
    pdf("CARonlydifferentially_present_SpeciesbyAATlevelDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(rank8, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different Species by AAT level, CAR only")
    dev.off()
    
    
    
 #### on Genus/rank7 level and country, not correcting for inflammation ####
    dfcleanfeces_g<-tax_glom(dfcleanfeces, "rank7")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g, ~  read_count + stunted + anemie2 + pays)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + stunted + anemie2)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    res_p_ordered_filt_2
    
    #### on Genus/rank7 level and country, correcting for calprotectine ####
    dfcleanfeces_g<-tax_glom(dfcleanfeces, "rank7")
    dfcleanfeces_g2<-subset_samples(dfcleanfeces_g, calprotectinelevel!="")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g2, ~  read_count + stunted + anemie2 + calprotectinelevel + pays)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + stunted + anemie2 + calprotectinelevel)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    res_p_ordered_filt_2
    
    # show results preserving taxa table
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_g2)[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    # Make an abundance and prevalence table at Rank7/Genus level
    #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
    prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
      x[x >= 2] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
      x[x >= 0] <- 1
      return(x)
    }
    
    prev_counts.pays <- dfcleanfeces_g2 %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.pays <- dfcleanfeces_g2 %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.pays/prev_possible.pays)*100
    
    tax_table.pays =  as.data.frame(tax_table(dfcleanfeces_g2))
    
    merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
    
    Afribiota_feces_abundance_OTU <- dfcleanfeces_g2 %>%
      transform_sample_counts(function(x) {
        x/sum(x)} ) 
    
    Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "pays")
    
    Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
    
    Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
    
    Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
      otu_table() %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    
    row.names(merge.pays.prev) <-merge.pays.prev$Row.names
    
    merge.abunanceprev.feces.OTU= merge(merge.pays.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
    merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
    row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
    
    # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
    payscorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrcountrygenderagecalpro_reduced_prevabdeseqDAL.csv")
    
    #plot
    
    # Order order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank4, function(x) max(x))
    x = sort(x, TRUE)
    # Class order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank5), levels=names(x))
    # Family order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank6), levels=names(x))
    # Genus order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank7), levels=names(x))
    # Species order
    x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8, function(x) max(x))
    x = sort(x, TRUE)
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8 = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$rank8), levels=names(x))
    
    
    #make the actual plot
    
    pdf("differentially_present_taxa_Genus_pays_corr_country_gender_age_calpro_feces_LRTDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 8)
    ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=rank7, y=log2FoldChange, color=rank4)) + geom_point(size=2) + 
      theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & calprotectine levels, LRT model")
    dev.off()
    
    # make  graph plot
       pdf("differentially_present_GenusbypayscontrolledcalproDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x = reorder(rank7, log2FoldChange), y=log2FoldChange)) +
      geom_bar(position="dodge",stat="identity", color="black") +
      coord_flip() + 
      theme(axis.text.x = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("log2 fold change")+
      ggtitle("Significantly different Genus by country of origin, controlled by calprotectine levels")
    dev.off()
    
    
    #### on Genus/rank7 level and country,  correcting for alphaantitrypsin ####
    dfcleanfeces_g<-tax_glom(dfcleanfeces, "rank7")
    dfcleanfeces_g2<-subset_samples(dfcleanfeces_g, alphaantitrypsinlevel!="")
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g2, ~ read_count +  stunted + anemie2 + alphaantitrypsinlevel + pays)
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ stunted + anemie2 + alphaantitrypsinlevel)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    # show results preserving taxa table
    res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_g2)[rownames(res_p_ordered_filt_2), ], "matrix"))
    
    
    
    #### on Genus/rank7 level and stunting, not correcting for inflammation ####
    dfcleanfeces_g<-tax_glom(dfcleanfeces, "rank7")
    dfcleanfeces_g<-subset_samples(dfcleanfeces_g, stunted!="")
    dfcleanfeces_g<-subset_samples(dfcleanfeces_g, pays!="")
    dfcleanfeces_g<-subset_samples(dfcleanfeces_g, ageyears!="")
    dfcleanfeces_g<-subset_samples(dfcleanfeces_g, sexe!="")
    dfcleanfeces_g = filter_taxa(dfcleanfeces_g, function(x) sum(x) > 0, TRUE)
    dfcleanfeces_g<-subset_samples(dfcleanfeces_g, (rowSums(otu_table(dfcleanfeces))!=0))
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g, ~  read_count + pays + anemie2 + stunted)
    
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + pays  + anemie2)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) # 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    View(res_p_ordered_filt_2)
    
    
    #### on Genus/rank7 level and stunting, not correcting for inflammation, non-filtered dataset  ####
    dfcleanfeces_g<-tax_glom(dfcleanfeces, "rank7")
    dfcleanfeces_g<-subset_samples(dfcleanfeces_g, stunted!="")
    dfcleanfeces_g<-subset_samples(dfcleanfeces_g, pays!="")
    dfcleanfeces_g<-subset_samples(dfcleanfeces_g, ageyears!="")
    dfcleanfeces_g<-subset_samples(dfcleanfeces_g, sexe!="")
    dfcleanfeces_g<-subset_samples(dfcleanfeces_g, anemie2!="")
    dfcleanfeces_g = filter_taxa(dfcleanfeces_g, function(x) sum(x) > 0, TRUE)
    dfcleanfeces_g<-subset_samples(dfcleanfeces_g, (rowSums(otu_table(dfcleanfeces_g))!=0))
    dfcleanfeces_g<-subset_samples(dfcleanfeces_g, (colSums(otu_table(dfcleanfeces_g))!=0))
    
    diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g, ~  pays + anemie2 + stunted)
    
    gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
    geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
    diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
    
    diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ pays + anemie2)
    
    resultsNames(diagddsfecesDM)
    res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") #log2( non_stunted / stunted ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
    
    # view a summary of the results table with a padj value < 0.01
    summary(res, alpha = 0.01) 
    
    #filtering the results table #
    # reorder the results table by adjusted p-value and remove any "NA" entries
    res_p_ordered <- res[order(res$padj, na.last = NA), ]
    # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
    res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
    res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
    View(res_p_ordered_filt_2)
    
    
  #### Make chi2 test on the presence/absence of given bacteria not correcting for confounding factors ####    
    #### On cluster level no correction for country ####
    dffiltered_feces_cluster<-dfcleanfeces
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -11] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    fecesprevalence <- dffiltered_feces_cluster%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    
    fecesprevalence2<-subset_samples(fecesprevalence, stunted!="") 
    fecesprevalence2<-subset_samples(fecesprevalence2, pays!="") # 252 samples remaining, no additional ones taken out
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Blastocystis")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$pays<-as.factor(meta_chi2$pays)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence) #take taxtable
    clusternames<-taxa_chi2[, 12]
    dim(clusternames)
    
    chi2results = data.frame(clusternames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig

    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ClusterpaysDAL.csv")
    
    
    #### On cluster level Blastocytis Bangui- stunted ####
    dfcleanfecesBangui<-subset_samples(dfcleanfeces, pays=="RCA")
    dffiltered_feces_cluster<-dfcleanfecesBangui
    
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -11] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    fecesprevalence <- dffiltered_feces_cluster%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    
    fecesprevalence2<-subset_samples(fecesprevalence, stunted!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Blastocystis")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$stunted<-as.factor(meta_chi2$stunted)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    clusternames<-taxa_chi2[, 12]
    dim(clusternames)
    
    chi2results = data.frame(clusternames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ClusterstuntedDALBangui.csv")
    
    
    #### On cluster level Blastocytis Bangui- AAT ####
    dfcleanfecesBangui<-subset_samples(dfcleanfeces, pays=="RCA")
    dffiltered_feces_cluster<-dfcleanfecesBangui
    
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -11] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    fecesprevalence <- dffiltered_feces_cluster%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, alphaantitrypsinlevel!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Blastocystis")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$alphaantitrypsinlevel<-as.factor(meta_chi2$alphaantitrypsinlevel)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    clusternames<-taxa_chi2[, 12]
    dim(clusternames)
    
    chi2results = data.frame(clusternames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ClusterAATDALBangui.csv")
    
    
    #### On cluster level Blastocytis Bangui-  Calpro ####
    dfcleanfecesBangui<-subset_samples(dfcleanfeces, pays=="RCA")
    dffiltered_feces_cluster<-dfcleanfecesBangui
    
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -11] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    fecesprevalence <- dffiltered_feces_cluster%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, calprotectinelevel!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Blastocystis")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    clusternames<-taxa_chi2[, 12]
    dim(clusternames)
    
    chi2results = data.frame(clusternames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ClustercalproDALBangui.csv")
    
    
    #### On cluster level Blastocytis Bangui-  Anemia ####
    dfcleanfecesBangui<-subset_samples(dfcleanfeces, pays=="RCA")
    dffiltered_feces_cluster<-dfcleanfecesBangui
    
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -11] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    fecesprevalence <- dffiltered_feces_cluster%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, anemie2!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Blastocystis")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    clusternames<-taxa_chi2[, 12]
    dim(clusternames)
    
    chi2results = data.frame(clusternames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ClustercalproanemieBangui.csv")
    
    
    #### On cluster level Blastocytis Bangui-  Age ####
    dfcleanfecesBangui<-subset_samples(dfcleanfeces, pays=="RCA")
    dffiltered_feces_cluster<-dfcleanfecesBangui
    
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -11] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    fecesprevalence <- dffiltered_feces_cluster%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, ageyears!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Blastocystis")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    clusternames<-taxa_chi2[, 12]
    dim(clusternames)
    
    chi2results = data.frame(clusternames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ClustercalproageBangui.csv")
    
    
    #### On cluster level Blastocytis Tana- stunted ####
    dfcleanfecesTana<-subset_samples(dfcleanfeces, pays=="Madagascar")
    dffiltered_feces_cluster<-dfcleanfecesTana
    
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -11] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    fecesprevalence <- dffiltered_feces_cluster%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    
    fecesprevalence2<-subset_samples(fecesprevalence, stunted!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Blastocystis")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$stunted<-as.factor(meta_chi2$stunted)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    clusternames<-taxa_chi2[, 12]
    dim(clusternames)
    
    chi2results = data.frame(clusternames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ClusterstuntedDALTana.csv")
    
    
    #### On cluster level Blastocytis Tana- AAT ####
    dfcleanfecesTana<-subset_samples(dfcleanfeces, pays=="Madagascar")
    dffiltered_feces_cluster<-dfcleanfecesTana
    
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -11] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    fecesprevalence <- dffiltered_feces_cluster%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, alphaantitrypsinlevel!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Blastocystis")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$alphaantitrypsinlevel<-as.factor(meta_chi2$alphaantitrypsinlevel)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    clusternames<-taxa_chi2[, 12]
    dim(clusternames)
    
    chi2results = data.frame(clusternames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ClusterAATDALTana.csv")
    
    
    #### On cluster level Blastocytis Tana-  Calpro ####
    dfcleanfecesTana<-subset_samples(dfcleanfeces, pays=="Madagascar")
    dffiltered_feces_cluster<-dfcleanfecesTana
    
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -11] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    fecesprevalence <- dffiltered_feces_cluster%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, calprotectinelevel!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Blastocystis")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    clusternames<-taxa_chi2[, 12]
    dim(clusternames)
    
    chi2results = data.frame(clusternames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ClustercalproDALTana.csv")
    
    
    #### On cluster level Blastocytis Tana-  Anemia ####
    dfcleanfecesTana<-subset_samples(dfcleanfeces, pays=="Madagascar")
    dffiltered_feces_cluster<-dfcleanfecesTana
    
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -11] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    fecesprevalence <- dffiltered_feces_cluster%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, anemie2!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Blastocystis")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    clusternames<-taxa_chi2[, 12]
    dim(clusternames)
    
    chi2results = data.frame(clusternames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ClustercalproanemieTana.csv")
    
    
    #### On cluster level Blastocytis Tana-  Age ####
    dfcleanfecesTana<-subset_samples(dfcleanfeces, pays=="Madagascar")
    dffiltered_feces_cluster<-dfcleanfecesTana
    
    tax_table(dffiltered_feces_cluster)= tax_table(dffiltered_feces_cluster)[, -11] # need to take away accession number as this does not let collapse correctly
    dffiltered_feces_cluster= tax_glom(dffiltered_feces_cluster, "cluster")
    dffiltered_feces_cluster<-subset_taxa(dffiltered_feces_cluster, cluster!="")
    
    fecesprevalence <- dffiltered_feces_cluster%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, ageyears!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Blastocystis")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    clusternames<-taxa_chi2[, 12]
    dim(clusternames)
    
    chi2results = data.frame(clusternames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ClustercalproageTana.csv")
    
    
    #### On cluster level Entamoeba Bangui- stunted ####
    dfcleanfecesBangui<-subset_samples(dfcleanfeces, pays=="RCA")
    dfcleanfecesBanguispecies<-tax_glom(dfcleanfecesBangui, "rank8")
    dfcleanfecesBanguispeciesEntamoeba<-subset_taxa(dfcleanfecesBanguispecies, rank7=="Entamoeba")
    
    
    tax_table(dfcleanfecesBanguispeciesEntamoeba)= tax_table(dfcleanfecesBanguispeciesEntamoeba)[, -11] # need to take away accession number as this does not let collapse correctly
    dfcleanfecesBanguispeciesEntamoeba= tax_glom(dfcleanfecesBanguispeciesEntamoeba, "cluster")
    dfcleanfecesBanguispeciesEntamoeba<-subset_taxa(dfcleanfecesBanguispeciesEntamoeba, cluster!="")
    
    fecesprevalence <- dfcleanfecesBanguispeciesEntamoeba%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    
    fecesprevalence2<-subset_samples(fecesprevalence, stunted!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Entamoeba")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$stunted<-as.factor(meta_chi2$stunted)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    rank8names<-taxa_chi2[, 8]
    dim(rank8names)
    
    chi2results = data.frame(rank8names,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebastuntedDALBangui.csv")
    
    
    #### On cluster level Entamoeba Bangui- AAT ####
    
    fecesprevalence <- dfcleanfecesBanguispeciesEntamoeba%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, alphaantitrypsinlevel!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Entamoeba")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$alphaantitrypsinlevel<-as.factor(meta_chi2$alphaantitrypsinlevel)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    rank8names<-taxa_chi2[, 8]
    dim(rank8names)
    
    chi2results = data.frame(rank8names,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebaAATDALBangui.csv")
    
    
    #### On cluster level Entamoeba Bangui-  Calpro ####
    
    fecesprevalence <- dfcleanfecesBanguispeciesEntamoeba%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, calprotectinelevel!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Entamoeba")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    rank8names<-taxa_chi2[, 8]
    dim(rank8names)
    
    chi2results = data.frame(rank8names,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebacalproDALBangui.csv")
    
    
    #### On cluster level Entamoeba Bangui-  Anemia ####
    
    fecesprevalence <- dfcleanfecesBanguispeciesEntamoeba%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, anemie2!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Entamoeba")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    rank8names<-taxa_chi2[, 8]
    dim(rank8names)
    
    chi2results = data.frame(rank8names,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebacalproanemieBangui.csv")
    
    
    #### On cluster level Entamoeba Bangui-  Age ####
    
    fecesprevalence <- dfcleanfecesBanguispeciesEntamoeba%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, ageyears!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Entamoeba")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    rank8names<-taxa_chi2[, 12]
    dim(rank8names)
    
    chi2results = data.frame(rank8names,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebacalproageBangui.csv")
    
    
    
    
    #### On cluster level Entamoeba Tana- stunted ####
    dfcleanfecesTana<-subset_samples(dfcleanfeces, pays=="Madagascar")
    dfcleanfecesTanaspecies<-tax_glom(dfcleanfecesTana, "rank8")
    dfcleanfecesTanaspeciesEntamoeba<-subset_taxa(dfcleanfecesTanaspecies, rank7=="Entamoeba")
    
    
    tax_table(dfcleanfecesTanaspeciesEntamoeba)= tax_table(dfcleanfecesTanaspeciesEntamoeba)[, -11] # need to take away accession number as this does not let collapse correctly
    dfcleanfecesTanaspeciesEntamoeba= tax_glom(dfcleanfecesTanaspeciesEntamoeba, "cluster")
    dfcleanfecesTanaspeciesEntamoeba<-subset_taxa(dfcleanfecesTanaspeciesEntamoeba, cluster!="")
    
    fecesprevalence <- dfcleanfecesTanaspeciesEntamoeba%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    
    fecesprevalence2<-subset_samples(fecesprevalence, stunted!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Entamoeba")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$stunted<-as.factor(meta_chi2$stunted)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    rank8names<-taxa_chi2[, 8]
    dim(rank8names)
    
    chi2results = data.frame(rank8names,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebastuntedDALTana.csv")
    
    
    #### On cluster level Entamoeba Tana- AAT ####
    
    fecesprevalence <- dfcleanfecesTanaspeciesEntamoeba%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, alphaantitrypsinlevel!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Entamoeba")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$alphaantitrypsinlevel<-as.factor(meta_chi2$alphaantitrypsinlevel)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    rank8names<-taxa_chi2[, 8]
    dim(rank8names)
    
    chi2results = data.frame(rank8names,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebaAATDALTana.csv")
    
    
    #### On cluster level Entamoeba Tana-  Calpro ####
    
    fecesprevalence <- dfcleanfecesTanaspeciesEntamoeba%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, calprotectinelevel!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Entamoeba")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    rank8names<-taxa_chi2[, 8]
    dim(rank8names)
    
    chi2results = data.frame(rank8names,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebacalproDALTana.csv")
    
    
    #### On cluster level Entamoeba Tana-  Anemia ####
    
    fecesprevalence <- dfcleanfecesTanaspeciesEntamoeba%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, anemie2!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Entamoeba")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    rank8names<-taxa_chi2[, 8]
    dim(rank8names)
    
    chi2results = data.frame(rank8names,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebacalproanemieTana.csv")
    
    
    #### On cluster level Entamoeba Tana-  Age ####
    
    fecesprevalence <- dfcleanfecesTanaspeciesEntamoeba%>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    colnames(sample_data(fecesprevalence))
    
    fecesprevalence2<-subset_samples(fecesprevalence, ageyears!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2 = subset_taxa(fecesprevalence2, rank7=="Entamoeba")
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence2) #take taxtable
    rank8names<-taxa_chi2[, 12]
    dim(rank8names)
    
    chi2results = data.frame(rank8names,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    chi2results_sig
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.EntamoebacalproageTana.csv")
    
    
    #### On ASV level no correction for stunting ####
    fecesprevalence <- dfcleanfeces %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    
    fecesprevalence2<-subset_samples(fecesprevalence, stunted!="") 
    fecesprevalence2<-subset_samples(fecesprevalence2, pays!="") 
    fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
    fecesprevalence2
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
    meta_chi2$stunted<-as.factor(meta_chi2$stunted)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
    ASVnames<-colnames(df_chi2)
    chi2results = data.frame(ASVnames,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence) #take taxtable
    ASVnames<-taxa_chi2[, 8]
    dim(ASVnames)
    
    chi2results = data.frame(ASVnames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.stuntedASVnoncorrectedchi2DAL.csv")
    
    
    
    
    #### On ASV level no correction for country ####
    meta_chi2$pays<-as.factor(meta_chi2$pays)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value)
    ASVnames<-colnames(df_chi2)
    chi2results = data.frame(ASVnames,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence) #take taxtable
    ASVnames<-taxa_chi2[, 8]
    dim(ASVnames)
    
    chi2results = data.frame(ASVnames,chi2results)
    
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.paysASVnoncorrectedchi2DAL.csv")
    
    #### On ASV level no correction for ageyears in three levels ####
    meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
    ASVnames<-colnames(df_chi2)
    chi2results = data.frame(ASVnames,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence) #take taxtable
    ASVnames<-taxa_chi2[, 8]
    dim(ASVnames)
    
    chi2results = data.frame(ASVnames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsASVnoncorrectedchi2DAL.csv")
    
    #### On ASV level no correction for calprotectinelevel ####
    meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
    ASVnames<-colnames(df_chi2)
    chi2results = data.frame(ASVnames,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence) #take taxtable
    ASVnames<-taxa_chi2[, 8]
    dim(ASVnames)
    
    chi2results = data.frame(ASVnames,chi2results)
    
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.calprotectinelevelASVnoncorrectedchi2DAL.csv")
    
    #### On ASV level no correction for alphaantitrypsinlevel ####
    meta_chi2$alphaantitrypsinlevel<-as.factor(meta_chi2$alphaantitrypsinlevel)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value)
    ASVnames<-colnames(df_chi2)
    chi2results = data.frame(ASVnames,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence) #take taxtable
    ASVnames<-taxa_chi2[, 8]
    dim(ASVnames)
    
    chi2results = data.frame(ASVnames,chi2results)
    
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.alphaantitrypsinlevelASVnoncorrectedchi2DAL.csv")
    
    #### On ASV level no correction for anemie2 ####
    meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
    ASVnames<-colnames(df_chi2)
    chi2results = data.frame(ASVnames,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(fecesprevalence) #take taxtable
    ASVnames<-taxa_chi2[, 8]
    dim(ASVnames)
    
    chi2results = data.frame(ASVnames,chi2results)

    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.anemie2ASVnoncorrectedchi2DAL.csv")
    
    #### On rank8 level no correction for stunting ####
    dfcleanfeces_s<- tax_glom(fecesprevalence, "rank8")
    fecesprevalence_s <- dfcleanfeces_s %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    
    fecesprevalence_s2<-subset_samples(fecesprevalence_s, stunted!="") # 575 samples remaining
    fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
    fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence_s2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
    meta_chi2$stunted<-as.factor(meta_chi2$stunted)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
      
    #### On rank8 level no correction for country ####
    
    meta_chi2$pays<-as.factor(meta_chi2$pays)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(dfcleanfeces_s) #take taxtable
    ASVnames<-taxa_chi2[, 8]
    dim(ASVnames)
    
    chi2results = data.frame(ASVnames,chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.paysrank8noncorrectedchi2DAL.csv")
    
    #### On rank8 level no correction for ageyears ####
    meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(dfcleanfeces_s) #take taxtable
    ASVnames<-taxa_chi2[, 8]
    dim(ASVnames)
    
    chi2results = data.frame(ASVnames,chi2results)

    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsrank8noncorrectedchi2DAL.csv")
    
    #### On rank8 level no correction for calprotectine levels ####
    meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(dfcleanfeces_s) #take taxtable
    ASVnames<-taxa_chi2[, 8]
    dim(ASVnames)
    
    chi2results = data.frame(ASVnames,chi2results)
    
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.calprotectinelevelrank8noncorrectedchi2DAL.csv")
    
    
    #### On rank8 level no correction for alphaantitrypsin levels ####
    meta_chi2$alphaantitrypsinlevel<-as.factor(meta_chi2$alphaantitrypsinlevel)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(dfcleanfeces_s) #take taxtable
    ASVnames<-taxa_chi2[, 8]
    dim(ASVnames)
    
    chi2results = data.frame(ASVnames,chi2results)

    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.alphaantitrypsinlevelrank8noncorrectedchi2DAL.csv")
    
    
    #### On rank8 level no correction for anemie ####
    meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
    
    Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
    rank8names<-colnames(df_chi2)
    chi2results = data.frame(rank8names,Chi2.p)
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    chi2results$name<-row.names(chi2results)
    
    # Merge with tax info
    taxa_chi2<- tax_table(dfcleanfeces_s) #take taxtable
    ASVnames<-taxa_chi2[, 8]
    dim(ASVnames)
    
    chi2results = data.frame(ASVnames,chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.anemie2rank8noncorrectedchi2DAL.csv")
    
    #### On rank7 level no correction for stunting ####
    dfcleanfeces_g<- tax_glom(dfcleanfeces, "rank7")
    fecesprevalence_g <- dfcleanfeces_g %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence_g))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    
    fecesprevalence_g2<-subset_samples(fecesprevalence_g, stunted!="") # 190 samples remaining
    fecesprevalence_g2<-subset_samples(fecesprevalence_g2, pays!="") # 190 samples remaining, no additional ones taken out
    fecesprevalence_g2 = filter_taxa(fecesprevalence_g2, function(x) sum(x) > 0, TRUE)
    
    df_chi2 <- as.matrix(t(otu_table(fecesprevalence_g2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_g2)) #take metadata
    meta_chi2$stunted<-as.factor(meta_chi2$stunted)
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
   
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.stuntedrank7noncorrectedchi2DAL.csv")

    
    #### On rank7 level no correction for country ####
    meta_chi2$pays<-as.factor(meta_chi2$pays)
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.paysrank7noncorrectedchi2DAL.csv")
    
    #### On rank7 level no correction for ageyears ####
    meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
    
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    
    
    write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsrank7noncorrectedchi2DAL.csv")
    

    #### On rank7 level no correction for calprotectinelevel ####
    meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
    
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    
    
    write.csv(chi2results_sig,"Chi2resultsFeces.calprotectinelevelrank7noncorrectedchi2DAL.csv")
    
    #### On rank7 level no correction for alphaantitrypsinlevel ####
    meta_chi2$alphaantitrypsinlevel<-as.factor(meta_chi2$alphaantitrypsinlevel)
    
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)

    write.csv(chi2results_sig,"Chi2resultsFeces.alphaantitrypsinlevelrank7noncorrectedchi2DAL.csv")
    
    
    #### On rank7 level no correction for anemie2 ####
    meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
    
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    
    
    write.csv(chi2results_sig,"Chi2resultsFeces.anemie2rank7noncorrectedchi2DAL.csv")
    
    #### On rank7 level no correction for stunting- BANGUI ####
    dfcleanfecesBangui<-subset_samples(dfcleanfeces, pays=="RCA")
    dfcleanfeces_g<- tax_glom(dfcleanfecesBangui, "rank7")
    fecesprevalence_g <- dfcleanfeces_g %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence_g))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    
    fecesprevalence_g2<-subset_samples(fecesprevalence_g, stunted!="") # 91 samples remaining
    fecesprevalence_g2<-subset_samples(fecesprevalence_g2, pays!="") 
    fecesprevalence_g2 = filter_taxa(fecesprevalence_g2, function(x) sum(x) > 0, TRUE)
    
    df_chi2 <- as.matrix(t(otu_table(fecesprevalence_g2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_g2)) #take metadata
    meta_chi2$stunted<-as.factor(meta_chi2$stunted)
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.stuntedrank7noncorrectedchi2DALBangui.csv")
    
    
    #### On rank7 level no correction for ageyears- BANGUI ####
    meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
    
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    
    
    write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsrank7noncorrectedchi2DALBangui.csv")
    
    
    #### On rank7 level no correction for calprotectinelevel- BANGUI  ####
    meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
    
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    
    
    write.csv(chi2results_sig,"Chi2resultsFeces.calprotectinelevelrank7noncorrectedchi2DALBangui.csv")
    
    #### On rank7 level no correction for alphaantitrypsinlevel- BANGUI  ####
    meta_chi2$alphaantitrypsinlevel<-as.factor(meta_chi2$alphaantitrypsinlevel)
    
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    
    write.csv(chi2results_sig,"Chi2resultsFeces.alphaantitrypsinlevelrank7noncorrectedchi2DALBangui.csv")
    
    
    #### On rank7 level no correction for anemie2- BANGUI  ####
    meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
    
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    
    
    write.csv(chi2results_sig,"Chi2resultsFeces.anemie2rank7noncorrectedchi2DALBangui.csv")
    
    
    #### On rank7 level no correction for stunting- MADAGASCAR ####
    dfcleanfecesBangui<-subset_samples(dfcleanfeces, pays=="Madagascar")
    dfcleanfeces_g<- tax_glom(dfcleanfecesBangui, "rank7")
    fecesprevalence_g <- dfcleanfeces_g %>% #this produces prevalence "counts" for each species, but not percentages
      transform_sample_counts(fun = prevalence) 
    
    View(head(otu_table(fecesprevalence_g))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
    
    fecesprevalence_g2<-subset_samples(fecesprevalence_g, stunted!="") # 91 samples remaining
    fecesprevalence_g2<-subset_samples(fecesprevalence_g2, pays!="") 
    fecesprevalence_g2 = filter_taxa(fecesprevalence_g2, function(x) sum(x) > 0, TRUE)
    
    df_chi2 <- as.matrix(t(otu_table(fecesprevalence_g2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_g2)) #take metadata
    meta_chi2$stunted<-as.factor(meta_chi2$stunted)
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    write.csv(chi2results_sig,"Chi2resultsFeces.stuntedrank7noncorrectedchi2DALMada.csv")
    
    
    #### On rank7 level no correction for ageyears- MADAGASCAR ####
    meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
    
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    
    
    write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsrank7noncorrectedchi2DALMada.csv")
    
    
    #### On rank7 level no correction for calprotectinelevel- MADAGASCAR  ####
    meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
    
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    
    
    write.csv(chi2results_sig,"Chi2resultsFeces.calprotectinelevelrank7noncorrectedchi2DALMada.csv")
    
    #### On rank7 level no correction for alphaantitrypsinlevel- MADAGASCAR  ####
    meta_chi2$alphaantitrypsinlevel<-as.factor(meta_chi2$alphaantitrypsinlevel)
    
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$alphaantitrypsinlevel), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    
    write.csv(chi2results_sig,"Chi2resultsFeces.alphaantitrypsinlevelrank7noncorrectedchi2DALMada.csv")
    
    
    #### On rank7 level no correction for anemie2- MADAGASCAR  ####
    meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
    
    taxa_chi2<- tax_table(fecesprevalence_g2) #take taxtable
    
    chi2results<-as.data.frame(apply(df_chi2, 1, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value))
    colnames(chi2results)<-"Chi2.p"
    chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
    
    # Merge with tax info
    Rank7names<-taxa_chi2[, 7]
    dim(Rank7names)
    dim(chi2results)
    
    chi2results = data.frame(Rank7names, chi2results)
    chi2results_sig <- dplyr::filter(chi2results, chi2results$rel.fdr<0.05)
    
    View(chi2results_sig)
    
    
    write.csv(chi2results_sig,"Chi2resultsFeces.anemie2rank7noncorrectedchi2DALMada.csv")
    
    
#### logistic regression on the presence/absence of given bacteria correcting for confounding factors ####       
    #### rank8 level: logistic model for stunting correcting for covariables None associated with multiple correction #####
    fecesprevalence_s2<-transform(dfcleanfeces_s, )
    df_chi2 <- as.matrix((otu_table(fecesprevalence_s2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
    data1 <- df_chi2
    taxa_chi2<-tax_table(fecesprevalence_s2)
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_chi2$age
    df$sexe <- meta_chi2$sexe
    df$Country <- meta_chi2$pays
    df$stunted <- meta_chi2$stunted
    df$calpro <- meta_chi2$calprotectinelevel
    df$aat <- meta_chi2$alphaantitrypsinlevel
    df$anemie <- meta_chi2$anemie2
   
    df$totalreads<-meta_chi2$read_count
    
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted",  "age", "sexe", "Country", "totalreads", "aat", "calpro", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value+  age  + totalreads + calpro + anemie , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    row.names(logresults)<-logresults$variable
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")

    Rank8names<-data.frame(taxa_chi2[, 8])
    Rank8names$variable<-row.names(Rank8names)
    logresults = merge(Rank8names, logresults)
    View(logresults)

    write.csv(logresults,"LogresultsFecespresabs.stuntedSpecieswithinflaDAL.csv")
    
    #### rank8 level: logistic model for country correcting for covariables#####
    df_chi2 <- as.matrix((otu_table(fecesprevalence_s2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
    data1 <- df_chi2
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$run <- meta_chi2$run
    df$age <- meta_chi2$age
    df$sexe <- meta_chi2$sexe
    df$Country <- meta_chi2$pays
    df$stunted <- meta_chi2$stunted
    df$calpro <- meta_chi2$calprotectinelevel
    df$aat <- meta_chi2$alphaantitrypsinlevel
    df$anemie<-meta_chi2$anemie2
       
    df$totalreads<-meta_chi2$read_count
    
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "Country", "totalreads", "microeukreads", "aat", "calpro", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, Country!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(Country ~ value+ stunted + totalreads + calpro +  anemie , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    logresults = merge(Rank8names, logresults)
    
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    write.csv(logresults,"LogresultsFecespresabs.paysSpecieswithinflaDAL.csv")
    
    prev_counts.pays <- dfcleanfeces_s %>% #this produces prevalence "counts" for each OTU, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_possible.pays <- dfcleanfeces_s %>% #this produces a maximum possible prevalence count per OTU
      transform_sample_counts(fun = allones) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
    test.prev = (prev_counts.pays/prev_possible.pays)*100
    
    # filter out any results which have a adjust p-value less than alpha (0.05)
    logresults_sig <- logresults[which(logresults$rel.fdr < 0.05), ]
    
    # filter out from the prevalence and absence table the unneeded lines
    test.prev2<-filter (test.prev, row.names(test.prev) %in% logresults_sig$variable)
    
    logisticregression_prevabdeseq <- cbind(as(logresults_sig, "data.frame"), as(test.prev2, "matrix"))
    logisticregression_prevabdeseq$taxonomy<- paste0(logisticregression_prevabdeseq$rank6, "/", logisticregression_prevabdeseq$rank7, "/", logisticregression_prevabdeseq$rank8)
    
    write.csv(logisticregression_prevabdeseq, "logisticregression_prevabdeseq_countrySpeciesDAL.csv")
    
    # Grouped Bar Plot
    data<-data_frame(logisticregression_prevabdeseq$taxonomy, logisticregression_prevabdeseq$prevalence.Madagascar, logisticregression_prevabdeseq$prevalence.RCA)
    View(data)
    colnames(data)<-c("taxonomy", "Madagascar","CAR")
    df2 = melt(data, id.vars = c('taxonomy'), 
               variable.name ='country', value.name = "prevalence")
    
    pdf("differentially_present_SpeciesbypayspresabsDAL.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    ggplot(data=df2, aes(x = taxonomy, y=prevalence, fill=factor(Country)), group = country) +
      geom_bar(aes(fill = country), position = "dodge", stat="identity") +
      theme(axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1 ))+ 
      theme(axis.text.y = element_text(size=12))+ 
      theme(axis.text.y = element_text(size=12))+ 
      xlab("")+ 
      theme(axis.title.y = element_text(size=14))+
      theme(title = element_text(size=16, face="bold"))+
      ylab("Prevalence")+
      ggtitle("Prevalence of Species by country of origin")
    dev.off()
    
    
    
    #### rank8 level: logistic model for stunting correcting for covariables only in Bangui  #####
    fecesprevalence_s2b<-subset_samples(fecesprevalence_s2, pays=="RCA")
    df_chi2 <- as.matrix((otu_table(fecesprevalence_s2b))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_s2b)) #take metadata
    data1 <- df_chi2
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_chi2$age
    df$sexe <- meta_chi2$sexe
    df$stunted <- meta_chi2$stunted
    df$calpro <- meta_chi2$calprotectinelevel
    df$aat <- meta_chi2$alphaantitrypsinlevel
    df$anemie<-meta_chi2$anemie2
    
    df$totalreads<-meta_chi2$read_count
    
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "totalreads", "aat", "calpro", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    long=filter(long, calpro!="") ## keep only the ones with valid data 
    long=filter(long, aat!="") ## keep only the ones with valid data 
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value+ totalreads + calpro + anemie , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    taxa_chi2<-tax_table(fecesprevalence_s2b)
    row.names(logresults)<-logresults$variable
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    Rank8names<-data.frame(taxa_chi2[, 8])
    Rank8names$variable<-row.names(Rank8names)
    logresults = merge(Rank8names, logresults)
    View(logresults)
    
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    write.csv(logresults,"LogresultsFecespresabs.stuntedSpecieswithinflaBanguiDAL.csv")
    
    
    #### rank8 level: logistic model for stunting correcting for covariables only in Tana NOTHING ASSOCIATED WITH MULTIPLE CORRECTION #####
    fecesprevalence_s2a<-subset_samples(fecesprevalence_s2, pays=="Madagascar")
    df_chi2 <- as.matrix((otu_table(fecesprevalence_s2a))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_s2a)) #take metadata
    data1 <- df_chi2
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_chi2$age
    df$run <- meta_chi2$run
    df$sexe <- meta_chi2$sexe
    df$stunted <- meta_chi2$stunted
    df$calpro <- meta_chi2$calprotectinelevel
    df$aat <- meta_chi2$alphaantitrypsinlevel
    df$anemie <- meta_chi2$anemie2
    
    df$totalreads<-meta_chi2$read_count
    
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "totalreads",  "aat", "calpro", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value+  totalreads + calpro + anemie , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    taxa_chi2<-tax_table(fecesprevalence_s2a)
    row.names(logresults)<-logresults$variable
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    Rank8names<-data.frame(taxa_chi2[, 8])
    Rank8names$variable<-row.names(Rank8names)
    logresults = merge(Rank8names, logresults)
    View(logresults)
    
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    write.csv(logresults,"LogresultsFecespresabs.stuntedSpecieswithinflaTanaDAL.csv")
    
    
    #### rank8 level: logistic model for age (more or less than three years) correcting for covariables  #####
    #create a variable age in years which has only two modalities
    sample_data(fecesprevalence_s2)$ageyears2<-cut(as.numeric(as.character(sample_data(fecesprevalence_s2)$age)), c(24,36,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
    which(is.na(sample_data(fecesprevalence_s2)$ageyears2)) # the controls have no assocaited age year
    levels(sample_data(fecesprevalence_s2)$ageyears2) <- c("2-3 years", "3-5 years")
    levels(sample_data(fecesprevalence_s2)$ageyears2)
    
    df_chi2 <- as.matrix((otu_table(fecesprevalence_s2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
    data1 <- df_chi2
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_chi2$ageyears2
    df$sexe <- meta_chi2$sexe
    df$run<-meta_chi2$run
    df$Country <- meta_chi2$pays
    df$stunted <- meta_chi2$stunted
    df$calpro <- meta_chi2$calprotectinelevel
    df$aat <- meta_chi2$alphaantitrypsinlevel
    df$anemie<-meta_chi2$anemie2
    
    df$totalreads<-meta_chi2$read_count
    
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted",  "age", "sexe", "Country", "totalreads", "aat", "calpro", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, age!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(age ~ value+ Country +  stunted + totalreads + calpro + anemie , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    taxa_chi2<-tax_table(fecesprevalence_s2a)
    row.names(logresults)<-logresults$variable
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    Rank8names<-data.frame(taxa_chi2[, 8])
    Rank8names$variable<-row.names(Rank8names)
    logresults = merge(Rank8names, logresults)
    View(logresults)
    
   
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    write.csv(logresults,"LogresultsFecespresabs.agetwocatSpecieswithinflaDAL.csv")
    
    #### rank7 level: logistic model for stunting correcting for covariables None associated when correcting for multiple testing #####
    df_chi2 <- as.matrix((otu_table(fecesprevalence_g2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_g2)) #take metadata
    data1 <- df_chi2
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_chi2$age
    df$sexe <- meta_chi2$sexe
    df$Country <- meta_chi2$pays
    df$stunted <- meta_chi2$stunted
    df$calpro <- meta_chi2$calprotectinelevel
    df$aat <- meta_chi2$alphaantitrypsinlevel
    df$anemie <- meta_chi2$anemie2
    
    df$totalreads<-meta_chi2$read_count
    
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted",  "age", "sexe", "Country", "totalreads", "aat", "calpro", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value+  age  + totalreads + calpro , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    taxa_chi2<-tax_table(fecesprevalence_g2)
    row.names(logresults)<-logresults$variable
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    Rank7names<-data.frame(taxa_chi2[, 7])
    Rank7names$variable<-row.names(Rank7names)
    logresults = merge(Rank7names, logresults)
    View(logresults)
    
    write.csv(logresults,"LogresultsFecespresabs.stuntedGenuswithinflaDAL.csv")
    
    
    
    #### rank7 level: logistic model for country correcting for covariables#####
    df_chi2 <- as.matrix((otu_table(fecesprevalence_g2))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_g2)) #take metadata
    data1 <- df_chi2
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$run <- meta_chi2$run
    df$age <- meta_chi2$age
    df$sexe <- meta_chi2$sexe
    df$Country <- meta_chi2$pays
    df$stunted <- meta_chi2$stunted
    df$calpro <- meta_chi2$calprotectinelevel
    df$aat <- meta_chi2$alphaantitrypsinlevel
    df$totalreads<-meta_chi2$read_count
   
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "Country", "totalreads",  "aat", "calpro")) ## use here the variables that showed to be associated in dispersion tes
    long=filter(long, Country!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(Country ~ value+ age  + stunted + totalreads + calpro , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    taxa_chi2<-tax_table(fecesprevalence_g2)
    row.names(logresults)<-logresults$variable
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    Rank7names<-data.frame(taxa_chi2[, 7])
    Rank7names$variable<-row.names(Rank7names)
    logresults = merge(Rank7names, logresults)
    View(logresults)
    
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    write.csv(logresults,"LogresultsFecespresabs.paysGenuswithinflaDAL.csv")
    
  
    #### rank7 level: logistic model for stunting correcting for covariables only in Bangui NOTHING ASSOCIATED #####
    fecesprevalence_g2b<-subset_samples(fecesprevalence_g2, pays=="RCA")
    df_chi2 <- as.matrix((otu_table(fecesprevalence_g2b))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_g2b)) #take metadata
    data1 <- df_chi2
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_chi2$age
    df$sexe <- meta_chi2$sexe
    df$stunted <- meta_chi2$stunted
    df$calpro <- meta_chi2$calprotectinelevel
    df$aat <- meta_chi2$alphaantitrypsinlevel
    df$anemie <- meta_chi2$anemie2
    
    df$totalreads<-meta_chi2$read_count
    
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "totalreads", "aat", "calpro", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    long=filter(long, calpro!="") ## keep only the ones with valid data 
    long=filter(long, aat!="") ## keep only the ones with valid data 
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value+ age  + totalreads + calpro + anemie , .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    taxa_chi2<-tax_table(fecesprevalence_g2b)
    row.names(logresults)<-logresults$variable
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    Rank7names<-data.frame(taxa_chi2[, 7])
    Rank7names$variable<-row.names(Rank7names)
    logresults = merge(Rank7names, logresults)
    View(logresults)
    
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    write.csv(logresults,"LogresultsFecespresabs.stuntedGenuswithinflaBanguiDAL.csv")
    
    
    #### rank7 level: logistic model for stunting correcting for covariables only in Tana NOTHING ASSOCIATED #####
    fecesprevalence_g2a<-subset_samples(fecesprevalence_g2, pays=="Madagascar")
    df_chi2 <- as.matrix((otu_table(fecesprevalence_g2a))) # take presence absence table
    meta_chi2 <- data.frame(sample_data(fecesprevalence_g2a)) #take metadata
    data1 <- df_chi2
    
    #add other categorical factors, etc. 
    df <- data.frame(data1)
    dim(df)
    df$age <- meta_chi2$age
    df$run <- meta_chi2$run
    df$sexe <- meta_chi2$sexe
    df$stunted <- meta_chi2$stunted
    df$calpro <- meta_chi2$calprotectinelevel
    df$aat <- meta_chi2$alphaantitrypsinlevel
    df$anemie<- meta_chi2$anemie2
    
    df$totalreads<-meta_chi2$read_count
    
    library(broom)
    library(dplyr)
    long = melt(df, id.vars = c("stunted", "age", "sexe", "totalreads", "aat", "calpro", "anemie")) ## use here the variables that showed to be associated in dispersion test
    long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
    
    logresults<- long %>% 
      group_by(variable) %>% 
      do(tidy(glm(stunted ~ value+ age  + totalreads +  calpro + anemie, .,  family=binomial))) %>% 
      filter(term == "value") %>% 
      mutate(Beta = as.character(round(estimate, 3)), "p.value" = round(p.value, 5), SE = round(std.error, 3)) %>% 
      ungroup()%>% 
      dplyr::select(variable, Beta, SE, "p.value") %>% 
      as.data.frame()
    
    taxa_chi2<-tax_table(fecesprevalence_g2a)
    row.names(logresults)<-logresults$variable
    View(logresults)
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    Rank7names<-data.frame(taxa_chi2[, 7])
    Rank7names$variable<-row.names(Rank7names)
    logresults = merge(Rank7names, logresults)
    View(logresults)
    
    logresults <- logresults[order(logresults$p.value), ]
    logresults$rel.fdr <- p.adjust(logresults$p.value, method="fdr")
    
    write.csv(logresults,"LogresultsFecespresabs.stuntedGenuswithinflaTanaDAL.csv")
    
    
    