  ### Analysis of ITS2 sample set Afribiota- Vonaesch et al., October 2022  #########
  
  library(tidyr)
  library(DESeq2)
  library(qualpalr)
  library(ggplot2)
  library(biomformat)
  library(plyr)
  library(dplyr)
  library(lme4)
  library(microbiome)
  library(data.table)
  library(nlme)
  library(SparseM)
  library(reshape2)
  library(readstata13)
  require(devtools)
  library(FactoMineR)
  library(factoextra)
  library(vegan)
  library(phyloseq)
  packageDescription("phyloseq")$Version 
  packageDescription("vegan")$Version
  packageDescription("ggplot2")$Version 
  packageDescription("microbiome")$Version
  packageDescription("DeSeq2")$Version
  version # to get the version from R
  library(metagMisc)
  #install_github("zdk123/SpiecEasi")
  library(SpiecEasi)
  library(GGally)
  
  ##### set the work directory and set working environment ####
  
  setwd("~/Users/admin/Desktop/Papers and manuscripts/under revision/Eukaryome paper/Revised version October 22 including last datasets/Test")
  rm(list=ls(all=TRUE)) # go get a clean workspace 
  
  ##### read in the data #####
  metadata<-read.csv("Table S10.csv", sep=",", header=TRUE)
  metadata<-sample_data(metadata)
  row.names(metadata)<-metadata$samplename
  View(metadata)
  
  taxtable<-read.csv("Table S12.csv", sep=",", header=FALSE)
  View(taxtable)
  colnames(taxtable)=as.character(unlist(taxtable[1, ]))
  taxtable=taxtable[-1, ]
  dim(taxtable)
  row.names(taxtable)<-taxtable[, 1]
  taxtable<-taxtable[, -1]
  
  taxtable=as.matrix(taxtable)
  taxonomy=tax_table(taxtable)
  #View(taxonomy)
  taxa_names(taxonomy)
  
  otutable<-read.csv("Table S11.csv", sep=";", header=TRUE)
  otutable$X<-as.character(otutable$X)
  row.names(otutable)<-otutable$X
  otutable<-otutable[, -1]
  otutable<-as.matrix(otutable)
  class(otutable)
  otutable<-otu_table(otutable, taxa_are_rows = TRUE)
  
  head(sample_names(otutable))
  head(sample_names(metadata))
  
  head(taxa_names(otutable))
  head(taxa_names(taxonomy))
  
  sample_names(otutable)<-gsub("AG.", "AG-", sample_names(otutable))
  sample_names(otutable)<-gsub("AD.", "AD-", sample_names(otutable))
  sample_names(otutable)<-gsub("SE.", "SE-", sample_names(otutable))
  sample_names(otutable)<-gsub("SE_", "SE-", sample_names(otutable))
  sample_names(otutable)<-gsub("ER.", "ER-", sample_names(otutable))
  sample_names(otutable)<-gsub("DNA.", "DNA-", sample_names(otutable))
  sample_names(otutable)<-gsub("S0079.", "S0079-", sample_names(otutable))
  sample_names(otutable)<-gsub("S003K.", "S003K-", sample_names(otutable))
  
  physeq= phyloseq(otutable, taxonomy)
  sample_names(physeq)
  
  dim(metadata)
  physeq
  
  setdiff(sample_names(metadata), sample_names(physeq))
  setdiff(sample_names(physeq), sample_names(metadata))
  
  df<-merge_phyloseq(physeq, metadata)
  df #886
  
  #### Prune singlets  ####      
        prune_singlets= function(x) {x[x <= 1] <- 0
        return(x)}
        
        df2 <- df %>% #this produces a maximum possible prevalence count per species per ASV
          transform_sample_counts(fun = prune_singlets)
        
             
  #### Define sample type  ####
        
        duodenalsamples<-grep("AD-", sample_data(df)$id_afri, value=TRUE)
        gastricsamples<-grep("AG-", sample_data(df)$id_afri, value=TRUE)
        fecalsamples<-grep("SE", sample_data(df)$id_afri, value=TRUE)
       
        sample_data(df)$SampleType<-""
        sample_data(df)$SampleType[sample_data(df)$id_afri %in% duodenalsamples] <- "duodenal"
        sample_data(df)$SampleType[sample_data(df)$id_afri %in% gastricsamples] <- "gastric"
        sample_data(df)$SampleType[sample_data(df)$id_afri %in% fecalsamples] <- "feces"
        levels(as.factor(sample_data(df)$SampleType))
        df<-subset_samples(df, SampleType=="feces" | SampleType=="gastric" | SampleType=="duodenal")
      
        duodenalsamples<-grep("AD-", sample_data(df2)$id_afri, value=TRUE)
        gastricsamples<-grep("AG-", sample_data(df2)$id_afri, value=TRUE)
        fecalsamples<-grep("SE", sample_data(df2)$id_afri, value=TRUE)
       
        sample_data(df2)$SampleType<-""
        sample_data(df2)$SampleType[sample_data(df2)$id_afri %in% duodenalsamples] <- "duodenal"
        sample_data(df2)$SampleType[sample_data(df2)$id_afri %in% gastricsamples] <- "gastric"
        sample_data(df2)$SampleType[sample_data(df2)$id_afri %in% fecalsamples] <- "feces"
        df2<-subset_samples(df2, SampleType=="feces" | SampleType=="gastric" | SampleType=="duodenal")
        
    
        length(which(sample_data(df2)$SampleType=="gastric")) # we start with 218 samples
        length(which(sample_data(df2)$SampleType=="duodenal")) # we start with 231 samples
        length(which(sample_data(df2)$SampleType=="feces")) # we start with 436 samples
      
        
  #### Tell the program how to interprete the names ####
        
        colnames(tax_table(df2))
        
  #### clean metadata and tax-table, filter your samples according to the blanks and analyze what leads to differences in total sample read count: Data here on how many have sequences! #####
        sample_data(df2)$SampleType<-as.factor(sample_data(df2)$SampleType)
        levels(as.factor(sample_data(df2)$SampleType))
        
        
        ##clean up your metadata for full afribiota dataset
        #create a variable age in years
        sample_data(df2)$ageyears<-cut(as.numeric(as.character(sample_data(df2)$age)), c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
        which(is.na(sample_data(df2)$ageyears)) # the controls have no assocaited age year
        levels(sample_data(df2)$ageyears) <- c("2-3 years", "3-4 years", "4-5 years")
        levels(sample_data(df2)$ageyears)
        
        sampledata=as.data.frame(sample_data(df2))
        sum(sample_sums(df2)) ##33458189
        
        mean=mean(sample_sums(df2)) ##37805.86 this is the mean sample sums
        mean
        
        # kick out samples that do not match the criteria
        sample_data(df2)$ph_estomac<-as.numeric(sample_data(df2)$ph_estomac)
        sample_data(df2)$ph_intestin<-as.numeric(sample_data(df2)$ph_intestin)
        
        df3 <- df2 %>%
          subset_samples(raison_hospi=="Recrutementcommunautaire") # to filter out samples which are not recruited in the community
        
        df4a <- df3 %>%
          subset_samples(whz_cont <=2)
                         
        df4b <- df4a %>%
        subset_samples(whz_cont >=-2)  # to keep only samples with no acute undernutrition or obesity
        
        df4f <- df4b %>%
          subset_samples(SampleType=="feces")
        
        df4g <- df4b %>%
          subset_samples(SampleType == "gastric" & ph_estomac<=4)
        
        df4d <- df4b %>%
          subset_samples(SampleType == "duodenal" & ph_intestin>=5)
        
         df5<-merge_phyloseq(df4f, df4g, df4d)
          
        df5 # we are left with 754 samples
        
        
        table(sample_data(df2)$SampleType)
        table(sample_data(df5)$SampleType)
            
        # now filter out low sample sums samples
        plot(sort(sample_data(df5)$read_count), ylim=c(0, 100000))
        abline(h=1000)
        
        plot(sort(sample_data(df5)$read_count), ylim=c(0, 100000))
        abline(h=10000)
        
        plot(sort(sample_data(df5)$read_count), ylim=c(0, 100000))
        abline(h=5000)

        dff = prune_samples((sample_data(df5)$read_count>=5000)==TRUE, df5) # only keep samples with more than 5000 reads
        dim(sample_data(dff)) ## we are left with 767 samples 
        dim(tax_table(dff)) 
        table(sample_data(dff)$SampleType)
        table(sample_data(df5)$SampleType)
        table(sample_data(dff)$SampleType, sample_data(dff)$pays)
        
        # filter out low abundance taxa, as they are likely contaminants. Here, we decided that an ASV needs to have at least 50 sequences in at least 1% of the samples
        dff=filter_taxa(dff, function(x) sum(x > 50) > (0.01*length(x)), TRUE)
        dim(tax_table(dff)) #240 ASV's. Note: we lost a lot of taxa, mainly due to the fact that they are spuriously distributed in between the samples!
        
        #track the evolution of your sample numbers
        table(sample_data(dff)$SampleType)
           
        # now make subsets of samples according to SampleType
        df_feces=subset_samples(dff, sample_data(dff)$SampleType=="feces") 
        df_duodenal=subset_samples(dff, sample_data(dff)$SampleType=="duodenal")
        df_gastric=subset_samples(dff, sample_data(dff)$SampleType=="gastric")
        
        sum_feces=sum(sample_sums(df_feces))
        sum_feces 
        
        median_feces=median(sample_sums(df_feces))
        median_feces 
        
        sum_gastric=sum(sample_sums(df_gastric))
        sum_gastric 
        
        median_gastric=median(sample_sums(df_gastric))
        median_gastric 
        
        sum_duodenal=sum(sample_sums(df_duodenal))
        sum_duodenal 
        
        median_duodenal=median(sample_sums(df_duodenal))
        median_duodenal 
        
  ####  Make a subset of samples which do have all three compartments covered ####  
        test<-sample_data(dff)
        test_gastric<-sample_data(df_gastric)
        test_duodenal<-sample_data(df_duodenal)
        test_feces<-sample_data(df_feces)
        dim(test_feces)
        which(duplicated(test_feces$id)==TRUE) # none is duplicated.

        test_gastric2<-test_gastric[test_gastric$id %in% test_duodenal$id, ]
        test_gastric3<-test_gastric2 [test_gastric2$id %in% test_feces$id, ]
        
        test_gastric3<-test_gastric2[test_gastric2$id %in% test_feces$id, ]
        dim(test_gastric3) # 27 samples remaining!
        
        test_duodenal2<-test_duodenal[test_duodenal$id %in% test_gastric3$id, ]
        dim(test_duodenal2) # 27 samples remaining!
        
        test_feces2<-test_feces[test_feces$id %in% test_gastric3$id, ]
        dim(test_feces2) # 27 samples remaining!
        
        dfclean_rest<-subset_samples(dff, sample_data(dff)$id %in% test_feces2$id)
        dfclean_rest # 81 samples are remaining!
        
        table(sample_data(dfclean_rest)$SampleType) # ok!
        table(sample_data(dfclean_rest)$SampleType, sample_data(dfclean_rest)$pays) # ok! 
        
        table(sample_data(dfclean_rest)$SampleType) # ok! 
        
  #### Look what is influecing the overall sample count for feces: run, Sample Type, age, stunting, country ####
        levels(as.factor(sample_data(dff)$SampleType))
        dff<-subset_samples(dff, SampleType!="")
        levels(as.factor(sample_data(dff)$ageyears))
        dff<-subset_samples(dff, ageyears!="")

        sampledata= as.data.frame(sample_data(dff))
        sampledata$SampleType<-as.factor(sampledata$SampleType)
        levels(sampledata$SampleType)
        sampledata$ageyears<-as.factor(sampledata$ageyears)
      
        boxplot(sampledata$read_count~sampledata$pays, main= "Sample Sums", xlab="pays of origin", ylab= "")
        kruskal.test(sampledata$read_count~sampledata$pays) # ns
        
        boxplot(sampledata$read_count~sampledata$stunted + sampledata$pays, main= "Sample Sums", xlab="Stunting", ylab= "")
        kruskal.test(sampledata$read_count~sampledata$stunted) # p-value =3.66e-05
        
        boxplot(sampledata$read_count~sampledata$run, main= "Sample Sums", xlab="Different runs", ylab= "") ## quite a big effect of the run!!
        kruskal.test(sampledata$read_count~sampledata$run) #p-value < 2.2e-16
        
        boxplot(sampledata$read_count~sampledata$SampleType, main= "Sample Sums", xlab="Sample Type", ylab= "")
        kruskal.test(sampledata$read_count~sampledata$SampleType) # p-value = 7.586e-13
        
        boxplot(sampledata$read_count~sampledata$ageyears, main= "Sample Sums", xlab="Age in years", ylab= "")
        kruskal.test(sampledata$read_count~sampledata$ageyears) # ns
        
        sampledata$row.names=row.names(sampledata)
        sampledata=as.data.frame(sampledata)
        
        # now make a multivariate model to see who is contributing independently 
        sampledata= as.data.frame(sample_data(dff))
        
        samplesums_total <- lm(sampledata$read_count ~ sampledata$run+ sampledata$haz + sampledata$ageyears + sampledata$pays + sampledata$SampleType, na.action=na.omit)
        summary(samplesums_total) # to see the results
  
        # redo the analysis on the reduced dataset
        sampledata_rest<-sampledata[sampledata$id %in% sample_data(dfclean_rest)$id]
        dim(sampledata_rest) # 150 samples are remaining!
        
        table(sampledata_rest$SampleType) 
        
        boxplot(sampledata_rest$read_count~sampledata_rest$SampleType, main= "Sample Sums", xlab="Sample Type", ylab= "")
        kruskal.test(sampledata_rest$read_count~sampledata_rest$SampleType) # p-value = ns
        
        samplesums_total <- lm(read_count ~ run+ haz + ageyears + pays + SampleType, data= sampledata_rest, na.action=na.omit)
        summary(samplesums_total) # to see the results
        
        
        # redo the analysis just for the feces
        sampledata= sample_data(df_feces)
        
        boxplot(sampledata$read_count~sampledata$pays, main= "Sample Sums", xlab="pays of origin", ylab= "")
        kruskal.test(sampledata$read_count~sampledata$pays) # 0.011
        
        boxplot(sampledata$read_count~sampledata$stunted + sampledata$pays, main= "Sample Sums", xlab="Stunting", ylab= "")
        kruskal.test(sampledata$read_count~sampledata$stunted) # ns
        
        boxplot(sampledata$read_count~sampledata$run, main= "Sample Sums", xlab="Different runs", ylab= "") ## quite a big effect of the run!!
        kruskal.test(sampledata$read_count~sampledata$run)  ##1.451e-13
        
        boxplot(sampledata$read_count~sampledata$ageyears, main= "Sample Sums", xlab="Age in years", ylab= "")
        kruskal.test(sampledata$read_count~sampledata$ageyears) # ns
        
        sampledata$row.names=row.names(sampledata)
        sampledata=as.data.frame(sampledata)
        View(sampledata$read_count)
        Fungalreads<-cbind(sampledata$id, sampledata$read_count)
        View(Fungalreads)
        colnames(Fungalreads)<-c("id", "fungalreads")
        write.csv(Fungalreads, "fungalreadsITS.csv")

        
        #now make a multivariate model to see who is contributing independently
        
        samplesums_total <- lm(read_count ~  run+ haz + ageyears + pays, data=sampledata, na.action=na.omit)
        summary(samplesums_total) # to see the results
        
        
  #### describe your dataset ####
  
        df_clean<-dff 
        df_clean## 618 samples, 219 taxa 
        colnames(tax_table(df_clean))
        
        colnames(tax_table(df_clean))<- colnames(tax_table(df2))
        View(tax_table(df_clean))
        
        get_taxa_unique(df_clean, "Family") # 38 + NA
      
        get_taxa_unique(df_clean, "Genus") ## 54 + NA unique taxa at genus level in dataset
        get_taxa_unique(df_clean, "Species") ## 72 + NA unique taxa at species level in dataset plus 1 listed as NA and 1 as uncultured
        get_taxa_unique(df_clean, "Phylum") ## 2+ NA unique taxa at phylum level
        get_taxa_unique(df_clean, "Order") ## 22 + NA unique taxa at order level
        get_taxa_unique(df_clean, "Class") ## 10 + NA unique taxa at class level
  
  #### make subset for feces, gastric and duodenal samples #####                 
        
        dffiltered_feces= subset_samples(df_clean, SampleType=="feces")
        dffiltered_duodenal= subset_samples(df_clean, SampleType=="duodenal")
        dffiltered_gastric= subset_samples(df_clean, SampleType=="gastric")
        
        dffiltered_duodenal = filter_taxa(dffiltered_duodenal, function(x) mean(x) > 0, TRUE)
        get_taxa_unique(dffiltered_duodenal, "Kingdom") # to see how many kingdoms we have, we have 1 different kingdoms
        get_taxa_unique(dffiltered_duodenal, "Phylum") ## 2 + NA unique taxa at phylum level in dataset 
        get_taxa_unique(dffiltered_duodenal, "Genus") ## 42 + NA unique taxa at genus level in dataset, Plus 1 listed as NA and 1 listed as uncultured, 1 as "
        get_taxa_unique(dffiltered_duodenal, "Species") ## 57 + NA unique taxa at species level in dataset plus 1 listed as NA; Beta_vulgaris_subsp._vulgaris is a plant -> look up if wrongly labelled in database
        
        dffiltered_gastric = filter_taxa(dffiltered_gastric, function(x) mean(x) > 0, TRUE)
        get_taxa_unique(dffiltered_gastric, "Kingdom") # to see how many kingdoms we have, we have 1 different kingdoms
        get_taxa_unique(dffiltered_gastric, "Phylum") ## 2+ NA  unique taxa at phylum level in dataset 
        get_taxa_unique(dffiltered_gastric, "Genus") ## 48 + NA unique taxa at genus level in dataset, Plus 1 listed as NA and 1 listed as uncultured, 1 as "
        get_taxa_unique(dffiltered_gastric, "Species") ## 61 + NA unique taxa at species level in dataset plus 1 listed as NA; Beta_vulgaris_subsp._vulgaris is a plant -> look up if wrongly labelled in database
        
        dffiltered_feces = filter_taxa(dffiltered_feces, function(x) mean(x) > 0, TRUE)
        get_taxa_unique(dffiltered_feces, "Kingdom") # to see how many kingdoms we have, we have 1 different kingdoms
        get_taxa_unique(dffiltered_feces, "Phylum") ## 2 + NA unique taxa at phylum level in dataset 
        get_taxa_unique(dffiltered_feces, "Genus") ## 53 + NA unique taxa at genus level in dataset, Plus 1 listed as NA and 1 listed as uncultured, 1 as "
        get_taxa_unique(dffiltered_feces, "Species") ## 72 + NA unique taxa at species level in dataset plus 1 listed as NA; Beta_vulgaris_subsp._vulgaris is a plant -> look up if wrongly labelled in database
        
        dfcleanfeces<-dffiltered_feces
        dfcleangastric<-dffiltered_gastric
        dfcleanduodenal<-dffiltered_duodenal
        
        dfclean<-df_clean
        
               
        
  #### make further subsets and assess for number of samples/group ####      

        dfcleanfecesA<-subset_samples(dfcleanfeces, pays=="Madagascar")
        dfcleanfecesB<-subset_samples(dfcleanfeces, pays=="RCA")
        
        dfcleanduodenalA<-subset_samples(dfcleanduodenal, pays=="Madagascar")
        dfcleanduodenalB<-subset_samples(dfcleanduodenal, pays=="RCA")
        
        dfcleangastricA<-subset_samples(dfcleangastric, pays=="Madagascar")
        dfcleangastricB<-subset_samples(dfcleangastric, pays=="RCA")
        
        
        #how many fecal samples do we have (how many subjects, 1 to 1 correspondence) 
        length(which(sample_data(dfclean)$SampleType == "feces"))  ## we have 315 fecal samples/subjects
        length(which(sample_data(dfcleanfeces)$pays == "RCA"))  ## we have 104 fecal samples/subjects from RCA
        length(which(sample_data(dfcleanfeces)$pays == "Madagascar"))  ## we have 211 fecal samples/subjects from Mada
        
        # how many from non-stunted or stunted children in feces and how are they divided according to pays
        length(which(sample_data(dfcleanfeces)$stunted == "stunted"))  ## we have 154 fecal samples/subjects
        length(which(sample_data(dfcleanfeces)$stunted == "non-stunted"))  ## we have 161 fecal samples/subjects
        
        sample_data(dfcleanfeces)$haz<-as.factor(sample_data(dfcleanfeces)$haz)
        length(which(sample_data(dfcleanfeces)$haz == "normonutris"))  ## we have 161 fecal samples/subjects
        length(which(sample_data(dfcleanfeces)$haz == "malnutris chronique modere"))  ## we have 70 fecal samples/subjects
        length(which(sample_data(dfcleanfeces)$haz == "malnutris chronique severe"))  ## we have 84 fecal samples/subjects
        
        length(which(sample_data(dfcleanfecesB)$stunted == "stunted"))  ## we have 54 fecal samples/subjects
        length(which(sample_data(dfcleanfecesB)$stunted == "non-stunted"))  ## we have 50 fecal samples/subjects
        length(which(sample_data(dfcleanfecesB)$haz == "malnutris chronique modere"))  ## we have 17 fecal samples/subjects
        length(which(sample_data(dfcleanfecesB)$haz == "malnutris chronique severe"))  ## we have 37 fecal samples/subjects
        
        length(which(sample_data(dfcleanfecesA)$stunted == "stunted"))  ## we have 100 fecal samples/subjects
        length(which(sample_data(dfcleanfecesA)$stunted == "non-stunted"))  ## we have 111 fecal samples/subjects
        length(which(sample_data(dfcleanfecesA)$haz == "malnutris chronique modere"))  ## we have 53 fecal samples/subjects
        length(which(sample_data(dfcleanfecesA)$haz == "malnutris chronique severe"))  ## we have 47 fecal samples/subjects
        
        #how many duodenal samples do we have and how do they distribute?
        length(which(sample_data(dfclean)$SampleType == "duodenal"))  ## we have 145 duodenal samples
        length(which(sample_data(dfcleanduodenal)$haz == "malnutris chronique modere"))  ## we have 69 duodenal samples
        length(which(sample_data(dfcleanduodenal)$haz == "malnutris chronique severe"))  ## we have 76 duodenal samples
        
        length(which(sample_data(dfcleanduodenal)$pays == "RCA"))  ## we have 80 duodenal samples
        length(which(sample_data(dfcleanduodenal)$pays == "Madagascar"))  ## we have 65 duodenal samples
        
        length(which(sample_data(dfcleanduodenalB)$haz == "malnutris chronique modere"))  ## we have 40 duodenal samples
        length(which(sample_data(dfcleanduodenalB)$haz == "malnutris chronique severe"))  ## we have 40 duodenal samples
        
        length(which(sample_data(dfcleanduodenalA)$haz == "malnutris chronique modere"))  ## we have 29 duodenal samples
        length(which(sample_data(dfcleanduodenalA)$haz == "malnutris chronique severe"))  ## we have 36 duodenal samples
        
 
        
  #### get information on the overall characteristics of the study group####
        table(sample_data(dffiltered_feces)$pays)
        table(sample_data(dffiltered_feces)$stunted, sample_data(dffiltered_feces)$pays)
        table(sample_data(dffiltered_feces)$haz, sample_data(dffiltered_feces)$pays)
        table(sample_data(dffiltered_feces)$sexe, sample_data(dffiltered_feces)$pays)
        table(sample_data(dffiltered_feces)$ageyears, sample_data(dffiltered_feces)$pays)
        table(sample_data(dffiltered_feces)$calprotectinelevel, sample_data(dffiltered_feces)$pays)
        table(sample_data(dffiltered_feces)$alphaantitrypsinlevel, sample_data(dffiltered_feces)$pays)
        table(sample_data(dffiltered_feces)$anemie2, sample_data(dffiltered_feces)$pays)
        
        
   ##### Get information on sequencing depths ####
  dfclean<-prune_samples(sample_sums(dfclean)>=1, dfclean)
  minimum= min(sample_sums(dfclean))
  minimum # minimum number of sequences is 32
  maximum=max(sample_sums(dfclean))
  maximum # maximum is 179260 sequences
  median=median(sample_sums(dfclean))
  median # median is 24047
  mean=mean(sample_sums(dfclean))
  mean # mean is 32347.03 sequences
  
  sample_sums(dfcleanfeces)
  minimumfeces= min(sample_sums(dfcleanfeces))
  minimumfeces  # minimum number of sequences is 999
  maximumfeces=max(sample_sums(dfcleanfeces))
  maximumfeces # maximum is 166411 sequences
  median=median(sample_sums(dfcleanfeces))
  median # median is 22075 sequences
  mean=mean(sample_sums(dfcleanfeces))
  mean # mean is 33487.64 sequences
  
  sample_sums(dfcleanfecesB)
  minimumfeces= min(sample_sums(dfcleanfecesB))
  minimumfeces  # minimum number of sequences is 2050
  maximumfeces=max(sample_sums(dfcleanfecesB))
  maximumfeces # maximum is 166411 sequences
  median=median(sample_sums(dfcleanfecesB))
  median # maximum is 16813 sequences
  mean=mean(sample_sums(dfcleanfecesB))
  mean # maximum is 29551.37 sequences
  
  sample_sums(dfcleanfecesA)
  minimumfeces= min(sample_sums(dfcleanfecesA))
  minimumfeces  # minimum number of sequences is 999
  maximumfeces=max(sample_sums(dfcleanfecesA))
  maximumfeces # maximum is 161912 sequences
  median=median(sample_sums(dfcleanfecesA))
  median # maximum is 25300 sequences
  mean=mean(sample_sums(dfcleanfecesA))
  mean # maximum is 35427.8 sequences
  
  sample_sums(dfcleanduodenal)
  minimumduodenal= min(sample_sums(dfcleanduodenal))
  minimumduodenal # minimum number of sequences is 0
  maximumduodenal=max(sample_sums(dfcleanduodenal))
  maximumduodenal # maximum is 167970 sequences
  median=median(sample_sums(dfcleanduodenal))
  median # maximum is 28332 sequences
  mean=mean(sample_sums(dfcleanduodenal))
  mean # maximum is 33408.63 sequences
  
  sample_sums(dfcleangastric)
  minimumduodenal= min(sample_sums(dfcleangastric))
  minimumduodenal # minimum number of sequences is 32
  maximumduodenal=max(sample_sums(dfcleangastric))
  maximumduodenal # maximum is 179260 sequences
  median=median(sample_sums(dfcleangastric))
  median # median is 23282 sequences
  mean=mean(sample_sums(dfcleangastric))
  mean # mean is 28894.06 sequences
  
  
      #### Get a visual idea of the sampling depth of each sample before rarefaction ####
    
    pdf("samplesequeningdepthfulldatasetnonfilteredallITSSMicrobiomeInsightsnofiltering.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 15, # define plot width and height. completely up to user.
        height = 8)
    sample_sum_df_clean <- data.frame(sum = sample_sums(df2))
    ggplot(sample_sum_df_clean, aes(x = sum)) + 
      geom_histogram(color = "black", fill = "blue", binwidth = 50) +
      ggtitle("Distribution of sample sequencing depth before rarefaction, all samples") + 
      xlab("Read counts") + theme(axis.text = element_text(size = 16), axis.title = element_text(size=18, face="bold"), plot.title=element_text(size=22, face="bold")) +
      theme(axis.title.y = element_blank())
    dev.off()  
    
  
    pdf("samplesequeningdepthfecesMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    sample_sum_dffilteredfeces <- data.frame(sum = sample_sums(dfcleanfeces))
    ggplot(sample_sum_dffilteredfeces, aes(x = sum)) + 
      geom_histogram(color = "black", fill = "red", binwidth = 50) +
      ggtitle("Distribution of fecal sample sequencing depth before rarefaction") + 
      xlab("Read counts") + theme(axis.text = element_text(size = 16), axis.title = element_text(size=18, face="bold"), plot.title=element_text(size=22, face="bold")) +
      theme(axis.title.y = element_blank())
    dev.off()  
    
    pdf("samplesequeningdepthfecesBMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    sample_sum_dffilteredfecesB <- data.frame(sum = sample_sums(dfcleanfecesB))
    ggplot(sample_sum_dffilteredfecesB, aes(x = sum)) + 
      geom_histogram(color = "red", fill = "red", binwidth = 50) +
      ggtitle("Distribution of fecal sample sequencing depth before rarefaction for Bangui") + 
      xlab("Read counts") + theme(axis.text = element_text(size = 16), axis.title = element_text(size=18, face="bold"), plot.title=element_text(size=22, face="bold")) +
      theme(axis.title.y = element_blank())
    dev.off()  
    
    pdf("samplesequeningdepthfecesAMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 20, # define plot width and height. completely up to user.
        height = 8)
    sample_sum_dfcleanfecesA <- data.frame(sum = sample_sums(dfcleanfecesA))
    ggplot(sample_sum_dfcleanfecesA, aes(x = sum)) + 
      geom_histogram(color = "green", fill = "green", binwidth = 50) +
      ggtitle("Distribution of fecal sample sequencing depth before rarefaction for Antananarivo") + 
      xlab("Read counts") + theme(axis.text = element_text(size = 16), axis.title = element_text(size=18, face="bold"), plot.title=element_text(size=22, face="bold")) +
      theme(axis.title.y = element_blank())
    dev.off()  
    
    
    
    pdf("samplesequeningdepthduodenalMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 18, # define plot width and height. completely up to user.
        height = 8)
    sample_sum_dfcleanduodenal <- data.frame(sum = sample_sums(dfcleanduodenal))
    ggplot(sample_sum_dfcleanduodenal, aes(x = sum)) + 
      geom_histogram(color = "black", fill = "green", binwidth = 10) +
      ggtitle("Distribution of fecal duodenal sequencing depth before rarefaction") + 
      xlab("Read counts") +theme(axis.text = element_text(size = 16), axis.title = element_text(size=18, face="bold"), plot.title=element_text(size=22, face="bold")) +
      theme(axis.title.y = element_blank())
    dev.off()  
    
    pdf("samplesequeningdepthgastricMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 18, # define plot width and height. completely up to user.
        height = 8)
    sample_sum_dfcleangastric <- data.frame(sum = sample_sums(dfcleangastric))
    ggplot(sample_sum_dfcleangastric, aes(x = sum)) + 
      geom_histogram(color = "black", fill = "green", binwidth = 10) +
      ggtitle("Distribution of fecal duodenal sequencing depth before rarefaction") + 
      xlab("Read counts") +theme(axis.text = element_text(size = 16), axis.title = element_text(size=18, face="bold"), plot.title=element_text(size=22, face="bold")) +
      theme(axis.title.y = element_blank())
    dev.off()  
    
  
    
    
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
   
    
      #### rarefy your data to 5000 sequences ####
    
    # Now make a loop and rarefy to even depth in repetitive manner to avoid inducing bias. It is always important to set a seed when you subsample so your result is replicable 
    
    set.seed(3)
    for (i in 1:100) {
      # Subsample
      dfrar5000 <- rarefy_even_depth(dfclean, sample.size = 5000, verbose = FALSE, replace = TRUE)}
    
    length(which(sample_data(dfrar5000)$pays=="RCA")) #227 are from RCA
    length(which(sample_data(dfrar5000)$pays=="Madagascar")) #352 are from Mada
    
    length(which(sample_data(dfrar5000)$SampleType=="gastric")) #148 are gastric
    length(which(sample_data(dfrar5000)$SampleType=="duodenal")) #132 are duodenal
    length(which(sample_data(dfrar5000)$SampleType=="feces")) #299 are feces
    
    
    set.seed(3)
    for (i in 1:100) {
      # Subsample
      dffecesrar5000 <- rarefy_even_depth(dfcleanfeces, sample.size = 5000, verbose = FALSE, replace = TRUE)}
    
    length(which(sample_data(dffecesrar5000)$pays=="RCA")) #100 are from RCA
    length(which(sample_data(dffecesrar5000)$pays=="Madagascar")) #202 are from Mada
    
    set.seed(3)
    for (i in 1:100) {
      # Subsample
      dfduodenalrar5000 <- rarefy_even_depth(dfcleanduodenal, sample.size = 5000, verbose = FALSE, replace = TRUE)}
    
    length(which(sample_data(dfduodenalrar5000)$pays=="RCA")) #71 are from RCA
    length(which(sample_data(dfduodenalrar5000)$pays=="Madagascar")) #61 are from Mada
    
    set.seed(3)
    for (i in 1:100) {
      # Subsample
      dfgastricrar5000 <- rarefy_even_depth(dfcleangastric, sample.size = 5000, verbose = FALSE, replace = TRUE)}
    
    length(which(sample_data(dfgastricrar5000)$pays=="RCA")) #59 are from RCA
    length(which(sample_data(dfgastricrar5000)$pays=="Madagascar")) #89 are from Mada
    
    
    set.seed(3)
    for (i in 1:100) {
      # Subsample
      dfcleanred5000 <- rarefy_even_depth(dfclean_rest, sample.size = 5000, verbose = FALSE, replace = TRUE)}
    
    length(which(sample_data(dfcleanred5000)$SampleType=="gastric")) #28 are from RCA
    length(which(sample_data(dfcleanred5000)$SampleType=="duodenal")) #27 are from Mada
    length(which(sample_data(dfcleanred5000)$SampleType=="feces")) #28 are from Mada
    
    
  #### Alpha diversity measures on untrimmed, rarefied datasets ####
    ## on sample type
    # comment: do this on untrimmed dataset
    sample_data(df)$ph_estomac<-as.numeric(sample_data(df)$ph_estomac)
    sample_data(df)$ph_intestin<-as.numeric(sample_data(df)$ph_intestin)
    sample_data(df)$whz_cont<-as.numeric(sample_data(df)$whz_cont)
 
    df <- df %>%
      subset_samples(raison_hospi=="Recrutement communautaire") # to filter out samples which are not recruited in the community
    
    df <- df %>%
      subset_samples((SampleType == "feces") | ph_intestin>=5 | ph_estomac<=4) # to keep only samples with good pH range (<=4 for gastric, >=5 for duodenal))
    
    df <- df %>%
      subset_samples((whz_cont >2) |(whz_cont <-2) ) # to keep only samples with no acute undernutrition or obesity
    
    df # we are left with 827 samples
    
    
    set.seed(3)
    for (i in 1:100) {
      dfrar5000untrimmed<-rarefy_even_depth(df, sample.size = 5000, verbose = FALSE, replace = TRUE)}
    
    sample_data(dfrar5000untrimmed)$ageyears<-cut(as.numeric(as.character(sample_data(dfrar5000untrimmed)$age)), c(24,36,48,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
    levels(sample_data(dfrar5000untrimmed)$ageyears) <- c("2-3 years", "3-4 years", "4-5 years")
    levels(sample_data(dfrar5000untrimmed)$ageyears)
    
    
    df_SampleType<-subset_samples(dfrar5000untrimmed, (SampleType=="gastric" | SampleType=="duodenal" | SampleType=="feces"))
    sample_data(df_SampleType)$SampleType<- factor(sample_data(df_SampleType)$SampleType,
                                 levels = c("gastric", "duodenal", "feces"))
    
    pdf("alphadiversityASVlevelSampleType.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_SampleType, x= "SampleType", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to SampleType")
    p + geom_boxplot(data = p$data, aes(x = SampleType, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Sample Type")
    dev.off()
    
    # calculate if there is a significant difference in Shannon Index
    results = estimate_richness(df_SampleType, measures = 'Shannon')
    d = sample_data(df_SampleType)
    feces = results[d[,'SampleType'] == 'feces',]
    duodenal = results[d[,'SampleType'] == 'duodenal',]
    gastric = results[d[,'SampleType'] == 'gastric',]
    wilcox.test(gastric, duodenal) #p-value = ns
    wilcox.test(gastric, feces) #p-value = 1.205e-08
    wilcox.test(duodenal, feces) #p-value = 1.282e-09
    
    
    df_feces<-subset_samples(df_SampleType, SampleType=="feces")
    df_feces<-subset_samples(df_feces, pays!="")
    
    # on country of origin
    pdf("alphadiversityASVlevelpaysfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces, x= "pays", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to Country of origin")
    p + geom_boxplot(data = p$data, aes(x = pays, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Country of origin")
    dev.off()
    
    # calculate if there is a significant difference in Shannon Index
    results = estimate_richness(df_feces, measures = 'Shannon')
    d = sample_data(df_feces)
    CAR = results[d[,'pays'] == 'RCA',]
    Mada = results[d[,'pays'] == 'Madagascar',]
    wilcox.test(Mada, CAR) # p-value = ns
    
    # on age in years
    pdf("alphadiversityASVlevelageyearsfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces, x= "ageyears", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to age in years")
    p + geom_boxplot(data = p$data, aes(x = ageyears, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Age in years")
    dev.off()
    
    # calculate if there is a significant difference in Shannon Index
    results = estimate_richness(df_feces, measures = 'Shannon')
    d = sample_data(df_feces)
    lowest = results[d[,'ageyears'] == '2-3 years',]
    middle = results[d[,'ageyears'] == '3-4 years',]
    highest= results[d[,'ageyears'] == '4-5 years',]
    wilcox.test(lowest, highest) # p-value = ns
    wilcox.test(lowest, middle) # p-value = ns
    wilcox.test(middle, highest) # p-value = ns
    
    
    df_duodenal<-subset_samples(df_SampleType, SampleType=="duodenal")
    df_duodenal<-subset_samples(df_duodenal, pays!="")
    
    pdf("alphadiversityASVlevelpaysduodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_duodenal, x= "pays", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to Country of origin")
    p + geom_boxplot(data = p$data, aes(x = pays, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Country of origin")
    dev.off()
    
    # calculate if there is a significant difference in Shannon Index
    results = estimate_richness(df_duodenal, measures = 'Shannon')
    d = sample_data(df_duodenal)
    CAR = results[d[,'pays'] == 'RCA',]
    Mada = results[d[,'pays'] == 'Madagascar',]
    wilcox.test(Mada, CAR) # p-value = 0.0004154
    
    
    df_gastric<-subset_samples(df_SampleType, SampleType=="gastric")
    df_gastric<-subset_samples(df_gastric, pays!="")
    
    pdf("alphadiversityASVlevelpaysgastric.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_gastric, x= "pays", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to Country of origin")
    p + geom_boxplot(data = p$data, aes(x = pays, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Country of origin")
    dev.off()
    
    # calculate if there is a significant difference in Shannon Index
    results = estimate_richness(df_gastric, measures = 'Shannon')
    d = sample_data(df_gastric)
    CAR = results[d[,'pays'] == 'RCA',]
    Mada = results[d[,'pays'] == 'Madagascar',]
    wilcox.test(Mada, CAR) # p-value = 0.01312
    
    
    # on stunting status
    pdf("alphadiversityASVlevelhazfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces, x= "haz", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to height-for-age z-score")
    p + geom_boxplot(data = p$data, aes(x = haz, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Height-for-age z-score")
    dev.off()
    
    # calculate if there is a significant difference in Shannon Index
    results = estimate_richness(df_gastric, measures = 'Shannon')
    d = sample_data(df_gastric)
    MCM = results[d[,'haz'] == 'malnutris chronique modere',]
    MCS = results[d[,'haz'] == 'malnutris chronique severe',]
    wilcox.test(MCM, MCS) # p-value = ns
    
    df_feces_B<-subset_samples(df_feces, pays=="RCA")
    pdf("alphadiversityASVlevelhazfecesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces_B, x= "haz", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to height-for-age z-score")
    p + geom_boxplot(data = p$data, aes(x = haz, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Height-for-age z-score")
    dev.off()
    
    # calculate if there is a significant difference in Shannon Index
    results = estimate_richness(df_feces_B, measures = 'Shannon')
    d = sample_data(df_feces_B)
    NN = results[d[,'haz'] == 'normonutris',]
    MCM = results[d[,'haz'] == 'malnutris chronique modere',]
    MCS = results[d[,'haz'] == 'malnutris chronique severe',]
    wilcox.test(NN, MCM) # p-value = 0.04314
    wilcox.test(NN, MCS) # p-value = ns
    wilcox.test(MCM, MCS) # p-value = ns
    
    df_feces_A<-subset_samples(df_feces, pays=="Madagascar")
    pdf("alphadiversityASVlevelhazfecesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces_A, x= "haz", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to height-for-age z-score")
    p + geom_boxplot(data = p$data, aes(x = haz, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Height-for-age z-score")
    dev.off()
    
    # calculate if there is a significant difference in Shannon Index
    results = estimate_richness(df_feces_A, measures = 'Shannon')
    d = sample_data(df_feces_A)
    NN = results[d[,'haz'] == 'normonutris',]
    MCM = results[d[,'haz'] == 'malnutris chronique modere',]
    MCS = results[d[,'haz'] == 'malnutris chronique severe',]
    wilcox.test(NN, MCM) # p-value = ns
    wilcox.test(NN, MCS) # p-value = ns
    wilcox.test(MCM, MCS) # p-value = ns
    

    pdf("alphadiversityASVlevelhazduodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_duodenal, x= "haz", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to height-for-age z-score")
    p + geom_boxplot(data = p$data, aes(x = haz, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Height-for-age z-score")
    dev.off()
    
    # calculate if there is a significant difference in Shannon Index
    results = estimate_richness(df_duodenal, measures = 'Shannon')
    d = sample_data(df_duodenal)
    MCM = results[d[,'haz'] == 'malnutris chronique modere',]
    MCS = results[d[,'haz'] == 'malnutris chronique severe',]
    wilcox.test(MCM, MCS) # p-value = ns
    
    
       pdf("alphadiversityASVlevelhazgastric.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_gastric, x= "haz", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to height-for-age z-score")
    p + geom_boxplot(data = p$data, aes(x = haz, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Height-for-age z-score")
    dev.off()
    
    # calculate if there is a significant difference in Shannon Index
    results = estimate_richness(df_gastric, measures = 'Shannon')
    d = sample_data(df_gastric)
    MCM = results[d[,'haz'] == 'malnutris chronique modere',]
    MCS = results[d[,'haz'] == 'malnutris chronique severe',]
    wilcox.test(MCM, MCS) # p-value = ns
    
    
    df_feces_AAT<-subset_samples(df_feces, alphaantitrypsinlevel!="")
    df_feces_AAT_A<-subset_samples(df_feces_AAT, pays=="Madagascar")
    df_feces_AAT_B<-subset_samples(df_feces_AAT, pays=="RCA")
   
     # on inflammation
    pdf("alphadiversityASVlevelAATfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces_AAT, x= "alphaantitrypsinlevel", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to AAT level")
    p + geom_boxplot(data = p$data, aes(x = alphaantitrypsinlevel, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Alphaantitrypsin Level")
    dev.off()
    
       
    pdf("alphadiversityASVlevelAATfecesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces_AAT_B, x= "alphaantitrypsinlevel", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to AAT level")
    p + geom_boxplot(data = p$data, aes(x = alphaantitrypsinlevel, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Alphaantitrypsin Level")
    dev.off()
    
    pdf("alphadiversityASVlevelAATfecesTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces_AAT_A, x= "alphaantitrypsinlevel", measures=c("Chao1", "Shannon"), title = "Alpha Diversity according to AAT level")
    p + geom_boxplot(data = p$data, aes(x = alphaantitrypsinlevel, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Alphaantitrypsin Level")
    dev.off()
    

    df_feces_calpro<-subset_samples(df_feces, calprotectinelevel!="")
    df_feces_calpro_A<-subset_samples(df_feces_calpro, pays=="Madagascar")
    df_feces_calpro_B<-subset_samples(df_feces_calpro, pays=="RCA")
    
    
    pdf("alphadiversityASVlevelCalprofeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces_calpro, x= "calprotectinelevel", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to calprotectine level")
    p + geom_boxplot(data = p$data, aes(x = calprotectinelevel, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Calprotectine Level")
    dev.off()
    
    pdf("alphadiversityASVlevelCalprofecesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces_calpro_A, x= "calprotectinelevel", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to SampleType")
    p + geom_boxplot(data = p$data, aes(x = calprotectinelevel, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Calprotectine Level")
    dev.off()
    
    pdf("alphadiversityASVlevelCalprofecesTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces_calpro_B, x= "calprotectinelevel", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to SampleType")
    p + geom_boxplot(data = p$data, aes(x = calprotectinelevel, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Calprotectine Level")
    dev.off()
    
    df_feces_anemie<-subset_samples(df_feces, anemie2!="")
    df_feces_anemie_A<-subset_samples(df_feces_anemie, pays=="Madagascar")
    df_feces_anemie_B<-subset_samples(df_feces_anemie, pays=="RCA")
    
    # on anemia
    pdf("alphadiversityASVlevelanemiefeces.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces_anemie, x= "anemie2", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to anemia status")
    p + geom_boxplot(data = p$data, aes(x = anemie2, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Anemia yes/no")
    dev.off()
    
    
    pdf("alphadiversityASVlevelanemiefecesBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces_anemie_B, x= "anemie2", measures=c("Chao1", "Shannon", "Observed", "InvSimpson"), title = "Alpha Diversity according to anemia status")
    p + geom_boxplot(data = p$data, aes(x = anemie2, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Anemia yes/no")
    dev.off()
    
    pdf("alphadiversityASVlevelanemiefecesTana.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 17, # define plot width and height. completely up to user.
        height = 8)
    p = plot_richness(df_feces_anemie_A, x= "anemie2", measures=c("Chao1", "Shannon"), title = "Alpha Diversity according to anemia status")
    p + geom_boxplot(data = p$data, aes(x = anemie2, y = value, color = NULL), alpha = 0.1) +   
      theme(axis.text=element_text(size=18), axis.title = element_text(size=18, face="bold"), plot.title =element_text(size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=18, face="bold" )) +
      xlab("Anemia yes/no")
    dev.off()
    

 #### Ordination plots ####
  
  ####  Make a PERMANOVA to look for factors influencing  dispersion of whole sample set ####
  dffiltered2<-subset_samples(dfrar5000, calprotectinelevel!="")
  dffiltered2<-subset_samples(dffiltered2, anemie2!="")
  dffiltered2<-subset_samples(dffiltered2, haz!="")
  dffiltered2<-subset_samples(dffiltered2, ageyears!="")
  dffiltered2<-subset_samples(dffiltered2, pays!="")
  dffiltered2<-subset_samples(dffiltered2, ageyears!="")
  
  dffiltered2<-subset_samples(dffiltered2, (rowSums(otu_table(dffiltered2))!=0))
  dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
  dffiltered2
  
  transformlog10= function(x) {log10(x)
  return(x)}
  
  prevalence = function(x){ #this only returns prevalence counts per phylum
    x[x >= 1] <- 1
    return(x)
  }
  
  
  dffiltered2transf<-transform_sample_counts(dffiltered2, fun=transformlog10)
  
  project_bray <- phyloseq::distance(dffiltered2transf, method = "bray")
  sample_df <- data.frame(sample_data(dffiltered2transf))
  table(sample_df$SampleType)
  
  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_bray ~  run + haz + SampleType+ anemie2 + calprotectinelevel + ageyears + read_count + pays, data=sample_df, method="bray")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "Disperionstestallsamples.csv")
  
  results$variable<-row.names(results)
  results$Contribution<-(results$R2*100)
  results_filt<-filter(results, results$Pr..F.<0.05)
  
  # make  graph plot
  pdf("Contributionfulldataset.pdf", #name of file to print. can also include relative or absolute path before filename.
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
  
  
  ## now on presence absence
  
  dffiltered2transf<-dffiltered2 %>% #this produces prevalence "counts" for each species, but not percentages
    transform_sample_counts(fun = prevalence) 
  
  project_jacc <- phyloseq::distance(dffiltered2transf, method = "jaccard")
  sample_df <- data.frame(sample_data(dffiltered2transf))
  table(sample_df$SampleType)
  
  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_jacc ~  run + haz + SampleType+ anemie2 + calprotectinelevel + ageyears + read_count + pays + sexe, data=sample_df, method="jaccard")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "Disperionstestallfecalsamplespresabs.csv")
  
  results$variable<-row.names(results)
  results$Contribution<-(results$R2*100)
  results_filt<-filter(results, results$Pr..F.<0.05)
  
  # make  graph plot
  pdf("Contributionfulldatasetpresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
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
  
  
  ####  Make a PERMANOVA to look for factors influencing  dispersion of reduced sample set ####
  
  test<-sample_data(dfrar5000)
  test_gastric<-subset(test, SampleType=="gastric")
  test_duodenal<-subset(test, SampleType=="duodenal")
  test_feces<-subset(test, SampleType=="feces")
  
  test_gastric2<-test_gastric[test_gastric$id %in% test_duodenal$id, ]
  test_gastric3<-test_gastric2 [test_gastric2$id %in% test_feces$id, ]
  
  test_gastric3<-test_gastric2[test_gastric2$id %in% test_feces$id, ]
  dim(test_gastric3) # 25 samples remaining!
  
  test_duodenal2<-test_duodenal[test_duodenal$id %in% test_gastric3$id, ]
  dim(test_duodenal2) # 25 samples remaining!
  
  test_feces2<-test_feces[test_feces$id %in% test_gastric3$id, ]
  dim(test_feces2) # 25 samples remaining!
  
  dfrar5000_rest<-subset_samples(dfrar5000, sample_data(dfrar5000)$id %in% test_feces2$id)
  dfrar5000_rest # 75 samples are remaining!
  
  dffiltered2<-subset_samples(dfrar5000_rest, haz!="")
  dffiltered2<-subset_samples(dffiltered2, ageyears!="")
  dffiltered2<-subset_samples(dffiltered2, pays!="")
  
  dffiltered2<-subset_samples(dffiltered2, (rowSums(otu_table(dffiltered2))!=0))
  dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
  dffiltered2transf<-transform(dffiltered2, "log10") # we are left with 94 samples
  
  project_bray <- phyloseq::distance(dffiltered2transf, method = "bray")
  sample_df <- data.frame(sample_data(dffiltered2transf))
  table(sample_df$SampleType)
  
  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_bray ~  run + haz + SampleType+ + ageyears + read_count + pays, data=sample_df, method="bray")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "Disperionstestallsamplesreducedset.csv")
  
  results$variable<-row.names(results)
  results$Contribution<-(results$R2*100)
  results_filt<-filter(results, results$Pr..F.<0.05)
  
  # make  graph plot
  pdf("Contributionfulldataset_reducedforsharedcompartments.pdf", #name of file to print. can also include relative or absolute path before filename.
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
  
  ####  Make a PERMANOVA to look for factors influencing  dispersion of fecal samples ####
  dffiltered2<-subset_samples(dffecesrar5000, calprotectinelevel!="")
  dffiltered2<-subset_samples(dffiltered2, anemie2!="")
  dffiltered2<-subset_samples(dffiltered2, haz!="")
  dffiltered2<-subset_samples(dffiltered2, ageyears!="")
  dffiltered2<-subset_samples(dffiltered2, pays!="")
  dffiltered2<-subset_samples(dffiltered2, ageyears!="")
  
  dffiltered2<-subset_samples(dffiltered2, (rowSums(otu_table(dffiltered2))!=0))
  dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
  dffiltered2transf<-transform(dffiltered2, "log10")
  
  project_bray <- phyloseq::distance(dffiltered2transf, method = "bray")
  sample_df <- data.frame(sample_data(dffiltered2transf))
  table(sample_df$pays)
 
  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_bray ~  run + haz + anemie2 + calprotectinelevel + ageyears + read_count + pays, data=sample_df, method="bray")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "Disperionstestallfecalsamples.csv")
  
  results$variable<-row.names(results)
  results$Contribution<-(results$R2*100)
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
  
  # now in Bangui dataset only
  dfrar5000B<-subset_samples(dffiltered2transf, pays=="RCA")
  dffiltered2<-subset_samples(dfrar5000B, (rowSums(otu_table(dffiltered2))!=0))
  dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
  
  project_bray <- phyloseq::distance(dffiltered2, method = "bray")
  sample_df <- data.frame(sample_data(dffiltered2))
  
  
  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_bray ~  run + haz  + calprotectinelevel + ageyears + read_count, data=sample_df, method="bray")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "DisperionstestallfecalsamplesBangui.csv")
  
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
  dfrar5000A<-subset_samples(dffiltered2transf, pays=="Madagascar")
  dffiltered2<-subset_samples(dfrar5000A, (rowSums(otu_table(dffiltered2))!=0))
  dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
  
  
  project_bray <- phyloseq::distance(dffiltered2, method = "bray")
  sample_df <- data.frame(sample_data(dffiltered2))

  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_bray ~ run + haz + anemie2 + calprotectinelevel + read_count + ageyears, data=sample_df, method="bray")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "DisperionstestallfecalsamplesTana.csv")
  
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
  
  
  
  ####  Make a PERMANOVA to look for factors influencing  dispersion of duodenal samples ####
  dffiltered2<-subset_samples(dfduodenalrar5000, calprotectinelevel!="")
  dffiltered2<-subset_samples(dffiltered2, anemie2!="")
  dffiltered2<-subset_samples(dffiltered2, haz!="")
  dffiltered2<-subset_samples(dffiltered2, ageyears!="")
  dffiltered2<-subset_samples(dffiltered2, pays!="")
  dffiltered2<-subset_samples(dffiltered2, ageyears!="")
  
  dffiltered2<-subset_samples(dffiltered2, (rowSums(otu_table(dffiltered2))!=0))
  dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
  dffiltered2transf<-transform(dffiltered2, "log10")
  
  project_bray <- phyloseq::distance(dffiltered2transf, method = "bray")
  sample_df <- data.frame(sample_data(dffiltered2transf))
  table(sample_df$pays)
  
  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_bray ~  run + haz + anemie2 + calprotectinelevel + ageyears + read_count + pays, data=sample_df, method="bray")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "Disperionstestallduodenalsamples.csv")
  
  results$variable<-row.names(results)
  results$Contribution<-(results$R2*100)
  results_filt<-filter(results, results$Pr..F.<0.05)
  
  # make  graph plot
  pdf("Contributionfulldatasetduodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
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
  
  ## now on presence/absence
  
  dfduodenalrar5000presabs<-dfduodenalrar5000 %>% #this produces prevalence "counts" for each species, but not percentages
    transform_sample_counts(fun = prevalence) 
  
  dffiltered2<-subset_samples(dfduodenalrar5000presabs, calprotectinelevel!="")
  dffiltered2<-subset_samples(dffiltered2, anemie2!="")
  dffiltered2<-subset_samples(dffiltered2, haz!="")
  dffiltered2<-subset_samples(dffiltered2, ageyears!="")
  dffiltered2<-subset_samples(dffiltered2, pays!="")
  dffiltered2<-subset_samples(dffiltered2, ageyears!="")
  
  dffiltered2<-subset_samples(dffiltered2, (rowSums(otu_table(dffiltered2))!=0))
  dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
  
  project_jacc <- phyloseq::distance(dffiltered2transf, method = "jaccard")
  sample_df <- data.frame(sample_data(dffiltered2transf))
  table(sample_df$pays)
  
  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_jacc ~  run + haz + anemie2 + calprotectinelevel + ageyears + read_count + pays, data=sample_df, method="jaccard")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "Disperionstestallduodenalsamplespresabs.csv")
  
  results$variable<-row.names(results)
  results$Contribution<-(results$R2*100)
  results_filt<-filter(results, results$Pr..F.<0.05)
  
  # make  graph plot
  pdf("Contributionfulldatasetduodenalpresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
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
  dfrar5000B<-subset_samples(dffiltered2transf, pays=="RCA")
  dffiltered2<-subset_samples(dfrar5000B, (rowSums(otu_table(dffiltered2))!=0))
  dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
  
  project_bray <- phyloseq::distance(dffiltered2, method = "bray")
  sample_df <- data.frame(sample_data(dffiltered2))
  
  
  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_bray ~  run + haz  + calprotectinelevel + ageyears + anemie2 + read_count, data=sample_df, method="bray")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "DisperionstestallduodenalsamplesBangui.csv")
  
  results$variable<-row.names(results)
  results$Contribution<-(results$R2*100)
  results_filt<-filter(results, results$Pr..F.<0.05)
  
  # make  graph plot
  pdf("ContributionfulldatasetduodenalBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
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
  dfrar5000A<-subset_samples(dffiltered2transf, pays=="Madagascar")
  dffiltered2<-subset_samples(dfrar5000A, (rowSums(otu_table(dffiltered2))!=0))
  dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
  
  
  project_bray <- phyloseq::distance(dffiltered2, method = "bray")
  sample_df <- data.frame(sample_data(dffiltered2))
  
  
  
  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_bray ~ run + haz + anemie2 + calprotectinelevel + read_count + ageyears, data=sample_df, method="bray")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "DisperionstestallduodenalsamplesTana.csv")
  
  results$variable<-row.names(results)
  results$Contribution<-(results$R2*100)
  results_filt<-filter(results, results$Pr..F.<0.05)
  
  # make  graph plot
  pdf("ContributionfulldatasetduodenalTana.pdf", #name of file to print. can also include relative or absolute path before filename.
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
  
  
  
  ####  Make a PERMANOVA to look for factors influencing  dispersion of fecal samples- presence/absence ####
  dffiltered2<-subset_samples(dffecesrar5000, calprotectinelevel!="")
  dffiltered2<-subset_samples(dffiltered2, anemie2!="")
  dffiltered2<-subset_samples(dffiltered2, haz!="")
  dffiltered2<-subset_samples(dffiltered2, ageyears!="")
  dffiltered2<-subset_samples(dffiltered2, pays!="")
  dffiltered2<-subset_samples(dffiltered2, ageyears!="")
  
  dffiltered2<-subset_samples(dffiltered2, (rowSums(otu_table(dffiltered2))!=0))
  dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
  
  dffiltered2transf<-dffiltered2 %>% #this produces prevalence "counts" for each species, but not percentages
    transform_sample_counts(fun = prevalence) 
  
  project_jacc <- phyloseq::distance(dffiltered2transf, method = "jaccard")
  sample_df <- data.frame(sample_data(dffiltered2transf))
  table(sample_df$pays)
  
  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_jacc ~  run + haz + anemie2 + calprotectinelevel + ageyears + read_count + pays, data=sample_df, method="jaccard")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "Disperionstestallfecalsamples.csv")
  
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
  dfrar5000B<-subset_samples(dffiltered2transf, pays=="RCA")
  dffiltered2<-subset_samples(dfrar5000B, (rowSums(otu_table(dffiltered2))!=0))
  dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
  
  dffiltered2transf<-dffiltered2 %>% #this produces prevalence "counts" for each species, but not percentages
    transform_sample_counts(fun = prevalence) 
  
  
  project_jacc <- phyloseq::distance(dffiltered2transf, method = "jaccard")
  sample_df <- data.frame(sample_data(dffiltered2transf))
  
  
  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_jacc ~  run + haz  + calprotectinelevel + ageyears + read_count, data=sample_df, method="jaccard")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "DisperionstestallfecalsamplesBanguipresabs.csv")
  
  results$variable<-row.names(results)
  results$Contribution<-(results$R2*100)
  results_filt<-filter(results, results$Pr..F.<0.05)
  
  # make  graph plot
  pdf("ContributionfulldatasetfecesBanguipresabs.pdf", #name of file to print. can also include relative or absolute path before filename.
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
  dfrar5000A<-subset_samples(dffiltered2transf, pays=="Madagascar")
  dffiltered2<-subset_samples(dfrar5000A, (rowSums(otu_table(dffiltered2))!=0))
  dffiltered2<-subset_samples(dffiltered2, (colSums(otu_table(dffiltered2))!=0))
  
  dffiltered2transf<-dffiltered2 %>% #this produces prevalence "counts" for each species, but not percentages
    transform_sample_counts(fun = prevalence) 
  
  
  project_jacc <- phyloseq::distance(dffiltered2transf, method = "jaccard")
  sample_df <- data.frame(sample_data(dffiltered2transf))
  
  
  
  #now the adonis test to see if there is a signficant difference according to different variables -> agemonths, country, nutstatus: agemonths and country indepdentely associated with diversity
  res.adonis <- adonis(project_jacc ~ run + haz + anemie2 + calprotectinelevel + read_count  + ageyears, data=sample_df, method="jaccard")
  res.adonis 
  results<-data.frame(as.data.frame(res.adonis$aov.tab))
  View(results)
  write.csv(results, "DisperionstestallfecalsamplesTanapresabs.csv")
  
  results$variable<-row.names(results)
  results$Contribution<-(results$R2*100)
  results_filt<-filter(results, results$Pr..F.<0.05)
  
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
  
  
  
  
##### Plot Ordination with significant taxa as explanatory arrows on Genus level in all sample types ####

  dfrar2genus<-tax_glom(dfrar5000, "Genus")
  OTU1 = as(otu_table(dfrar2genus), "matrix")
  ft.df.c<-as.data.frame(OTU1)
  
  TAX1 = as(tax_table(dfrar2genus), "matrix")
  ft.df.taxa<-as.data.frame(TAX1)
  
  META1= as(sample_data(dfrar2genus), "matrix")
  ft.df.meta<-as.data.frame(META1)
  
  nrow(ft.df.taxa)
  nrow(ft.df.c)
  
  test<-cbind(ft.df.taxa, ft.df.c)
  colnames(test)
  row.names(test)
  
  ##### Make genus-level OTU table with Sample IDs as row names 
  test2<-test[, -(1:5)]
  test2<-test2[, -(2:4)]
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
  
  test4<-test3[rowSums(test3[, 1:54]!=0), ]
  genus.dat <- colnames(test4[, 1:54])
  
  which(rowSums(test3[, 1:54])==0) #AG-CPB334_S139  AG-CPB358_S71  AD-CPB343_S82 
  
  rownames(test3) 
  
  rownames(ft.df.meta)
  
  ft.df.meta2<-ft.df.meta[-468, ]
  ft.df.meta2<-ft.df.meta2[-333, ]
  ft.df.meta2<-ft.df.meta2[-328, ]
  
  dim(test3)
  dim(test4)
  dim(ft.df.meta2)
  
  # First step is to calculate a distance matrix. 
  # Here we use Bray-Curtis distance metric
  test4_log<-transform(test4[, 1:54], 'log10')
  dist <- vegdist(test4_log[, 1:54],  method = "bray")
  
  # PCoA is not included in vegan. 
  # We will use the ape package instead
  library(ape)
  PCOA <- pcoa(dist)
  
  # plot the eigenvalues and interpret
  barplot(PCOA$values$Relative_eig[1:10])
  biplot.pcoa(PCOA)
  biplot.pcoa(PCOA, test4[, 1:54], scale.=F, center=T)
  
  vectorPCOA<-cbind(PCOA$vectors[, 1:2])
  
  #### Include genus level Abundances in Ordination 
  env.fit <- envfit(vectorPCOA, test4[, 1:54], perm = 999) 
  
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
  df_envfit <- df_envfit[df_envfit$genus %in% Cred$genus, ]
  
  #### Get PCOA scores for axis 1 and 2
  nrow(ft.df.meta2)
  nrow(vectorPCOA) # we kicked-out 0 samples in the process. need to filter metadata file
  
  ft.scores <- as.data.frame(vectorPCOA)
  ft.scores$sampleid <- rownames(test4)
  
  # Paste metadata to PCoA scores for plotting
  ft.scores$pays <- ft.df.meta2$pays
  ft.scores$SampleType <- ft.df.meta2$SampleType
  
  # get distinct rows
  ft.scores.d <- distinct(ft.scores)
  
  #### Make ordination plot with genera as explanatory variables
  ft.shapes <- c(0,1,2,5,16)
  
  pdf(file="~/Desktop/GenusSampleTypeITS_explantory_arrowsrelabun.pdf",
      width = 8,
      height = 5)
  p<-ggplot() 
  p+geom_point(aes(ft.scores.d$Axis.1,ft.scores.d$Axis.2, colour=ft.scores.d$SampleType)) + 
    scale_color_manual(values = c("red","green", "blue")) +
    scale_shape_manual(values = ft.shapes) +
    geom_segment(aes(x = 0, y = 0, xend = df_envfit$Axis.1*0.0005, yend = df_envfit$Axis.2*0.0005), arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)+
    geom_text(aes(df_envfit$Axis.1*0.0005, df_envfit$Axis.2*0.0005, label = df_envfit$genus),color="#808080",alpha=0.5, vjust=2.5, hjust=.5) +
    theme_classic()
  dev.off()
  
  
  
  ##### Plot Ordination with significant taxa as explanatory arrows on Genus level in feces ####
  dfrar2genus<-tax_glom(dffecesrar5000, "Genus")
  OTU2 = as(otu_table(dfrar2genus), "matrix")
  ft.df.c<-as.data.frame(OTU2)
  
  TAX2 = as(tax_table(dfrar2genus), "matrix")
  ft.df.taxa<-as.data.frame(TAX2)
  
  META2= as(sample_data(dfrar2genus), "matrix")
  ft.df.meta<-as.data.frame(META2)
  
  nrow(ft.df.taxa)
  nrow(ft.df.c)
  
  test<-cbind(ft.df.taxa, ft.df.c)
  colnames(test)
  row.names(test)
  
  ##### Make genus-level OTU table with Sample IDs as row names 
  test2<-test[, -(1:5)]
  test2<-test2[, -(2:4)]
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
  
  test4<-filter(test3, (rowSums(test3[, 1:53])!=0))
  genus.dat <- colnames(test4[, 1:53])
  
  dim(test3)
  dim(test4)
  dim(ft.df.meta)
  
  
  #### use PCoA instead
  # First step is to calculate a distance matrix. 
  # Here we use Bray-Curtis distance metric
  test4_log<-transform(test4[, 1:53], 'log10')
  dist <- vegdist(test4_log[, 1:53],  method = "bray")
  
  # PCoA is not included in vegan. 
  # We will use the ape package instead
  library(ape)
  PCOA <- pcoa(dist)
  
  # plot the eigenvalues and interpret
  barplot(PCOA$values$Relative_eig[1:10])
  biplot.pcoa(PCOA)
  biplot.pcoa(PCOA, test4[, 1:53], scale.=F, center=T)
  
  vectorPCOA<-cbind(PCOA$vectors[, 1:2])
  
  #### Include genus level Abundances in Ordination 
  env.fit <- envfit(vectorPCOA, test4[, 1:53], perm = 999) 
  
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
  
  pdf(file="~/Desktop/GenuspaysITS_explantory_arrowsrelabun.pdf",
      width = 8,
      height = 5)
  p<-ggplot() 
  p+geom_point(aes(ft.scores.d$Axis.1,ft.scores.d$Axis.2, colour=ft.scores.d$pays)) + 
    scale_color_manual(values = c("red","blue")) +
    scale_shape_manual(values = ft.shapes) +
    geom_segment(aes(x = 0, y = 0, xend = df_envfit$Axis.1*0.0001, yend = df_envfit$Axis.2*0.0001), arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)+
    geom_text(aes(df_envfit$Axis.1*0.0001, df_envfit$Axis.2*0.0001, label = df_envfit$genus),color="#808080",alpha=0.5, vjust=2.5, hjust=.5) +
    theme_classic()
  dev.off()
  
  
  
  ##### Plot Ordination with significant taxa as explanatory arrows on Genus level pres abs ####
  dfrar2genus<-tax_glom(dffecesrar5000, "Genus")
  dfrar2genuspresabs<-transform_sample_counts(dfrar2genus, fun = prevalence)
  
  dfrar2genus<-tax_glom(dfrar2genuspresabs, "Genus")
  OTU2 = as(otu_table(dfrar2genus), "matrix")
  ft.df.c<-as.data.frame(OTU2)
  
  TAX2 = as(tax_table(dfrar2genus), "matrix")
  ft.df.taxa<-as.data.frame(TAX2)
  
  META2= as(sample_data(dfrar2genus), "matrix")
  ft.df.meta<-as.data.frame(META2)
  
  nrow(ft.df.taxa)
  nrow(ft.df.c)
  
  test<-cbind(ft.df.taxa, ft.df.c)
  colnames(test)
  row.names(test)
  
  ##### Make genus-level OTU table with Sample IDs as row names 
  test2<-test[, -(1:5)]
  test2<-test2[, -(2:4)]
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
  
  test4<-filter(test3, (rowSums(test3[, 1:53])!=0))
  genus.dat <- colnames(test4[, 1:53])
  
  dim(test3)
  dim(test4)
  dim(ft.df.meta)
  
  #### use PCoA 
  # First step is to calculate a distance matrix. 
  # Here we use Jaccard distance metric
  
  dist <- vegdist(test4[, 1:53],  method = "jaccard")
  
  # PCoA is not included in vegan. 
  # We will use the ape package instead
  library(ape)
  PCOA <- pcoa(dist)
  
  # plot the eigenvalues and interpret
  barplot(PCOA$values$Relative_eig[1:10])
  biplot.pcoa(PCOA)
  biplot.pcoa(PCOA, test4[, 1:53], scale.=F, center=T)
  
  vectorPCOA<-cbind(PCOA$vectors[, 1:2])
  
  #### Include genus level Abundances in Ordination 
  env.fit <- envfit(vectorPCOA, test4[, 1:53], perm = 999) 
  
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
  nrow(vectorPCOA) # we kicked-out several samples in the process. need to filter metadata file
  
  ft.scores <- as.data.frame(vectorPCOA)
  ft.scores$sampleid <- rownames(test4)
  
  # Paste metadata to PCoA scores for plotting
  ft.scores$pays <- ft.df.meta$pays
  ft.scores$stunted <- ft.df.meta$stunted
  
  # get distinct rows
  ft.scores.d <- distinct(ft.scores)
  
  
  #### Make ordination plot with Families as explanatory variables
  ft.shapes <- c(0,1,2,5,16)
  
  pdf(file="~/Desktop/GenuspaysMI_explantory_arrowspresabs.pdf",
      width = 8,
      height = 5)
  p<-ggplot() 
  p+geom_point(aes(ft.scores.d$Axis.1,ft.scores.d$Axis.2, colour=ft.scores.d$pays)) + 
    scale_color_manual(values = c("red","blue")) +
    scale_shape_manual(values = ft.shapes) +
    geom_segment(aes(x = 0, y = 0, xend = df_envfit$Axis.1*0.35, yend = df_envfit$Axis.2*0.35), arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)+
    geom_text(aes(df_envfit$Axis.1*0.35, df_envfit$Axis.2*0.35, label = df_envfit$genus),color="#808080",alpha=0.5, vjust=2.5, hjust=.5) +
    theme_classic()
  dev.off()
  
  
  
  #####   Plot Ordination with significant taxa as explanatory arrows on Species level ####
  dfrar2species<-tax_glom(dffecesrar5000, "Species")
  
  OTU2 = as(otu_table(dfrar2species), "matrix")
  ft.df.c<-as.data.frame(OTU2)
  
  TAX2 = as(tax_table(dfrar2species), "matrix")
  ft.df.taxa<-as.data.frame(TAX2)
  
  META2= as(sample_data(dfrar2species), "matrix")
  ft.df.meta<-as.data.frame(META2)
  
  nrow(ft.df.taxa)
  nrow(ft.df.c)
  
  test<-cbind(ft.df.taxa, ft.df.c)
  colnames(test)
  row.names(test)
  
  ##### Make species-level OTU table with Sample IDs as row names (need to tweek this to get the genus level in there as well)
  test$taxonomy<-paste0(test$Genus, test$Species)
  
  colnames(test)
  
  test2<-test[, -(1:9)]
  colnames(test2)
  
  test2[, 300]<-as.character(test2[, 300])
  
  test2[, 300]<-make.unique(test2[, 300])
  
  test3<-t(test2)
  colnames(test3)<-test3[300, ]
  test3<-test3[-300, ]
  class(test3) <- "numeric"
  test3<-as.data.frame(test3)
  test3[is.na(test3)] <- 0 #convert NAs to zeros
  which(is.na(test3))==TRUE
  
  test3$sampleid<-row.names(test3)
  dim(test3)
  
  test4<-test3[rowSums(test3[, 1:72]!=0), ] # we loose here three samples, take them out of the metadata
  
  which(rowSums(test3[, 1:72])==0) #S0079-0449_S221 , SE-HJRA028_S31 SE-HMET028_S87 
  
  rownames(test3) # we got rid of lines 113 and 193
  
  rownames(ft.df.meta)
  
  ft.df.meta2<-ft.df.meta[-193, ]
  ft.df.meta2<-ft.df.meta2[-113, ]
  
  species.dat <- colnames(test4[, 1:72])
  
  dim(test3)
  dim(test4)
  dim(ft.df.meta2)
  
  
  # First step is to calculate a distance matrix. 
  # Here we use Bray-Curtis distance metric
  test4_log<-transform(test4[, 1:72], 'log10')
  dist <- vegdist(test4_log[, 1:72],  method = "bray")
  
  # PCoA is not included in vegan. 
  # We will use the ape package instead
  library(ape)
  PCOA <- pcoa(dist)
  
  # plot the eigenvalues and interpret
  barplot(PCOA$values$Relative_eig[1:10])
  biplot.pcoa(PCOA)
  biplot.pcoa(PCOA, test4[, 1:72], scale.=F, center=T)
  
  vectorPCOA<-cbind(PCOA$vectors[, 1:2])
  
  #### Include species level Abundances in Ordination 
  env.fit <- envfit(vectorPCOA, test4[, 1:72], perm = 999) 
  
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
  df_envfit <- df_envfit[df_envfit$species %in% Cred$species, ]
  
  #### Get PCOA scores for axis 1 and 2
  nrow(ft.df.meta2)
  nrow(vectorPCOA) # we kicked-out 0 samples in the process. need to filter metadata file
  
  ft.scores <- as.data.frame(vectorPCOA)
  ft.scores$sampleid <- rownames(test4)
  
  # Paste metadata to PCoA scores for plotting
  
  ft.scores$pays <- ft.df.meta2$pays
  
  # get distinct rows
  ft.scores.d <- distinct(ft.scores)

  #### Make ordination plot with Families as explanatory variables
  ft.shapes <- c(0,1,2,5,16)
  
  pdf(file="~/Desktop/speciespaysMI_explantory_arrowsrelabun.pdf",
      width = 8,
      height = 5)
  p<-ggplot() 
  p+geom_point(aes(ft.scores.d$Axis.1,ft.scores.d$Axis.2, colour=ft.scores.d$pays)) + 
    scale_color_manual(values = c("red","blue")) +
    scale_shape_manual(values = ft.shapes) +
    geom_segment(aes(x = 0, y = 0, xend = df_envfit$Axis.1*0.0001, yend = df_envfit$Axis.2*0.0001), arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)+
    geom_text(aes(df_envfit$Axis.1*0.0001, df_envfit$Axis.2*0.0001, label = df_envfit$species),color="#808080",alpha=0.5, vjust=2.5, hjust=.5) +
    theme_classic()
  dev.off()
  
  
  ##### Plot Ordination with significant taxa as explanatory arrows on species level pres abs ####
  dfrar2species<-tax_glom(dffecesrar5000, "Species")
  dfrar2speciespresabs<-transform_sample_counts(dfrar2species, fun = prevalence)
  
  OTU1 = as(otu_table(dfrar2speciespresabs), "matrix")
  ft.df.c<-as.data.frame(OTU1)
  
  TAX1 = as(tax_table(dfrar2speciespresabs), "matrix")
  ft.df.taxa<-as.data.frame(TAX1)
  
  META1= as(sample_data(dfrar2speciespresabs), "matrix")
  ft.df.meta<-as.data.frame(META1)
  
  nrow(ft.df.taxa)
  nrow(ft.df.c)
  
  test<-cbind(ft.df.taxa, ft.df.c)
  colnames(test)
  row.names(test)
  
  ##### Make species-level OTU table with Sample IDs as row names 
  test$taxonomy<-paste0(test$Genus, test$Species)
  
  colnames(test)
  
  test2<-test[, -(1:9)]
  colnames(test2)
  
  test2[, 300]<-as.character(test2[, 300])
  
  test2[, 300]<-make.unique(test2[, 300])
  
  test3<-t(test2)
  colnames(test3)<-test3[300, ]
  test3<-test3[-300, ]
  
  class(test3) <- "numeric"
  test3<-as.data.frame(test3)
  test3[is.na(test3)] <- 0 #convert NAs to zeros
  which(is.na(test3))==TRUE
  
  test3$sampleid<-row.names(test3)
  dim(test3)
  
  test4<-test3[rowSums(test3[, 1:72])!=0, ]
  
  which(rowSums(test3[, 1:72])==0) #S0079-0449_S221 , SE-HJRA028_S31 
  
  rownames(test3) 
  
  rownames(ft.df.meta)
  
  ft.df.meta2<-ft.df.meta[-193, ]
  ft.df.meta2<-ft.df.meta2[-113, ]
  
  species.dat <- colnames(test4[, 1:72])
  
  dim(test3)
  dim(test4)
  dim(ft.df.meta)
  
   # First step is to calculate a distance matrix. 
  # Here we use Bray-Curtis distance metric
  test4_log<-transform(test4[, 1:72], 'log10')
  dist <- vegdist(test4_log[, 1:72],  method = "bray")
  
  # PCoA is not included in vegan. 
  # We will use the ape package instead
  library(ape)
  PCOA <- pcoa(dist)
  
  # plot the eigenvalues and interpret
  barplot(PCOA$values$Relative_eig[1:10])
  biplot.pcoa(PCOA)
  biplot.pcoa(PCOA, test4[, 1:72], scale.=F, center=T)
  
  vectorPCOA<-cbind(PCOA$vectors[, 1:2])
  
  #### Include species level Abundances in Ordination 
  env.fit <- envfit(vectorPCOA, test4[, 1:72], perm = 999) 
  
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
  df_envfit <- df_envfit[species %in% Cred$species, ]
  
  #### Get PCOA scores for axis 1 and 2
  nrow(ft.df.meta2)
  nrow(vectorPCOA) # we kicked-out 0 samples in the process. need to filter metadata file
  
  ft.scores <- as.data.frame(vectorPCOA)
  ft.scores$sampleid <- rownames(test4)
  
  # Paste metadata to PCoA scores for plotting
  ft.scores$pays <- ft.df.meta2$pays
  ft.scores$stunted <- ft.df.meta2$stunted
  
  # get distinct rows
  ft.scores.d <- distinct(ft.scores)
  
  
  #### Make ordination plot with Families as explanatory variables
  ft.shapes <- c(0,1,2,5,16)
  
  pdf(file="~/Desktop/speciespaysITS_explantory_arrowspresabs.pdf",
      width = 8,
      height = 5)
  p<-ggplot() 
  p+geom_point(aes(ft.scores.d$Axis.1,ft.scores.d$Axis.2, colour=ft.scores.d$pays)) + 
    scale_color_manual(values = c("red","blue")) +
    scale_shape_manual(values = ft.shapes) +
    geom_segment(aes(x = 0, y = 0, xend = df_envfit$Axis.1*0.5, yend = df_envfit$Axis.2*0.5), arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)+
    geom_text(aes(df_envfit$Axis.1*0.5, df_envfit$Axis.2*0.5, label = df_envfit$species),color="#808080",alpha=0.5, vjust=2.5, hjust=.5) +
    theme_classic()
  dev.off()
  
#### calculate prevalence for Genus for filtered dataset in  for SampleType only for samples from stunted children ####
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
  dffiltered= tax_glom(dfcleanstunted, "Genus")
  
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
  
  tax_table.SampleType =  as.data.frame(tax_table(dffiltered))
  
  merge.prev_counts.dffiltered_s_SampleType= merge(tax_table.SampleType, test.prev, by="row.names")
  
  write.csv(merge.prev_counts.dffiltered_s_SampleType, "merge.prev_counts.dffilteredGenus_SampleType.csv")
  
    #### export abundance tables for Genus and SampleType of origin and fuse them with prevalence tables to get a single one  ####
  Afribiota__abundance_Genus <- dffiltered %>% #this produces prevalence "counts" for each phylum, but not percentages
    merge_samples("SampleType") %>%
    t() %>%
    otu_table() %>%
    transform_sample_counts(function(x) x / sum(x) * 100) %>%
    as.data.frame()
  
  
  colnames(Afribiota__abundance_Genus) <- paste("abundance", colnames(Afribiota__abundance_Genus), sep = ".") #add something to distinguish between relative abundance and prevalence
  
  row.names(merge.prev_counts.dffiltered_s_SampleType) <- merge.prev_counts.dffiltered_s_SampleType$Row.names
  
  merge.abunanceprevGenus.SampleType= merge(merge.prev_counts.dffiltered_s_SampleType, Afribiota__abundance_Genus, by="row.names")
  
  test=merge.abunanceprevGenus.SampleType[,1]==merge.abunanceprevGenus.SampleType[,2] # to see if the two colums are actually similar
  which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
  
  merge.abunanceprevGenus.SampleType <- merge.abunanceprevGenus.SampleType[,-1:-2]
  
  write.csv(merge.abunanceprevGenus.SampleType, "merge.abunanceprevGenus.SampleType.csv")
  
  
  
  
    #### now make a plot of the prevalence of genus by SampleType ####
  ## first add genus level data of taxa to metadata ##
  dffiltered_prev<-transform_sample_counts(dffiltered, fun = prevalence) 
  otu_table_g=otu_table(dffiltered_prev) %>%
    unclass() %>%
    as.data.frame()
  
  View(otu_table_g)
  
  otu_table_g$taxonomy= tax_table(dffiltered_prev)[, 6]
  otu_table_g$taxonomy=gsub("g__", "", otu_table_g$taxonomy)
  
  row.names(otu_table_g)=otu_table_g$taxonomy
  otu_table_g=t(otu_table_g)
  
    GenusdataITS = GenusdataITS[-457,] # to take away the line with taxonomy text
  
  GenusdataITS$Row.names=as.character(GenusdataITS$Row.names)
  write.csv(GenusdataITS, "GenusdataITS.csv")
sample_data_g=as.data.frame(sample_data(dffiltered_prev))
  GenusdataITS=merge(otu_table_g, sample_data_g, by="row.names", all=TRUE)
  
  ## now make the actual graph
  library(binom)
  
  percentage= function(x) {per <- (sum(as.numeric(as.character(x)), na.rm=T)/length(x))*100 
  return (per)}
  
  
  GenusdataITSgastric <-as.data.frame(GenusdataITS[which(GenusdataITS$SampleType=="gastric"), ]) # to generate a list of all the ones that are from gastric
  GenusdataITSduodenal <-as.data.frame(GenusdataITS[which(GenusdataITS$SampleType=="duodenal"), ]) # to generate a list of all the ones that are from duodenal
  GenusdataITSfeces <-as.data.frame(GenusdataITS[which(GenusdataITS$SampleType=="feces"), ]) # to generate a list of all the ones that are from feces
  
  nrow(GenusdataITSgastric) # 158
  nrow(GenusdataITSduodenal) # 144
  nrow(GenusdataITSfeces) # 154
  colnames(GenusdataITSfeces)
  
  GenusdataITSgastric2<-GenusdataITSgastric[, 2:55] 
  GenusdataITSgastric2<- as.matrix(sapply(GenusdataITSgastric2, as.numeric))
  measured<-colSums(GenusdataITSgastric2)
  length(measured)
  samples<-numeric( length = 54)
  samples[1:54] <- 158
  
  test<-as.data.frame(rbind(measured, samples))
  test<-t(test)
  
  CIgastric <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
  
  GenusdataITSduodenal2<-GenusdataITSduodenal[, 2:55] 
  GenusdataITSduodenal2<- as.matrix(sapply(GenusdataITSduodenal2, as.numeric))
  measured<-colSums(GenusdataITSduodenal2)
  length(measured)
  samples<-numeric( length = 54)
  samples[1:54] <- 144
  
  test<-as.data.frame(rbind(measured, samples))
  test<-t(test)
  
  CIduodenal <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
  
  GenusdataITSfeces2<-GenusdataITSfeces[, 2:55] 
  GenusdataITSfeces2<- as.matrix(sapply(GenusdataITSfeces2, as.numeric))
  measured<-colSums(GenusdataITSfeces2)
  length(measured)
  samples<-numeric( length = 54)
  samples[1:54] <- 154
  
  test<-as.data.frame(rbind(measured, samples))
  test<-t(test)
  
  CIfeces <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
  
  Lower=rbind(CIgastric$lower*100, CIduodenal$lower*100, CIfeces$lower*100)
  row.names(Lower)=c("gastric", "duodenal", "feces")
  colnames(Lower)<-colnames(GenusdataITSfeces2)
  
  Upper=rbind(CIgastric$upper*100, CIduodenal$upper*100, CIfeces$upper*100)
  row.names(Upper)=c("gastric", "duodenal", "feces")
  colnames(Upper)<-colnames(GenusdataITSfeces2)
  
  GenusdataITSgastriccollapsed=colwise(percentage)(GenusdataITSgastric)
  GenusdataITSduodenalcollapsed=colwise(percentage)(GenusdataITSduodenal)
  GenusdataITSfecescollapsed=colwise(percentage)(GenusdataITSfeces)
  
  PercentageGenus=rbind(GenusdataITSgastriccollapsed, GenusdataITSduodenalcollapsed, GenusdataITSfecescollapsed)
  
  row.names(PercentageGenus)=c("gastric", "duodenal", "feces")
  View(PercentageGenus)
  ncol(PercentageGenus)
  
  PercentageGenus <- PercentageGenus[,-(55:98), drop=FALSE] # kick-out the metadata
  PercentageGenus <- PercentageGenus[,-(1),drop=FALSE] # kick-out the metadata
  
  write.csv(PercentageGenus, "PercentageGenusSampleType.csv")
  
  PercentageGenus2=t(PercentageGenus)
  Upper2=t(Upper)
  Lower2=t(Lower)
  
  PercentageGenus2toplot=melt(PercentageGenus2, value.name="PercentageGenus2", varnames=c("Taxon", "SampleType"))
  Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "SampleType"))
  Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "SampleType"))
  PercentageGenus2toplot<-merge(PercentageGenus2toplot, Uppertoplot, by=c("Taxon", "SampleType"))
  PercentageGenus2toplot<-merge(PercentageGenus2toplot, Lowertoplot, by=c("Taxon", "SampleType"))
  
  PercentageGenus2toplot<-PercentageGenus2toplot[order(-PercentageGenus2),]
  
  pdf("generaaccordingtoSampleType.pdf", #name of file to print. can also include relative or absolute path before filename.
      width = 10, # define plot width and height. completely up to user.
      height = 5)
  ggplot(data=PercentageGenus2toplot, aes(x=reorder(Taxon, -PercentageGenus2), y=PercentageGenus2,  color=SampleType, fill=SampleType)) +
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
    labs(title="Percentage of ITS2 genera in gastric, duodenal and fecal samples")
  dev.off()
  
  
  # now filter out the low abundance taxa
  PercentageGenus2 <- PercentageGenus[, (PercentageGenus[1, ] > 10) | (PercentageGenus[2, ] > 10) | (PercentageGenus[3, ] > 10)] # now kick-out taxa with less than 10% prevalence in any of the two countries
  
  PercentageGenus2=t(PercentageGenus2)
  
  PercentageGenus2toplot=melt(PercentageGenus2, value.name="PercentageGenus2", varnames=c("Taxon", "SampleType"))
  Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "SampleType"))
  Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "SampleType"))
  PercentageGenus2toplot<-merge(PercentageGenus2toplot, Uppertoplot, by=c("Taxon", "SampleType"))
  PercentageGenus2toplot<-merge(PercentageGenus2toplot, Lowertoplot, by=c("Taxon", "SampleType"))
  
  PercentageGenus2toplot<-PercentageGenus2toplot[order(-PercentageGenus2),]
  
  pdf("generaaccordingtoSampleTypefiltered.pdf", #name of file to print. can also include relative or absolute path before filename.
      width = 10, # define plot width and height. completely up to user.
      height = 5)
  ggplot(data=PercentageGenus2toplot, aes(x=reorder(Taxon, -PercentageGenus2), y=PercentageGenus2,  color=SampleType, fill=SampleType)) +
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
    labs(title="Percentage of ITS2 genera in gastric, duodenal and fecal samples")
  dev.off()
  
  
    #### calculate prevalence for Genus for filtered dataset in  for SampleType on reduced dataset ####
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
    dffiltered= tax_glom(dfcleanstunted, "Genus")
    
    test_gastric<-subset_samples(dfcleanstunted, SampleType=="gastric")
    test_duodenal<-subset_samples(dfcleanstunted, SampleType=="duodenal")
    test_feces<-subset_samples(dfcleanstunted, SampleType=="feces")
    
    test_gastric2<-subset_samples(test_gastric, sample_data(test_gastric)$id %in% sample_data(test_duodenal)$id)
    test_gastric3<-subset_samples(test_gastric2, sample_data(test_gastric2)$id %in% sample_data(test_feces)$id)
    
    test_gastric3<-subset_samples(test_gastric2, sample_data(test_gastric2)$id %in% sample_data(test_feces)$id)
    dim(sample_data(test_gastric3)) # 27 samples remaining!
    
    test_duodenal2<-subset_samples(test_duodenal, sample_data(test_duodenal)$id %in% sample_data(test_gastric3)$id)
    dim(sample_data(test_duodenal2)) # 27 samples remaining!
    
    test_feces2<-subset_samples(test_feces, sample_data(test_feces)$id %in% sample_data(test_gastric3)$id)
    dim(sample_data(test_feces2)) # 27 samples remaining!

    dffiltered_red<-subset_samples(dffiltered, sample_data(dffiltered)$id %in% sample_data(test_feces2)$id)
    dffiltered_red # 81 samples are remaining!
    
    prev_counts.dffiltered_s <- dffiltered_red %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("SampleType") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_s) <- paste("prevalence", colnames(prev_counts.dffiltered_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_possible <- dffiltered_red %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("SampleType") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_s) <- paste("prevalence", colnames(prev_counts.dffiltered_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_s/prev_counts.dffiltered_possible)*100
    
    tax_table.SampleType =  as.data.frame(tax_table(dffiltered))
    
    merge.prev_counts.dffiltered_s_SampleType= merge(tax_table.SampleType, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_s_SampleType, "merge.prev_counts.dffilteredGenus_SampleType_reduced.csv")
    
    
    #### export abundance tables for Genus and SampleType of origin and fuse them with prevalence tables to get a single one  ####
    Afribiota__abundance_Genus <- dffiltered_red %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("SampleType") %>%
      t() %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      as.data.frame()
    
    colnames(Afribiota__abundance_Genus) <- paste("abundance", colnames(Afribiota__abundance_Genus), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    row.names(merge.prev_counts.dffiltered_s_SampleType) <- merge.prev_counts.dffiltered_s_SampleType$Row.names
    
    merge.abunanceprevGenus.SampleType= merge(merge.prev_counts.dffiltered_s_SampleType, Afribiota__abundance_Genus, by="row.names")
    
    test=merge.abunanceprevGenus.SampleType[,1]==merge.abunanceprevGenus.SampleType[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprevGenus.SampleType <- merge.abunanceprevGenus.SampleType[,-1:-2]
    
    write.csv(merge.abunanceprevGenus.SampleType, "merge.abunanceprevGenus.SampleType_reduced.csv")

    
    #### now make a plot of the prevalence of genus by SampleType ####
    ## first add genus level data of taxa to metadata ##
    dffiltered_prev<-transform_sample_counts(dffiltered_red, fun = prevalence) 
    otu_table_g=otu_table(dffiltered_prev) %>%
      unclass() %>%
      as.data.frame()
    
    View(otu_table_g)
    
    otu_table_g$taxonomy= tax_table(dffiltered_prev)[, 6]
    otu_table_g$taxonomy=gsub("g__", "", otu_table_g$taxonomy)
    
    row.names(otu_table_g)=otu_table_g$taxonomy
    otu_table_g=t(otu_table_g)
    
    sample_data_g=as.data.frame(sample_data(dffiltered_prev))
    GenusdataITS=merge(otu_table_g, sample_data_g, by="row.names", all=TRUE)
    View(GenusdataITS)
    
    GenusdataITS = GenusdataITS[-82,] # to take away the line with taxonomy text
    
    GenusdataITS$Row.names=as.character(GenusdataITS$Row.names)
    write.csv(GenusdataITS, "GenusdataITS_red.csv")
    
    ## now make the actual graph
    library(binom)
    
    percentage= function(x) {per <- (sum(as.numeric(as.character(x)), na.rm=T)/length(x))*100 
    return (per)}
    
    GenusdataITSgastric <-as.data.frame(GenusdataITS[which(GenusdataITS$SampleType=="gastric"), ]) # to generate a list of all the ones that are from gastric
    GenusdataITSduodenal <-as.data.frame(GenusdataITS[which(GenusdataITS$SampleType=="duodenal"), ]) # to generate a list of all the ones that are from duodenal
    GenusdataITSfeces <-as.data.frame(GenusdataITS[which(GenusdataITS$SampleType=="feces"), ]) # to generate a list of all the ones that are from feces
    
    nrow(GenusdataITSgastric) # 27
    nrow(GenusdataITSduodenal) # 27
    nrow(GenusdataITSfeces) # 27
    colnames(GenusdataITSfeces)
    
    GenusdataITSgastric2<-GenusdataITSgastric[, 2:55] 
    GenusdataITSgastric2<- as.matrix(sapply(GenusdataITSgastric2, as.numeric))
    measured<-colSums(GenusdataITSgastric2)
    length(measured)
    samples<-numeric( length = 54)
    samples[1:54] <- 27
    
    test<-as.data.frame(rbind(measured, samples))
    test<-t(test)
    
    CIgastric <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    GenusdataITSduodenal2<-GenusdataITSduodenal[, 2:55] 
    GenusdataITSduodenal2<- as.matrix(sapply(GenusdataITSduodenal2, as.numeric))
    measured<-colSums(GenusdataITSduodenal2)
    length(measured)
    samples<-numeric( length = 54)
    samples[1:54] <- 27
    
    test<-as.data.frame(rbind(measured, samples))
    test<-t(test)
    
    CIduodenal <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    GenusdataITSfeces2<-GenusdataITSfeces[, 2:55] 
    GenusdataITSfeces2<- as.matrix(sapply(GenusdataITSfeces2, as.numeric))
    measured<-colSums(GenusdataITSfeces2)
    length(measured)
    samples<-numeric( length = 54)
    samples[1:54] <- 27
    
    test<-as.data.frame(rbind(measured, samples))
    test<-t(test)
    
    CIfeces <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Lower=rbind(CIgastric$lower*100, CIduodenal$lower*100, CIfeces$lower*100)
    row.names(Lower)=c("gastric", "duodenal", "feces")
    colnames(Lower)<-colnames(GenusdataITSfeces2)
    
    Upper=rbind(CIgastric$upper*100, CIduodenal$upper*100, CIfeces$upper*100)
    row.names(Upper)=c("gastric", "duodenal", "feces")
    colnames(Upper)<-colnames(GenusdataITSfeces2)
    
    GenusdataITSgastriccollapsed=colwise(percentage)(GenusdataITSgastric)
    GenusdataITSduodenalcollapsed=colwise(percentage)(GenusdataITSduodenal)
    GenusdataITSfecescollapsed=colwise(percentage)(GenusdataITSfeces)
    
    PercentageGenus=rbind(GenusdataITSgastriccollapsed, GenusdataITSduodenalcollapsed, GenusdataITSfecescollapsed)
    
    row.names(PercentageGenus)=c("gastric", "duodenal", "feces")
    View(PercentageGenus)
    ncol(PercentageGenus)
    
    PercentageGenus <- PercentageGenus[,-(55:98), drop=FALSE] # kick-out the metadata
    PercentageGenus <- PercentageGenus[,-(1),drop=FALSE] # kick-out the metadata
    
    write.csv(PercentageGenus, "PercentageGenusSampleType_red.csv")
    
    PercentageGenus2=t(PercentageGenus)
    Upper2=t(Upper)
    Lower2=t(Lower)
    
    PercentageGenus2toplot=melt(PercentageGenus2, value.name="PercentageGenus2", varnames=c("Taxon", "SampleType"))
    Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "SampleType"))
    Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "SampleType"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Uppertoplot, by=c("Taxon", "SampleType"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Lowertoplot, by=c("Taxon", "SampleType"))
    
    PercentageGenus2toplot<-PercentageGenus2toplot[order(-PercentageGenus2),]
    
    pdf("generaaccordingtoSampleType_red.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 5)
    ggplot(data=PercentageGenus2toplot, aes(x=reorder(Taxon, -PercentageGenus2), y=PercentageGenus2,  color=SampleType, fill=SampleType)) +
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
      labs(title="Percentage of ITS2 genera in gastric, duodenal and fecal samples")
    dev.off()
    
    
    # now filter out the low abundance taxa
    PercentageGenus2 <- PercentageGenus[, (PercentageGenus[1, ] > 10) | (PercentageGenus[2, ] > 10) | (PercentageGenus[3, ] > 10)] # now kick-out taxa with less than 10% prevalence in any of the two countries
    
    PercentageGenus2=t(PercentageGenus2)
    
    PercentageGenus2toplot=melt(PercentageGenus2, value.name="PercentageGenus2", varnames=c("Taxon", "SampleType"))
    Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "SampleType"))
    Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "SampleType"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Uppertoplot, by=c("Taxon", "SampleType"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Lowertoplot, by=c("Taxon", "SampleType"))
    
    PercentageGenus2toplot<-PercentageGenus2toplot[order(-PercentageGenus2),]
    
    pdf("generaaccordingtoSampleTypefiltered_red.pdf", #name of file to print. can also include relative or absolute path before filename.
        width = 10, # define plot width and height. completely up to user.
        height = 5)
    ggplot(data=PercentageGenus2toplot, aes(x=reorder(Taxon, -PercentageGenus2), y=PercentageGenus2,  color=SampleType, fill=SampleType)) +
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
      labs(title="Percentage of ITS2 genera in gastric, duodenal and fecal samples")
    dev.off()
    
    
    
    
    
#### calculate prevalence for Genus for filtered dataset in  for stunting for fecal samples ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered= tax_glom(dfcleanfeces, "Genus")
    
    prev_counts.dffiltered_s <- dffiltered %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_s) <- paste("prevalence", colnames(prev_counts.dffiltered_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_possible <- dffiltered %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_s) <- paste("prevalence", colnames(prev_counts.dffiltered_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_s/prev_counts.dffiltered_possible)*100
    
    tax_table.stunted =  as.data.frame(tax_table(dffiltered))
    
    merge.prev_counts.dffiltered_s_SampleType= merge(tax_table.stunted, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_s_SampleType, "merge.prev_counts.dffilteredGenus_SampleType.csv")
    
    #### export abundance tables for Genus and stunted of origin and fuse them with prevalence tables to get a single one  ####
    Afribiota__abundance_Genus <- dffiltered %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      as.data.frame()
    
    
    colnames(Afribiota__abundance_Genus) <- paste("abundance", colnames(Afribiota__abundance_Genus), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    row.names(merge.prev_counts.dffiltered_s_SampleType) <- merge.prev_counts.dffiltered_s_SampleType$Row.names
    
    merge.abunanceprevGenus.stunted= merge(merge.prev_counts.dffiltered_s_SampleType, Afribiota__abundance_Genus, by="row.names")
    
    test=merge.abunanceprevGenus.stunted[,1]==merge.abunanceprevGenus.stunted[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprevGenus.stunted <- merge.abunanceprevGenus.stunted[,-1:-2]
    
    write.csv(merge.abunanceprevGenus.stunted, "merge.abunanceprevGenus.stunted.csv")
    
    
    
    
    #### now make a plot of the prevalence of genus by stunted ####
    ## first add genus level data of taxa to metadata ##
    dffiltered_prev<-transform_sample_counts(dffiltered, fun = prevalence) 
    otu_table_g=otu_table(dffiltered_prev) %>%
      unclass() %>%
      as.data.frame()
    
    View(otu_table_g)
    
    otu_table_g$taxonomy= tax_table(dffiltered_prev)[, 6]
    otu_table_g$taxonomy=gsub("g__", "", otu_table_g$taxonomy)
    
    row.names(otu_table_g)=otu_table_g$taxonomy
    otu_table_g=t(otu_table_g)
    
    otu_table_g = otu_table_g[-316,] # to take away the line with taxonomy text
    
    write.csv(GenusdataITS, "GenusdataITS.csv")
    sample_data_g=as.data.frame(sample_data(dffiltered_prev))
    GenusdataITS=merge(otu_table_g, sample_data_g, by="row.names", all=TRUE)
    View(GenusdataITS)
    
    ## now make the actual graph
    library(binom)
    
    percentage= function(x) {per <- (sum(as.numeric(as.character(x)), na.rm=T)/length(x))*100 
    return (per)}
    
    
    GenusdataITSstunted <-as.data.frame(GenusdataITS[which(GenusdataITS$stunted=="stunted"), ]) # to generate a list of all the ones that are from gastric
    GenusdataITSnonstunted <-as.data.frame(GenusdataITS[which(GenusdataITS$stunted=="non-stunted"), ]) # to generate a list of all the ones that are from duodenal
    
    nrow(GenusdataITSstunted) # 154
    nrow(GenusdataITSnonstunted) # 161
    
    colnames(GenusdataITSstunted)
    
    GenusdataITSstunted2<-GenusdataITSstunted[, 2:54] 
    GenusdataITSstunted2<- as.matrix(sapply(GenusdataITSstunted2, as.numeric))
    measured<-colSums(GenusdataITSstunted2)
    length(measured)
    samples<-numeric( length = 53)
    samples[1:53] <- 154
    
    test<-as.data.frame(rbind(measured, samples))
    test<-t(test)
    
    CIstunted <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    GenusdataITSnonstunted2<-GenusdataITSnonstunted[, 2:54] 
    GenusdataITSnonstunted2<- as.matrix(sapply(GenusdataITSnonstunted2, as.numeric))
    measured<-colSums(GenusdataITSnonstunted2)
    length(measured)
    samples<-numeric( length = 53)
    samples[1:53] <- 154
    
    test<-as.data.frame(rbind(measured, samples))
    test<-t(test)
    
    CInonstunted <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    GenusdataITSfeces2<-GenusdataITSfeces[, 2:54] 
    GenusdataITSfeces2<- as.matrix(sapply(GenusdataITSfeces2, as.numeric))
    measured<-colSums(GenusdataITSfeces2)
    length(measured)
    samples<-numeric( length = 53)
    samples[1:53] <- 161
    
    test<-as.data.frame(rbind(measured, samples))
    test<-t(test)
    
    CInonstunted <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Lower=rbind(CIstunted$lower*100, CInonstunted$lower*100)
    row.names(Lower)=c("stunted", "non-stunted")
    colnames(Lower)<-colnames(GenusdataITSstunted2)
    
    Upper=rbind(CIstunted$upper*100, CInonstunted$upper*100)
    row.names(Upper)=c("stunted", "non-stunted")
    colnames(Upper)<-colnames(GenusdataITSstunted2)
    
    GenusdataITSstuntedcollapsed=colwise(percentage)(GenusdataITSstunted)
    GenusdataITSnonstuntedcollapsed=colwise(percentage)(GenusdataITSnonstunted)
    
    PercentageGenus=rbind(GenusdataITSstuntedcollapsed, GenusdataITSnonstuntedcollapsed)
    
    row.names(PercentageGenus)=c("stunted", "non-stunted")
    View(PercentageGenus)
    colnames(PercentageGenus)
    ncol(PercentageGenus)
    
    PercentageGenus <- PercentageGenus[,-(55:96), drop=FALSE] # kick-out the metadata
    PercentageGenus <- PercentageGenus[,-(1),drop=FALSE] # kick-out the metadata
    
    write.csv(PercentageGenus, "PercentageGenusstunted.csv")
    
    PercentageGenus2=t(PercentageGenus)
    Upper2=t(Upper)
    Lower2=t(Lower)
    
    PercentageGenus2toplot=melt(PercentageGenus2, value.name="PercentageGenus2", varnames=c("Taxon", "stunted"))
    Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "stunted"))
    Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "stunted"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Uppertoplot, by=c("Taxon", "stunted"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Lowertoplot, by=c("Taxon", "stunted"))
    
    PercentageGenus2toplot<-PercentageGenus2toplot[order(-PercentageGenus2),]
    
    pdf("generaaccordingtostunting.pdf", #name of file to print. can also include relative or absolute path before filename.
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
      labs(title="Percentage of ITS2 genera according to stunting status, feces")
    dev.off()
    
    
    # now filter out the low abundance taxa
    PercentageGenus2 <- PercentageGenus[, (PercentageGenus[1, ] > 10) | (PercentageGenus[2, ] > 10) ] # now kick-out taxa with less than 10% prevalence in any of the two countries
    
    PercentageGenus2=t(PercentageGenus2)
    
    PercentageGenus2toplot=melt(PercentageGenus2, value.name="PercentageGenus2", varnames=c("Taxon", "stunted"))
    Uppertoplot=melt(Upper2, value.name="Upper", varnames=c("Taxon", "stunted"))
    Lowertoplot=melt(Lower2, value.name="Lower", varnames=c("Taxon", "stunted"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Uppertoplot, by=c("Taxon", "stunted"))
    PercentageGenus2toplot<-merge(PercentageGenus2toplot, Lowertoplot, by=c("Taxon", "stunted"))
    
    PercentageGenus2toplot<-PercentageGenus2toplot[order(-PercentageGenus2),]
    
    pdf("generaaccordingtoStuntingfiltered.pdf", #name of file to print. can also include relative or absolute path before filename.
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
      labs(title="Percentage of ITS2 genera according to stunting status, fecal samples")
    dev.off()
    
    
#### calculate prevalence for species for filtered dataset in feces for whole dataset merged ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_feces= tax_glom(dfcleanfeces, "Species")
    
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
    tax_table.fecespays<-as.data.frame(TAX1)

    merge.prev_counts.dffiltered_feces_s_pays= merge(tax_table.fecespays, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_feces_s_pays, "merge.prev_counts.dffiltered_feces_s.csv")
    
    #### export abundance tables for species overall and fuse them with prevalence tables to get a single one  ####
    
    Afribiota_feces_abundance_species <- dffiltered_feces %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("collapse") %>%
      t() %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      as.data.frame()
    
    colnames(Afribiota_feces_abundance_species) <- paste("abundance", colnames(Afribiota_feces_abundance_species), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    row.names(merge.prev_counts.dffiltered_feces_s_pays) <- merge.prev_counts.dffiltered_feces_s_pays$Row.names
    
    merge.abunanceprev.feces.species.pays= merge(merge.prev_counts.dffiltered_feces_s_pays, Afribiota_feces_abundance_species, by="row.names")
    
    test=merge.abunanceprev.feces.species.pays[,1]==merge.abunanceprev.feces.species.pays[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.feces.species.pays <- merge.abunanceprev.feces.species.pays[,-1:-2]
    
    write.csv(merge.abunanceprev.feces.species.pays, "merge.abunanceprev.feces.species.csv")
    
    
    #### calculate prevalence for species for filtered dataset in duodenal for whole dataset merged ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_duodenal= tax_glom(dfcleanduodenal, "Species")
    
    sample_data(dffiltered_duodenal)$collapse<-"collapse"
    
    prev_counts.dffiltered_duodenal_s <- dffiltered_duodenal %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("collapse") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_duodenal_s) <- paste("prevalence", colnames(prev_counts.dffiltered_duodenal_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_duodenal_possible <- dffiltered_duodenal %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("collapse") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_duodenal_s) <- paste("prevalence", colnames(prev_counts.dffiltered_duodenal_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_duodenal_s/prev_counts.dffiltered_duodenal_possible)*100
    
    TAX1 = as(tax_table(dfcleanduodenal), "matrix")
    tax_table.duodenalpays<-as.data.frame(TAX1)
    
    
    merge.prev_counts.dffiltered_duodenal_s_pays= merge(tax_table.duodenalpays, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_duodenal_s_pays, "merge.prev_counts.dffiltered_duodenal_s.csv")
    
    #### export abundance tables for species overall and fuse them with prevalence tables to get a single one  ####
    
    Afribiota_duodenal_abundance_species <- dffiltered_duodenal %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("collapse") %>%
      t() %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      as.data.frame()
    
    colnames(Afribiota_duodenal_abundance_species) <- paste("abundance", colnames(Afribiota_duodenal_abundance_species), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    row.names(merge.prev_counts.dffiltered_duodenal_s_pays) <- merge.prev_counts.dffiltered_duodenal_s_pays$Row.names
    
    merge.abunanceprev.duodenal.species.pays= merge(merge.prev_counts.dffiltered_duodenal_s_pays, Afribiota_duodenal_abundance_species, by="row.names")
    
    test=merge.abunanceprev.duodenal.species.pays[,1]==merge.abunanceprev.duodenal.species.pays[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.duodenal.species.pays <- merge.abunanceprev.duodenal.species.pays[,-1:-2]
    
    write.csv(merge.abunanceprev.duodenal.species.pays, "merge.abunanceprev.duodenal.species.csv")
    
    
    
    
    
    
    #### calculate prevalence for species for filtered dataset in gastric for whole dataset merged ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_gastric= tax_glom(dfcleangastric, "Species")
    
    sample_data(dffiltered_gastric)$collapse<-"collapse"
    
    prev_counts.dffiltered_gastric_s <- dffiltered_gastric %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("collapse") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_gastric_s) <- paste("prevalence", colnames(prev_counts.dffiltered_gastric_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_gastric_possible <- dffiltered_gastric %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("collapse") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_gastric_s) <- paste("prevalence", colnames(prev_counts.dffiltered_gastric_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_gastric_s/prev_counts.dffiltered_gastric_possible)*100
    
    TAX1 = as(tax_table(dfcleangastric), "matrix")
    tax_table.gastricpays<-as.data.frame(TAX1)
    
    
    merge.prev_counts.dffiltered_gastric_s_pays= merge(tax_table.gastricpays, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_gastric_s_pays, "merge.prev_counts.dffiltered_gastric_s.csv")
    
    #### export abundance tables for species overall and fuse them with prevalence tables to get a single one  ####
    
    Afribiota_gastric_abundance_species <- dffiltered_gastric %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("collapse") %>%
      t() %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      as.data.frame()
    
    colnames(Afribiota_gastric_abundance_species) <- paste("abundance", colnames(Afribiota_gastric_abundance_species), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    row.names(merge.prev_counts.dffiltered_gastric_s_pays) <- merge.prev_counts.dffiltered_gastric_s_pays$Row.names
    
    merge.abunanceprev.gastric.species.pays= merge(merge.prev_counts.dffiltered_gastric_s_pays, Afribiota_gastric_abundance_species, by="row.names")
    
    test=merge.abunanceprev.gastric.species.pays[,1]==merge.abunanceprev.gastric.species.pays[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.gastric.species.pays <- merge.abunanceprev.gastric.species.pays[,-1:-2]
    
    write.csv(merge.abunanceprev.gastric.species.pays, "merge.abunanceprev.gastric.species.csv")
    
    
    
    
    
    
    
    
    #### calculate prevalence for species for filtered dataset in feces for pays of origin ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_feces= tax_glom(dfcleanfeces, "Species")
    
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
    tax_table.fecespays<-as.data.frame(TAX1)
    merge.prev_counts.dffiltered_feces_s_pays= merge(tax_table.fecespays, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_feces_s_pays, "merge.prev_counts.dffiltered_feces_s_pays.csv")
    
 #### export abundance tables for species and pays of origin and fuse them with prevalence tables to get a single one  ####
    Afribiota_feces_abundance_species <- dffiltered_feces %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      as.data.frame()
    
    
    colnames(Afribiota_feces_abundance_species) <- paste("abundance", colnames(Afribiota_feces_abundance_species), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    row.names(merge.prev_counts.dffiltered_feces_s_pays) <- merge.prev_counts.dffiltered_feces_s_pays$Row.names
    
    merge.abunanceprev.feces.species.pays= merge(merge.prev_counts.dffiltered_feces_s_pays, Afribiota_feces_abundance_species, by="row.names")
    
    test=merge.abunanceprev.feces.species.pays[,1]==merge.abunanceprev.feces.species.pays[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.feces.species.pays <- merge.abunanceprev.feces.species.pays[,-1:-2]
    
    write.csv(merge.abunanceprev.feces.species.pays, "merge.abunanceprev.feces.species.pays.csv")
    
    
    
    #### calculate prevalence for Genus for filtered dataset in feces for pays of origin ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_feces= tax_glom(dfcleanfeces, "Genus")
    
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
    tax_table.fecespays<-as.data.frame(TAX1)
    
    
    merge.prev_counts.dffiltered_feces_s_pays= merge(tax_table.fecespays, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_feces_s_pays, "merge.prev_counts.dffiltered_feces_Genus_pays.csv")
    
    
    
    
    
    #### export abundance tables for Genus and pays of origin and fuse them with prevalence tables to get a single one  ####
    Afribiota_feces_abundance_Genus <- dffiltered_feces %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      as.data.frame()
    
    
    colnames(Afribiota_feces_abundance_Genus) <- paste("abundance", colnames(Afribiota_feces_abundance_Genus), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    row.names(merge.prev_counts.dffiltered_feces_s_pays) <- merge.prev_counts.dffiltered_feces_s_pays$Row.names
    
    merge.abunanceprev.feces.Genus.pays= merge(merge.prev_counts.dffiltered_feces_s_pays, Afribiota_feces_abundance_Genus, by="row.names")
    
    test=merge.abunanceprev.feces.Genus.pays[,1]==merge.abunanceprev.feces.Genus.pays[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.feces.Genus.pays <- merge.abunanceprev.feces.Genus.pays[,-1:-2]
    
    write.csv(merge.abunanceprev.feces.Genus.pays, "merge.abunanceprev.feces.Genus.pays.csv")
    
    
    
    
    #### now make a plot of the prevalence of genus by country of origin ####
    ## first add genus level data of taxa to metadata ##
    dffiltered_prev<-transform_sample_counts(dffiltered_feces, fun = prevalence) 
    otu_table_g=otu_table(dffiltered_prev) %>%
      unclass() %>%
      as.data.frame()
    
    View(otu_table_g)
    
    otu_table_g$taxonomy= tax_table(dffiltered_prev)[, 6]
    otu_table_g$taxonomy=gsub("g__", "", otu_table_g$taxonomy)
    
    row.names(otu_table_g)=otu_table_g$taxonomy
    otu_table_g=t(otu_table_g)
    
    sample_data_g=as.data.frame(sample_data(dffiltered_prev))
    GenusdataITS=merge(otu_table_g, sample_data_g, by="row.names", all=TRUE)
    dim(GenusdataITS)
    GenusdataITS = GenusdataITS[-296,] # to take away the line with taxonomy text
    
    GenusdataITS$Row.names=as.character(GenusdataITS$Row.names)
    write.csv(GenusdataITS, "GenusdataITSfeces.csv")
    View(GenusdataITS)
    
    ## now make the actual graph
    library(binom)
    
    percentage= function(x) {per <- (sum(as.numeric(as.character(x)), na.rm=T)/length(x))*100 
    return (per)}
    
    
    GenusdataITSCAR <-as.data.frame(GenusdataITS[which(GenusdataITS$pays=="RCA"), ]) # to generate a list of all the ones that are from CAR
    GenusdataITSMada <-as.data.frame(GenusdataITS[which(GenusdataITS$pays=="Madagascar"), ]) # to generate a list of all the ones that are from Madagascar
    
    nrow(GenusdataITSCAR) # 100
    nrow(GenusdataITSMada) # 195
    
    GenusdataITSCAR2<-GenusdataITSCAR[, 2:54] 
    GenusdataITSCAR2<- as.matrix(sapply(GenusdataITSCAR2, as.numeric))
    measured<-colSums(GenusdataITSCAR2)
    length(measured)
    samples<-numeric( length = 53)
    samples[1:53] <- 100
    
    test<-as.data.frame(rbind(measured, samples))
    test<-t(test)
    
    CICAR <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    GenusdataITSMada2<-GenusdataITSMada[, 2:54] 
    GenusdataITSMada2<- as.matrix(sapply(GenusdataITSMada2, as.numeric))
    measured<-colSums(GenusdataITSMada2)
    length(measured)
    samples<-numeric( length = 53)
    samples[1:53] <- 195
    
    test<-as.data.frame(rbind(measured, samples))
    test<-t(test)
    
    CIMada <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Lower=rbind(CICAR$lower*100, CIMada$lower*100)
    row.names(Lower)=c("CAR", "Madagascar")
    colnames(Lower)<-colnames(GenusdataITSCAR2)
    
    Upper=rbind(CICAR$upper*100, CIMada$upper*100)
    row.names(Upper)=c("CAR", "Madagascar")
    colnames(Upper)<-colnames(GenusdataITSMada2)
    
    GenusdataITSCARcollapsed=colwise(percentage)(GenusdataITSCAR)
    GenusdataITSMadacollapsed=colwise(percentage)(GenusdataITSMada)
    
    PercentageGenus=rbind(GenusdataITSCARcollapsed, GenusdataITSMadacollapsed)
    
    row.names(PercentageGenus)=c("CAR", "Madagascar")
    View(PercentageGenus)
    ncol(PercentageGenus)
    colnames(PercentageGenus)
    
    PercentageGenus <- PercentageGenus[,-(55:97), drop=FALSE] # kick-out the metadata
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
      labs(title="Percentage of ITS2 genera in fecal samples according to country of origin")
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
      labs(title="Percentage of ITS2 genera in fecal samples according to country of origin")
    dev.off()  
    
    #### calculate prevalence for Genus for filtered dataset in feces for stunting status ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_feces= tax_glom(dfcleanfeces, "Genus")
    
    prev_counts.dffiltered_feces_s <- dffiltered_feces %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_s) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_feces_possible <- dffiltered_feces %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_feces_s) <- paste("prevalence", colnames(prev_counts.dffiltered_feces_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_feces_s/prev_counts.dffiltered_feces_possible)*100
    
    TAX1 = as(tax_table(dffiltered_feces), "matrix")
    tax_table.fecesstunted<-as.data.frame(TAX1)
    
    
    merge.prev_counts.dffiltered_feces_s_stunted= merge(tax_table.fecesstunted, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_feces_s_stunted, "merge.prev_counts.dffiltered_feces_Genus_stunted.csv")
    
    
    
    
    
    #### export abundance tables for Genus and stunting status and fuse them with prevalence tables to get a single one  ####
    Afribiota_feces_abundance_Genus <- dffiltered_feces %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("stunted") %>%
      t() %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      as.data.frame()
    
    
    colnames(Afribiota_feces_abundance_Genus) <- paste("abundance", colnames(Afribiota_feces_abundance_Genus), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    row.names(merge.prev_counts.dffiltered_feces_s_stunted) <- merge.prev_counts.dffiltered_feces_s_stunted$Row.names
    
    merge.abunanceprev.feces.Genus.stunted= merge(merge.prev_counts.dffiltered_feces_s_stunted, Afribiota_feces_abundance_Genus, by="row.names")
    
    test=merge.abunanceprev.feces.Genus.stunted[,1]==merge.abunanceprev.feces.Genus.stunted[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.feces.Genus.stunted <- merge.abunanceprev.feces.Genus.stunted[,-1:-2]
    
    write.csv(merge.abunanceprev.feces.Genus.stunted, "merge.abunanceprev.feces.Genus.stunted.csv")
    
    
    
    
    #### now make a plot of the prevalence of genus by stunting status ####
    ## first add genus level data of taxa to metadata ##
    dffiltered_prev<-transform_sample_counts(dffiltered_feces, fun = prevalence) 
    otu_table_g=otu_table(dffiltered_prev) %>%
      unclass() %>%
      as.data.frame()
    
    View(otu_table_g)
    
    otu_table_g$taxonomy= tax_table(dffiltered_prev)[, 6]
    otu_table_g$taxonomy=gsub("g__", "", otu_table_g$taxonomy)
    
    row.names(otu_table_g)=otu_table_g$taxonomy
    otu_table_g=t(otu_table_g)
    
    sample_data_g=as.data.frame(sample_data(dffiltered_prev))
    GenusdataITS=merge(otu_table_g, sample_data_g, by="row.names", all=TRUE)
    GenusdataITS = GenusdataITS[-296,] # to take away the line with taxonomy text
    
    GenusdataITS$Row.names=as.character(GenusdataITS$Row.names)
    write.csv(GenusdataITS, "GenusdataITSfeces.csv")
    View(GenusdataITS)
    
    ## now make the actual graph
    library(binom)
    
    percentage= function(x) {per <- (sum(as.numeric(as.character(x)), na.rm=T)/length(x))*100 
    return (per)}
    
    
    GenusdataITSstunted <-as.data.frame(GenusdataITS[which(GenusdataITS$stunted=="stunted"), ]) # to generate a list of all the ones that are from stunted
    GenusdataITSnonstunted <-as.data.frame(GenusdataITS[which(GenusdataITS$stunted=="non-stunted"), ]) # to generate a list of all the ones that are from non-stunted
    
    nrow(GenusdataITSstunted) # 142
    nrow(GenusdataITSnonstunted) # 153
    
    GenusdataITSstunted2<-GenusdataITSstunted[, 2:54] 
    GenusdataITSstunted2<- as.matrix(sapply(GenusdataITSstunted2, as.numeric))
    measured<-colSums(GenusdataITSstunted2)
    length(measured)
    samples<-numeric( length = 53)
    samples[1:53] <- 142
    
    test<-as.data.frame(rbind(measured, samples))
    test<-t(test)
    
    CIstunted <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    GenusdataITSnonstunted2<-GenusdataITSnonstunted[, 2:54] 
    GenusdataITSnonstunted2<- as.matrix(sapply(GenusdataITSnonstunted2, as.numeric))
    measured<-colSums(GenusdataITSnonstunted2)
    length(measured)
    samples<-numeric( length = 53)
    samples[1:53] <- 153
    
    test<-as.data.frame(rbind(measured, samples))
    test<-t(test)
    
    CInonstunted <- binom.confint(x=test[, 1], n=test[, 2], methods="wilson")
    
    Lower=rbind(CIstunted$lower*100, CInonstunted$lower*100)
    row.names(Lower)=c("stunted", "non-stunted")
    colnames(Lower)<-colnames(GenusdataITSstunted2)
    
    Upper=rbind(CIstunted$upper*100, CInonstunted$upper*100)
    row.names(Upper)=c("stunted", "non-stunted")
    colnames(Upper)<-colnames(GenusdataITSnonstunted2)
    
    GenusdataITSstuntedcollapsed=colwise(percentage)(GenusdataITSstunted)
    GenusdataITSnonstuntedcollapsed=colwise(percentage)(GenusdataITSnonstunted)
    
    PercentageGenus=rbind(GenusdataITSstuntedcollapsed, GenusdataITSnonstuntedcollapsed)
    
    row.names(PercentageGenus)=c("stunted", "non-stunted")
    View(PercentageGenus)
    ncol(PercentageGenus)
    colnames(PercentageGenus)
    
    PercentageGenus <- PercentageGenus[,-(55:97), drop=FALSE] # kick-out the metadata
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
    
    pdf("generaaccordingtostuntedfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
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
      labs(title="Percentage of ITS2 genera in fecal samples according to stunting status")
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
    
    pdf("generaaccordingtostuntedfecesfiltered.pdf", #name of file to print. can also include relative or absolute path before filename.
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
      labs(title="Percentage of ITS2 genera in fecal samples according to stunting status")
    dev.off()  
    
    #### calculate prevalence for species for filtered dataset in duodenal for pays of origin ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_duodenal= tax_glom(dfcleanduodenal, "Species")
    
    prev_counts.dffiltered_duodenal_s <- dffiltered_duodenal %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_duodenal_s) <- paste("prevalence", colnames(prev_counts.dffiltered_duodenal_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_duodenal_possible <- dffiltered_duodenal %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_duodenal_s) <- paste("prevalence", colnames(prev_counts.dffiltered_duodenal_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_duodenal_s/prev_counts.dffiltered_duodenal_possible)*100
    
    TAX1 = as(tax_table(dfcleanduodenal), "matrix")
    tax_table.duodenalpays<-as.data.frame(TAX1)
    
    
    merge.prev_counts.dffiltered_duodenal_s_pays= merge(tax_table.duodenalpays, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_duodenal_s_pays, "merge.prev_counts.dffiltered_duodenal_s_pays.csv")
    
    #### export abundance tables for species and pays of origin and fuse them with prevalence tables to get a single one  ####
    Afribiota_duodenal_abundance_species <- dffiltered_duodenal %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      as.data.frame()
    
    
    colnames(Afribiota_duodenal_abundance_species) <- paste("abundance", colnames(Afribiota_duodenal_abundance_species), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    row.names(merge.prev_counts.dffiltered_duodenal_s_pays) <- merge.prev_counts.dffiltered_duodenal_s_pays$Row.names
    
    merge.abunanceprev.duodenal.species.pays= merge(merge.prev_counts.dffiltered_duodenal_s_pays, Afribiota_duodenal_abundance_species, by="row.names")
    
    test=merge.abunanceprev.duodenal.species.pays[,1]==merge.abunanceprev.duodenal.species.pays[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.duodenal.species.pays <- merge.abunanceprev.duodenal.species.pays[,-1:-2]
    
    write.csv(merge.abunanceprev.duodenal.species.pays, "merge.abunanceprev.duodenal.species.pays.csv")
    
    
    
    
    
    #### calculate prevalence for species for filtered dataset in gastric for pays of origin ####
    #prevalence is the percentage of samples that an phylum shows up in (compared to the total number of samples).
    prevalence = function(x){ #this only returns prevalence counts per phylum
      x[x >= 1] <- 1
      return(x)
    }
    
    allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix
      x[x >= 0] <- 1
      return(x)
    }
    
    dffiltered_gastric= tax_glom(dfcleangastric, "Species")
    
    prev_counts.dffiltered_gastric_s <- dffiltered_gastric %>% #this produces prevalence "counts" for each phylum, but not percentages
      transform_sample_counts(fun = prevalence) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_gastric_s) <- paste("prevalence", colnames(prev_counts.dffiltered_gastric_s), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
    
    prev_counts.dffiltered_gastric_possible <- dffiltered_gastric %>% #this produces a maximum possible prevalence count per phylum
      transform_sample_counts(fun = allones) %>%
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      as.data.frame()
    colnames(prev_counts.dffiltered_gastric_s) <- paste("prevalence", colnames(prev_counts.dffiltered_gastric_possible), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-phylum basis
    test.prev = (prev_counts.dffiltered_gastric_s/prev_counts.dffiltered_gastric_possible)*100
    
    TAX1 = as(tax_table(dfcleangastric), "matrix")
    tax_table.gastricpays<-as.data.frame(TAX1)

    merge.prev_counts.dffiltered_gastric_s_pays= merge(tax_table.gastricpays, test.prev, by="row.names")
    
    write.csv(merge.prev_counts.dffiltered_gastric_s_pays, "merge.prev_counts.dffiltered_gastric_s_pays.csv")
    
    #### export abundance tables for species and pays of origin and fuse them with prevalence tables to get a single one  ####
    Afribiota_gastric_abundance_species <- dffiltered_gastric %>% #this produces prevalence "counts" for each phylum, but not percentages
      merge_samples("pays") %>%
      t() %>%
      otu_table() %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      as.data.frame()
    
    
    colnames(Afribiota_gastric_abundance_species) <- paste("abundance", colnames(Afribiota_gastric_abundance_species), sep = ".") #add something to distinguish between relative abundance and prevalence
    
    row.names(merge.prev_counts.dffiltered_gastric_s_pays) <- merge.prev_counts.dffiltered_gastric_s_pays$Row.names
    
    merge.abunanceprev.gastric.species.pays= merge(merge.prev_counts.dffiltered_gastric_s_pays, Afribiota_gastric_abundance_species, by="row.names")
    
    test=merge.abunanceprev.gastric.species.pays[,1]==merge.abunanceprev.gastric.species.pays[,2] # to see if the two colums are actually similar
    which(test==FALSE) # they are actually similar, so we can the two first rows once we actually set the rownames
    
    merge.abunanceprev.gastric.species.pays <- merge.abunanceprev.gastric.species.pays[,-1:-2]
    
    write.csv(merge.abunanceprev.gastric.species.pays, "merge.abunanceprev.gastric.species.pays.csv")
   
  #### Generate core fungom for given taxon level (ASV level) in feces in merged dataset and depending on coutnry and rel. abundance and also using presence/absence and at least prevalence of 90%, filtered dataset, non rarefied ####    
      ## Anything conserved on ASV level?
      eukaryome.rel <- microbiome::transform(dfcleanfeces, "compositional")
      eukaryome.rel
      View(head(otu_table(eukaryome.rel)))
      
      core.taxa <- core(eukaryome.rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      # no taxon is present in 90% of all samples, nor anything above 50
      
      core.taxa <- core(eukaryome.rel, detection = 0.00001, prevalence = 50/100, include.lowest = TRUE) 
      core.taxa # 2 taxa are conserved
      View(tax_table(core.taxa))
      # Humicola grisea, Saccharomyces
      
      core.taxa <- core(eukaryome.rel, detection = 0.00000001, prevalence = 50/100, include.lowest = TRUE) 
      core.taxa ##2 taxa
     
      core.taxa <- core(eukaryome.rel, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE) 
      core.taxa #9 taxa!
      View(tax_table(core.taxa)) # 
      
      core.taxa <- core(eukaryome.rel, detection = 0.00000001, prevalence = 10/100, include.lowest = TRUE)
      core.taxa # 49 taxa 
      View(tax_table(core.taxa))
      
      ## ASV level stratified by country?
      eukaryome.rel.B=subset_samples(eukaryome.rel, pays=="RCA")
      
      core.taxa <- core(eukaryome.rel.B, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      # no taxon is present in 90% of all samples
      
      core.taxa <- core(eukaryome.rel.B, detection = 0.00001, prevalence = 50/100, include.lowest = TRUE) 
      core.taxa # 3 are conserved
      tax_table(core.taxa) # 
      
      core.taxa <- core(eukaryome.rel.B, detection = 0.00000001, prevalence = 50/100, include.lowest = TRUE) 
      core.taxa ##3 taxa  
      View(tax_table(core.taxa))
      
      core.taxa <- core(eukaryome.rel.B, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE)
      core.taxa # 17 taxa
      View(tax_table(core.taxa)) ## several possible gut and a few environmental/diet
      
      core.taxa <- core(eukaryome.rel.B, detection = 0.00000001, prevalence = 10/100, include.lowest = TRUE)
      core.taxa # 51 taxa 
      View(tax_table(core.taxa)) # especially Cladisporium and several others
      
      eukaryome.rel.A=subset_samples(eukaryome.rel, pays=="Madagascar")
      
      core.taxa <- core(eukaryome.rel.A, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      # no taxon is present in 90% of all samples
      
      core.taxa <- core(eukaryome.rel.A, detection = 0.00001, prevalence = 50/100, include.lowest = TRUE) 
      core.taxa
      # 3 are conserved
      tax_table(core.taxa) #Saccharomyces, Cladisporium, Humicola
      
      core.taxa <- core(eukaryome.rel.A, detection = 0.00000001, prevalence = 50/100, include.lowest = TRUE) 
      core.taxa ##3 taxa
      
      core.taxa <- core(eukaryome.rel.A, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE)
      core.taxa # 14 taxa   
      View(tax_table(core.taxa)) # almost all gut
      
      core.taxa <- core(eukaryome.rel.A, detection = 0.00000001, prevalence = 10/100, include.lowest = TRUE)
      core.taxa # 51 taxa 
      View(tax_table(core.taxa)) # now also several enviromental
      
      ## Anything conserved on Species level?
      dfcleanfeces_Species<-tax_glom(dfcleanfeces, "Species")
      eukaryomeSpecies.rel <- microbiome::transform(dfcleanfeces_Species, "compositional")
      core.taxa <- core(eukaryomeSpecies.rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeSpecies.rel, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeSpecies.rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
      core.taxa <- core(eukaryomeSpecies.rel, detection = 0.0001, prevalence = 50/100, include.lowest = TRUE)
      core.taxa # 2 taxa #Humicola grisea, Malassezia restricta
      tax_table(core.taxa)
      
      eukaryomeSpecies.rel.A=subset_samples(eukaryomeSpecies.rel, pays=="Madagascar")
      eukaryomeSpecies.rel.B=subset_samples(eukaryomeSpecies.rel, pays=="RCA")
      
      core.taxa <- core(eukaryomeSpecies.rel.A, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeSpecies.rel.A, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeSpecies.rel.A, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeSpecies.rel.A, detection = 0.0001, prevalence = 65/100, include.lowest = TRUE)
      core.taxa
      tax_table(core.taxa) #Humicola grisea, Malassezia restricta
      
      core.taxa <- core(eukaryomeSpecies.rel.A, detection = 0.0001, prevalence = 10/100, include.lowest = TRUE) 
      core.taxa
      View(tax_table(core.taxa)) #27 taxa
      
      core.taxa <- core(eukaryomeSpecies.rel.B, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeSpecies.rel.B, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeSpecies.rel.B, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeSpecies.rel.B, detection = 0.0001, prevalence = 65/100, include.lowest = TRUE) 
      core.taxa
      View(tax_table(core.taxa)) #, Trichosporon asahii
      
      core.taxa <- core(eukaryomeSpecies.rel.B, detection = 0.0001, prevalence = 10/100, include.lowest = TRUE) 
      core.taxa
      View(tax_table(core.taxa)) #20 taxa
      
      core.taxa <- core(eukaryomeSpecies.rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeSpecies.rel, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeSpecies.rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
      core.taxa <- core(eukaryomeSpecies.rel, detection = 0.0001, prevalence = 50/100, include.lowest = TRUE)
      core.taxa <- core(eukaryomeSpecies.rel, detection = 0.0001, prevalence = 25/100, include.lowest = TRUE)
      core.taxa
      View(tax_table(core.taxa)) #5 species are conserved in 25 percent of samples
      
      core.taxa <- core(eukaryomeSpecies.rel, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE)
      core.taxa
      tax_table(core.taxa) #5 species are conserved in 25 percent of samples
      
      core.taxa <- core(eukaryomeSpecies.rel, detection = 0.0001, prevalence = 10/100, include.lowest = TRUE) 
      core.taxa
      tax_table(core.taxa) #23 taxa
      

      ## Anything conserved on Family rank?
      dfcleanfeces_Family<-tax_glom(dfcleanfeces, "Family")
      eukaryomeFamily.rel <- microbiome::transform(dfcleanfeces_Family, "compositional")
      core.taxa <- core(eukaryomeFamily.rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeFamily.rel, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeFamily.rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
      core.taxa
      View(tax_table(core.taxa)) #only 1
      
      eukaryomeFamily.rel.A=subset_samples(eukaryomeFamily.rel, pays=="Madagascar")
      eukaryomeFamily.rel.B=subset_samples(eukaryomeFamily.rel, pays=="RCA")
      
      core.taxa <- core(eukaryomeFamily.rel.A, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeFamily.rel.A, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeFamily.rel.A, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE) 
      core.taxa
      View(tax_table(core.taxa)) #3, va. Saccharomycetes and Malasseziaceae
      
      core.taxa <- core(eukaryomeFamily.rel.B, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeFamily.rel.B, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeFamily.rel.B, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE) 
      core.taxa
      View(tax_table(core.taxa)) #2 families,two conserved with Tana
      
      
      ## Anything conserved on Genus level?
      dfcleanfeces_Genus<-tax_glom(dfcleanfeces, "Genus")
      eukaryomeGenus.rel <- microbiome::transform(dfcleanfeces_Genus, "compositional")
      core.taxa <- core(eukaryomeGenus.rel, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeGenus.rel, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeGenus.rel, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE)
      core.taxa <- core(eukaryomeGenus.rel, detection = 0.00000001, prevalence = 75/100, include.lowest = TRUE)
      core.taxa # 1 taxa are shared
      View(tax_table(core.taxa)) # 
      core.taxa <- core(eukaryomeGenus.rel, detection = 0.00000001, prevalence = 50/100, include.lowest = TRUE)
      core.taxa # 5 taxa are shared
      View(tax_table(core.taxa)) # 
      
      core.taxa <- core(eukaryomeGenus.rel, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE)
      core.taxa # 12 taxa are shared
      View(tax_table(core.taxa)) # 
      
      core.taxa <- core(eukaryomeGenus.rel, detection = 0.00000001, prevalence = 10/100, include.lowest = TRUE)
      core.taxa # 22 taxa are shared
      View(tax_table(core.taxa)) # 
      
      eukaryomeGenus.rel.A=subset_samples(eukaryomeGenus.rel, pays=="Madagascar")
      eukaryomeGenus.rel.B=subset_samples(eukaryomeGenus.rel, pays=="RCA")
      
      core.taxa <- core(eukaryomeGenus.rel.A, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeGenus.rel.A, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeGenus.rel.A, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeGenus.rel.A, detection = 0.0001, prevalence = 65/100, include.lowest = TRUE)
      core.taxa # 4 taxon
      core.taxa <- core(eukaryomeGenus.rel.A, detection = 0.00000001, prevalence = 50/100, include.lowest = TRUE)
      core.taxa # 5 taxa
      View(tax_table(core.taxa)) # 
      
      core.taxa <- core(eukaryomeGenus.rel.A, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE)
      core.taxa # 14 taxa
      View(tax_table(core.taxa))
      
      core.taxa <- core(eukaryomeGenus.rel.A, detection = 0.00000001, prevalence = 10/100, include.lowest = TRUE)
      core.taxa # 22 taxa
      View(tax_table(core.taxa))
      
      core.taxa <- core(eukaryomeGenus.rel.B, detection = 0.0001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeGenus.rel.B, detection = 0.00001, prevalence = 90/100, include.lowest = TRUE) 
      core.taxa <- core(eukaryomeGenus.rel.B, detection = 0.0001, prevalence = 75/100, include.lowest = TRUE) 
      core.taxa
      tax_table(core.taxa) #only Trichosporon
      
      core.taxa <- core(eukaryomeGenus.rel.B, detection = 0.00000001, prevalence = 50/100, include.lowest = TRUE)
      core.taxa # 4 taxa
      View(tax_table(core.taxa)) # two uncultured Saccharomycetales, Pichia, Blastocystis, Saccharomyces
      
      core.taxa <- core(eukaryomeGenus.rel.B, detection = 0.00000001, prevalence = 25/100, include.lowest = TRUE)
      core.taxa # 15 taxa
      View(tax_table(core.taxa)) # 
      
      core.taxa <- core(eukaryomeGenus.rel.B, detection = 0.00000001, prevalence = 10/100, include.lowest = TRUE)
      core.taxa # 20 taxa
      View(tax_table(core.taxa))
      
  #### Generate tables of rel. abundance by different categories ####    
      
      
      #### feces ####
      dfcleanfeces_rel<-subset_samples(dfcleanfeces_rel, pays!="")
      # Plot Phylum content
      
      dfcleanfeces_rel.p <- tax_glom(dfcleanfeces_rel, "Phylum")
      test<-filter_taxa(dfcleanfeces_rel.p, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      
      dfcleanfeces_rel.p= test2  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Phylum) 
      
      pdf("AbundancetablePhylumpersampleMicrobiomeInsightsfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleanfeces_rel.p, aes(x = Sample, y = Abundance, fill = Phylum)) + 
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
      
      pdf("AbundancetablePhylumpersamplepaysstuntedMicrobiomeInsightsfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleanfeces_rel.p, aes(x = Sample, y = Abundance, fill = Phylum)) + 
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
      
      # Plot Order content
      
      dfcleanfeces_rel.4 <- tax_glom(dfcleanfeces_rel, "Order")
      test<-filter_taxa(dfcleanfeces_rel.4, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      
      dfcleanfeces_rel.4 = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      dfcleanfeces_rel.4_melt= dfcleanfeces_rel.4  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      
      pdf("AbundancetableOrderpersampleMicrobiomeInsightsfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleanfeces_rel.4_melt, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples") 
      dev.off()
      
      pdf("AbundancetableOrderpersamplepaysstuntedMicrobiomeInsightsfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleanfeces_rel.4_melt, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples") 
      dev.off()
      
      # now restricted to top 10
      TopNOTUs <- names(sort(taxa_sums(dfcleanfeces_rel.4), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleanfeces_rel.4))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      feces10 <- prune_taxa(TopNOTUs, dfcleanfeces_rel.4)
      print(feces10)
      
      feces10= feces10  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      pdf("AbundancetableOrderpersampleMicrobiomeInsightstop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(feces10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples, top 10 orders only") 
      dev.off()
      
      pdf("AbundancetableOrderpersamplepaysstuntedMicrobiomeInsightstop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(feces10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples top 10 orders") 
      dev.off()
      
      # only restricted on potential gut fungi
      dfcleanfeces_rel_gut<-subset_taxa(dfcleanfeces_rel, gut_environmental_diet=="gut" | gut_environmental_diet=="possible_gut")
      dfcleanfeces_rel.4_gut <- tax_glom(dfcleanfeces_rel_gut, "Order")
      test<-filter_taxa(dfcleanfeces_rel.4_gut, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleanfeces_rel.4_gut = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      TopNOTUs <- names(sort(speciesSums(dfcleanfeces_rel.4_gut), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleanfeces_rel.4_gut))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      feces10 <- prune_taxa(TopNOTUs, dfcleanfeces_rel.4_gut)
      print(feces10)
      
      feces10= feces10  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      
      
      pdf("AbundancetableOrderpersampleMicrobiomeInsightsguttop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(feces10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples, top 10 orders only from gut fungi") 
      dev.off()
      
      pdf("AbundancetableOrderpersamplepaysstuntedMicrobiomeInsightsguttop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(feces10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples top 10 orders only from gut fungi") 
      dev.off()
      
      # only restricted on potential environmental fungi
      dfcleanfeces_rel_env<-subset_taxa(dfcleanfeces_rel, gut_environmental_diet=="environmental_diet")
      dfcleanfeces_rel.4_env <- tax_glom(dfcleanfeces_rel_env, "Order")
      test<-filter_taxa(dfcleanfeces_rel.4_env, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleanfeces_rel.4_env = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      TopNOTUs <- names(sort(taxa_sums(dfcleanfeces_rel.4_env), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleanfeces_rel.4_env))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      feces10env <- prune_taxa(TopNOTUs, dfcleanfeces_rel.4_env)
      print(feces10env)
      
      feces10env= feces10env  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      pdf("AbundancetableOrderpersampleMicrobiomeInsightsenvtop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(feces10env, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples, top 10 orders only from environmental fungi") 
      dev.off()
      
      pdf("AbundancetableOrderpersamplepaysstuntedMicrobiomeInsightsenvtop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(feces10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples top 10 orders only from environmental fungi") 
      dev.off()
      
      # Plot Family content
      dfcleanfeces_rel.4 <- tax_glom(dfcleanfeces_rel, "Family")
      test<-filter_taxa(dfcleanfeces_rel.4, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleanfeces_rel.4 = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      dfcleanfeces_rel.4_melt= dfcleanfeces_rel.4  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      
      pdf("AbundancetableFamilypersampleMicrobiomeInsightsfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleanfeces_rel.4_melt, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples") 
      dev.off()
      
      pdf("AbundancetableFamilypersamplepaysstuntedMicrobiomeInsightsfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleanfeces_rel.4_melt, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples") 
      dev.off()
      
      # now restricted to top 10
      TopNOTUs <- names(sort(speciesSums(dfcleanfeces_rel.4), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleanfeces_rel.4))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      feces10 <- prune_taxa(TopNOTUs, dfcleanfeces_rel.4)
      print(feces10)
      
      feces10= feces10  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      pdf("AbundancetableFamilypersampleMicrobiomeInsightstop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(feces10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples, top 10 Familys only") 
      dev.off()
      
      pdf("AbundancetableFamilypersamplepaysstuntedMicrobiomeInsightstop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(feces10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples top 10 Familys") 
      dev.off()
      
      # only restricted on potential gut fungi
      dfcleanfeces_rel_gut<-subset_taxa(dfcleanfeces_rel, gut_environmental_diet=="gut" | gut_environmental_diet=="possible_gut")
      dfcleanfeces_rel.4_gut <- tax_glom(dfcleanfeces_rel_gut, "Family")
      test<-filter_taxa(dfcleanfeces_rel.4_gut, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleanfeces_rel.4_gut = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      
      TopNOTUs <- names(sort(speciesSums(dfcleanfeces_rel.4_gut), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleanfeces_rel.4_gut))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      feces10 <- prune_taxa(TopNOTUs, dfcleanfeces_rel.4_gut)
      print(feces10)
      
      feces10= feces10  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      
      
      pdf("AbundancetableFamilypersampleMicrobiomeInsightsguttop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(feces10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples, top 10 Familys only from gut fungi") 
      dev.off()
      
      pdf("AbundancetableFamilypersamplepaysstuntedMicrobiomeInsightsguttop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(feces10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples top 10 Familys only from gut fungi") 
      dev.off()
      
      # only restricted on potential environmental fungi
      dfcleanfeces_rel_env<-subset_taxa(dfcleanfeces_rel, gut_environmental_diet=="environmental_diet")
      dfcleanfeces_rel.4_env <- tax_glom(dfcleanfeces_rel_env, "Family")
      test<-filter_taxa(dfcleanfeces_rel.4_env, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleanfeces_rel.4_env = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      TopNOTUs <- names(sort(taxa_sums(dfcleanfeces_rel.4_env), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleanfeces_rel.4_env))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      feces10env <- prune_taxa(TopNOTUs, dfcleanfeces_rel.4_env)
      print(feces10env)
      
      feces10env= feces10env  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      pdf("AbundancetableFamilypersampleMicrobiomeInsightsenvtop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(feces10env, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples, top 10 Familys only from environmental fungi") 
      dev.off()
      
      pdf("AbundancetableFamilypersamplepaysstuntedMicrobiomeInsightsenvtop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(feces10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples top 10 Familys only from environmental fungi") 
      dev.off()
      
   
      
      #### duodenal ####
      ### now on duodenal samples and country
           
      # make ordination using PCoA log 10 transformed
      GP.ord <- ordinate(dfduodenalsrar5000transf, "PCoA", "bray")
      p1 = plot_ordination(dfduodenalsrar5000transf, GP.ord, type="taxa", color="Class", title="taxa")
      print(p1)
      
      p1 = plot_ordination(dfduodenalsrar5000transf, GP.ord, type="taxa", color="gut_environmental_diet", title="taxa")
      print(p1)
      
      p2 = plot_ordination(dfduodenalsrar5000transf, GP.ord, type="samples", color="pays", shape="stunted") 
      p2 + geom_point(size=1) + ggtitle("samples")
      
      p2 = plot_ordination(dfduodenalsrar5000transf, GP.ord, type="samples", color="pays", shape="stunted") 
      p2 + geom_polygon(aes(fill=stunted)) + geom_point(size=5) + ggtitle("samples")
      
      
         # Plot Phylum content
      
      dfcleanduodenal_rel.p <- tax_glom(dfcleanduodenal_rel, "Phylum")
      test<-filter_taxa(dfcleanduodenal_rel.p, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      
      dfcleanduodenal_rel.p= test2  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Phylum) 
      
      pdf("AbundancetablePhylumpersampleMicrobiomeInsightsduodenum.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleanduodenal_rel.p, aes(x = Sample, y = Abundance, fill = Phylum)) + 
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
      
      pdf("AbundancetablePhylumpersamplepaysstuntedMicrobiomeInsightsduodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleanduodenal_rel.p, aes(x = Sample, y = Abundance, fill = Phylum)) + 
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
      
      # Plot Order content
      
      dfcleanduodenal_rel.4 <- tax_glom(dfcleanduodenal_rel, "Order")
      test<-filter_taxa(dfcleanduodenal_rel.4, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      
      dfcleanduodenal_rel.4 = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      dfcleanduodenal_rel.4_melt= dfcleanduodenal_rel.4  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      
      pdf("AbundancetableOrderpersampleMicrobiomeInsightsduodenum.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleanduodenal_rel.4_melt, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples") 
      dev.off()
      
      pdf("AbundancetableOrderpersamplepaysstuntedMicrobiomeInsightsduodenum.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleanduodenal_rel.4_melt, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples") 
      dev.off()
      
      # now restricted to top 10
      TopNOTUs <- names(sort(taxa_sums(dfcleanduodenal_rel.4), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleanduodenal_rel.4))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      duodenal10 <- prune_taxa(TopNOTUs, dfcleanduodenal_rel.4)
      print(duodenal10)
      
      duodenal10= duodenal10  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      pdf("AbundancetableOrderpersampleMicrobiomeInsightstop10duodenum.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(duodenal10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples, top 10 orders only") 
      dev.off()
      
      pdf("AbundancetableOrderpersamplepaysstuntedMicrobiomeInsightstop10duodenum.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(duodenal10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples top 10 orders") 
      dev.off()
      
      # only restricted on potential gut fungi
      dfcleanduodenal_rel_gut<-subset_taxa(dfcleanduodenal_rel, gut_environmental_diet=="gut" | gut_environmental_diet=="possible_gut")
      dfcleanduodenal_rel.4_gut <- tax_glom(dfcleanduodenal_rel_gut, "Order")
      test<-filter_taxa(dfcleanduodenal_rel.4_gut, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleanduodenal_rel.4_gut = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      TopNOTUs <- names(sort(speciesSums(dfcleanduodenal_rel.4_gut), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleanduodenal_rel.4_gut))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      duodenal10 <- prune_taxa(TopNOTUs, dfcleanduodenal_rel.4_gut)
      print(duodenal10)
      
      duodenal10= duodenal10  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      
      
      pdf("AbundancetableOrderpersampleMicrobiomeInsightsguttop10duodenum.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(duodenal10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples, top 10 orders only from gut fungi") 
      dev.off()
      
      pdf("AbundancetableOrderpersamplepaysstuntedMicrobiomeInsightsguttop10duodenum.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(duodenal10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples top 10 orders only from gut fungi") 
      dev.off()
      
      # only restricted on potential environmental fungi
      dfcleanduodenal_rel_env<-subset_taxa(dfcleanduodenal_rel, gut_environmental_diet=="environmental_diet")
      dfcleanduodenal_rel.4_env <- tax_glom(dfcleanduodenal_rel_env, "Order")
      test<-filter_taxa(dfcleanduodenal_rel.4_env, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleanduodenal_rel.4_env = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      TopNOTUs <- names(sort(taxa_sums(dfcleanduodenal_rel.4_env), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleanduodenal_rel.4_env))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      duodenal10env <- prune_taxa(TopNOTUs, dfcleanduodenal_rel.4_env)
      print(duodenal10env)
      
      duodenal10env= duodenal10env  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      pdf("AbundancetableOrderpersampleMicrobiomeInsightsenvtop10duodenum.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(duodenal10env, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples, top 10 orders only from environmental fungi") 
      dev.off()
      
      pdf("AbundancetableOrderpersamplepaysstuntedMicrobiomeInsightsenvtop10duodenum.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(duodenal10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples top 10 orders only from environmental fungi") 
      dev.off()
      
      # Plot Family content
      dfcleanduodenal_rel.4 <- tax_glom(dfcleanduodenal_rel, "Family")
      test<-filter_taxa(dfcleanduodenal_rel.4, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleanduodenal_rel.4 = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      dfcleanduodenal_rel.4_melt= dfcleanduodenal_rel.4  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      
      pdf("AbundancetableFamilypersampleMicrobiomeInsightsduodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleanduodenal_rel.4_melt, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples") 
      dev.off()
      
      pdf("AbundancetableFamilypersamplepaysstuntedMicrobiomeInsightsduodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleanduodenal_rel.4_melt, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples") 
      dev.off()
      
      # now restricted to top 10
      TopNOTUs <- names(sort(speciesSums(dfcleanduodenal_rel.4), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleanduodenal_rel.4))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      duodenal10 <- prune_taxa(TopNOTUs, dfcleanduodenal_rel.4)
      print(duodenal10)
      
      duodenal10= duodenal10  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      pdf("AbundancetableFamilypersampleMicrobiomeInsightstop10duodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(duodenal10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples, top 10 Familys only") 
      dev.off()
      
      pdf("AbundancetableFamilypersamplepaysstuntedMicrobiomeInsightstop10duodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(duodenal10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples top 10 Familys") 
      dev.off()
      
      # only restricted on potential gut fungi
      dfcleanduodenal_rel_gut<-subset_taxa(dfcleanduodenal_rel, gut_environmental_diet=="gut" | gut_environmental_diet=="possible_gut")
      dfcleanduodenal_rel.4_gut <- tax_glom(dfcleanduodenal_rel_gut, "Family")
      test<-filter_taxa(dfcleanduodenal_rel.4_gut, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleanduodenal_rel.4_gut = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      
      TopNOTUs <- names(sort(speciesSums(dfcleanduodenal_rel.4_gut), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleanduodenal_rel.4_gut))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      duodenal10 <- prune_taxa(TopNOTUs, dfcleanduodenal_rel.4_gut)
      print(duodenal10)
      
      duodenal10= duodenal10  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      
      
      pdf("AbundancetableFamilypersampleMicrobiomeInsightsguttop10duodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(duodenal10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples, top 10 Familys only from gut fungi") 
      dev.off()
      
      pdf("AbundancetableFamilypersamplepaysstuntedMicrobiomeInsightsguttop10duodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(duodenal10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples top 10 Familys only from gut fungi") 
      dev.off()
      
      # only restricted on potential environmental fungi
      dfcleanduodenal_rel_env<-subset_taxa(dfcleanduodenal_rel, gut_environmental_diet=="environmental_diet")
      dfcleanduodenal_rel.4_env <- tax_glom(dfcleanduodenal_rel_env, "Family")
      test<-filter_taxa(dfcleanduodenal_rel.4_env, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleanduodenal_rel.4_env = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      TopNOTUs <- names(sort(taxa_sums(dfcleanduodenal_rel.4_env), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleanduodenal_rel.4_env))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      duodenal10env <- prune_taxa(TopNOTUs, dfcleanduodenal_rel.4_env)
      print(duodenal10env)
      
      duodenal10env= duodenal10env  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      pdf("AbundancetableFamilypersampleMicrobiomeInsightsenvtop10duodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(duodenal10env, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples, top 10 Familys only from environmental fungi") 
      dev.off()
      
      pdf("AbundancetableFamilypersamplepaysstuntedMicrobiomeInsightsenvtop10duodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(duodenal10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples top 10 Familys only from environmental fungi") 
      dev.off()
      
      
      
      
      
      
      
      #### gastric ####
      dfcleangastric_rel<-subset_samples(dfcleangastric_rel, pays!="")
      # Plot Phylum content
      
      dfcleangastric_rel.p <- tax_glom(dfcleangastric_rel, "Phylum")
      test<-filter_taxa(dfcleangastric_rel.p, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      
      dfcleangastric_rel.p= test2  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Phylum) 
      
      pdf("AbundancetablePhylumpersampleMicrobiomeInsightsgastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleangastric_rel.p, aes(x = Sample, y = Abundance, fill = Phylum)) + 
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
      
      pdf("AbundancetablePhylumpersamplepaysstuntedMicrobiomeInsightsgastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleangastric_rel.p, aes(x = Sample, y = Abundance, fill = Phylum)) + 
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
      
      # Plot Order content
      
      dfcleangastric_rel.4 <- tax_glom(dfcleangastric_rel, "Order")
      test<-filter_taxa(dfcleangastric_rel.4, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      
      dfcleangastric_rel.4 = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      dfcleangastric_rel.4_melt= dfcleangastric_rel.4  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      
      pdf("AbundancetableOrderpersampleMicrobiomeInsightsgastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleangastric_rel.4_melt, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples") 
      dev.off()
      
      pdf("AbundancetableOrderpersamplepaysstuntedMicrobiomeInsightsgastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleangastric_rel.4_melt, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples") 
      dev.off()
      
      # now restricted to top 10
      TopNOTUs <- names(sort(taxa_sums(dfcleangastric_rel.4), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleangastric_rel.4))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      gastric10 <- prune_taxa(TopNOTUs, dfcleangastric_rel.4)
      print(gastric10)
      
      gastric10= gastric10  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      pdf("AbundancetableOrderpersampleMicrobiomeInsightstop10gastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(gastric10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples, top 10 orders only") 
      dev.off()
      
      pdf("AbundancetableOrderpersamplepaysstuntedMicrobiomeInsightstop10gastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(gastric10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples top 10 orders") 
      dev.off()
      
      # only restricted on potential gut fungi
      dfcleangastric_rel_gut<-subset_taxa(dfcleangastric_rel, gut_environmental_diet=="gut" | gut_environmental_diet=="possible_gut")
      dfcleangastric_rel.4_gut <- tax_glom(dfcleangastric_rel_gut, "Order")
      test<-filter_taxa(dfcleangastric_rel.4_gut, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleangastric_rel.4_gut = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      TopNOTUs <- names(sort(speciesSums(dfcleangastric_rel.4_gut), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleangastric_rel.4_gut))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      gastric10 <- prune_taxa(TopNOTUs, dfcleangastric_rel.4_gut)
      print(gastric10)
      
      gastric10= gastric10  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      
      
      pdf("AbundancetableOrderpersampleMicrobiomeInsightsguttop10gastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(gastric10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples, top 10 orders only from gut fungi") 
      dev.off()
      
      pdf("AbundancetableOrderpersamplepaysstuntedMicrobiomeInsightsguttop10gastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(gastric10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples top 10 orders only from gut fungi") 
      dev.off()
      
      # only restricted on potential environmental fungi
      dfcleangastric_rel_env<-subset_taxa(dfcleangastric_rel, gut_environmental_diet=="environmental_diet")
      dfcleangastric_rel.4_env <- tax_glom(dfcleangastric_rel_env, "Order")
      test<-filter_taxa(dfcleangastric_rel.4_env, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleangastric_rel.4_env = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      TopNOTUs <- names(sort(taxa_sums(dfcleangastric_rel.4_env), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleangastric_rel.4_env))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      gastric10env <- prune_taxa(TopNOTUs, dfcleangastric_rel.4_env)
      print(gastric10env)
      
      gastric10env= gastric10env  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      pdf("AbundancetableOrderpersampleMicrobiomeInsightsenvtop10gastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(gastric10env, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples, top 10 orders only from environmental fungi") 
      dev.off()
      
      pdf("AbundancetableOrderpersamplepaysstuntedMicrobiomeInsightsenvtop10gastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(gastric10, aes(x = Sample, y = Abundance, fill = Order)) + 
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
        ggtitle("Order Composition of Afribiota samples top 10 orders only from environmental fungi") 
      dev.off()
      
      # Plot Family content
      dfcleangastric_rel.4 <- tax_glom(dfcleangastric_rel, "Family")
      test<-filter_taxa(dfcleangastric_rel.4, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleangastric_rel.4 = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      dfcleangastric_rel.4_melt= dfcleangastric_rel.4  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      
      pdf("AbundancetableFamilypersampleMicrobiomeInsightsgastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleangastric_rel.4_melt, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples") 
      dev.off()
      
      pdf("AbundancetableFamilypersamplepaysstuntedMicrobiomeInsightsgastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(dfcleangastric_rel.4_melt, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples") 
      dev.off()
      
      # now restricted to top 10
      TopNOTUs <- names(sort(speciesSums(dfcleangastric_rel.4), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleangastric_rel.4))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      gastric10 <- prune_taxa(TopNOTUs, dfcleangastric_rel.4)
      print(gastric10)
      
      gastric10= gastric10  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      pdf("AbundancetableFamilypersampleMicrobiomeInsightstop10gastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(gastric10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples, top 10 Familys only") 
      dev.off()
      
      pdf("AbundancetableFamilypersamplepaysstuntedMicrobiomeInsightstop10gastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(gastric10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples top 10 Familys") 
      dev.off()
      
      # only restricted on potential gut fungi
      dfcleangastric_rel_gut<-subset_taxa(dfcleangastric_rel, gut_environmental_diet=="gut" | gut_environmental_diet=="possible_gut")
      dfcleangastric_rel.4_gut <- tax_glom(dfcleangastric_rel_gut, "Family")
      test<-filter_taxa(dfcleangastric_rel.4_gut, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleangastric_rel.4_gut = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      
      TopNOTUs <- names(sort(speciesSums(dfcleangastric_rel.4_gut), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleangastric_rel.4_gut))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      gastric10 <- prune_taxa(TopNOTUs, dfcleangastric_rel.4_gut)
      print(gastric10)
      
      gastric10= gastric10  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      
      
      pdf("AbundancetableFamilypersampleMicrobiomeInsightsguttop10gastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(gastric10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples, top 10 Familys only from gut fungi") 
      dev.off()
      
      pdf("AbundancetableFamilypersamplepaysstuntedMicrobiomeInsightsguttop10gastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(gastric10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples top 10 Familys only from gut fungi") 
      dev.off()
      
      # only restricted on potential environmental fungi
      dfcleangastric_rel_env<-subset_taxa(dfcleangastric_rel, gut_environmental_diet=="environmental_diet")
      dfcleangastric_rel.4_env <- tax_glom(dfcleangastric_rel_env, "Family")
      test<-filter_taxa(dfcleangastric_rel.4_env, function(x) sum(x)>0, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      dfcleangastric_rel.4_env = transform_sample_counts(test2, function(x) 100 * x/sum(x))
      
      TopNOTUs <- names(sort(taxa_sums(dfcleangastric_rel.4_env), TRUE)[1:10])
      taxtable<-as.data.frame(tax_table(dfcleangastric_rel.4_env))
      TopNOTUannotated<-filter(taxtable, (row.names(taxtable) %in% TopNOTUs))
      View(TopNOTUannotated)
      
      gastric10env <- prune_taxa(TopNOTUs, dfcleangastric_rel.4_env)
      print(gastric10env)
      
      gastric10env= gastric10env  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      pdf("AbundancetableFamilypersampleMicrobiomeInsightsenvtop10gastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(gastric10env, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples, top 10 Familys only from environmental fungi") 
      dev.off()
      
      pdf("AbundancetableFamilypersamplepaysstuntedMicrobiomeInsightsenvtop10gastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 40, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(gastric10, aes(x = Sample, y = Abundance, fill = Family)) + 
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
        ggtitle("Family Composition of Afribiota samples top 10 Familys only from environmental fungi") 
      dev.off()
      
      
      
      
      
      
      
      #### by Sample Type   ####
      sample_data(dfcleanstunted)$SampleType <- factor(sample_data(dfcleanstunted)$SampleType, levels = c("gastric", "duodenal", "feces"))
      dfcleansample = merge_samples(dfcleanstunted, "SampleType")
      
      # phylum level
      dfclean_p <- tax_glom(dfcleansample, "Phylum")
      
      dfcleanrel_p = transform_sample_counts(dfclean_p, function(x) 100 * x/sum(x))
      sample_data(dfcleanrel_p)$SampleType <- factor(sample_data(dfcleanrel_p)$SampleType, levels = c("gastric", "duodenal", "feces"))
      
      dfcleanrel_p_s= dfcleanrel_p  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Phylum) 
      
      class_colors <- c(
        "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
        "#AD6F3B", "#673770","black", "#652926", "#C84248", 
        "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
      )
      
      dfcleanrel_p_s$Sample <- factor(dfcleanrel_p_s$Sample, levels=c("gastric", "duodenal", "feces"))
      
      pdf("AbundancetablePhylumSampleTypenopruning.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 6, # define plot width and height. completely up to user.
          height = 8) 
      ggplot(dfcleanrel_p_s, aes(x =Sample, y = Abundance, fill = Phylum)) + 
        theme_bw() +
        theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
        theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
        geom_bar(stat = "identity", width = 1.0) +
        scale_y_continuous(expand = c(0.01,0.01)) +
        theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
        ylab("Relative Abundance") +
        xlab("SampleType") + ggtitle ("Relative Abundance of Phyla according to sample type")
      dev.off()
      
      # order level
      dfclean_o <- tax_glom(dfcleansample, "Order")
      dfcleanrel_o = transform_sample_counts(dfclean_o, function(x) 100 * x/sum(x))
      
      dfcleanrel_o_s= dfcleanrel_o  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      dfcleanrel_o_s$Sample <- factor(dfcleanrel_o_s$Sample, levels=c("gastric", "duodenal", "feces"))
      
      
      class_colors <- c(
        "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
        "#AD6F3B", "#673770","black", "#652926", "#C84248", 
        "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
      )
      
      pdf("AbundancetableOrderSampleTypenopruning.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 8, # define plot width and height. completely up to user.
          height = 6) 
      ggplot(dfcleanrel_o_s, aes(x =Sample, y = Abundance, fill = Order)) + 
        theme_bw() +
        theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
        theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
        geom_bar(stat = "identity", width = 1.0) +
        scale_y_continuous(expand = c(0.01,0.01)) +
        theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
        ylab("Relative Abundance") +
        xlab("SampleType") + ggtitle ("Relative Abundance of orders according to sample type")
      dev.off()
      
      #### by Sample Type reduced dataset  ####
      sample_data(dfcleanstunted)$SampleType <- factor(sample_data(dfcleanstunted)$SampleType, levels = c("gastric", "duodenal", "feces"))
      
      dfcleanstunted_red<-subset_samples(dfcleanstunted, sample_data(dfcleanstunted)$id %in% test_feces2$id)
      dfcleanstunted_red # 130 samples are remaining!
      dfcleanstunted_red<-subset_samples(dfcleanstunted_red, row.names(sample_data(dfcleanstunted_red))!="AG-CPB180_S260")
      table(sample_data(dfcleanstunted_red)$SampleType) # ok!
      
      dfcleansample = merge_samples(dfcleanstunted_red, "SampleType")
      
      # phylum level
      dfclean_p <- tax_glom(dfcleansample, "Phylum")
      
      dfcleanrel_p = transform_sample_counts(dfclean_p, function(x) 100 * x/sum(x))
      sample_data(dfcleanrel_p)$SampleType <- factor(sample_data(dfcleanrel_p)$SampleType, levels = c("gastric", "duodenal", "feces"))
      
      dfcleanrel_p_s= dfcleanrel_p  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Phylum) 
      
      class_colors <- c(
        "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
        "#AD6F3B", "#673770","black", "#652926", "#C84248", 
        "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
      )
      
      dfcleanrel_p_s$Sample <- factor(dfcleanrel_p_s$Sample, levels=c("gastric", "duodenal", "feces"))
      
      pdf("AbundancetablePhylumSampleTypenopruning_red.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 8, # define plot width and height. completely up to user.
          height = 8) 
      ggplot(dfcleanrel_p_s, aes(x =Sample, y = Abundance, fill = Phylum)) + 
        theme_bw() +
        theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
        theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
        geom_bar(stat = "identity", width = 1.0) +
        scale_y_continuous(expand = c(0.01,0.01)) +
        theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
        ylab("Relative Abundance") +
        xlab("SampleType") + ggtitle ("Relative Abundance of Phyla according to sample type")
      dev.off()
      
      # order level
      dfclean_o <- tax_glom(dfcleanstunted_red, "Order")
      dfcleanrel_o = transform_sample_counts(dfclean_o, function(x) 100 * x/sum(x))
      
      dfcleanrel_o_s= dfcleanrel_o  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      dfcleanrel_o_s$Sample <- factor(dfcleanrel_o_s$Sample, levels=c("gastric", "duodenal", "feces"))
      
      
      class_colors <- c(
        "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
        "#AD6F3B", "#673770","black", "#652926", "#C84248", 
        "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
      )
      
      pdf("AbundancetableOrderSampleTypenopruning_red.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 11, # define plot width and height. completely up to user.
          height = 8) 
      ggplot(dfcleanrel_o_s, aes(x =Sample, y = Abundance, fill = Order)) + 
        theme_bw() +
        theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
        theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
        geom_bar(stat = "identity", width = 1.0) +
        scale_y_continuous(expand = c(0.01,0.01)) +
        theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
        ylab("Relative Abundance") +
        xlab("SampleType") + ggtitle ("Relative Abundance of orders according to sample type")
      dev.off()
      
      #### on different sampling sites stratified: country of origin ####
      ### feces
      # phylum level
      dfcleanfecespays = merge_samples(dfcleanfeces, "pays")
      dfcleanfecespays_p <- tax_glom(dfcleanfecespays, "Phylum")
      
      dfcleanfecespaysrel_p = transform_sample_counts(dfcleanfecespays_p, function(x) 100 * x/sum(x))
      
      
      dfcleanfecespaysrel_p_s= dfcleanfecespaysrel_p  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Phylum) 
      
      class_colors <- c(
        "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
        "#AD6F3B", "#673770","black", "#652926", "#C84248", 
        "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
      )
      
      
      pdf("AbundancetablePhylumpaysnopruningfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 11, # define plot width and height. completely up to user.
          height = 8) 
      ggplot(dfcleanfecespaysrel_p_s, aes(x =Sample, y = Abundance, fill = Phylum)) + 
        theme_bw() +
        theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
        theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
        geom_bar(stat = "identity", width = 1.0) +
        scale_y_continuous(expand = c(0.01,0.01)) +
        theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
        ylab("Relative Abundance") +
        xlab("pays") + ggtitle ("Relative Abundance of feces by country of origin")
      dev.off()
      
      # order level
      dfclean_o <- tax_glom(dfcleanfecespays, "Order")
      dfcleanrel_o = transform_sample_counts(dfclean_o, function(x) 100 * x/sum(x))
      
      dfcleanrel_o_s= dfcleanrel_o  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      class_colors <- c(
        "pink", "white", "#CBD588", "yellow", "red", #5F7FC7", "orange", "#508578", "green", #CD9BCD",
        "#AD6F3B", "#673770","black", "#652926", "#C84248", "cyan", "magenta", "turquoise",
        "#8569D5","blue", "#5E738F","#D1A33D", "purple", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
      )
      
      pdf("AbundancetableOrderpaysnopruningfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 8, # define plot width and height. completely up to user.
          height = 8) 
      ggplot(dfcleanrel_o_s, aes(x =Sample, y = Abundance, fill = Order)) + 
        theme_bw() +
        theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
        theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
        geom_bar(stat = "identity", width = 1.0) +
        scale_fill_manual(values = class_colors) +
        scale_y_continuous(expand = c(0.01,0.01)) +
        theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
        ylab("Relative Abundance") +
        xlab("pays") + ggtitle ("Relative Abundance of orders according to country of origin")
      dev.off()
      
      # order level
      dfclean_o <- tax_glom(dfcleanfecespays, "Order")
      dfcleanrel_o = transform_sample_counts(dfclean_o, function(x) 100 * x/sum(x))
      
      dfcleanrel_o_s= dfcleanrel_o  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      class_colors <- c(
        "pink", "white", "#CBD588", "yellow", "red", #5F7FC7", "orange", "#508578", "green", #CD9BCD",
        "#AD6F3B", "#673770","black", "#652926", "#C84248", "cyan", "magenta", "turquoise",
        "#8569D5","blue", "#5E738F","#D1A33D", "purple", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
      )
      
      pdf("AbundancetableOrderpaysnopruningfeces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 8, # define plot width and height. completely up to user.
          height = 8) 
      ggplot(dfcleanrel_o_s, aes(x =Sample, y = Abundance, fill = Order)) + 
        theme_bw() +
        theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
        theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
        geom_bar(stat = "identity", width = 1.0) +
        scale_fill_manual(values = class_colors) +
        scale_y_continuous(expand = c(0.01,0.01)) +
        theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
        ylab("Relative Abundance") +
        xlab("pays") + ggtitle ("Relative Abundance of orders according to country of origin")
      dev.off()
      
      #Sort the Order by abundance and pick the top 10
      top10Order.names = names(sort(taxa_sums(dfcleanrel_o), TRUE)[1:10])
      top10OrderfOrdercountry = prune_taxa(top10Order.names, dfcleanrel_o)
      
      top10Order.names
      top10OrderfOrdercountry
      
      sample_data(top10OrderfOrdercountry)
      
      top10OrderfOrdercountry_s= top10OrderfOrdercountry  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      View(sample_data(top10OrderfOrdercountry_s))
      
      pdf("AbundancetableOrdercountrynopruningtop10feces.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 10) 
      ggplot(top10OrderfOrdercountry_s, aes(x =Sample, y = Abundance, fill = Order)) + 
        theme_bw() +
        theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
        theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=16, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=16, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0, size=24, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
        geom_bar(stat = "identity", width = 1.0) +
        scale_y_continuous(expand = c(0.01,0.01)) +
        scale_fill_manual(values = class_colors) +
        theme(legend.text=element_text(size=18, ), legend.title = element_text(size=20, face="bold")) + 
        ylab("Relative Abundance") +
        xlab("Country") + ggtitle ("Relative Abundance of fungal orders according to country, 10 most abundant in feces")
      dev.off()  
      
      
      ### duodenal
      # phylum level
      dfcleanduodenalpays = merge_samples(dfcleanduodenal, "pays")
      dfcleanduodenalpays_p <- tax_glom(dfcleanduodenalpays, "Phylum")
      
      dfcleanduodenalpaysrel_p = transform_sample_counts(dfcleanduodenalpays_p, function(x) 100 * x/sum(x))
      
      
      dfcleanduodenalpaysrel_p_s= dfcleanduodenalpaysrel_p  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Phylum) 
      
      class_colors <- c(
        "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
        "#AD6F3B", "#673770","black", "#652926", "#C84248", 
        "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
      )
      
      
      pdf("AbundancetablePhylumpaysnopruningduodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 16, # define plot width and height. completely up to user.
          height = 8) 
      ggplot(dfcleanduodenalpaysrel_p_s, aes(x =Sample, y = Abundance, fill = Phylum)) + 
        theme_bw() +
        theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
        theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
        geom_bar(stat = "identity", width = 1.0) +
        scale_y_continuous(expand = c(0.01,0.01)) +
        theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
        ylab("Relative Abundance") +
        xlab("pays") + ggtitle ("Relative Abundance of duodenal by country of origin")
      dev.off()
      
      # order level
      dfclean_o <- tax_glom(dfcleanduodenalpays, "Order")
      dfcleanrel_o = transform_sample_counts(dfclean_o, function(x) 100 * x/sum(x))
      
      dfcleanrel_o_s= dfcleanrel_o  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      class_colors <- c(
        "pink", "white", "#CBD588", "yellow", "red", #5F7FC7", "orange", "#508578", "green", #CD9BCD",
        "#AD6F3B", "#673770","black", "#652926", "#C84248", "cyan", "magenta", "turquoise",
        "#8569D5","blue", "#5E738F","#D1A33D", "purple", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
      )
      
      pdf("AbundancetableOrderpaysnopruningduodenal.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 8, # define plot width and height. completely up to user.
          height = 8) 
      ggplot(dfcleanrel_o_s, aes(x =Sample, y = Abundance, fill = Order)) + 
        theme_bw() +
        theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
        theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
        geom_bar(stat = "identity", width = 1.0) +
        scale_fill_manual(values = class_colors) +
        scale_y_continuous(expand = c(0.01,0.01)) +
        theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
        ylab("Relative Abundance") +
        xlab("pays") + ggtitle ("Relative Abundance of orders according to country of origin")
      dev.off()
      
      
      ### gastric
      # phylum level
      dfcleangastricpays = merge_samples(dfcleangastric, "pays")
       dfcleangastricpays_p <- tax_glom(dfcleangastricpays, "Phylum")
      
      dfcleangastricpaysrel_p = transform_sample_counts(dfcleangastricpays_p, function(x) 100 * x/sum(x))
      
      
      dfcleangastricpaysrel_p_s= dfcleangastricpaysrel_p  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Phylum) 
      
      class_colors <- c(
        "#CBD588", "yellow", "#5F7FC7", "orange", "#508578", "#CD9BCD",
        "#AD6F3B", "#673770","black", "#652926", "#C84248", 
        "#8569D5","blue", "#5E738F","#D1A33D", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
      )
      
      
      pdf("AbundancetablePhylumpaysnopruninggastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 16, # define plot width and height. completely up to user.
          height = 8) 
      ggplot(dfcleangastricpaysrel_p_s, aes(x =Sample, y = Abundance, fill = Phylum)) + 
        theme_bw() +
        theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
        theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
        geom_bar(stat = "identity", width = 1.0) +
        scale_y_continuous(expand = c(0.01,0.01)) +
        theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
        ylab("Relative Abundance") +
        xlab("pays") + ggtitle ("Relative Abundance of gastric by country of origin")
      dev.off()
      
      # order level
      dfclean_o <- tax_glom(dfcleangastricpays, "Order")
      dfcleanrel_o = transform_sample_counts(dfclean_o, function(x) 100 * x/sum(x))
      
      dfcleanrel_o_s= dfcleanrel_o  %>%
        psmelt() %>%                                         # Melt to long format
        arrange(Order) 
      
      class_colors <- c(
        "pink", "white", "#CBD588", "yellow", "red", #5F7FC7", "orange", "#508578", "green", #CD9BCD",
        "#AD6F3B", "#673770","black", "#652926", "#C84248", "cyan", "magenta", "turquoise",
        "#8569D5","blue", "#5E738F","#D1A33D", "purple", "#8A7C64", "#599861", "grey", "#D14285", "#DA5724"
      )
      
      pdf("AbundancetableOrderpaysnopruninggastric.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 8, # define plot width and height. completely up to user.
          height = 8) 
      ggplot(dfcleanrel_o_s, aes(x =Sample, y = Abundance, fill = Order)) + 
        theme_bw() +
        theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
        theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.text.y = element_text(angle=90, hjust = 1, size=18), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=18, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.5, size=18, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title =element_text(angle = 0, hjust = 0.5, size=22, face="bold", margin = margin(t = 0, r = 0, b = 20, l = 0)) ) +
        geom_bar(stat = "identity", width = 1.0) +
        scale_fill_manual(values = class_colors) +
        scale_y_continuous(expand = c(0.01,0.01)) +
        theme(legend.text=element_text(size=14), legend.title = element_text(size=18, face="bold")) + 
        ylab("Relative Abundance") +
        xlab("pays") + ggtitle ("Relative Abundance of orders according to country of origin")
      dev.off()
      
 #### Are given ASV associated with the nutritional status/stunting?  ####
 #### all feces, stunted  on Species level ####      
      #### Species level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. Do not filter the species ####
      eukaryome<-subset_samples(df_clean, SampleType=="feces")
      eukaryome.s <- tax_glom(eukaryome, "Species")
      eukaryome.filt.s = filter_taxa(eukaryome.s, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.s<-microbiome::transform(eukaryome.filt.s, "compositional")
      View(tax_table(eukaryome.rel.filt.s)) 
      eukaryome.rel.filt.s=subset_samples(eukaryome.rel.filt.s, sample_sums(eukaryome.rel.filt.s)!="0")
      eukaryome.rel.filt.s = filter_taxa(eukaryome.rel.filt.s, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.s))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.s)) #take metadata
      
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
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.s))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"Feces.stunted.SpeciesWilcoxMicrobiomeInsights.csv")
      View(p.res) ## signficiant: none after multiple testing!
    
      #### Species level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
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
      write.csv(corres,"Feces.stuntedSpeciesCorrelationMicrobiomeInsights.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Species level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Species level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$stunted <- meta_wilcox$stunted
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel.filt.s))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      
      #### Species level: make a loop for logistic regressions for stunting with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ age + run+  Country + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.stuntedSpecieswithinflaMicrobiomeInsights.csv")
      
      #### Species level: make a loop for logistic regressions for country with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "age", "run", "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, Country!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(Country ~ value+ run+ age + stunted + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.paysSpecieswithinflaMicrobiomeInsights.csv")
      
      #### Species level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run+ age + Country + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.stuntedSpecieswithoutinflaMicrobiomeInsights.csv")
      
  
  #### Species level: Step 1: Calprotectine Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. Do not filter the species ####
      eukaryome<-subset_samples(df_clean, SampleType=="feces")
      eukaryome.s <- tax_glom(eukaryome, "Species")
      eukaryome.filt.s = filter_taxa(eukaryome.s, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.s<-microbiome::transform(eukaryome.filt.s, "compositional")
      View(tax_table(eukaryome.rel.filt.s)) 
      eukaryome.rel.filt.s=subset_samples(eukaryome.rel.filt.s, sample_sums(eukaryome.rel.filt.s)!="0")
      eukaryome.rel.filt.s = filter_taxa(eukaryome.rel.filt.s, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.s))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.s)) #take metadata
      
      dim(df_wilcox)
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$calprotectinelevel)$p.value)
      Speciesnames<-colnames(df_wilcox)
      p.res = data.frame(Speciesnames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.s))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"Feces.calprotectinelevel.SpeciesWilcoxMicrobiomeInsights.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Species level: Step 2: Calprotectine continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$calprotectineggdeps), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$calprotectineggdeps), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"Feces.calprotectinelevelSpeciesCorrelationMicrobiomeInsights.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
 
  #### Species level: Step 1: anemie Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. Do not filter the species ####
      eukaryome<-subset_samples(df_clean, SampleType=="feces")
      eukaryome.s <- tax_glom(eukaryome, "Species")
      eukaryome.filt.s = filter_taxa(eukaryome.s, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.s<-microbiome::transform(eukaryome.filt.s, "compositional")
      View(tax_table(eukaryome.rel.filt.s)) 
      eukaryome.rel.filt.s=subset_samples(eukaryome.rel.filt.s, sample_sums(eukaryome.rel.filt.s)!="0")
      eukaryome.rel.filt.s = filter_taxa(eukaryome.rel.filt.s, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.s))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.s)) #take metadata
      
      dim(df_wilcox)
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$anemie2)$p.value)
      Speciesnames<-colnames(df_wilcox)
      p.res = data.frame(Speciesnames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.s))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"Feces.anemie2.SpeciesWilcoxMicrobiomeInsights.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Species level: Step 2: anemie continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$hemoglobine2), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$hemoglobine2), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"Feces.anemie2SpeciesCorrelationMicrobiomeInsights.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
           
  #### Species level: Step 1: alphaantitrypsinlevel Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. Do not filter the species ####
      eukaryome<-subset_samples(df_clean, SampleType=="feces")
      eukaryome.s <- tax_glom(eukaryome, "Species")
      eukaryome.filt.s = filter_taxa(eukaryome.s, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.s<-microbiome::transform(eukaryome.filt.s, "compositional")
      View(tax_table(eukaryome.rel.filt.s)) 
      eukaryome.rel.filt.s=subset_samples(eukaryome.rel.filt.s, sample_sums(eukaryome.rel.filt.s)!="0")
      eukaryome.rel.filt.s = filter_taxa(eukaryome.rel.filt.s, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.s))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.s)) #take metadata
      
      dim(df_wilcox)
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$alphaantitrypsinlevel)$p.value)
      Speciesnames<-colnames(df_wilcox)
      p.res = data.frame(Speciesnames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.s))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"Feces.alphaantitrypsinlevel.SpeciesWilcoxMicrobiomeInsights.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Species level: Step 2: alphaantitrypsinlevel continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$aatmggdeps), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$aatmggdeps), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"Feces.alphaantitrypsinlevelSpeciesCorrelationMicrobiomeInsights.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
    
      #### Species level: prepare data for without inflammation ####
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$stunted <- meta_wilcox$stunted
      df$run <- meta_wilcox$run
      df$totalreads<-sampledata2$read_count
      df$alphaantitrypsinlevel <- meta_wilcox$alphaantitrypsinlevel
      
      
      #### Species level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c( "stunted", "run", "age",  "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(alphaantitrypsinlevel ~ value+ run+ age + Country + totalreads + stunted , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.aatSpecieswithoutinflaMicrobiomeInsights.csv")
      
      
      
#### all feces, stunted  on Genus level ####
     
      #### Genus level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
      eukaryome<-subset_samples(df_clean, SampleType=="feces")
      eukaryome.g <- tax_glom(eukaryome, "Genus")
      eukaryome.filt.g = filter_taxa(eukaryome.g, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.g<-microbiome::transform(eukaryome.filt.g, "compositional")
      View(tax_table(eukaryome.rel.filt.g)) 
      eukaryome.rel.filt.g=subset_samples(eukaryome.rel.filt.g, sample_sums(eukaryome.rel.filt.g)!="0")
      eukaryome.rel.filt.g = filter_taxa(eukaryome.rel.filt.g, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.g))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.g)) #take metadata
      
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
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.g))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"Feces.stunted.GenusWilcoxMicrobiomeInsights.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Genus level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"Feces.stuntedGenusCorrelationMicrobiomeInsights.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Genus level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Genus level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$stunted <- meta_wilcox$stunted
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel.filt.g))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      #### Genus level: make a loop for logistic regressions with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ age + run +  Country + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.stuntedGenuswitinflaMicrobiomeInsights.csv")
      
      #### Genus level: make a loop for logistic regressions for country with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, Country!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(Country ~ value+ run+ age +  stunted + totalreads + calpro, .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.paysGenuswithoutinflaMicrobiomeInsights.csv")
      
      #### Genus level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age",  "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run+ age + Country + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.stuntedGenuswithoutinflaMicrobiomeInsights.csv")
  
  
 #### Are given Order associated with the nutritional status/stunting? logistic regression in loop first only on feces and on dataset filtered for 0.01% rel abundance ####
      ### all feces, stunted  on Order level
      #### Order level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
      df_clean_order<-tax_glom(df_clean, "Order")
      df_clean_order_filt = filter_taxa(df_clean_order, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel<-microbiome::transform(df_clean_order_filt, "compositional")
      eukaryome.rel<-subset_samples(eukaryome.rel, SampleType=="feces")
      eukaryome.rel.filt=subset_samples(eukaryome.rel, sample_sums(eukaryome.rel)!="0")
      eukaryome.rel.filt = filter_taxa(eukaryome.rel, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt)) #take metadata
      
      dim(df_wilcox)
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$stunted)$p.value)
      Ordernames<-colnames(df_wilcox)
      p.res = data.frame(Ordernames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"Feces.stunted.OrderWilcoxMicrobiomeInsights.csv")
      View(p.res) ## signficiant: Hypocreales
      
      # now for country of origin
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$pays)$p.value)
      Ordernames<-colnames(df_wilcox)
      p.res = data.frame(Ordernames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"Feces.pays.OrderWilcoxMicrobiomeInsights.csv")
      View(p.res) ## signficiant: none after multiple correction
      
      # now for calpro levels
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$calprotectinelevel)$p.value)
      Ordernames<-colnames(df_wilcox)
      p.res = data.frame(Ordernames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"Feces.calpro.OrderWilcoxMicrobiomeInsights.csv")
      View(p.res) ## signficiant: none after multiple correction
      
      
      #### Order level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      meta_wilcox$haz_cont<-as.numeric(meta_wilcox$haz_cont)
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
      write.csv(corres,"Feces.stuntedCorrelationMicrobiomeInsightsOrder.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Order level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Order level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$stunted <- meta_wilcox$stunted
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      View(df) ## numbers look ok now!
      
      #### Order level: make a loop for logistic regressions for stunting with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run+ totalreads + age +  Country + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.stuntedOrderwithinflaMicrobiomeInsights.csv")
      
      #### Order level: make a loop for logistic regressions for country with and without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, Country!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(Country ~ value + run + totalreads + age + stunted , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.paysOrderwithoutinflaMicrobiomeInsights.csv")
      
      # now with inflammation (use calpro)
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(Country ~ value + run + totalreads + age + stunted + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.paysOrderwithinflaMicrobiomeInsights.csv")
      
      
      
      #### Order level: make a loop for logistic regressions for calpro levels  ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age",  "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, calpro!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ run+ totalreads + age +  Country + stunted , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.calproOrderMicrobiomeInsights.csv")
      
      
  #### Species level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. Do not filter the species ####
      eukaryome<-subset_samples(df_clean, SampleType=="feces")
      eukaryome.s <- tax_glom(eukaryome, "Species")
      eukaryome.s <- subset_samples(eukaryome.s, pays=="Madagascar")
      eukaryome.filt.s = filter_taxa(eukaryome.s, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.s<-microbiome::transform(eukaryome.filt.s, "compositional")
      View(tax_table(eukaryome.rel.filt.s)) #We have 45 taxa 
      eukaryome.rel.filt.s=subset_samples(eukaryome.rel.filt.s, sample_sums(eukaryome.rel.filt.s)!="0")
      eukaryome.rel.filt.s = filter_taxa(eukaryome.rel.filt.s, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.s))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.s)) #take metadata
      
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
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.s))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"Feces.stunted.SpeciesWilcoxMicrobiomeInsightsMada.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Species level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"Feces.stuntedSpeciesCorrelationMicrobiomeInsightsMada.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Species level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Species level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$stunted <- meta_wilcox$stunted
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel.filt.s))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      
      #### Species level: make a loop for logistic regressions for stunting with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ age + run  + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.stuntedSpecieswithinflaMicrobiomeInsightsMada.csv")
      
      #### Species level: make a loop for logistic regressions for calpro ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, calpro!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ run+ age  + totalreads + stunted , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.calproSpecieswithoutinflaMicrobiomeInsightsMada.csv")
      
      
      #### Species level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run+ age  + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.stuntedSpecieswithoutinflaMicrobiomeInsightsMada.csv")
      

 #### all feces, stunted  on Genus level ####
      #### Genus level: make a loop for logistic regressions with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ age + run   + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.stuntedGenuswitinflaMicrobiomeInsightsMada.csv")
      
      #### Genus level: make a loop for logistic regressions with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age",  "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, calpro!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ age + run   + totalreads + stunted , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.calproGenuswitinflaMicrobiomeInsightsMada.csv")
      
      #### Genus level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run+ age   + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.stuntedGenuswithoutinflaMicrobiomeInsightsMada.csv")
      
      
      
      
 #### CAR feces, stunted  on ASV level ####
      #### Species level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. Do not filter the species ####
      eukaryome<-subset_samples(df_clean, SampleType=="feces")
      eukaryome.s <- tax_glom(eukaryome, "Species")
      eukaryome.s <- subset_samples(eukaryome.s, pays=="RCA")
      eukaryome.filt.s = filter_taxa(eukaryome.s, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.s<-microbiome::transform(eukaryome.filt.s, "compositional")
      View(tax_table(eukaryome.rel.filt.s)) #We have 46 taxa (78 if it would be unfiltered)
      eukaryome.rel.filt.s=subset_samples(eukaryome.rel.filt.s, sample_sums(eukaryome.rel.filt.s)!="0")
      eukaryome.rel.filt.s = filter_taxa(eukaryome.rel.filt.s, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.s))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.s)) #take metadata
      
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
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.s))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"Feces.stunted.SpeciesWilcoxMicrobiomeInsightsCAR.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Species level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"Feces.stuntedSpeciesCorrelationMicrobiomeInsightsCAR.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Species level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Species level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$stunted <- meta_wilcox$stunted
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel.filt.s))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      
      #### Species level: make a loop for logistic regressions for stunting with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age",  "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ age + run  + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.stuntedSpecieswithinflaMicrobiomeInsightsCAR.csv")
      
      #### Species level: make a loop for logistic regressions for calpro ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, calpro!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ run+ age  + totalreads + stunted , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.calproSpecieswithoutinflaMicrobiomeInsightsCAR.csv")
      
      
      #### Species level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age",  "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run+ age  + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"LogresultsFeces.stuntedSpecieswithoutinflaMicrobiomeInsightsCAR.csv")
      
      
      
      
 #### all feces, stunted  on Genus level ####
      #### Genus level: make a loop for logistic regressions with inflammation ####
      eukaryome<-subset_samples(df_clean, SampleType=="feces")
      eukaryome.g <- tax_glom(eukaryome, "Species")
      eukaryome.g <- subset_samples(eukaryome.g, pays=="RCA")
      eukaryome.filt.g = filter_taxa(eukaryome.g, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.g<-microbiome::transform(eukaryome.filt.g, "compositional")
      View(tax_table(eukaryome.rel.filt.g)) 
      eukaryome.rel.filt.g=subset_samples(eukaryome.rel.filt.g, sample_sums(eukaryome.rel.filt.g)!="0")
      eukaryome.rel.filt.g = filter_taxa(eukaryome.rel.filt.g, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.g))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.g)) #take metadata
      
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$stunted <- meta_wilcox$stunted
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ age + run   + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.stuntedGenuswitinflaMicrobiomeInsightsCAR.csv")
      
      #### Genus level: make a loop for logistic regressions with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, calpro!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ age + run   + totalreads + stunted , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.calproGenuswitinflaMicrobiomeInsightsCAR.csv")
      
      #### Genus level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age",  "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run+ age   + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.stuntedGenuswithoutinflaMicrobiomeInsightsCAR.csv")

      
#### CAR duodenal, stunted  on ASV level ####
      #### ASV level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
      df_clean.M<-subset_samples(df_clean, pays=="RCA")
      df_clean.M_filt = filter_taxa(df_clean.M, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 71 taxa!
      eukaryome.rel<-microbiome::transform(df_clean.M_filt, "compositional")
      eukaryome.rel<-subset_samples(eukaryome.rel, SampleType=="duodenal")
      eukaryome.rel.filt=subset_samples(eukaryome.rel, sample_sums(eukaryome.rel)!="0")
      eukaryome.rel.filt = filter_taxa(eukaryome.rel, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt)) #take metadata
      
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
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"duodenal.stunted.ASVWilcoxMicrobiomeInsightsCAR.csv")
      View(p.res) ## signficiant: none after multiple correction
      
      #### ASV level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      meta_wilcox$haz_cont<-as.numeric(meta_wilcox$haz_cont)
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
      write.csv(corres,"duodenal.stuntedCorrelationMicrobiomeInsightsCAR.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### ASV level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### ASV level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$stunted <- meta_wilcox$stunted
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      View(df) ## numbers look ok now!
      
      #### ASV level: make a loop for logistic regressions for stunting without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age",  "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run+ totalreads + age  , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.stuntedASVwithoutinflaMicrobiomeInsightsCAR.csv")
      
      #### ASV level: make a loop for logistic regressions for stunting with inflammation ####
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run+ totalreads + age  + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.stuntedASVwithinflaMicrobiomeInsightsCAR.csv")
      
      
      #### ASV level: make a loop for logistic regressions for calpro  ####
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ run+ totalreads + age  + stunted , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.calproASVMicrobiomeInsightsCAR.csv")
      
      
      
      #### Species level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. Do not filter the species ####
      eukaryome<-subset_samples(df_clean, SampleType=="duodenal")
      eukaryome.s <- tax_glom(eukaryome, "Species")
      eukaryome.s <- subset_samples(eukaryome.s, pays=="RCA")
      eukaryome.filt.s = filter_taxa(eukaryome.s, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.s<-microbiome::transform(eukaryome.filt.s, "compositional")
      View(tax_table(eukaryome.rel.filt.s)) #We have 11 taxa 
      eukaryome.rel.filt.s=subset_samples(eukaryome.rel.filt.s, sample_sums(eukaryome.rel.filt.s)!="0")
      eukaryome.rel.filt.s = filter_taxa(eukaryome.rel.filt.s, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.s))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.s)) #take metadata
      
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
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.s))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"duodenal.stunted.SpeciesWilcoxMicrobiomeInsightsCAR.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Species level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"duodenal.stuntedSpeciesCorrelationMicrobiomeInsightsCAR.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Species level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Species level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$stunted <- meta_wilcox$stunted
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel.filt.s))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      
      #### Species level: make a loop for logistic regressions for stunting with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age",  "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ age + run  + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.stuntedSpecieswithinflaMicrobiomeInsightsCAR.csv")
      
      #### Species level: make a loop for logistic regressions for calpro ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, calpro!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ run+ age   + totalreads + stunted , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.calproSpecieswithoutinflaMicrobiomeInsightsCAR.csv")
      
      
      #### Species level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run+ age   + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.stuntedSpecieswithoutinflaMicrobiomeInsightsCAR.csv")
      
      #### Genus level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
      eukaryome<-subset_samples(df_clean, SampleType=="duodenal")
      eukaryome.g <- tax_glom(eukaryome, "Genus")
      eukaryome.g<-subset_samples(eukaryome.g, pays=="RCA")
      eukaryome.filt.g = filter_taxa(eukaryome.g, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.g<-microbiome::transform(eukaryome.filt.g, "compositional")
      View(tax_table(eukaryome.rel.filt.g)) #We have 35 taxa (59 if it would be unfiltered)
      eukaryome.rel.filt.g=subset_samples(eukaryome.rel.filt.g, sample_sums(eukaryome.rel.filt.g)!="0")
      eukaryome.rel.filt.g = filter_taxa(eukaryome.rel.filt.g, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.g))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.g)) #take metadata
      
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
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.g))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"duodenal.stunted.GenusWilcoxMicrobiomeInsightsCAR.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Genus level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"duodenal.stuntedGenusCorrelationMicrobiomeInsightsCAR.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Genus level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Genus level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$stunted <- meta_wilcox$stunted
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel.filt.g))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      #### Genus level: make a loop for logistic regressions with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ age + run   + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.stuntedGenuswitinflaMicrobiomeInsightsCAR.csv")
      
      #### Genus level: make a loop for logistic regressions with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, calpro!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ age + run   + totalreads + stunted , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.calproGenuswitinflaMicrobiomeInsightsCAR.csv")
      
      #### Genus level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age",  "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run+ age   + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.stuntedGenuswithoutinflaMicrobiomeInsightsCAR.csv")
      
      
      
 #### Are given ASV associated with the nutritional status/stunting?  ####
 ### all duodenal, stunted  on ASV level ####
      ####  ASV level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
      df_clean_filt = filter_taxa(df_clean, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 3 seqs in at least 5% of samples, we keep 67 taxa
      eukaryome.rel<-microbiome::transform(df_clean_filt, "compositional")
      eukaryome.rel<-subset_samples(eukaryome.rel, SampleType=="duodenal")
      eukaryome.rel<-subset_samples(eukaryome.rel, stunted!="")
      eukaryome.rel.filt=subset_samples(eukaryome.rel, sample_sums(eukaryome.rel)!="0")
      eukaryome.rel.filt = filter_taxa(eukaryome.rel.filt, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt)) #take metadata
      
      dim((df_wilcox))
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$haz)$p.value)
      ASVnames<-colnames(df_wilcox)
      p.res = data.frame(ASVnames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"duodenal.stunted.ASVWilcoxMicrobiomeInsights.csv")
      View(p.res) ## signficiant: none after multiple correction
      
      #### ASV level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      meta_wilcox$haz_cont<-as.numeric(meta_wilcox$haz_cont)
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
      write.csv(corres,"duodenal.stuntedCorrelationMicrobiomeInsights.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### ASV level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### ASV level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$sexe <- meta_wilcox$sexe
      df$Country <- meta_wilcox$pays
      df$haz <- meta_wilcox$haz
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      View(df) ## numbers look ok now!
      
      #### ASV level: make a loop for logistic regressions for stunting without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run+ totalreads + age  + Country , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazASVwithoutinflaMicrobiomeInsights.csv")
      
      #### ASV level: make a loop for logistic regressions on haz without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "age", "run",  "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run + age  + Country , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazASVwithoutinflaMicrobiomeInsights.csv")
      
      
 ### all duodenal, haz  on Species level ####
      
      #### Species level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. Do not filter the species ####
      eukaryome<-subset_samples(df_clean, SampleType=="duodenal")
      eukaryome.s <- tax_glom(eukaryome, "Species")
      eukaryome.filt.s = filter_taxa(eukaryome.s, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.s<-microbiome::transform(eukaryome.filt.s, "compositional")
      View(tax_table(eukaryome.rel.filt.s)) #We have 40 taxa (78 if it would be unfiltered)
      eukaryome.rel.filt.s=subset_samples(eukaryome.rel.filt.s, sample_sums(eukaryome.rel.filt.s)!="0")
      eukaryome.rel.filt.s = filter_taxa(eukaryome.rel.filt.s, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.s))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.s)) #take metadata
      
      dim(df_wilcox)
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$haz)$p.value)
      Speciesnames<-colnames(df_wilcox)
      p.res = data.frame(Speciesnames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.s))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"duodenal.haz.SpeciesWilcoxMicrobiomeInsights.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Species level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
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
      write.csv(corres,"duodenal.hazSpeciesCorrelationMicrobiomeInsights.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Species level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Species level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$haz <- meta_wilcox$haz
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel.filt.s))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      
      #### Species level: make a loop for logistic regressions for haz with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age", "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ age + run + Country + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazSpecieswithinflaMicrobiomeInsights.csv")
      
      #### Species level: make a loop for logistic regressions for country with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "age", "run",  "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, Country!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(Country ~ value+ run+ age  + haz + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.paysSpecieswithinflaMicrobiomeInsights.csv")
      
      #### Species level: make a loop for logistic regressions for haz without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run+ age  + Country + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazSpecieswithoutinflaMicrobiomeInsights.csv")
      
      
      
      
 #### all duodenal, haz  on Genus level ####
      
      #### Genus level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
      eukaryome<-subset_samples(df_clean, SampleType=="duodenal")
      eukaryome.g <- tax_glom(eukaryome, "Genus")
      eukaryome.filt.g = filter_taxa(eukaryome.g, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.g<-microbiome::transform(eukaryome.filt.g, "compositional")
      View(tax_table(eukaryome.rel.filt.g)) #We have 30 taxa (59 if it would be unfiltered)
      eukaryome.rel.filt.g=subset_samples(eukaryome.rel.filt.g, sample_sums(eukaryome.rel.filt.g)!="0")
      eukaryome.rel.filt.g = filter_taxa(eukaryome.rel.filt.g, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.g))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.g)) #take metadata
      
      dim(df_wilcox)
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$haz)$p.value)
      Genusnames<-colnames(df_wilcox)
      p.res = data.frame(Genusnames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.g))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"duodenal.haz.GenusWilcoxMicrobiomeInsights.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Genus level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"duodenal.hazGenusCorrelationMicrobiomeInsights.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Genus level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Genus level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$haz <- meta_wilcox$haz
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel.filt.g))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      #### Genus level: make a loop for logistic regressions for haz with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ age + run  + Country + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazGenuswitinflaMicrobiomeInsights.csv")
      
      #### Genus level: make a loop for logistic regressions for country with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, Country!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(Country ~ value+ run+ age  + haz + totalreads + calpro, .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.paysGenuswithoutinflaMicrobiomeInsights.csv")
      
      #### Genus level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run+ age  + Country + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazGenuswithoutinflaMicrobiomeInsights.csv")
      
      
 #### Mada duodenal, haz  on ASV level ####
      #### ASV level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
      df_clean.M<-subset_samples(df_clean, pays=="Madagascar")
      df_clean.M_filt = filter_taxa(df_clean.M, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 71 taxa!
      eukaryome.rel<-microbiome::transform(df_clean.M_filt, "compositional")
      eukaryome.rel<-subset_samples(eukaryome.rel, SampleType=="duodenal")
      eukaryome.rel.filt=subset_samples(eukaryome.rel, sample_sums(eukaryome.rel)!="0")
      eukaryome.rel.filt = filter_taxa(eukaryome.rel, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt)) #take metadata
      
      dim(df_wilcox)
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$haz)$p.value)
      ASVnames<-colnames(df_wilcox)
      p.res = data.frame(ASVnames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"duodenal.haz.ASVWilcoxMicrobiomeInsightsMada.csv")
      View(p.res) ## signficiant: none after multiple correction
      
      #### ASV level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      meta_wilcox$haz_cont<-as.numeric(meta_wilcox$haz_cont)
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
      write.csv(corres,"duodenal.hazCorrelationMicrobiomeInsightsMada.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### ASV level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### ASV level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$haz <- meta_wilcox$haz
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      View(df) ## numbers look ok now!
      
      #### ASV level: make a loop for logistic regressions for stunting without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run+ totalreads + age  , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazASVwithoutinflaMicrobiomeInsightsMada.csv")
      
      #### ASV level: make a loop for logistic regressions for stunting with inflammation ####
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run+ totalreads + age  + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazASVwithinflaMicrobiomeInsightsMada.csv")
      
      
      #### ASV level: make a loop for logistic regressions for calpro  ####
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ run+ totalreads + age  + haz , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.calproASVMicrobiomeInsightsMada.csv")
      
      
      
 #### Species level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. Do not filter the species ####
      eukaryome<-subset_samples(df_clean, SampleType=="duodenal")
      eukaryome.s <- tax_glom(eukaryome, "Species")
      eukaryome.s <- subset_samples(eukaryome.s, pays=="Madagascar")
      eukaryome.filt.s = filter_taxa(eukaryome.s, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.s<-microbiome::transform(eukaryome.filt.s, "compositional")
      View(tax_table(eukaryome.rel.filt.s)) #We have 14 taxa 
      eukaryome.rel.filt.s=subset_samples(eukaryome.rel.filt.s, sample_sums(eukaryome.rel.filt.s)!="0")
      eukaryome.rel.filt.s = filter_taxa(eukaryome.rel.filt.s, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.s))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.s)) #take metadata
      
      dim(df_wilcox)
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$haz)$p.value)
      Speciesnames<-colnames(df_wilcox)
      p.res = data.frame(Speciesnames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.s))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"duodenal.haz.SpeciesWilcoxMicrobiomeInsightsMada.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Species level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"duodenal.hazSpeciesCorrelationMicrobiomeInsightsMada.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Species level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Species level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$haz <- meta_wilcox$haz
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel.filt.s))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      
      #### Species level: make a loop for logistic regressions for stunting with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ age + run  + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazSpecieswithinflaMicrobiomeInsightsMada.csv")
      
      #### Species level: make a loop for logistic regressions for calpro ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, calpro!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ run+ age   + totalreads + haz , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.calproSpecieswithoutinflaMicrobiomeInsightsMada.csv")
      
      
      #### Species level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run+ age   + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazSpecieswithoutinflaMicrobiomeInsightsMada.csv")
      
      
      
      
      ### all duodenal, haz  on Genus level
      
 #### Genus level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
      eukaryome<-subset_samples(df_clean, SampleType=="duodenal")
      eukaryome.g <- tax_glom(eukaryome, "Genus")
      eukaryome.g<-subset_samples(eukaryome.g, pays=="Madagascar")
      eukaryome.filt.g = filter_taxa(eukaryome.g, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.g<-microbiome::transform(eukaryome.filt.g, "compositional")
      View(tax_table(eukaryome.rel.filt.g)) #We have 17 taxa 
      eukaryome.rel.filt.g=subset_samples(eukaryome.rel.filt.g, sample_sums(eukaryome.rel.filt.g)!="0")
      eukaryome.rel.filt.g = filter_taxa(eukaryome.rel.filt.g, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.g))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.g)) #take metadata
      
      dim(df_wilcox)
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$haz)$p.value)
      Genusnames<-colnames(df_wilcox)
      p.res = data.frame(Genusnames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.g))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"duodenal.haz.GenusWilcoxMicrobiomeInsightsMada.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Genus level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"duodenal.hazGenusCorrelationMicrobiomeInsightsMada.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Genus level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Genus level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$haz <- meta_wilcox$haz
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel.filt.g))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      #### Genus level: make a loop for logistic regressions with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age", "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ age + run   + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazGenuswitinflaMicrobiomeInsightsMada.csv")
      
      #### Genus level: make a loop for logistic regressions with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, calpro!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ age + run   + totalreads + haz , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.calproGenuswitinflaMicrobiomeInsightsMada.csv")
      
      #### Genus level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run+ age   + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazGenuswithoutinflaMicrobiomeInsightsMada.csv")
      
      
      
      
 #### CAR duodenal, haz  on ASV level ####
   #### ASV level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
      df_clean.M<-subset_samples(df_clean, pays=="RCA")
      df_clean.M_filt = filter_taxa(df_clean.M, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 71 taxa!
      eukaryome.rel<-microbiome::transform(df_clean.M_filt, "compositional")
      eukaryome.rel<-subset_samples(eukaryome.rel, SampleType=="duodenal")
      eukaryome.rel.filt=subset_samples(eukaryome.rel, sample_sums(eukaryome.rel)!="0")
      eukaryome.rel.filt = filter_taxa(eukaryome.rel, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt)) #take metadata
      
      dim(df_wilcox)
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$haz)$p.value)
      ASVnames<-colnames(df_wilcox)
      p.res = data.frame(ASVnames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"duodenal.haz.ASVWilcoxMicrobiomeInsightsCAR.csv")
      View(p.res) ## signficiant: none after multiple correction
      
      #### ASV level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      meta_wilcox$haz_cont<-as.numeric(meta_wilcox$haz_cont)
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
      write.csv(corres,"duodenal.hazCorrelationMicrobiomeInsightsCAR.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### ASV level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### ASV level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$haz <- meta_wilcox$haz
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      View(df) ## numbers look ok now!
      
      #### ASV level: make a loop for logistic regressions for stunting without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age", "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run+ totalreads + age  , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazASVwithoutinflaMicrobiomeInsightsCAR.csv")
      
      #### ASV level: make a loop for logistic regressions for stunting with inflammation ####
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run+ totalreads + age  + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazASVwithinflaMicrobiomeInsightsCAR.csv")
      
      
      #### ASV level: make a loop for logistic regressions for calpro  ####
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ run+ totalreads + age  + haz , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.calproASVMicrobiomeInsightsCAR.csv")
      
      
      
   #### Species level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. Do not filter the species ####
      eukaryome<-subset_samples(df_clean, SampleType=="duodenal")
      eukaryome.s <- tax_glom(eukaryome, "Species")
      eukaryome.s <- subset_samples(eukaryome.s, pays=="RCA")
      eukaryome.filt.s = filter_taxa(eukaryome.s, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.s<-microbiome::transform(eukaryome.filt.s, "compositional")
      View(tax_table(eukaryome.rel.filt.s)) #We have 46 taxa (78 if it would be unfiltered)
      eukaryome.rel.filt.s=subset_samples(eukaryome.rel.filt.s, sample_sums(eukaryome.rel.filt.s)!="0")
      eukaryome.rel.filt.s = filter_taxa(eukaryome.rel.filt.s, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.s))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.s)) #take metadata
      
      dim(df_wilcox)
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$haz)$p.value)
      Speciesnames<-colnames(df_wilcox)
      p.res = data.frame(Speciesnames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.s))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"duodenal.haz.SpeciesWilcoxMicrobiomeInsightsCAR.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Species level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"duodenal.hazSpeciesCorrelationMicrobiomeInsightsCAR.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Species level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Species level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$haz <- meta_wilcox$haz
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel.filt.s))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      
      #### Species level: make a loop for logistic regressions for stunting with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age", "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ age + run  + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazSpecieswithinflaMicrobiomeInsightsCAR.csv")
      
      #### Species level: make a loop for logistic regressions for calpro ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, calpro!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ run+ age   + totalreads + haz , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.calproSpecieswithoutinflaMicrobiomeInsightsCAR.csv")
      
      
      #### Species level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run+ age   + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazSpecieswithoutinflaMicrobiomeInsightsCAR.csv")
      
      
      
      
      
      
      
      #### Species level: make a loop for logistic regressions CAR only ####
      library(broom)
      library(dplyr)
      df<-filter(df, Country=="RCA")
      long = melt(df, id.vars = c("haz", "run", "age",  "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run+ age + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsfeces.hazSpecieswithoutinflaMicrobiomeInsightsCAR.csv")
      
      
      
      
      
      
      
   #### Genus level: Step 1: Pick Taxa on Individual Relative Signifiance with no correction for cofactors using Wilcox Test. ####
      eukaryome<-subset_samples(df_clean, SampleType=="duodenal")
      eukaryome.g <- tax_glom(eukaryome, "Genus")
      eukaryome.g<-subset_samples(eukaryome.g, pays=="RCA")
      eukaryome.filt.g = filter_taxa(eukaryome.g, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # only keep taxa with at least 10 seqs in at least 5% of samples; we keep 64 taxa!
      eukaryome.rel.filt.g<-microbiome::transform(eukaryome.filt.g, "compositional")
      View(tax_table(eukaryome.rel.filt.g)) #We have 35 taxa (59 if it would be unfiltered)
      eukaryome.rel.filt.g=subset_samples(eukaryome.rel.filt.g, sample_sums(eukaryome.rel.filt.g)!="0")
      eukaryome.rel.filt.g = filter_taxa(eukaryome.rel.filt.g, function(x) sum(x) > 0, TRUE)
      
      df_wilcox <- as.matrix(t(otu_table(eukaryome.rel.filt.g))) #take rel abund and Wilcoxon rank-sum
      meta_wilcox <- as.data.frame(sample_data(eukaryome.rel.filt.g)) #take metadata
      
      dim(df_wilcox)
      dim(meta_wilcox)
      
      MW.p = apply(df_wilcox,2,
                   function(x) wilcox.test(c(x)~meta_wilcox$haz)$p.value)
      Genusnames<-colnames(df_wilcox)
      p.res = data.frame(Genusnames,MW.p)
      # Perform multiple comparison correction using a given method of choice
      p.res$rel.fdr <- p.adjust(p.res$MW.p, method="fdr")
      p.res$bonferroni <- p.adjust(p.res$MW.p, method="bonferroni") 
      # Merge with tax info
      Tax_corr<-as.data.frame(tax_table(eukaryome.rel.filt.g))
      p.res = merge(p.res,Tax_corr, by="row.names")
      #export results
      write.csv(p.res,"duodenal.haz.GenusWilcoxMicrobiomeInsightsCAR.csv")
      View(p.res) ## signficiant: none after multiple testing!
      
      #### Genus level: Step 2: continuous variable, ex. haz or  age multiple correlation tests, careful: you need to subset first if you have missing variables ####
      mcorr <- apply(df_wilcox, 2,
                     function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$p.value)
      estcorr <- apply(df_wilcox,2,
                       function(x) cor.test(c(x), as.numeric(meta_wilcox$haz_cont), method="spearman")$estimate)
      corres = data.frame(colnames(df_wilcox),mcorr,estcorr)
      # Perform multiple comparison correction using a given method of choice
      corres$rel.fdr <- p.adjust(corres$mcorr, method="fdr")
      corres$bonferroni <- p.adjust(corres$mcorr, method="bonferroni")
      # Merge with tax info
      corres = merge(corres,Tax_corr, by="row.names")
      
      #export results
      write.csv(corres,"duodenal.hazGenusCorrelationMicrobiomeInsightsCAR.csv")
      View(corres) ## several are associated, but they do not survive multiple testing!
      
      #### Genus level: Step 3: LOGISTIC MODELS CORRECTING FOR COVARIABLES#####
      #### Genus level: prepare your data ####
      dim(df_wilcox)
      dim(meta_wilcox)
      data1 <- df_wilcox
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_wilcox$ageyears
      df$Country <- meta_wilcox$pays
      df$haz <- meta_wilcox$haz
      df$calpro <- meta_wilcox$calprotectinelevel
      df$aat <- meta_wilcox$alphaantitrypsinlevel
      df$run <- meta_wilcox$run
      
      sampledata= as.data.frame(sample_data(eukaryome.rel.filt.g))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      #### Genus level: make a loop for logistic regressions with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ age + run   + totalreads + calpro , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazGenuswitinflaMicrobiomeInsightsCAR.csv")
      
      #### Genus level: make a loop for logistic regressions with inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age",  "calpro", "aat", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, calpro!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(calpro ~ value+ age + run   + totalreads + haz , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.calproGenuswitinflaMicrobiomeInsightsCAR.csv")
      
      #### Genus level: make a loop for logistic regressions without inflammation ####
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("haz", "run", "age", "Country", "totalreads")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, haz!="") ## keep only the ones with valid data for haz
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(haz ~ value+ run+ age   + totalreads , .,  family=binomial))) %>% 
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
      write.csv(logresults,"Logresultsduodenal.hazGenuswithoutinflaMicrobiomeInsightsCAR.csv")
    
      
#### Ascomycota to Basidiomycota levels ####
      Asco<-subset_taxa(dfclean, Phylum=="p__Ascomycota")
      Basidio<-subset_taxa(dfclean, Phylum=="p__Basidiomycota")
      AscoBasidio<-sample_sums(Asco)/sample_sums(Basidio)
      sample_data(dfclean)$AscotoBasidio<-AscoBasidio
     
      boxplot(sample_data(dfclean)$AscotoBasidio~sample_data(dfclean)$SampleType)
      boxplot(log10~sample_data(dfclean)$SampleType)
      kruskal.test(sample_data(dfclean)$AscotoBasidio~sample_data(dfclean)$SampleType) ## p-value =2.835e-11
      
      Asco<-subset_taxa(dfcleanfeces, Phylum=="p__Ascomycota")
      Basidio<-subset_taxa(dfcleanfeces, Phylum=="p__Basidiomycota")
      AscoBasidio<-sample_sums(Asco)/sample_sums(Basidio)
      sample_data(dfcleanfeces)$AscotoBasidio<-AscoBasidio
      log10<-transform(sample_data(dfcleanfeces)$AscotoBasidio,  'log10')
      
      boxplot(sample_data(dfcleanfeces)$AscotoBasidio~sample_data(dfcleanfeces)$pays)
      boxplot(log10~sample_data(dfcleanfeces)$pays)
     
      wilcox.test(sample_data(dfcleanfeces)$AscotoBasidio~sample_data(dfcleanfeces)$pays) ##p-value = 0.003
      
      boxplot(sample_data(dfcleanfeces)$AscotoBasidio~sample_data(dfcleanfeces)$stunted)
      wilcox.test(sample_data(dfcleanfeces)$AscotoBasidio~sample_data(dfcleanfeces)$stunted) ##ns
      kruskal.test(sample_data(dfcleanfeces)$AscotoBasidio~sample_data(dfcleanfeces)$haz) ##ns
      
      boxplot(sample_data(dfcleanfeces)$AscotoBasidio~sample_data(dfcleanfeces)$ageyears)
      kruskal.test(sample_data(dfcleanfeces)$AscotoBasidio~sample_data(dfcleanfeces)$ageyears) ##ns
      
      dfcleanfecesB<-subset_samples(dfcleanfeces, pays=="RCA")
      boxplot(sample_data(dfcleanfecesB)$AscotoBasidio~sample_data(dfcleanfecesB)$calprotectinelevel)
      wilcox.test(sample_data(dfcleanfecesB)$AscotoBasidio~sample_data(dfcleanfecesB)$calprotectinelevel) ##ns
      
      boxplot(sample_data(dfcleanfecesB)$AscotoBasidio~sample_data(dfcleanfecesB)$alphaantitrypsinlevel)
      wilcox.test(sample_data(dfcleanfecesB)$AscotoBasidio~sample_data(dfcleanfecesB)$alphaantitrypsinlevel) ##ns
      
      boxplot(sample_data(dfcleanfecesB)$AscotoBasidio~sample_data(dfcleanfecesB)$stunted)
      wilcox.test(sample_data(dfcleanfecesB)$AscotoBasidio~sample_data(dfcleanfecesB)$stunted) ##ns
      
      dfcleanfecesM<-subset_samples(dfcleanfeces, pays=="Madagascar")
      boxplot(sample_data(dfcleanfecesM)$AscotoBasidio~sample_data(dfcleanfecesM)$calprotectinelevel)
      wilcox.test(sample_data(dfcleanfecesM)$AscotoBasidio~sample_data(dfcleanfecesM)$calprotectinelevel) ##ns
      
      boxplot(sample_data(dfcleanfecesM)$AscotoBasidio~sample_data(dfcleanfecesM)$alphaantitrypsinlevel)
      wilcox.test(sample_data(dfcleanfecesM)$AscotoBasidio~sample_data(dfcleanfecesM)$alphaantitrypsinlevel) ##ns
    
      boxplot(sample_data(dfcleanfecesM)$AscotoBasidio~sample_data(dfcleanfecesM)$stunted)
      wilcox.test(sample_data(dfcleanfecesM)$AscotoBasidio~sample_data(dfcleanfecesM)$stunted) ##ns
      
      Asco<-subset_taxa(dfcleanduodenal, Phylum=="p__Ascomycota")
      Basidio<-subset_taxa(dfcleanduodenal, Phylum=="p__Basidiomycota")
      AscoBasidio<-sample_sums(Asco)/sample_sums(Basidio)
      sample_data(dfcleanduodenal)$AscotoBasidio<-AscoBasidio
      
      boxplot(sample_data(dfcleanduodenal)$AscotoBasidio~sample_data(dfcleanduodenal)$pays)
      wilcox.test(sample_data(dfcleanduodenal)$AscotoBasidio~sample_data(dfcleanduodenal)$pays) ##ns
      
      boxplot(sample_data(dfcleanduodenal)$AscotoBasidio~sample_data(dfcleanduodenal)$haz)
      kruskal.test(sample_data(dfcleanduodenal)$AscotoBasidio~sample_data(dfcleanduodenal)$haz) ##ns
      
      boxplot(sample_data(dfcleanduodenal)$AscotoBasidio~sample_data(dfcleanduodenal)$ageyears)
      kruskal.test(sample_data(dfcleanduodenal)$AscotoBasidio~sample_data(dfcleanduodenal)$ageyears) ##ns
      
      dfcleanduodenalB<-subset_samples(dfcleanduodenal, pays=="RCA")
      boxplot(sample_data(dfcleanduodenalB)$AscotoBasidio~sample_data(dfcleanduodenalB)$haz)
      wilcox.test(sample_data(dfcleanduodenalB)$AscotoBasidio~sample_data(dfcleanduodenalB)$haz) ##ns
      boxplot(sample_data(dfcleanduodenalB)$AscotoBasidio~sample_data(dfcleanduodenalB)$anemie2)
      wilcox.test(sample_data(dfcleanduodenalB)$AscotoBasidio~sample_data(dfcleanduodenalB)$anemie2) ##ns
      
      
      dfcleanduodenalM<-subset_samples(dfcleanduodenal, pays=="Madagascar")
      boxplot(sample_data(dfcleanduodenalM)$AscotoBasidio~sample_data(dfcleanduodenalM)$haz)
      wilcox.test(sample_data(dfcleanduodenalM)$AscotoBasidio~sample_data(dfcleanduodenalM)$haz) ##ns
      boxplot(sample_data(dfcleanduodenalM)$AscotoBasidio~sample_data(dfcleanduodenalM)$anemie2)
      wilcox.test(sample_data(dfcleanduodenalM)$AscotoBasidio~sample_data(dfcleanduodenalM)$anemie2) ##ns
      
      Asco<-subset_taxa(dfcleanduodenal, Phylum=="p__Ascomycota")
      Basidio<-subset_taxa(dfcleanduodenal, Phylum=="p__Basidiomycota")
      AscoBasidio<-sample_sums(Asco)/sample_sums(Basidio)
      sample_data(dfcleanduodenal)$AscotoBasidio<-AscoBasidio
      
      boxplot(sample_data(dfcleanduodenal)$AscotoBasidio~sample_data(dfcleanduodenal)$pays)
      wilcox.test(sample_data(dfcleanduodenal)$AscotoBasidio~sample_data(dfcleanduodenal)$pays) ##ns
      
      boxplot(sample_data(dfcleanduodenal)$AscotoBasidio~sample_data(dfcleanduodenal)$haz)
      kruskal.test(sample_data(dfcleanduodenal)$AscotoBasidio~sample_data(dfcleanduodenal)$haz) ##ns
      
      boxplot(sample_data(dfcleanduodenal)$AscotoBasidio~sample_data(dfcleanduodenal)$ageyears)
      kruskal.test(sample_data(dfcleanduodenal)$AscotoBasidio~sample_data(dfcleanduodenal)$ageyears) ##ns
      
      dfcleanduodenalB<-subset_samples(dfcleanduodenal, pays=="RCA")
      boxplot(sample_data(dfcleanduodenalB)$AscotoBasidio~sample_data(dfcleanduodenalB)$haz)
      wilcox.test(sample_data(dfcleanduodenalB)$AscotoBasidio~sample_data(dfcleanduodenalB)$haz) ##ns
      boxplot(sample_data(dfcleanduodenalB)$AscotoBasidio~sample_data(dfcleanduodenalB)$anemie2)
      wilcox.test(sample_data(dfcleanduodenalB)$AscotoBasidio~sample_data(dfcleanduodenalB)$anemie2) ##ns
      
      
      dfcleanduodenalM<-subset_samples(dfcleanduodenal, pays=="Madagascar")
      boxplot(sample_data(dfcleanduodenalM)$AscotoBasidio~sample_data(dfcleanduodenalM)$haz)
      wilcox.test(sample_data(dfcleanduodenalM)$AscotoBasidio~sample_data(dfcleanduodenalM)$haz) ##ns
      boxplot(sample_data(dfcleanduodenalM)$AscotoBasidio~sample_data(dfcleanduodenalM)$anemie2)
      wilcox.test(sample_data(dfcleanduodenalM)$AscotoBasidio~sample_data(dfcleanduodenalM)$anemie2) ##ns
      
      
 #### DeSeq2 Analysis on differential abundance of taxa ####  
      #### on Genus level and SampleType ####
      dfcleanduodfeces<-subset_samples(dfclean, SampleType!="gastric")
      dfcleanduodfeces<-subset_samples(dfcleanduodfeces, stunted=="stunted")
      
      dfcleanduodfeces_g<-tax_glom(dfcleanduodfeces, "Genus")
      test<-filter_taxa(dfcleanduodfeces_g, function(x) sum(x)>=1, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      
      diagddsfecesDM = phyloseq_to_deseq2(test2, ~  read_count + pays + SampleType)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count + pays)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "SampleType_duodenal_vs_feces") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanduodfeces_g)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at rank 7/ Genus level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.SampleType <- dfcleanduodfeces_g %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("SampleType") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.SampleType) <- paste("prevalence", colnames(prev_counts.SampleType), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.SampleType <- dfcleanduodfeces_g %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("SampleType") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.SampleType) <- paste("prevalence", colnames(prev_possible.SampleType), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.SampleType/prev_possible.SampleType)*100
      
      tax_table.SampleType =  as.data.frame(tax_table(dfcleanduodfeces_g))
      
      merge.SampleType.prev= merge(tax_table.SampleType, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanduodfeces_g %>%
        transform_sample_counts(function(x) {
          x/sum(x)} ) 
      
      Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "SampleType")
      
      Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
      
      Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
      
      Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
        otu_table() %>%
        as.data.frame()
      
      colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      
      row.names(merge.SampleType.prev) <-merge.SampleType.prev$Row.names
      
      merge.abunanceprev.feces.OTU= merge(merge.SampleType.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
      merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
      row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
      
      # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
      SampleTypecorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      write.csv(SampleTypecorrcountrygenderage_reduced_prevabdeseq, "GenusSampleTypeDeSeq2ITS.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Genus_fecesvsduod_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Feces vs. duodenal samples, corrected for country of origin, LRT model")
      dev.off()
      
      # make  graph plot
      pdf("differentially_present_GenusbySampleTypeMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 6, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=SampleTypecorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different genera feces vs. duodenal samples")
      dev.off()
      

      #### on Genus level and SampleType in Bangui ####
      dfcleanduodfecesB<-subset_samples(dfcleanduodfeces, pays=="RCA")
      
      dfcleanduodfecesB_g<-tax_glom(dfcleanduodfecesB, "Genus")
      test<-filter_taxa(dfcleanduodfecesB_g, function(x) sum(x)>=1, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      
      diagddsfecesDM = phyloseq_to_deseq2(test2, ~  read_count + SampleType)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "SampleType_duodenal_vs_feces") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanduodfecesB_g)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at rank 7 Genus level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.SampleType <- dfcleanduodfecesB_g %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("SampleType") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.SampleType) <- paste("prevalence", colnames(prev_counts.SampleType), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.SampleType <- dfcleanduodfecesB_g %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("SampleType") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.SampleType) <- paste("prevalence", colnames(prev_possible.SampleType), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.SampleType/prev_possible.SampleType)*100
      
      tax_table.SampleType =  as.data.frame(tax_table(dfcleanduodfecesB_g))
      
      merge.SampleType.prev= merge(tax_table.SampleType, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanduodfecesB_g %>%
        transform_sample_counts(function(x) {
          x/sum(x)} ) 
      
      Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "SampleType")
      
      Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
      
      Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
      
      Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
        otu_table() %>%
        as.data.frame()
      
      colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      
      row.names(merge.SampleType.prev) <-merge.SampleType.prev$Row.names
      
      merge.abunanceprev.feces.OTU= merge(merge.SampleType.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
      merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
      row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
      
      # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
      SampleTypecorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      write.csv(SampleTypecorrcountrygenderage_reduced_prevabdeseq, "GenusSampleTypeDeSeq2BanguiITS.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Genus_fecesvsduod_LRTMicrobiomeInsightsBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Feces vs. duodenal samples, corrected for country of origin, LRT model")
      dev.off()
      
      # make  graph plot
      pdf("differentially_present_GenusbySampleTypeMicrobiomeInsightsBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 6, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=SampleTypecorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different genera feces vs. duodenal samples in Bangui")
      dev.off()      
     
      #### on Genus level and SampleType in Antananarivo ####
      dfcleanduodfecesA<-subset_samples(dfcleanduodfeces, pays=="Madagascar")
      
      dfcleanduodfecesA_g<-tax_glom(dfcleanduodfecesA, "Genus")
      test<-filter_taxa(dfcleanduodfecesA_g, function(x) sum(x)>=1, TRUE)
      test2<-prune_samples(sample_sums(test)>0,test )
      
      diagddsfecesDM = phyloseq_to_deseq2(test2, ~  read_count + SampleType)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "SampleType_duodenal_vs_feces") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanduodfecesA_g)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Genus level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.SampleType <- dfcleanduodfecesA_g %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("SampleType") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.SampleType) <- paste("prevalence", colnames(prev_counts.SampleType), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.SampleType <- dfcleanduodfecesA_g %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("SampleType") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.SampleType) <- paste("prevalence", colnames(prev_possible.SampleType), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.SampleType/prev_possible.SampleType)*100
      
      tax_table.SampleType =  as.data.frame(tax_table(dfcleanduodfecesA_g))
      
      merge.SampleType.prev= merge(tax_table.SampleType, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanduodfecesA_g %>%
        transform_sample_counts(function(x) {
          x/sum(x)} ) 
      
      Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "SampleType")
      
      Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
      
      Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(OTU) OTU *100/sum(OTU))
      
      Afribiota_feces_abundance_OTU_st= Afribiota_feces_abundance_OTU_st %>%
        otu_table() %>%
        as.data.frame()
      
      colnames(Afribiota_feces_abundance_OTU_st) <- paste("abundance", colnames(Afribiota_feces_abundance_OTU_st), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      
      row.names(merge.SampleType.prev) <-merge.SampleType.prev$Row.names
      
      merge.abunanceprev.feces.OTU= merge(merge.SampleType.prev, Afribiota_feces_abundance_OTU_st, by="row.names")
      merge.abunanceprev.feces.OTU<-merge.abunanceprev.feces.OTU[, -1]
      row.names(merge.abunanceprev.feces.OTU)<-merge.abunanceprev.feces.OTU$Row.names
      
      # Merge abundance and prevalence table from above with DeSeq results (no need to make again binding with Tax table)
      SampleTypecorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      write.csv(SampleTypecorrcountrygenderage_reduced_prevabdeseq, "GenusSampleTypeDeSeq2TanaITS.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Genus_fecesvsduod_LRTMicrobiomeInsightsAntananarivo.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Feces vs. duodenal samples, corrected for country of origin, LRT model")
      dev.off()
      
      # make  graph plot
      pdf("differentially_present_GenusbySampleTypeMicrobiomeInsightsAntananarivo.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 6, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=SampleTypecorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different genera feces vs. duodenal samples in Antananarivo")
      dev.off()
      
   
      #### Feces: on ASV level and stunting corrected country ####
      dfcleanfeces<-subset_samples(dfclean, SampleType=="feces")
      sample_data(dfcleanfeces)$stunted<-sample_data(dfcleanfeces)$stunted
      dfcleanfeces<-subset_samples(dfcleanfeces, stunted!="")
      dfcleanfeces<-subset_samples(dfcleanfeces, pays!="")
      dfcleanfeces<-subset_samples(dfcleanfeces, sexe!="")
      dfcleanfeces = filter_taxa(dfcleanfeces, function(x) sum(x) > 0, TRUE)
      dfcleanfeces<-subset_samples(dfcleanfeces, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces, ~  run + pays + age + read_count+ stunted)
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run+ pays + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
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
      
      Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(x) x *100/sum(x))
      
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
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "StuntedASVlevelcorrcountry.csv")
      
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Species = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Species), levels=names(x))
      
      #make the actual plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      
      pdf("differentially_present_taxa_ASV_stunted_corr_country_sexe_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=taxonomy, y=log2FoldChange, color=Order)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunted vs. non-stunted in feces, controlled for gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      pdf("differentially_present_ASVbystuntingstatusMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
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
        ggtitle("Significantly different ASV by stunting status")
      dev.off()
      
      
      #### Feces: on ASV level and country, not correcting for inflammation ####
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces, ~  run + stunted + age + read_count + pays)
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") 
      
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
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "ASVlevelfecescountrynotcorrinflammationDeSeq2MicrobiomeInsights.csv")
      
      #plot
      
      # Order order
      x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$Order, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$Class, function(x) max(x))
      x = sort(x, TRUE)
      payscorrcountrygenderage_reduced_prevabdeseq$Family = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$Class), levels=names(x))
      # Family order
      x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      payscorrcountrygenderage_reduced_prevabdeseq$Order = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$Family), levels=names(x))
      # Genus order
      x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$Genus, function(x) max(x))
      x = sort(x, TRUE)
      payscorrcountrygenderage_reduced_prevabdeseq$Genus = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$Genus), levels=names(x))
      # Species order
      x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$Species, function(x) max(x))
      x = sort(x, TRUE)
      payscorrcountrygenderage_reduced_prevabdeseq$Species = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$Species), levels=names(x))
      
      payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$Family, "/", payscorrcountrygenderage_reduced_prevabdeseq$Genus, "/", payscorrcountrygenderage_reduced_prevabdeseq$Species)
      
      
      #make the actual plot
      
      pdf("differentially_present_taxa_ASV_pays_corr_country_sexe_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(payscorrcountrygenderage_reduced_prevabdeseq, aes(x=taxonomy, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$Order, "/", payscorrcountrygenderage_reduced_prevabdeseq$Genus, "/", payscorrcountrygenderage_reduced_prevabdeseq$Species)
      pdf("differentially_present_ASVbypaysMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
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
      
      
      #### Feces: on ASV level and country, correcting for inflammation ####
      dfcleanfeces<-subset_samples(dfcleanfeces, calprotectinelevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces, ~  run + stunted +  age + read_count + calprotectinelevel + pays)
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted +  age + read_count + calprotectinelevel)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") 
      
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
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "ASVlevelcountrywithcorrinflammationDeSeq2MicrobiomeInsights.csv")
      
      #plot
      
      # Order order
      x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$Order, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$Class, function(x) max(x))
      x = sort(x, TRUE)
      payscorrcountrygenderage_reduced_prevabdeseq$Family = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$Class), levels=names(x))
      # Family order
      x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      payscorrcountrygenderage_reduced_prevabdeseq$Order = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$Family), levels=names(x))
      # Genus order
      x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$Genus, function(x) max(x))
      x = sort(x, TRUE)
      payscorrcountrygenderage_reduced_prevabdeseq$Genus = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$Genus), levels=names(x))
      # Species order
      x = tapply(payscorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, payscorrcountrygenderage_reduced_prevabdeseq$Species, function(x) max(x))
      x = sort(x, TRUE)
      payscorrcountrygenderage_reduced_prevabdeseq$Species = factor(as.character(payscorrcountrygenderage_reduced_prevabdeseq$Species), levels=names(x))
      
      payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$Family, "/", payscorrcountrygenderage_reduced_prevabdeseq$Genus, "/", payscorrcountrygenderage_reduced_prevabdeseq$Species)
      
      
      #make the actual plot
      
      pdf("differentially_present_taxa_ASV_pays_corr_country_sexe_age_calpro_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(payscorrcountrygenderage_reduced_prevabdeseq, aes(x=taxonomy, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$Order, "/", payscorrcountrygenderage_reduced_prevabdeseq$Genus, "/", payscorrcountrygenderage_reduced_prevabdeseq$Species)
      pdf("differentially_present_ASVbypayscorrinflaMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
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
      
      
      
      #### Feces: on Species level and stunting, not correcting for inflammation  ####
      dfcleanfecesfilt_s<-tax_glom(dfcleanfeces, "Species")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, stunted!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, pays!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, ageyears!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, sexe!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, run!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, !is.na(read_count))
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, (rowSums(otu_table(dfcleanfecesfilt_s))!=0))
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, (colSums(otu_table(dfcleanfecesfilt_s))!=0))
      dfcleanfecesfilt_s<-  phyloseq_rm_na_tax(dfcleanfecesfilt_s)
      
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s, ~  run + pays + age + read_count + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + pays + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
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
      
      prev_counts.stunted <- dfcleanfecesfilt_s %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfecesfilt_s %>% #this produces a maximum possible prevalence count per OTU
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
      
      Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(x) x *100/sum(x))
      
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
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "StuntedSpecieslevelcorrcountry.csv")
      
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Species = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Species), levels=names(x))
      
      #make the actual plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      
      pdf("differentially_present_taxa_Species_stunted_corr_country_sexe_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=taxonomy, y=log2FoldChange, color=Order)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunted vs. non-stunted in feces, controlled for gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      pdf("differentially_present_SpeciesbystuntingstatusMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
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
        ggtitle("Significantly different species by stunting status")
      dev.off()
      
     
      #### Feces: on Species level and stunting, correcting for calprotectine  ####
      dfcleanfecesfilt_s<-tax_glom(dfcleanfeces, "Species")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, stunted!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, pays!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, ageyears!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, sexe!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, (rowSums(otu_table(dfcleanfecesfilt_s))!=0))
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, (colSums(otu_table(dfcleanfecesfilt_s))!=0))
      dfcleanfecesfilt_s<-  phyloseq_rm_na_tax(dfcleanfecesfilt_s)
      
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + pays + age + read_count + calprotectinelevel + stunted)
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + pays + age + read_count + calprotectinelevel)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
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
      
      prev_counts.stunted <- dfcleanfecesfilt_s %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfecesfilt_s %>% #this produces a maximum possible prevalence count per OTU
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
      
      Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(x) x *100/sum(x))
      
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
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "StuntedSpecieslevelcorrcountrycalpro.csv")
      
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Species = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Species), levels=names(x))
      
      #make the actual plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      
      pdf("differentially_present_taxa_Species_stunted_corr_country_sexe_age_calpro_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=taxonomy, y=log2FoldChange, color=Order)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunted vs. non-stunted in feces, controlled for gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      pdf("differentially_present_SpeciesbystuntingstatusMicrobiomeInsightscorrcalpro.pdf", #name of file to print. can also include relative or absolute path before filename.
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
        ggtitle("Significantly different species by stunting status")
      dev.off()
      
      
      #### Feces: on Species level and stunting, correcting for alphaantitrypsin  (age taken out) ####
      dfcleanfecesfilt_s<-tax_glom(dfcleanfeces, "Species")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, stunted!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, pays!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, ageyears!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, sexe!="")
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, (rowSums(otu_table(dfcleanfecesfilt_s))!=0))
      dfcleanfecesfilt_s<-subset_samples(dfcleanfecesfilt_s, (colSums(otu_table(dfcleanfecesfilt_s))!=0))
      dfcleanfecesfilt_s<-  phyloseq_rm_na_tax(dfcleanfecesfilt_s)
      
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + pays + alphaantitrypsinlevel + read_count + stunted)
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + pays + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2 # none associated with at least two fold change
      
      
      
     
      #### Feces: on Species level and stunting, not corrected for inflammation, Mada only  ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, pays!="RCA")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + age +  read_count + stunted )
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + age +  read_count )
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt
      
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2 # none associated with at least two fold change
      
      # Make an abundance and prevalence table at Species level 
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.stunted <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfecesfilt_s2))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s2 %>%
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
      stuntedcorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt), ], "matrix"))
      View(stuntedcorrcountrygenderage_reduced_prevabdeseq)
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "stuntedcorrcountrygenderage_reduced_prevabdeseqMicrobiomeInsightsSpeciesnoncorrMadagascar.csv")
      
      #plot
      
      # Order order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Order = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order), levels=names(x))
      # Family order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Family = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family), levels=names(x))
      # Genus order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus), levels=names(x))
      # Species order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Species, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Species = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Species), levels=names(x))
      
      #make the actual plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      
      pdf("differentially_present_taxa_Species_stunted_corr_country_sexe_age_feces_LRTMicrobiomeInsightsMada.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=taxonomy, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      #make the graph
      pdf("differentially_present_SpeciesbystuntedMicrobiomeInsightsMada.pdf", #name of file to print. can also include relative or absolute path before filename.
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
        ggtitle("Significantly different Species by stunting, Madagascar only")
      dev.off()
      
      
      
      #### Feces: on Species level and stunting, not corrected for inflammation, CAR only ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, pays!="Madagascar")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + read_count + stunted)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run +  read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt
      
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2 
      
      # Make an abundance and prevalence table at Species level 
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.stunted <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfecesfilt_s2))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s2 %>%
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
      stuntedcorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt), ], "matrix"))
      View(stuntedcorrcountrygenderage_reduced_prevabdeseq)
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "SpeciesFecesStuntedCARonly.csv")
      
      #plot
      
      # Order order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Order = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order), levels=names(x))
      # Family order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Family = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family), levels=names(x))
      # Genus order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus), levels=names(x))
      # Species order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Species, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Species = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Species), levels=names(x))
      
      #make the actual plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      
      pdf("differentially_present_taxa_Species_stunted_corr_country_sexe_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=taxonomy, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      #make the graph
      pdf("differentially_present_SpeciesbystuntedMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
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
        ggtitle("Significantly different Species by stunting")
      dev.off()
      
     
      #### Feces: on Species level and stunting, correcting for calprotectine, Mada only  ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="RCA")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run +  calprotectinelevel + read_count + stunted)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run +  calprotectinelevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
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
      
      prev_counts.stunted <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfecesfilt_s2))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces %>%
        transform_sample_counts(function(x) {
          x/sum(x)} ) 
      
      Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "stunted")
      
      Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
      
      Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(x) x *100/sum(x))
      
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
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "StuntedSpecieslevelcorrcalproMadaonly.csv")
      
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Species = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Species), levels=names(x))
      
      #make the actual plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      
      pdf("differentially_present_taxa_Species_stunted_corr_sexe_calpro_feces_LRTMicrobiomeInsightsMadaonly.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=taxonomy, y=log2FoldChange, color=Order)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunted vs. non-stunted in feces, controlled for gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      pdf("differentially_present_SpeciesbystuntingstatusMicrobiomeInsightscorrcalproMadaonly.pdf", #name of file to print. can also include relative or absolute path before filename.
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
        ggtitle("Significantly different species by stunting status")
      dev.off()
      
        
      #### Feces: on Species level and stunting, correcting for calprotectine, CAR only NONE REMAINS SIGNIFICANT (age taken out) ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="Madagascar")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, stunted!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, age!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + calprotectinelevel + read_count + stunted)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + calprotectinelevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt 
      
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2 # none associated with at least two fold change
      
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
      
      prev_counts.stunted <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfecesfilt_s2))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces %>%
        transform_sample_counts(function(x) {
          x/sum(x)} ) 
      
      Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "stunted")
      
      Afribiota_feces_abundance_OTU_st=t(Afribiota_feces_abundance_OTU_st)
      
      Afribiota_feces_abundance_OTU_st=  transform_sample_counts(Afribiota_feces_abundance_OTU_st, function(x) x *100/sum(x))
      
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
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "StuntedSpecieslevelcorrcalproCARonly.csv")
      
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Species = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Species), levels=names(x))
      
      #make the actual plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      
      pdf("differentially_present_taxa_Species_stunted_corr_sexe_calpro_feces_LRTMicrobiomeInsightsCARonly.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=taxonomy, y=log2FoldChange, color=Order)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunted vs. non-stunted in feces, controlled for gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      pdf("differentially_present_SpeciesbystuntingstatusMicrobiomeInsightscorrcalproCARonly.pdf", #name of file to print. can also include relative or absolute path before filename.
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
        ggtitle("Significantly different species by stunting status")
      dev.off()
      
      #### on Species level and stunting, correcting for alphaantitrypsin, Mada only (age had to be taken out)  ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="RCA")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + alphaantitrypsinlevel + read_count + stunted)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + alphaantitrypsinlevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt
      
      #### on Species level and stunting, correcting for alphaantitrypsin, CAR only  (age had to be taken out) ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="Madagascar")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + alphaantitrypsinlevel + read_count + stunted)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run +  alphaantitrypsinlevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt 
      
      # Make an abundance and prevalence table at Species level 
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.stunted <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfecesfilt_s2))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s2 %>%
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
      stuntedcorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt), ], "matrix"))
      View(stuntedcorrcountrygenderage_reduced_prevabdeseq)
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "stuntedcorrcalproCARonly.csv")
      
      #plot
      
      # Order order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Order = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order), levels=names(x))
      # Family order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Family = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family), levels=names(x))
      # Genus order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus), levels=names(x))
      # Species order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Species, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Species = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Species), levels=names(x))
      
      #make the actual plot
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Species)
      
      pdf("differentially_present_taxa_Species_stunted_corr_country_sexe_age_AAT_feces_LRTMicrobiomeInsightsCAR.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=taxonomy, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      #make the graph
      pdf("differentially_present_SpeciesbystuntedcorrAATMicrobiomeInsightsCAR.pdf", #name of file to print. can also include relative or absolute path before filename.
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
        ggtitle("Significantly different Species by stunting, CAR only")
      dev.off()
      
      
      
      #### on Species level and country, not correcting for inflammation ####
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s, ~  run + stunted + sexe + age + read_count + pays)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      
      # show results preserving taxa table
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfecesfilt_s)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Species level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfecesfilt_s %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfecesfilt_s %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfecesfilt_s))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s %>%
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrstuntedgenderageMicrobiomeInsightsSpecies.csv")
      
      #plot
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$taxonomy<- paste0(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, "/", res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, "/", res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species)
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      # Species order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species), levels=names(x))
      
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Species_pays_corr_stunted_gender_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=taxonomy, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
      dev.off()
      
      # make a boxplot of rel. abundance
      dfcleanfecesfilt_s.rel <- microbiome::transform(dfcleanfecesfilt_s, "compositional")
      data<-data.frame(otu_table(dfcleanfecesfilt_s.rel))
      data$pays<-sample_data(dfcleanfecesfilt_s.rel)$pays
      
      # make  graph plot
      payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$Order, "/", payscorrcountrygenderage_reduced_prevabdeseq$Genus, "/", payscorrcountrygenderage_reduced_prevabdeseq$Species)
      pdf("differentially_present_SpeciesbypaysMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
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
      
      
      #### on Species level and country, correcting for calprotectine ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + calprotectinelevel + read_count + pays )
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + calprotectinelevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01)
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      # show results preserving taxa table
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt, "data.frame"), as(tax_table(dfcleanfecesfilt_s2)[rownames(res_p_ordered_filt), ], "matrix"))
      
      # Make an abundance and prevalence table at Species level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfecesfilt_s2))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s2 %>%
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
      payscorrcountrygenderage_reduced_prevabdeseq <- cbind(as(res_p_ordered_filt, "data.frame"), as(merge.abunanceprev.feces.OTU[rownames(res_p_ordered_filt), ], "matrix"))
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrcountrygenderagecalprotectine_reduced_prevabdeseqMicrobiomeInsightsSpecies.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      # Species order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species), levels=names(x))
      
      
      #make the actual plot
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$taxonomy<- paste0(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, "/", res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, "/", res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species)
      
      pdf("differentially_present_taxa_Species_pays_corr_country_gender_age_calpro_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=taxonomy, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & calpro,  LRT model")
      dev.off()
      
      # make a boxplot of rel. abundance
         
      # make  graph plot
         pdf("differentially_present_SpeciesbypayscorrcalproMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Species by country of origin, corrected for age, stunting status, calprotectine level")
      dev.off()
      
      
      #### on Species level and country, correcting for alphaantitrypsin ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + alphaantitrypsinlevel + read_count + pays )
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + alphaantitrypsinlevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfecesfilt_s2)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Species level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfecesfilt_s2))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s2 %>%
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrcountrygenderagealphaantitrypsinMicrobiomeInsights.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      # Species order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species), levels=names(x))
      
      
      #make the actual plot
      payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$Order, "/", payscorrcountrygenderage_reduced_prevabdeseq$Genus, "/", payscorrcountrygenderage_reduced_prevabdeseq$Species)
      
      pdf("differentially_present_taxa_Species_pays_corr_country_gender_age_aat_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Species, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & calpro,  LRT model")
      dev.off()
      
         
      # make a boxplot of rel. abundance
      dfcleanfecesfilt_s2.rel <- microbiome::transform(dfcleanfecesfilt_s2, "compositional")
      data<-data.frame(otu_table(dfcleanfecesfilt_s2.rel))
      data$pays<-sample_data(dfcleanfecesfilt_s2.rel)$pays
      
      # make  graph plot
      payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$Order, "/", payscorrcountrygenderage_reduced_prevabdeseq$Genus, "/", payscorrcountrygenderage_reduced_prevabdeseq$Species)
      pdf("differentially_present_SpeciesbypayscorraatMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
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
        ggtitle("Significantly different Species by country of origin, corrected for age, stunting status, alphaantitrypsin level")
      dev.off()
      
      
      #### on Species level and calprotectinelevel  ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + pays + read_count + calprotectinelevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + pays + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "calprotectinelevel") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2
      
      #### on Species level and calprotectinelevel, Madagascar only####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="RCA")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + read_count + calprotectinelevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "calprotectinelevel") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2
      
      #### on Species level and calprotectinelevel, CAR only  ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="Madagascar")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + read_count + calprotectinelevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "calprotectinelevel") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfecesfilt_s)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Species level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("calprotectinelevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("calprotectinelevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfecesfilt_s))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s2 %>%
        transform_sample_counts(function(x) {
          x/sum(x)} ) 
      
      Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "calprotectinelevel")
      
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "CARonlycalprolevelcorrcountrygenderage_reduced_prevabdeseqMicrobiomeInsights.csv")
      
      #plot
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      # Species order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species), levels=names(x))
      
      #make the actual plot
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$taxonomy<- paste0(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, "/", res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, "/", res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species)
      
      pdf("CARonlydifferentially_present_taxa_Species_Calprolevel_corr_country_gender_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=taxonomy, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("High vs. low AAT levels in feces, controlled for stunting, gender, age & country of origin LRT model")
      dev.off()
      
      # make a boxplot of rel. abundance
      dfcleanfecesfilt_s.rel <- microbiome::transform(dfcleanfecesfilt_s2, "compositional")
      data<-data.frame(otu_table(dfcleanfecesfilt_s.rel))
      data$pays<-sample_data(dfcleanfecesfilt_s.rel)$pays
      
      # make  graph plot
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$taxonomy<- paste0(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, "/", res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, "/", res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species)
      pdf("CARonlydifferentially_present_SpeciesbyCalprolevelMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Species by Calprotectine level, CAR only")
      dev.off()
      
      
      
      #### on Species level and AAT level ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + pays + read_count + alphaantitrypsinlevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + pays + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "alphaantitrypsinlevel") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfecesfilt_s)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Species level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfecesfilt_s %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("alphaantitrypsinlevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfecesfilt_s %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("alphaantitrypsinlevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfecesfilt_s))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s %>%
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "AATlevelcorrcountrygenderage_reduced_prevabdeseqMicrobiomeInsights.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      # Species order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species), levels=names(x))
      
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Species_AATlevel_corr_country_gender_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Species, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("High vs. low AAT levels in feces, controlled for stunting, gender, age & country of origin LRT model")
      dev.off()
      
          # make  graph plot
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$taxonomy<- paste0(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, "/", res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, "/", res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species)
      
      pdf("differentially_present_SpeciesbyAATlevelMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Species by AAT level")
      dev.off()
      
      
      #### on Species level and AAT level, Madagascar only  ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="RCA")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + read_count + alphaantitrypsinlevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "alphaantitrypsinlevel") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfecesfilt_s2)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Species level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("alphaantitrypsinlevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("alphaantitrypsinlevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfecesfilt_s2))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s2 %>%
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "AATlevelcorrgenderage_reduced_prevabdeseqMicrobiomeInsightsMada.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      # Species order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species), levels=names(x))
      
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Species_AATlevel_corr_gender_age_feces_LRTMicrobiomeInsightsMada.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Species, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("High vs. low AAT levels in feces, controlled for stunting, gender, age & country of origin LRT model")
      dev.off()
      
      # make  graph plot
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$taxonomy<- paste0(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, "/", res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, "/", res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Species)
      
      pdf("differentially_present_SpeciesbyAATlevelMicrobiomeInsightsMada.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x = reorder(taxonomy, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Species by AAT level")
      dev.off()
      
      
      
      #### on Species level and AAT level, RCA only  ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="Madagascar")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, stunted!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + read_count + alphaantitrypsinlevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "alphaantitrypsinlevel") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2
     
      
      #### on Genus level and country, not correcting for inflammation ####
      dfcleanfeces_g<-tax_glom(dfcleanfeces, "Genus")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g, ~  run + stunted + sexe + age + read_count + pays)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_g)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Species level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfeces_g %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfeces_g %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfeces_g))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_g %>%
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrgenderageMicrobiomeInsightsGenus.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Genus_pays_corr_gender_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      pdf("differentially_present_GenusbypaysMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Genus by country of origin")
      dev.off()
      
      #### on Genus level and country,  correcting for calprotectine ####
      dfcleanfeces_g<-tax_glom(dfcleanfeces, "Genus")
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g, calprotectinelevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g2, ~  run + stunted + sexe + age + calprotectinelevel + read_count + pays)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + calprotectinelevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") 
      
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
      
      # Make an abundance and prevalence table at Genus/Genus level
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrcountrygenderagecalprodeseqMicrobiomeInsightsGenus.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Genus_pays_corr_country_gender_age_calpro_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & calprotectine levels, LRT model")
      dev.off()
      
      # make  graph plot
      payscorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(payscorrcountrygenderage_reduced_prevabdeseq$Order, "/", payscorrcountrygenderage_reduced_prevabdeseq$Genus)
      pdf("differentially_present_GenusbypayscontrolledcalproMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
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
      
      
      
      
      
      
      #### on Genus level and country,  correcting for alphaantitrypsin ####
      dfcleanfeces_g<-tax_glom(dfcleanfeces, "Genus")
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g, alphaantitrypsinlevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g2, ~  run + stunted + sexe + age + alphaantitrypsinlevel + read_count + pays)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + alphaantitrypsinlevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") 
      
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
      
      # Make an abundance and prevalence table at Genus/Genus level
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrcountrygenderageaatMicrobiomeInsightsGenus.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Genus_pays_corr_gender_age_aat_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & alphaantitrypsin levels, LRT model")
      dev.off()
      
      # make  graph plot
      pdf("differentially_present_GenusbypayscontrolledaatMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Genus by country of origin, controlled by alphaantitrypsin levels")
      dev.off()
      
      
      
      
      
      
      #### on Genus level and stunting, not correcting for inflammation ####
      dfcleanfeces_g<-tax_glom(dfcleanfeces, "Genus")
      dfcleanfeces_g<-subset_samples(dfcleanfeces_g, stunted!="")
      dfcleanfeces_g<-subset_samples(dfcleanfeces_g, pays!="")
      dfcleanfeces_g<-subset_samples(dfcleanfeces_g, ageyears!="")
      dfcleanfeces_g<-subset_samples(dfcleanfeces_g, sexe!="")
      dfcleanfeces_g = filter_taxa(dfcleanfeces_g, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_g<-subset_samples(dfcleanfeces_g, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g, ~  run + pays + sexe + age + read_count + stunted)
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + pays + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
      # Make an abundance and prevalence table at Genus level 
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.stunted <- dfcleanfeces_g %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfeces_g %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces_g))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_g %>%
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
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "GenusstuntedcorrcountrygenderageMicrobiomeInsights.csv")
      
      #plot
      
      # Order order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Family = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family), levels=names(x))
      # Family order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Order = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order), levels=names(x))
      # Genus order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus), levels=names(x))
      # Genus order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus), levels=names(x))
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Genus_stunted_corr_country_sexe_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunted vs. non stunted in feces, controlled for stunting, gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      #make the graph
      pdf("differentially_present_GenusbystuntedMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 12, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Genus by stunting status")
      dev.off()
      
      
      
      #### on Genus level and stunting, not corrected for inflammation, Mada only ####
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g, pays!="RCA")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g2, ~  run + sexe + age + read_count + stunted)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + sexe + age + read_count )
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
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
      
      # Make an abundance and prevalence table at Genus level
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
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfeces_g2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "MadaonlyGenusstuntedcorrgenderageMicrobiomeInsights.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      
      
      #make the actual plot
      
      pdf("Madaonlydifferentially_present_taxa_Genus_stunted_corr_gender_age__feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunting status in feces from Mada, controlled for country, gender, age & ,  LRT model")
      dev.off()
      
      
      # make  graph plot
      pdf("Madaonlydifferentially_present_GenusbystuntingcorrMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Genus by stunting status in Mada, corrected for age, country of origin,  level")
      dev.off()
      
      
      
      #### on Genus level and stunting, not corrected for inflammation, CAR only ####
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g, pays!="Madagascar")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g2, ~  run + sexe + age + read_count + stunted)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
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
      
      # Make an abundance and prevalence table at Genus level
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
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfeces_g2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
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
      
      Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "stunted")
      
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "CARonlystuntedcorrgenderageqMicrobiomeInsights.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      
      
      #make the actual plot
      
      pdf("CARonlydifferentially_present_taxa_Genus_stunted_corr_gender_age__feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunting statusi in feces from CAR, controlled for country, gender, age & ,  LRT model")
      dev.off()
      
      # make  graph plot
      pdf("CARonlydifferentially_present_GenusbystuntingcorrMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Genus by stunting status in CAR, corrected for age, country of origin,  level")
      dev.off()
      
      
      
      
      #### on Genus level and stunting, correcting for calprotectin levels and age ####
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g, calprotectinelevel!="")
      
      dfcleanfeces_g2 = filter_taxa(dfcleanfeces_g2, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g2, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g2, ~  run + age + pays + sexe + calprotectinelevel + read_count + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + age + pays + sexe + calprotectinelevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
     
      #### on Genus level and stunting, correcting for alphaantitrypsin levels ####
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g, alphaantitrypsinlevel!="")
      
      dfcleanfeces_g2 = filter_taxa(dfcleanfeces_g2, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g2, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g2, ~  run + pays + sexe + age + alphaantitrypsinlevel + read_count + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + pays + sexe + age + alphaantitrypsinlevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
     
      #### on Genus level and stunting, correcting for calprotectin levels Mada only ####
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g, calprotectinelevel!="")
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g2, pays!="RCA")
      
      dfcleanfeces_g2 = filter_taxa(dfcleanfeces_g2, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g2, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g2, ~  run + sexe + age + calprotectinelevel + read_count + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + sexe + age + calprotectinelevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
     
      #### on Genus level and stunting, correcting for alphaantitrypsin levels Mada only ####
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g, alphaantitrypsinlevel!="")
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g2, pays!="RCA")
      
      dfcleanfeces_g2 = filter_taxa(dfcleanfeces_g2, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g2, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g2, ~  run + sexe + age + alphaantitrypsinlevel + read_count + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + sexe + age + alphaantitrypsinlevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
     
      #### on Genus level and stunting, correcting for calprotectin levels CAR only  ####
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g, calprotectinelevel!="")
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g2, pays!="Madagascar")
      
      dfcleanfeces_g2 = filter_taxa(dfcleanfeces_g2, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g2, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g2, ~  run + sexe + read_count + calprotectinelevel + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + sexe + read_count + calprotectinelevel)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
      # Make an abundance and prevalence table at Genus level 
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.stunted <- dfcleanfeces_g2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfeces_g2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces_g2))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_g2 %>%
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
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "CARonlygenusstuntedcorrgenderagecalpro_reduced_prevabdeseqMicrobiomeInsights.csv")
      
      #plot
      
      # Order order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Family = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family), levels=names(x))
      # Family order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Order = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order), levels=names(x))
      # Genus order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus), levels=names(x))
      
      #make the actual plot
      
      pdf("CARdifferentially_present_taxa_Genus_stunted_corr_sexe_age_calpro_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunted vs. non stunted in feces, controlled for stunting, gender, age & calprotectine level,  CAR only, LRT model")
      dev.off()
      
      # make  graph plot
      #make the graph
      pdf("CARonlydifferentially_present_GenusbystuntedcorrcalprotectinMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 12, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Genus by stunting status, controlled for calprotectine levels, CAR only")
      dev.off()
      
      
      
      #### on Genus level and stunting, correcting for alphaantitrypsin levels CAR only ####
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g, alphaantitrypsinlevel!="")
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g2, pays!="Madagascar")
      
      dfcleanfeces_g2 = filter_taxa(dfcleanfeces_g2, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_g2<-subset_samples(dfcleanfeces_g2, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_g2, ~  run + sexe + age + alphaantitrypsinlevel + read_count + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + sexe + age + alphaantitrypsinlevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
      # Make an abundance and prevalence table at Genus level 
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.stunted <- dfcleanfeces_g2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfeces_g2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces_g2))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_g2 %>%
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
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "CARonlygenusstuntedcorrgenderageaat_reduced_prevabdeseqMicrobiomeInsights.csv")
      
      #plot
      
      # Order order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Family = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family), levels=names(x))
      # Family order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Order = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order), levels=names(x))
      # Genus order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus), levels=names(x))
      # Genus order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Genus), levels=names(x))
      # Accession order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$accession, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Accession = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$accession), levels=names(x))
      
      #make the actual plot
      
      pdf("CARonlydifferentially_present_taxa_Genus_stunted_corr_sexe_age_aat_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunted vs. non stunted in feces, controlled for stunting, gender, age & alphaantitrypsin level,  CAR only, LRT model")
      dev.off()
      
      # make  graph plot
      #make the graph
      pdf("CARonlydifferentially_present_GenusbystuntedcorralphaantitrypsinMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 12, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Genus by stunting status, controlled for alphaantitrypsin levels, CAR only")
      dev.off()
      
 
      #### on Genus level and calprotectinelevel  ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + pays + read_count + calprotectinelevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + pays + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "calprotectinelevel") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2
      
     
      
      #### on Genus level and calprotectinelevel, Madagascar only####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="RCA")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + read_count + calprotectinelevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "calprotectinelevel") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2
      
     
      #### on Genus level and calprotectinelevel, CAR only  ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="Madagascar")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + read_count + calprotectinelevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "calprotectinelevel") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfecesfilt_s)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Genus level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("calprotectinelevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("calprotectinelevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfecesfilt_s))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s2 %>%
        transform_sample_counts(function(x) {
          x/sum(x)} ) 
      
      Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "calprotectinelevel")
      
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "CARonlycalprolevelcorrgenderage_reduced_prevabdeseqMicrobiomeInsightsGenus.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      
      #make the actual plot
      
      pdf("CARonlydifferentially_present_taxa_Genus_Calprolevel_corr_gender_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("High vs. low AAT levels in feces, controlled for stunting, gender, age & country of origin LRT model")
      dev.off()
      
      # make  graph plot
      pdf("CARonlydifferentially_present_GenusbyCalprolevelMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Genus by Calprotectine level, CAR only")
      dev.off()
      
      
      
      #### on Genus level and AAT level ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + pays + read_count + alphaantitrypsinlevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + pays + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "alphaantitrypsinlevel") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2
      
     
      
      #### on Genus level and AAT level, Madagascar only  ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="RCA")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, stunted!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + read_count + alphaantitrypsinlevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "alphaantitrypsinlevel") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfecesfilt_s2)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Genus level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("alphaantitrypsinlevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("alphaantitrypsinlevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfecesfilt_s2))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s2 %>%
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "Madagascaronlyaatlevelcorrgenderage_reduced_prevabdeseqMicrobiomeInsightsGenus.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      # Family order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Order), levels=names(x))
      # Genus order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Genus), levels=names(x))
      
      #make the actual plot
      
      pdf("Madagascaronlydifferentially_present_taxa_Genus_aatlevel_corr_gender_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("High vs. low AAT levels in feces, controlled for stunting, gender, age & country of origin LRT model")
      dev.off()
      
      # make  graph plot
      pdf("Madagascaronlydifferentially_present_GenusbyaatlevelMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x = reorder(Genus, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Genus by alphaantitrypsin level, Madagascar only")
      dev.off()
      
      
      
      #### on Genus level and AAT level, CAR only NO RESULTS ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="Madagascar")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, stunted!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + read_count + alphaantitrypsinlevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "alphaantitrypsinlevel") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2
      
      
      
      
      #### on Order level and country, not correcting for inflammation ####
      dfcleanfeces_o<-tax_glom(dfcleanfeces, "Order")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_o, ~  run + stunted + sexe + age + read_count + pays)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_o)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Order level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfeces_o %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfeces_o %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfeces_o))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_o %>%
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrgenderageMicrobiomeInsightsOrder.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Order_pays_corr_oender_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Order, y=log2FoldChange, color=Class)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age, LRT model")
      dev.off()
      
      # make  graph plot
      pdf("differentially_present_OrderbypaysMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Class, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Order by country of origin")
      dev.off()
      
         
      #### on Order level and country,  correcting for calprotectine ####
      dfcleanfeces_o<-tax_glom(dfcleanfeces, "Order")
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o, calprotectinelevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_o2, ~  run + stunted + sexe + age + calprotectinelevel + read_count + pays)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + calprotectinelevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_o2)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Order level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfeces_o2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfeces_o2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfeces_o2))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_o2 %>%
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrgenderagecalprodeseqMicrobiomeInsightsOrder.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Order_pays_corr_gender_age_calpro_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Order, y=log2FoldChange, color=Class)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & calprotectine levels, LRT model")
      dev.off()
      
      # make  graph plot
      pdf("differentially_present_OrderbypayscontrolledcalproMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Order, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Order by country of origin, controlled by calprotectine levels")
      dev.off()
      
      
      
      
      
      
      #### on Order level and country,  correcting for alphaantitrypsin ####
      dfcleanfeces_o<-tax_glom(dfcleanfeces, "Order")
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o, alphaantitrypsinlevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_o2, ~  run + stunted + sexe + age + alphaantitrypsinlevel + read_count + pays)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + alphaantitrypsinlevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "pays_RCA_vs_Madagascar") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_o2)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Order level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfeces_o2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfeces_o2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("pays") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfeces_o2))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_o2 %>%
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "payscorrgenderageaatMicrobiomeInsightsOrder.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Order_pays_corr_gender_age_aat_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Order, y=log2FoldChange, color=Class)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & alphaantitrypsin levels, LRT model")
      dev.off()
      
      # make  graph plot
      pdf("differentially_present_OrderbypayscontrolledaatMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=payscorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Order, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Order by country of origin, controlled by alphaantitrypsin levels")
      dev.off()
      
      
      
      
      
      
      #### on Order level and stunting, not correcting for inflammation  ####
      dfcleanfeces_o<-tax_glom(dfcleanfeces, "Order")
      dfcleanfeces_o<-subset_samples(dfcleanfeces_o, stunted!="")
      dfcleanfeces_o<-subset_samples(dfcleanfeces_o, pays!="")
      dfcleanfeces_o<-subset_samples(dfcleanfeces_o, ageyears!="")
      dfcleanfeces_o<-subset_samples(dfcleanfeces_o, sexe!="")
      dfcleanfeces_o = filter_taxa(dfcleanfeces_o, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_o<-subset_samples(dfcleanfeces_o, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_o, ~  run + pays + sexe + age + read_count + stunted)
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + pays + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
 
      # show results preserving taxa table
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_o2)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Order level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.stunted <- dfcleanfeces_o2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfeces_o2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces_o2))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_o2 %>%
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
      
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "stuntedcorrgenderageMicrobiomeInsightsOrder.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Order_stunted_corr_country_gender_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Order, y=log2FoldChange, color=Class)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & alphaantitrypsin levels, LRT model")
      dev.off()
      
      # make  graph plot
      pdf("differentially_present_OrderbystuntedcontrolledpaysMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Order, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Order by stunting status, controlled by country of origin, gender, age ")
      dev.off()
      
      
      
      #### on Order level and stunting, not corrected for inflammation, Mada only  ####
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o, pays!="RCA")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_o2, ~  run + sexe + age + read_count + stunted)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + sexe + age + read_count )
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2
      
      
      
      
      #### on Order level and stunting, not corrected for inflammation, CAR only  ####
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o, pays!="Madagascar")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_o2, ~  run + sexe + age + read_count + stunted)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_o22)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Order level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.stunted <- dfcleanfeces_o22 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfeces_o22 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces_o22))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_o22 %>%
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
      
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "stuntedcorrgenderageMicrobiomeInsightsOrderBangui.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Order_stunted_corr_gender_age_feces_LRTMicrobiomeInsightsBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Order, y=log2FoldChange, color=Class)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & alphaantitrypsin levels, LRT model")
      dev.off()
      
      # make  graph plot
      pdf("differentially_present_OrderbystuntedcontrolledMicrobiomeInsightsBangui.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Order, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Order by stunting status, controlled by country of origin, gender, age ")
      dev.off()
      
      
      
      
      #### on Order level and stunting, correcting for calprotectin levels  ####
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o, calprotectinelevel!="")
      
      dfcleanfeces_o2 = filter_taxa(dfcleanfeces_o2, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o2, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_o2, ~  run + pays + sexe + age + calprotectinelevel + read_count + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + pays + sexe + age + calprotectinelevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
      # show results preserving taxa table
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfeces_o22)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Order level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.stunted <- dfcleanfeces_o22 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfeces_o22 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces_o22))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_o22 %>%
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
      
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "stuntedcorrgenderagecalproMicrobiomeInsightsOrder.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Order_stunted_corr_gender_age_calpro_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Order, y=log2FoldChange, color=Class)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("CAR vs Madagascar in feces, controlled for stunting, gender, age & alphaantitrypsin levels, LRT model")
      dev.off()
      
      # make  graph plot
      pdf("differentially_present_OrderbystuntedcontrolledcalproMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Order, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Order by stunting status, controlled by country of origin, gender, age, calprotectine level ")
      dev.off()
      
      
      
      #### on Order level and stunting, correcting for alphaantitrypsin levels  ####
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o, alphaantitrypsinlevel!="")
      
      dfcleanfeces_o2 = filter_taxa(dfcleanfeces_o2, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o2, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_o2, ~  run + pays + sexe + age + alphaantitrypsinlevel + read_count + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + pays + sexe + age + alphaantitrypsinlevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
      # Make an abundance and prevalence table at Order level 
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.stunted <- dfcleanfeces_o2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfeces_o2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces_o2))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_o2 %>%
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
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "OrderstuntedcorrcountrygenderageaatMicrobiomeInsights.csv")
      
      #plot
      
      # Order order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Family = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family), levels=names(x))
      
      #make the actual plot
      
      pdf("differentially_present_taxa_Order_stunted_corr_country_sexe_age_aat_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=Order, y=log2FoldChange, color=Class)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunted vs. non stunted in feces, controlled for stunting, gender, age & alphaantitrypsin level,  LRT model")
      dev.off()
      
      # make  graph plot
      #make the graph
      pdf("differentially_present_OrderbystuntedcorralphaantitrypsinMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 12, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Order, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Order by stunting status, controlled for alphaantitrypsin levels")
      dev.off()
      
      
      
      #### on Order level and stunting, correcting for calprotectin levels Mada only  ####
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o, calprotectinelevel!="")
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o2, pays!="RCA")
      
      dfcleanfeces_o2 = filter_taxa(dfcleanfeces_o2, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o2, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_o2, ~  run + sexe + age + calprotectinelevel + read_count + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + sexe + age + calprotectinelevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
      
      
      #### on Order level and stunting, correcting for alphaantitrypsin levels Mada only ####
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o, alphaantitrypsinlevel!="")
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o2, pays!="RCA")
      
      dfcleanfeces_o2 = filter_taxa(dfcleanfeces_o2, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o2, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_o2, ~  run + sexe + age + alphaantitrypsinlevel + read_count + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + sexe + age + alphaantitrypsinlevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
      
      #### on Order level and stunting, correcting for calprotectin levels CAR only  ####
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o, calprotectinelevel!="")
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o2, pays!="Madagascar")
      
      dfcleanfeces_o2 = filter_taxa(dfcleanfeces_o2, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o2, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_o2, ~  run + sexe + age + calprotectinelevel + read_count + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + sexe + age + calprotectinelevel)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
      # Make an abundance and prevalence table at Order level 
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.stunted <- dfcleanfeces_o2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfeces_o2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces_o2))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_o2 %>%
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
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "CARonlyOrderstuntedcorrcountrygenderagecalpro_reduced_prevabdeseqMicrobiomeInsights.csv")
      
      #plot
      
      # Order order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Family, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Family = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Family), levels=names(x))
      # Family order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Order = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order), levels=names(x))
      # Order order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Order = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order), levels=names(x))
      # Order order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Order = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order), levels=names(x))
      # Accession order
      x = tapply(stuntedcorrcountrygenderage_reduced_prevabdeseq$log2FoldChange, stuntedcorrcountrygenderage_reduced_prevabdeseq$accession, function(x) max(x))
      x = sort(x, TRUE)
      stuntedcorrcountrygenderage_reduced_prevabdeseq$Accession = factor(as.character(stuntedcorrcountrygenderage_reduced_prevabdeseq$accession), levels=names(x))
      
      #make the actual plot
      
      pdf("CARdifferentially_present_taxa_Order_stunted_corr_country_sexe_age_calpro_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=Order, y=log2FoldChange, color=Family)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunted vs. non stunted in feces, controlled for stunting, gender, age & calprotectine level,  CAR only, LRT model")
      dev.off()
      
      # make  graph plot
      #make the graph
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Order)
      pdf("CARonlydifferentially_present_OrderbystuntedcorrcalprotectinMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
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
        ggtitle("Significantly different Order by stunting status, controlled for calprotectine levels, CAR only")
      dev.off()
      
      
      
      #### on Order level and stunting, correcting for alphaantitrypsin levels CAR only  ####
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o, alphaantitrypsinlevel!="")
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o2, pays!="Madagascar")
      
      dfcleanfeces_o2 = filter_taxa(dfcleanfeces_o2, function(x) sum(x) > 0, TRUE)
      dfcleanfeces_o2<-subset_samples(dfcleanfeces_o2, (rowSums(otu_table(dfcleanfeces))!=0))
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfeces_o2, ~  run + sexe + age + alphaantitrypsinlevel + read_count + stunted )
      
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + sexe + age + alphaantitrypsinlevel + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "stunted_stunted_vs_non.stunted") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      View(res_p_ordered_filt_2)
      
      # Make an abundance and prevalence table at Order level 
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.stunted <- dfcleanfeces_o2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.stunted) <- paste("prevalence", colnames(prev_counts.stunted), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.stunted <- dfcleanfeces_o2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("stunted") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.stunted) <- paste("prevalence", colnames(prev_possible.stunted), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.stunted/prev_possible.stunted)*100
      
      tax_table.stunted =  as.data.frame(tax_table(dfcleanfeces_o2))
      
      merge.stunted.prev= merge(tax_table.stunted, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfeces_o2 %>%
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
      write.csv(stuntedcorrcountrygenderage_reduced_prevabdeseq, "CARonlyOrderstuntedcorrcountrygenderageaat_reduced_prevabdeseqMicrobiomeInsights.csv")
      
      #plot
      
      #make the actual plot
      
      pdf("CARdifferentially_present_taxa_Order_stunted_corr_country_sexe_age_aat_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x=Order, y=log2FoldChange, color=Order)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("Stunted vs. non stunted in feces, controlled for stunting, gender, age & calprotectine level,  CAR only, LRT model")
      dev.off()
      
      # make  graph plot
      #make the graph
      stuntedcorrcountrygenderage_reduced_prevabdeseq$taxonomy<- paste0(stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Order, "/", stuntedcorrcountrygenderage_reduced_prevabdeseq$Order)
      pdf("CARonlydifferentially_present_OrderbystuntedcorraatMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 12, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=stuntedcorrcountrygenderage_reduced_prevabdeseq, aes(x = reorder(Order, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Order by stunting status, controlled for AAT levels, CAR only")
      dev.off()
      
      
      #### on Order level and calprotectinelevel ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + pays + read_count + calprotectinelevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + pays + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "calprotectinelevel") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2
      

      
      #### on Order level and calprotectinelevel, Madagascar only####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="RCA")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + read_count + calprotectinelevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "calprotectinelevel") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2
      
      #### on Order level and calprotectinelevel, CAR only  ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, calprotectinelevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="Madagascar")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + read_count + calprotectinelevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "calprotectinelevel") #log2( too high / normal ) #so when we say things are "upregulated" we mean more prevalent in stunted individuals.)
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfecesfilt_s)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Order level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("calprotectinelevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("calprotectinelevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfecesfilt_s))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s2 %>%
        transform_sample_counts(function(x) {
          x/sum(x)} ) 
      
      Afribiota_feces_abundance_OTU_st <- merge_samples(Afribiota_feces_abundance_OTU, "calprotectinelevel")
      
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "CARonlycalprolevelcorrcountrygenderagedeseqMicrobiomeInsightsOrder.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      
      
      #make the actual plot
      
      pdf("CARonlydifferentially_present_taxa_Order_Calprolevel_corr_country_gender_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Order, y=log2FoldChange, color=Class)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("High vs. low AAT levels in feces, controlled for stunting, gender, age & country of origin LRT model")
      dev.off()
      
      
      # make  graph plot
      pdf("CARonlydifferentially_present_OrderbyCalprolevelMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x = reorder(Class, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Order by Calprotectine level, CAR only")
      dev.off()
      
      
      
      #### on Order level and AAT level ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + pays + read_count + alphaantitrypsinlevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + pays + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "alphaantitrypsinlevel") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2
      
     
      
      #### on Order level and AAT level, Madagascar only ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="RCA")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + read_count + alphaantitrypsinlevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "alphaantitrypsinlevel") 
      
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
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced <- cbind(as(res_p_ordered_filt_2, "data.frame"), as(tax_table(dfcleanfecesfilt_s)[rownames(res_p_ordered_filt_2), ], "matrix"))
      
      # Make an abundance and prevalence table at Order level
      #prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). 
      prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function
        x[x >= 2] <- 1
        return(x)
      }
      
      allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each OTU
        x[x >= 0] <- 1
        return(x)
      }
      
      prev_counts.pays <- dfcleanfecesfilt_s2 %>% #this produces prevalence "counts" for each OTU, but not percentages
        transform_sample_counts(fun = prevalence) %>%
        merge_samples("alphaantitrypsinlevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_counts.pays) <- paste("prevalence", colnames(prev_counts.pays), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables
      
      prev_possible.pays <- dfcleanfecesfilt_s2 %>% #this produces a maximum possible prevalence count per OTU
        transform_sample_counts(fun = allones) %>%
        merge_samples("alphaantitrypsinlevel") %>%
        t() %>%
        otu_table() %>%
        as.data.frame()
      colnames(prev_possible.pays) <- paste("prevalence", colnames(prev_possible.pays), sep = ".") #add something to distinguish between relative abundance and prevalence
      
      #dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-OTU basis
      test.prev = (prev_counts.pays/prev_possible.pays)*100
      
      tax_table.pays =  as.data.frame(tax_table(dfcleanfecesfilt_s))
      
      merge.pays.prev= merge(tax_table.pays, test.prev, by="row.names")
      
      Afribiota_feces_abundance_OTU <- dfcleanfecesfilt_s2 %>%
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
      
      write.csv(payscorrcountrygenderage_reduced_prevabdeseq, "CARonlyaatlevelcorrcountrygenderagedeseqMicrobiomeInsightsOrder.csv")
      
      #plot
      
      # Order order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      # Class order
      x = tapply(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$log2FoldChange, res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family, function(x) max(x))
      x = sort(x, TRUE)
      res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family = factor(as.character(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced$Family), levels=names(x))
      
      
      #make the actual plot
      
      pdf("CARonlydifferentially_present_taxa_Order_aatlevel_corr_country_gender_age_feces_LRTMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 10, # define plot width and height. completely up to user.
          height = 8)
      ggplot(res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x=Order, y=log2FoldChange, color=Class)) + geom_point(size=2) + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size=8,  margin = margin(t = 10, r = 10, b = 10, l = 10))) + ggtitle("High vs. low AAT levels in feces, controlled for stunting, gender, age & country of origin LRT model")
      dev.off()
      
      
      # make  graph plot
      pdf("CARonlydifferentially_present_OrderbyaatlevelMicrobiomeInsights.pdf", #name of file to print. can also include relative or absolute path before filename.
          width = 20, # define plot width and height. completely up to user.
          height = 8)
      ggplot(data=res_p_ordered_wtaxa.stuntedcorrcountrygenderage_reduced, aes(x = reorder(Class, log2FoldChange), y=log2FoldChange)) +
        geom_bar(position="dodge",stat="identity", color="black") +
        coord_flip() + 
        theme(axis.text.x = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        theme(axis.text.y = element_text(size=12))+ 
        xlab("")+ 
        theme(axis.title.y = element_text(size=14))+
        theme(title = element_text(size=16, face="bold"))+
        ylab("log2 fold change")+
        ggtitle("Significantly different Order by alphaantitrypsin level, CAR only")
      dev.off()
      
      
      
      #### on Order level and AAT level, RCA only NONE ASSOCIATED ####
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s, alphaantitrypsinlevel!="")
      dfcleanfecesfilt_s2<-subset_samples(dfcleanfecesfilt_s2, pays!="Madagascar")
      diagddsfecesDM = phyloseq_to_deseq2(dfcleanfecesfilt_s2, ~  run + stunted + sexe + age + read_count + alphaantitrypsinlevel)
      gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      geoMeans = apply(counts(diagddsfecesDM), 1, gm_mean)
      diagddsfecesDM = estimateSizeFactors(diagddsfecesDM, geoMeans=geoMeans)
      diagddsfecesDM = DESeq(diagddsfecesDM, test="LRT", fitType="parametric", reduced= ~ run + stunted + sexe + age + read_count)
      
      resultsNames(diagddsfecesDM)
      res = results(diagddsfecesDM, cooksCutoff = FALSE, alpha = 0.01, pAdjustMethod = "BH", name = "alphaantitrypsinlevel") 
      
      # view a summary of the results table with a padj value < 0.01
      summary(res, alpha = 0.01) 
      
      #filtering the results table #
      # reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered <- res[order(res$padj, na.last = NA), ]
      # filter out any results which have a adjust p-value less than alpha (0.01) and for changes of at least 2
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < 0.01), ]
      res_p_ordered_filt_2 <- res_p_ordered_filt[which(res_p_ordered_filt$log2FoldChange >= 1 | res_p_ordered_filt$log2FoldChange <= -1), ]
      res_p_ordered_filt_2
      
      
      
 #### Are given ASV associated with the nutritional status/stunting? presence/absence and Chi2 in loop####    
      #### Make chi2 test on the presence/absence of given fungi not correcting for confounding factors ####    
      #### On ASV level no correction for stunting ####
      fecesprevalence <- dfcleanfeces %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence2<-subset_samples(fecesprevalence, stunted!="") # 575 samples remaining
      fecesprevalence2<-subset_samples(fecesprevalence2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence2 = filter_taxa(fecesprevalence2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
      ASVnames<-colnames(df_chi2)
      chi2results = data.frame(ASVnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.stuntedASVnoncorrectedchi2MicrobiomeInsights.csv")
      
      #### On ASV level no correction for country ####
      meta_chi2$pays<-as.factor(meta_chi2$pays)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value)
      ASVnames<-colnames(df_chi2)
      chi2results = data.frame(ASVnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.paysASVnoncorrectedchi2MicrobiomeInsights.csv")
      
      #### On ASV level no correction for ageyears in three levels ####
      meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
      ASVnames<-colnames(df_chi2)
      chi2results = data.frame(ASVnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsASVnoncorrectedchi2MicrobiomeInsights.csv")
      
      #### On ASV level no correction for SampleType in three levels ####
      allprevalence <- dfclean %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(allprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      allprevalence2<-subset_samples(allprevalence, stunted!="") # 575 samples remaining
      allprevalence2<-subset_samples(allprevalence2, pays!="") # 575 samples remaining, no additional ones taken out
      allprevalence2 = filter_taxa(allprevalence2, function(x) sum(x) > 0, TRUE)
      
      
      df_chi2 <- as.matrix(t(otu_table(allprevalence2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(allprevalence2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$SampleType), simulate.p.value = TRUE)$p.value)
      ASVnames<-colnames(df_chi2)
      chi2results = data.frame(ASVnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(allprevalence2))
      Tax_corr$name<-row.names(tax_table(allprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsSampleTypeASVnoncorrectedchi2MicrobiomeInsightsfulldataset.csv")
      
      
      
      #### On ASV level no correction for SampleType in three levels, stunted only ####
      allprevalence <- dfclean %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(allprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      allprevalence2<-subset_samples(allprevalence, stunted!="") # 575 samples remaining
      allprevalence2<-subset_samples(allprevalence2, pays!="") # 575 samples remaining, no additional ones taken out
      allprevalence2 = filter_taxa(allprevalence2, function(x) sum(x) > 0, TRUE)
      allprevalence2<-subset_samples(allprevalence2, stunted=="stunted")
      
      df_chi2 <- as.matrix(t(otu_table(allprevalence2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(allprevalence2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$SampleType), simulate.p.value = TRUE)$p.value)
      ASVnames<-colnames(df_chi2)
      chi2results = data.frame(ASVnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(allprevalence2))
      Tax_corr$name<-row.names(tax_table(allprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsSampleTypeASVnoncorrectedchi2MicrobiomeInsightsstuntedonly.csv")
      
      
      #### On ASV level no correction for SampleType in three levels, shared samples only ####
      allprevalence <- dfclean %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(allprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      allprevalence2<-subset_samples(allprevalence, stunted!="") # 575 samples remaining
      allprevalence2<-subset_samples(allprevalence2, pays!="") # 575 samples remaining, no additional ones taken out
      allprevalence2 = filter_taxa(allprevalence2, function(x) sum(x) > 0, TRUE)
      allprevalence2<-subset_samples(allprevalence2, sample_names(allprevalence2) %in% sample_names(dfrar5000_rest2))
      table(sample_data(allprevalence2)$SampleType)
      
      
      df_chi2 <- as.matrix(t(otu_table(allprevalence2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(allprevalence2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$SampleType), simulate.p.value = TRUE)$p.value)
      ASVnames<-colnames(df_chi2)
      chi2results = data.frame(ASVnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(allprevalence2))
      Tax_corr$name<-row.names(tax_table(allprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsSampleTypeASVnoncorrectedchi2MicrobiomeInsightsreduceddataset.csv")
      
      
      #### On Species level no correction for stunting ####
      dfcleanfeces_s<- tax_glom(dfcleanfeces, "Species")
      fecesprevalence_s <- dfcleanfeces_s %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, stunted!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      table(df_chi2[, 13], meta_chi2$stunted)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.stuntedSpeciesnoncorrectedchi2MicrobiomeInsights.csv")
      
 
      #### On Species level no correction for stunting, Bangui only ####
      dfcleanfeces_s_B<- tax_glom(dfcleanfecesB, "Species")
      fecesprevalence_s <- dfcleanfeces_s_B %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, stunted!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.stuntedSpeciesnoncorrectedchi2MicrobiomeInsightsBangui.csv")
      
      table(df_chi2[, 12], meta_chi2$stunted)
      
      #### On Species level no correction for stunting, Tana only Nothing is significantly associated ####
      dfcleanfeces_s_A<- tax_glom(dfcleanfecesA, "Species")
      fecesprevalence_s <- dfcleanfeces_s_A %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, stunted!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig) 
      
      table(df_chi2[, 11], meta_chi2$stunted)
      
      
      #### On Species level no correction for calpro ####
      dfcleanfeces_s<- tax_glom(dfcleanfeces, "Species")
      fecesprevalence_s <- dfcleanfeces_s %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, calprotectinelevel!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.calprotectinelevelSpeciesnoncorrectedchi2MicrobiomeInsights.csv")
      
      
      #### On Species level no correction for calpro, Bangui only ####
      dfcleanfeces_s_B<- tax_glom(dfcleanfecesB, "Species")
      fecesprevalence_s <- dfcleanfeces_s_B %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, calprotectinelevel!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(dfcleanfeces_s))
      Tax_corr$name<-row.names(tax_table(dfcleanfeces_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.calprotectinelevelSpeciesnoncorrectedchi2MicrobiomeInsightsBangui.csv")
      
      
      #### On Species level no correction for calpro, Tana only  ####
      dfcleanfeces_s_A<- tax_glom(dfcleanfecesA, "Species")
      fecesprevalence_s <- dfcleanfeces_s_A %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, calprotectinelevel!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig) 
      
      
      
      #### On Species level no correction for country ####
      meta_chi2$pays<-as.factor(meta_chi2$pays)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.paysSpeciesnoncorrectedchi2MicrobiomeInsights.csv")
      
      #### On Species level no correction for ageyears ####
      meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsSpeciesnoncorrectedchi2MicrobiomeInsights.csv")
      
      #### On Species level no correction for ageyears, Bangui ####
      dfcleanfeces_s_B<- tax_glom(dfcleanfecesB, "Species")
      fecesprevalence_s <- dfcleanfeces_s_B %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, pays!="") 
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, ageyears!="") 
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      
      meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsSpeciesnoncorrectedchi2MicrobiomeInsightsBangui.csv")
      
      #### On Species level no correction for ageyears, Tana ####
      dfcleanfeces_s_A<- tax_glom(dfcleanfecesB, "Species")
      fecesprevalence_s <- dfcleanfeces_s_A %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, pays!="") 
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, ageyears!="") 
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      
      meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsSpeciesnoncorrectedchi2MicrobiomeInsightsTana.csv")
      
      #### On Species level no correction for SampleType in three levels ####
      allprevalence <- dfclean %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(allprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      allprevalence<-tax_glom(allprevalence, "Species")
      allprevalence2<-subset_samples(allprevalence, stunted!="") # 575 samples remaining
      allprevalence2<-subset_samples(allprevalence2, pays!="") # 575 samples remaining, no additional ones taken out
      allprevalence2 = filter_taxa(allprevalence2, function(x) sum(x) > 0, TRUE)
      
      
      df_chi2 <- as.matrix(t(otu_table(allprevalence2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(allprevalence2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$SampleType), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(allprevalence2))
      Tax_corr$name<-row.names(tax_table(allprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsSampleTypeSpeciesnoncorrectedchi2MicrobiomeInsightsfulldataset.csv")
      
      
      
  ####  no correction for SampleType in three levels, stunted only ####
      allprevalence <- dfclean %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(allprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      allprevalence<-tax_glom(allprevalence, "Species")
      
      allprevalence2<-subset_samples(allprevalence, stunted!="") # 575 samples remaining
      allprevalence2<-subset_samples(allprevalence2, pays!="") # 575 samples remaining, no additional ones taken out
      allprevalence2 = filter_taxa(allprevalence2, function(x) sum(x) > 0, TRUE)
      allprevalence2<-subset_samples(allprevalence2, stunted=="stunted")
      
      df_chi2 <- as.matrix(t(otu_table(allprevalence2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(allprevalence2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$SampleType), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(allprevalence2))
      Tax_corr$name<-row.names(tax_table(allprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsSampleTypeSpeciesnoncorrectedchi2MicrobiomeInsightsstuntedonly.csv")
      
      
      #### On Species level no correction for SampleType in three levels, shared samples only ####
      allprevalence <- dfclean %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(allprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      allprevalence<-tax_glom(allprevalence, "Species")
      
      allprevalence2<-subset_samples(allprevalence, stunted!="") # 575 samples remaining
      allprevalence2<-subset_samples(allprevalence2, pays!="") # 575 samples remaining, no additional ones taken out
      allprevalence2 = filter_taxa(allprevalence2, function(x) sum(x) > 0, TRUE)
      allprevalence2<-subset_samples(allprevalence2, sample_names(allprevalence2) %in% sample_names(dfrar5000_rest2))
      table(sample_data(allprevalence2)$SampleType)
      
      
      df_chi2 <- as.matrix(t(otu_table(allprevalence2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(allprevalence2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$SampleType), simulate.p.value = TRUE)$p.value)
      Speciesnames<-colnames(df_chi2)
      chi2results = data.frame(Speciesnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(allprevalence2))
      Tax_corr$name<-row.names(tax_table(allprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsSampleTypeSpeciesnoncorrectedchi2MicrobiomeInsightsreduceddataset.csv")
      
      
      
      #### On Genus level no correction for stunting ####
      dfcleanfeces_g<- tax_glom(dfcleanfeces, "Genus")
      fecesprevalence_g <- dfcleanfeces_g %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_g))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_g2<-subset_samples(fecesprevalence_g, stunted!="") # 575 samples remaining
      fecesprevalence_g2<-subset_samples(fecesprevalence_g2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_g2 = filter_taxa(fecesprevalence_g2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_g2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_g2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(dfcleanfeces_g))
      Tax_corr$name<-row.names(tax_table(dfcleanfeces_g))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.stuntedGenusnoncorrectedchi2MicrobiomeInsights.csv")
      
      
      
      
      #### On Genus level no correction for stunting, Bangui only ####
      dfcleanfeces_s_B<- tax_glom(dfcleanfecesB, "Genus")
      fecesprevalence_s <- dfcleanfeces_s_B %>% #this produces prevalence "counts" for each Genus, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, stunted!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.stuntedGenusnoncorrectedchi2MicrobiomeInsightsBangui.csv")
      
      
      #### On Genus level no correction for stunting, Tana only Nothing is significantly associated ####
      dfcleanfeces_s_A<- tax_glom(dfcleanfecesA, "Genus")
      fecesprevalence_s <- dfcleanfeces_s_A %>% #this produces prevalence "counts" for each Genus, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, stunted!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig) 
      
      
      
      #### On Genus level no correction for country ####
      dfcleanfeces_g<- tax_glom(dfcleanfeces, "Genus")
      fecesprevalence_g <- dfcleanfeces_g %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_g))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_g2<-subset_samples(fecesprevalence_g, stunted!="") # 575 samples remaining
      fecesprevalence_g2<-subset_samples(fecesprevalence_g2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_g2 = filter_taxa(fecesprevalence_g2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_g2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_g2)) #take metadatameta_chi2$pays<-as.factor(meta_chi2$pays)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(dfcleanfeces_g))
      Tax_corr$name<-row.names(tax_table(dfcleanfeces_g))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.paysGenusnoncorrectedchi2MicrobiomeInsights.csv")
      
      #### On Genus level no correction for ageyears ####
      meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(dfcleanfeces_g))
      Tax_corr$name<-row.names(tax_table(dfcleanfeces_g))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsGenusnoncorrectedchi2MicrobiomeInsights.csv")
      
      
      #### On Genus level no correction for ageyears, Bangui ####
      dfcleanfeces_s_B<- tax_glom(dfcleanfecesB, "Genus")
      fecesprevalence_s <- dfcleanfeces_s_B %>% #this produces prevalence "counts" for each Genus, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, pays!="") 
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, ageyears!="") 
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      
      meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsGenusnoncorrectedchi2MicrobiomeInsightsBangui.csv")
      
      #### On Genus level no correction for ageyears, Tana ####
      dfcleanfeces_s_A<- tax_glom(dfcleanfecesB, "Genus")
      fecesprevalence_s <- dfcleanfeces_s_A %>% #this produces prevalence "counts" for each Genus, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, pays!="") 
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, ageyears!="") 
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      
      meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsGenusnoncorrectedchi2MicrobiomeInsightsTana.csv")
      
      
      #### On Genus level no correction for calpro ####
      dfcleanfeces_s<- tax_glom(dfcleanfeces, "Genus")
      fecesprevalence_s <- dfcleanfeces_s %>% #this produces prevalence "counts" for each Genus, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, calprotectinelevel!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.calprotectinelevelGenusnoncorrectedchi2MicrobiomeInsights.csv")
      
      
      #### On Genus level no correction for calpro, Bangui only Nothing is significantly associated ####
      dfcleanfeces_s_B<- tax_glom(dfcleanfecesB, "Genus")
      fecesprevalence_s <- dfcleanfeces_s_B %>% #this produces prevalence "counts" for each Genus, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, calprotectinelevel!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      View(chi2results_sig)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.calprotectinelevelGenusnoncorrectedchi2MicrobiomeInsightsBangui.csv")
      
      
      #### On Genus level no correction for calpro, Tana only Nothing is significantly associated ####
      dfcleanfeces_s_A<- tax_glom(dfcleanfecesA, "Genus")
      fecesprevalence_s <- dfcleanfeces_s_A %>% #this produces prevalence "counts" for each Genus, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, calprotectinelevel!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$calprotectinelevel<-as.factor(meta_chi2$calprotectinelevel)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$calprotectinelevel), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s2))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig) 
      
      
      
      #### On Genus level no correction for anemie2 ####
      meta_chi2$anemia<-as.factor(meta_chi2$anemie2)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(dfcleanfeces_g))
      Tax_corr$name<-row.names(tax_table(dfcleanfeces_g))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsGenusnoncorrectedchi2ITS.csv")
      
      
      #### On Genus level no correction for anemie2, Bangui ####
      dfcleanfeces_s_B<- tax_glom(dfcleanfecesB, "Genus")
      fecesprevalence_s <- dfcleanfeces_s_B %>% #this produces prevalence "counts" for each Genus, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, pays!="") 
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, anemie2!="") 
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      
      meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsGenusnoncorrectedchi2ITS.csv")
      
      #### On Genus level no correction for anemie2, Tana ####
      dfcleanfeces_s_A<- tax_glom(dfcleanfecesB, "Genus")
      fecesprevalence_s <- dfcleanfeces_s_A %>% #this produces prevalence "counts" for each Genus, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, pays!="") 
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, anemie2!="") 
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      
      meta_chi2$anemie2<-as.factor(meta_chi2$anemie2)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$anemie2), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsGenusnoncorrectedchi2ITSTana.csv")
      
      
      #### On Genus level no correction for SampleType in three levels ####
      allprevalence <- dfclean %>% #this produces prevalence "counts" for each Genus, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(allprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      allprevalence<-tax_glom(allprevalence, "Genus")
      allprevalence2<-subset_samples(allprevalence, stunted!="") # 575 samples remaining
      allprevalence2<-subset_samples(allprevalence2, pays!="") # 575 samples remaining, no additional ones taken out
      allprevalence2 = filter_taxa(allprevalence2, function(x) sum(x) > 0, TRUE)
      
      
      df_chi2 <- as.matrix(t(otu_table(allprevalence2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(allprevalence2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$SampleType), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(allprevalence2))
      Tax_corr$name<-row.names(tax_table(allprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsSampleTypeGenusnoncorrectedchi2MicrobiomeInsightsfulldataset.csv")
      
      
      
      #### On Genus level no correction for SampleType in three levels, stunted only ####
      allprevalence <- dfclean %>% #this produces prevalence "counts" for each Genus, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(allprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      allprevalence<-tax_glom(allprevalence, "Genus")
      
      allprevalence2<-subset_samples(allprevalence, stunted!="") # 575 samples remaining
      allprevalence2<-subset_samples(allprevalence2, pays!="") # 575 samples remaining, no additional ones taken out
      allprevalence2 = filter_taxa(allprevalence2, function(x) sum(x) > 0, TRUE)
      allprevalence2<-subset_samples(allprevalence2, stunted=="stunted")
      
      df_chi2 <- as.matrix(t(otu_table(allprevalence2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(allprevalence2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$SampleType), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(allprevalence2))
      Tax_corr$name<-row.names(tax_table(allprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsSampleTypeGenusnoncorrectedchi2MicrobiomeInsightsstuntedonly.csv")
      
      
      #### On Genus level no correction for SampleType in three levels, shared samples only ####
      allprevalence <- dfclean %>% #this produces prevalence "counts" for each Genus, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(allprevalence))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      allprevalence<-tax_glom(allprevalence, "Genus")
      
      allprevalence2<-subset_samples(allprevalence, stunted!="") # 575 samples remaining
      allprevalence2<-subset_samples(allprevalence2, pays!="") # 575 samples remaining, no additional ones taken out
      allprevalence2 = filter_taxa(allprevalence2, function(x) sum(x) > 0, TRUE)
      allprevalence2<-subset_samples(allprevalence2, sample_names(allprevalence2) %in% sample_names(dfrar5000_rest2))
      table(sample_data(allprevalence2)$SampleType)
      
      
      df_chi2 <- as.matrix(t(otu_table(allprevalence2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(allprevalence2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$SampleType), simulate.p.value = TRUE)$p.value)
      Genusnames<-colnames(df_chi2)
      chi2results = data.frame(Genusnames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(allprevalence2))
      Tax_corr$name<-row.names(tax_table(allprevalence2))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsSampleTypeGenusnoncorrectedchi2MicrobiomeInsightsreduceddataset.csv")
      
      
      
      #### On Family level no correction for stunting ####
      dfcleanfeces_f<- tax_glom(dfcleanfeces, "Family")
      fecesprevalence_f <- dfcleanfeces_f %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_f))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_f2<-subset_samples(fecesprevalence_f, stunted!="") # 575 samples remaining
      fecesprevalence_f2<-subset_samples(fecesprevalence_f2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_f2 = filter_taxa(fecesprevalence_f2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_f2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_f2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
      Familynames<-colnames(df_chi2)
      chi2results = data.frame(Familynames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(dfcleanfeces_f))
      Tax_corr$name<-row.names(tax_table(dfcleanfeces_f))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.stuntedFamilynoncorrectedchi2MicrobiomeInsights.csv")
      
      
      #### On Family level no correction for stunting, Bangui only ####
      dfcleanfeces_s_B<- tax_glom(dfcleanfecesB, "Family")
      fecesprevalence_s <- dfcleanfeces_s_B %>% #this produces prevalence "counts" for each Family, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, stunted!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
      Familynames<-colnames(df_chi2)
      chi2results = data.frame(Familynames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.stuntedFamilynoncorrectedchi2MicrobiomeInsightsBangui.csv")
      
      
      #### On Family level no correction for stunting, Tana only Nothing is significantly associated ####
      dfcleanfeces_s_A<- tax_glom(dfcleanfecesA, "Family")
      fecesprevalence_s <- dfcleanfeces_s_A %>% #this produces prevalence "counts" for each Family, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, stunted!="") # 575 samples remaining
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
      Familynames<-colnames(df_chi2)
      chi2results = data.frame(Familynames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig) 
      
      
      
      #### On Family level no correction for country ####
      meta_chi2$pays<-as.factor(meta_chi2$pays)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value)
      Familynames<-colnames(df_chi2)
      chi2results = data.frame(Familynames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(dfcleanfeces_f))
      Tax_corr$name<-row.names(tax_table(dfcleanfeces_f))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.paysFamilynoncorrectedchi2MicrobiomeInsights.csv")
      
      #### On Family level no correction for ageyears ####
      meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
      Familynames<-colnames(df_chi2)
      chi2results = data.frame(Familynames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(dfcleanfeces_f))
      Tax_corr$name<-row.names(tax_table(dfcleanfeces_f))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsFamilynoncorrectedchi2MicrobiomeInsights.csv")
      
      
      
      #### On Family level no correction for ageyears, Bangui ####
      dfcleanfeces_s_B<- tax_glom(dfcleanfecesB, "Family")
      fecesprevalence_s <- dfcleanfeces_s_B %>% #this produces prevalence "counts" for each Family, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, pays!="") 
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, ageyears!="") 
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      
      meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
      Familynames<-colnames(df_chi2)
      chi2results = data.frame(Familynames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsFamilynoncorrectedchi2MicrobiomeInsightsBangui.csv")
      
      #### On Family level no correction for ageyears, Tana ####
      dfcleanfeces_s_A<- tax_glom(dfcleanfecesB, "Family")
      fecesprevalence_s <- dfcleanfeces_s_A %>% #this produces prevalence "counts" for each Family, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, pays!="") 
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, ageyears!="") 
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
      
      meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
      Familynames<-colnames(df_chi2)
      chi2results = data.frame(Familynames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(fecesprevalence_s))
      Tax_corr$name<-row.names(tax_table(fecesprevalence_s))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsFamilynoncorrectedchi2MicrobiomeInsightsTana.csv")
      
      
      #### On Order level no correction for stunting ####
      dfcleanfeces_o<- tax_glom(dfcleanfeces, "Order")
      fecesprevalence_o <- dfcleanfeces_o %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_o))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_o2<-subset_samples(fecesprevalence_o, stunted!="") # 575 samples remaining
      fecesprevalence_o2<-subset_samples(fecesprevalence_o2, pays!="") # 575 samples remaining, no additional ones taken out
      fecesprevalence_o2 = filter_taxa(fecesprevalence_o2, function(x) sum(x) > 0, TRUE)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_o2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_o2)) #take metadata
      meta_chi2$stunted<-as.factor(meta_chi2$stunted)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$stunted), simulate.p.value = TRUE)$p.value)
      Ordernames<-colnames(df_chi2)
      chi2results = data.frame(Ordernames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(dfcleanfeces_o))
      Tax_corr$name<-row.names(tax_table(dfcleanfeces_o))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.stuntedOrdernoncorrectedchi2MicrobiomeInsights.csv")
      
      
      #### On Order level no correction for country ####
      meta_chi2$pays<-as.factor(meta_chi2$pays)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$pays), simulate.p.value = TRUE)$p.value)
      Ordernames<-colnames(df_chi2)
      chi2results = data.frame(Ordernames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(dfcleanfeces_o))
      Tax_corr$name<-row.names(tax_table(dfcleanfeces_o))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.paysOrdernoncorrectedchi2MicrobiomeInsights.csv")
      
      #### On Order level no correction for ageyears ####
      meta_chi2$ageyears<-as.factor(meta_chi2$ageyears)
      
      Chi2.p<-apply(df_chi2, 2, function(x) chisq.test(table(x,meta_chi2$ageyears), simulate.p.value = TRUE)$p.value)
      Ordernames<-colnames(df_chi2)
      chi2results = data.frame(Ordernames,Chi2.p)
      chi2results$rel.fdr <- p.adjust(chi2results$Chi2.p, method="fdr")
      chi2results$name<-row.names(chi2results)
      
      chi2results_sig <- dplyr::filter(chi2results, Chi2.p<0.05)
      
      # Merge with tax info
      Tax_corr<-data.frame(tax_table(dfcleanfeces_o))
      Tax_corr$name<-row.names(tax_table(dfcleanfeces_o))
      Tax_corr_chi2<-filter(Tax_corr, Tax_corr$name %in% chi2results_sig$name)
      chi2results_sig = merge(chi2results_sig,Tax_corr_chi2, by="name")
      View(chi2results_sig)
      write.csv(chi2results_sig,"Chi2resultsFeces.ageyearsOrdernoncorrectedchi2MicrobiomeInsights.csv")
      
      
      
#### Make logistic regression on the presence/absence of given fungi correcting for confounding factors RECHECK ####       
      #### Species level: logistic model for stunting correcting for covariables NOTHING ASSOCIATED #####
      dfcleanfeces_s<- tax_glom(dfcleanfeces, "Species")
      fecesprevalence_s <- dfcleanfeces_s %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_s))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_s2<-subset_samples(fecesprevalence_s, stunted!="") 
      fecesprevalence_s2<-subset_samples(fecesprevalence_s2, pays!="") 
      fecesprevalence_s2 = filter_taxa(fecesprevalence_s2, function(x) sum(x) > 0, TRUE)
      
      fecesprevalence_s2
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2)) #take metadata
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
      df$run <- meta_chi2$run
      df$totalreads <- meta_chi2$read_count
      
      sampledata= as.data.frame(sample_data(fecesprevalence_s2))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run",  "age", "sexe", "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+  run + sexe + totalreads + calpro + Country , .,  family=binomial))) %>% 
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
      Tax_corr<-data.frame(tax_table(fecesprevalence_s2))
      Tax_corr$variable<-row.names(tax_table(fecesprevalence_s2))
      Tax_corr_log<-filter(Tax_corr, Tax_corr$variable %in% logresults$variable)
      logresults = merge(logresults,Tax_corr_log, by="variable")
      write.csv(logresults,"LogresultsFecespresabs.stuntedSpecieswithinflaMicrobiomeInsights.csv")
      
      #### Species level: logistic model for country correcting for covariables#####
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
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
      df$totalreads <- meta_chi2$read_count
      
      
      sampledata= as.data.frame(sample_data(fecesprevalence_s2))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "sexe", "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, Country!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(Country ~ value+ run+ sexe + stunted + totalreads + calpro , .,  family=binomial))) %>% 
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
      Tax_corr<-data.frame(tax_table(fecesprevalence_s2))
      Tax_corr$variable<-row.names(tax_table(fecesprevalence_s2))
      Tax_corr_log<-filter(Tax_corr, Tax_corr$variable %in% logresults$variable)
      logresults = merge(logresults,Tax_corr_log, by="variable")
      write.csv(logresults,"LogresultsFecespresabs.paysSpecieswithinflaMicrobiomeInsights.csv")
    
      #### Species level: logistic model for stunting correcting for covariables only in Bangui NOTHING ASSOCIATED #####
      fecesprevalence_s2b<-subset_samples(fecesprevalence_s2, pays=="RCA")
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2b))) # take presence absence table
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
      df$run <- meta_chi2$run
      df$totalreads <- meta_chi2$read_count
      
      
      dffb<-subset_samples(fecesprevalence_s2, pays=="RCA")
      sampledata= as.data.frame(sample_data(dffb))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "sexe", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      long=filter(long, calpro!="") ## keep only the ones with valid data 
      long=filter(long, aat!="") ## keep only the ones with valid data 
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run + sexe + totalreads + calpro , .,  family=binomial))) %>% 
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
      Tax_corr<-data.frame(tax_table(dffb))
      Tax_corr$variable<-row.names(tax_table(dffb))
      logresults = merge(logresults,Tax_corr_log, by="variable")
      write.csv(logresults,"LogresultsFecespresabs.stuntedSpecieswithinflaBanguiMicrobiomeInsights.csv")
      
      #### Species level: logistic model for stunting correcting for covariables only in Tana NOTHING ASSOCIATED #####
      fecesprevalence_s2a<-subset_samples(fecesprevalence_s2, pays=="Madagascar")
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2a))) # take presence absence table
      meta_chi2 <- data.frame(sample_data(fecesprevalence_s2a)) #take metadata
      data1 <- df_chi2
      
      dffa<-subset_samples(fecesprevalence_s2, pays=="Madagascar")
      sampledata= as.data.frame(sample_data(dffa))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      
      
      #add other categorical factors, etc. 
      df <- data.frame(data1)
      dim(df)
      df$age <- meta_chi2$age
      df$run <- meta_chi2$run
      df$sexe <- meta_chi2$sexe
      df$stunted <- meta_chi2$stunted
      df$calpro <- meta_chi2$calprotectinelevel
      df$aat <- meta_chi2$alphaantitrypsinlevel
      df$totalreads<-meta_chi2$read_count
      
           
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run", "age", "sexe", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+ run + sexe + totalreads + calpro , .,  family=binomial))) %>% 
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
      Tax_corr<-data.frame(tax_table(dffa))
      Tax_corr$variable<-row.names(tax_table(dffa))
      Tax_corr_log<-filter(Tax_corr, Tax_corr$variable %in% logresults$variable)
      Tax_corr <- Tax_corr[order(Tax_corr$variable), ]
      logresults <- logresults[order(logresults$variable), ]
      logresults = data.frame(logresults,Tax_corr_log)
      write.csv(logresults,"LogresultsFecespresabs.stuntedSpecieswithinflaTanaMicrobiomeInsights.csv")
      
      
      #### Species level: logistic model for age (more or less than three years) correcting for covariables Aspergillus is associated #####
      #create a variable age in years which has only two modalities
      sample_data(fecesprevalence_s2)$ageyears2<-cut(as.numeric(as.character(sample_data(fecesprevalence_s2)$age)), c(24,36,61), include.lowest = TRUE, right=TRUE, dig.lab=5, ordered_result = TRUE)
      which(is.na(sample_data(fecesprevalence_s2)$ageyears2)) # the controls have no assocaited age year
      levels(sample_data(fecesprevalence_s2)$ageyears2) <- c("2-3 years", "3-5 years")
      levels(sample_data(fecesprevalence_s2)$ageyears2)
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_s2))) # take presence absence table
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
      df$totalreads <- meta_chi2$read_count
      
      
      sampledata= as.data.frame(sample_data(fecesprevalence_s2))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run",  "age", "sexe", "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, age!="") ## keep only the ones with valid data for age
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(age ~ value+ Country + run + sexe + stunted + totalreads + calpro , .,  family=binomial))) %>% 
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
      Tax_corr<-data.frame(tax_table(fecesprevalence_s2))
      Tax_corr$variable<-row.names(tax_table(fecesprevalence_s2))
      Tax_corr_log<-filter(Tax_corr, Tax_corr$variable %in% logresults$variable)
      Tax_corr <- Tax_corr[order(Tax_corr$variable), ]
      logresults <- logresults[order(logresults$variable), ]
      logresults = data.frame(logresults,Tax_corr_log)
      write.csv(logresults,"LogresultsFecespresabs.agetwocatSpecieswithinflaMicrobiomeInsights.csv")
      
      #### Genus level: logistic model for stunting correcting for covariables None associated #####
      dfcleanfeces_g<- tax_glom(dfcleanfeces, "Genus")
      fecesprevalence_g <- dfcleanfeces_g %>% #this produces prevalence "counts" for each species, but not percentages
        transform_sample_counts(fun = prevalence) 
      
      View(head(otu_table(fecesprevalence_g))) ## ok, this is the sample table with prevalence as 1 if at least 10 seqs in sample
      
      fecesprevalence_g2<-subset_samples(fecesprevalence_g, stunted!="") 
      fecesprevalence_g2<-subset_samples(fecesprevalence_g2, pays!="") 
      fecesprevalence_g2 = filter_taxa(fecesprevalence_g2, function(x) sum(x) > 0, TRUE)
      
      fecesprevalence_g2
      
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_g2))) # take presence absence table
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
      df$run <- meta_chi2$run
      sampledata= as.data.frame(sample_data(fecesprevalence_g2))
      sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
      df$totalreads<-sampledata2$read_count
      
      library(broom)
      library(dplyr)
      long = melt(df, id.vars = c("stunted", "run",  "age", "sexe", "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
      long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
      
      logresults<- long %>% 
        group_by(variable) %>% 
        do(tidy(glm(stunted ~ value+  run + sexe + totalreads + calpro, .,  family=binomial))) %>% 
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
      Tax_corr<-data.frame(tax_table(fecesprevalence_g2))
      Tax_corr$variable<-row.names(tax_table(fecesprevalence_g2))
      Tax_corr_log<-filter(Tax_corr, Tax_corr$variable %in% logresults$variable)
      logresults = merge(logresults,Tax_corr_log, by="variable")
      write.csv(logresults,"LogresultsFecespresabs.stuntedGenuswithinflaMicrobiomeInsights.csv")
      
      
      
      #### Genus level: logistic model for country correcting for covariables#####
      df_chi2 <- as.matrix(t(otu_table(fecesprevalence_g2))) # take presence absence table
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
      
      sampledata= as.data.frame(sample_data(fecesprevalence_g2))
      sampledata2<-filter(sampledata, row.names(sampledata %in% row.names(df)))
      df$totalreads<-sampledata2$read_count
                          
    library(broom)
     library(dplyr)
     long = melt(df, id.vars = c("stunted", "run", "age", "sexe", "Country", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion tes
    long=filter(long, Country!="") ## keep only the ones with valid data for stunted
                          
  logresults<- long %>% 
   group_by(variable) %>% 
    do(tidy(glm(Country ~ value+ run+ sexe + stunted + totalreads + calpro, .,  family=binomial))) %>% 
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
  Tax_corr<-data.frame(tax_table(fecesprevalence_g2))
   Tax_corr$variable<-row.names(tax_table(fecesprevalence_g2))
   Tax_corr_log<-filter(Tax_corr, Tax_corr$variable %in% logresults$variable)
 logresults = merge(logresults,Tax_corr_log, by="variable")
    write.csv(logresults,"LogresultsFecespresabs.paysGenuswithinflaMicrobiomeInsights.csv")
                          
                       
                          
      #### Genus level: logistic model for stunting correcting for covariables only in Bangui  #####
                          fecesprevalence_g2b<-subset_samples(fecesprevalence_g2, pays=="RCA")
                          df_chi2 <- as.matrix(t(otu_table(fecesprevalence_g2b))) # take presence absence table
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
                          df$run <- meta_chi2$run
                          
                          dffb<-subset_samples(fecesprevalence_g2, pays=="RCA")
                          sampledata= as.data.frame(sample_data(dffb))
                          sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
                          df$totalreads<-sampledata2$read_count
                          
                          library(broom)
                          library(dplyr)
                          long = melt(df, id.vars = c("stunted", "run", "age", "sexe", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
                          long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
                          long=filter(long, calpro!="") ## keep only the ones with valid data 
                          long=filter(long, aat!="") ## keep only the ones with valid data 
                          
                          logresults<- long %>% 
                            group_by(variable) %>% 
                            do(tidy(glm(stunted ~ value+  run + sexe + totalreads + calpro , .,  family=binomial))) %>% 
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
                          Tax_corr<-data.frame(tax_table(fecesprevalence_g2))
                          Tax_corr$variable<-row.names(tax_table(fecesprevalence_g2))
                          logresults = merge(logresults,Tax_corr_log, by="variable")
                          write.csv(logresults,"LogresultsFecespresabs.stuntedGenuswithinflaBanguiMicrobiomeInsights.csv")
                          
                          
      #### Genus level: logistic model for stunting correcting for covariables only in Tana NOTHING ASSOCIATED #####
                          fecesprevalence_g2a<-subset_samples(fecesprevalence_g2, pays=="Madagascar")
                          df_chi2 <- as.matrix(t(otu_table(fecesprevalence_g2a))) # take presence absence table
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
                          
                          dffa<-subset_samples(fecesprevalence_g2, pays=="Madagascar")
                          sampledata= as.data.frame(sample_data(dffa))
                          sampledata2<-sampledata[row.names(sampledata) %in% row.names(df), ]
                          df$totalreads<-sampledata2$read_count
                          
                          library(broom)
                          library(dplyr)
                          long = melt(df, id.vars = c("stunted", "run", "age", "sexe", "totalreads", "aat", "calpro")) ## use here the variables that showed to be associated in dispersion test
                          long=filter(long, stunted!="") ## keep only the ones with valid data for stunted
                          
                          logresults<- long %>% 
                            group_by(variable) %>% 
                            do(tidy(glm(stunted ~ value+ run + sexe + totalreads + calpro, .,  family=binomial))) %>% 
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
                          Tax_corr<-data.frame(tax_table(fecesprevalence_g2))
                          Tax_corr$variable<-row.names(tax_table(fecesprevalence_g2))
                          Tax_corr_log<-filter(Tax_corr, Tax_corr$variable %in% logresults$variable)
                          Tax_corr <- Tax_corr[order(Tax_corr$variable), ]
                          logresults <- logresults[order(logresults$variable), ]
                          logresults = data.frame(logresults,Tax_corr_log)
                          write.csv(logresults,"LogresultsFecespresabs.stuntedGenuswithinflaTanaMicrobiomeInsights.csv")
                          
                          