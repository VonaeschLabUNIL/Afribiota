####DADA2 PROCESSING FOR AFRIBIOTA MI 16S SEQUENCES####
## adapted after the original pipeline from the Holmes lab, see here: https://benjjneb.github.io/dada2/tutorial.html
####Libraries####
library(dada2)
library(phyloseq)
library(ShortRead)
library(Biostrings)
library(vegan)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(ape)
library(qualpalr)
library(viridis)

####Environment Setup####
theme_set(theme_bw())
setwd("/home/wslnx/Documents/Pascale/16S") # # CHANGE ME to your working directory

####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "/home/wslnx/Documents/Pascale/16S" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
fnFs <- sort(list.files(path, pattern="R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_001.fastq.gz", full.names = TRUE))
  sample.namesF<- sapply(strsplit(basename(fnFs), "_L"), `[`, 1)
  sample.namesR <- sapply(strsplit(basename(fnRs), "_L"), `[`, 1)
  intersect<-intersect(sample.namesF, sample.namesR)
  lenght(intersect) # for xx samples we have both forward and reverse
  length(fnFs) 
  length(fnRs) 
  sample.namesR [! sample.namesR %in% sample.namesF] # these are in the reverse but not in the forward
  sample.namesF [! sample.namesF %in% sample.namesR] # these are in the forward but not in the reverse
  
sample.names <- sapply(strsplit(basename(fnFs), "_L"), `[`, 1)

####fastq Quality Plots####
plotQualityProfile(fnFs[1:25]) #this plots the quality profiles for each sample, if you have a lot of samples, it's best to look at just a few of them, the plots take a minute or two to generate even only showing 10-20 samples at reasonable sequencing depth.
plotQualityProfile(fnFs[1100:1121]) 
plotQualityProfile(fnRs[1:25])
plotQualityProfile(fnRs[1100:1121])

####trim & filter####
filtFs <- file.path(path, "dada2_filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "dada2_filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#IMPORTANT: if you want to use this object later on you can either re-run this after eliminating samples (time consuming), or just remove them from the data table
#using paramters from massana lab
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,220), maxN=0, maxEE=c(6,8), truncQ=c(2,2), rm.phix=TRUE, verbose=TRUE, matchIDs=TRUE, compress=TRUE, multithread=TRUE) #the more aggressive you are with truncation, the more reads will be retained, but you still need to merge, so consider how much read overlap remains after truncation
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
View(retained)
#the samples beginning with S003K are all about the same quality (with close to 50% of the reads retained), all other samples, with a few exceptions, have nearly 100% read retention. Ask Pascale about this.

####learn error rates####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

####dereplication####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]

####OPTIONAL: remove low-sequence samples before merging####
#a "subscript out of bounds" error at the next step (merging) may indicate that you aren't merging any reads in one or more samples.
#NB, NOT getting this error doesn't necessarily mean that all of your samples end up with more than 0 merged reads, as i found out while processing a large 18s dataset. your guess is as good as mine as to why this error does or does not appear, but filtering out the samples that cause it is necessary for completion of the pipeline.
#samples_to_keep <- as.numeric(out[,"reads.out"]) > 100 #example of simple method used above after the filter and trim step. if you already did this but still got an error when merging, try the steps below
getN <- function(x) sum(getUniques(x)) #keeping track of read retention, number of unique sequences after ASV inference
track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
samples_to_keep <- track[,4] > 100 #your threshold. try different ones to get the lowest one that will work. #this method accounts for dereplication/ASVs left after inference

####merge paired reads####
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####construct sequence table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

####View Sequence Length Distribution Post-Merging####
#most useful with merged data. this plot will not show you much for forward reads only, which should have a uniform length distribution.
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot

####remove low-count singleton ASVs####
#create phyloseq otu_table
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#some metrics from the sequence table
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV
length(otu_pres_abs_rowsums) #how many ASVs
length(which(otu_pres_abs_rowsums == 1)) #how many ASVs only present in one sample

#what are the counts of each ASV
otu_rowsums <- rowSums(otus) #raw counts per ASV
otu_singleton_rowsums <- as.data.frame(otu_rowsums[which(otu_pres_abs_rowsums == 1)]) #raw read counts in ASVs only presesnt in one sample
hist(otu_singleton_rowsums[,1], breaks=500, xlim = c(0,200), xlab="# Reads in ASV") #histogram plot of above
length(which(otu_singleton_rowsums <= 50)) #how many are there with N reads or fewer?

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
dim(otus) # (this should be the same as last command, but the dimensions reversed)
otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step (a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
dim(otus_filt) #how many ASVs did you retain?
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=TRUE, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset

dim(seqtab.nosingletons.nochim)
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low


####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras")
rownames(track) <- sample.names[samples_to_keep]

####save output from sequnce table construction steps####
write.table(track, "read_retention.merged_16s_all.txt", quote=F, row.names=T, col.names=T, sep="\t")
write.table(seqtab.nosingletons.nochim, "sequence_table.afribiota_16s_MI.txt", sep="\t", quote=F, row.names=T, col.names=T)

####assign taxonomy####
#note, this takes ages if you have a large dataset. saving the sequences as a fasta file (with writeFasta) and using QIIME's taxonomy assignment command will save you time, and is only slightly less accurate than the dada2 package's taxonomy assignment function.
taxa <- assignTaxonomy(seqtab.nosingletons.nochim, "/home/wslnx/Desktop/silva_128.16s.99_rep_set.dada2.fa.gz", multithread=TRUE) # CHANGE ME to where your reference taxonomy dataset is located

#NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
unique(taxa[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
NAs <- is.na(taxa[,1]) #test for NA
NAs <- which(NAs == TRUE) #get indices of NA values
taxa[NAs,1] <- "Unassigned" #apply new label to identified indices
colnames(taxa) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Accession")

write.table(taxa, "taxonomy_table.afribiota_16s_MI_all.txt", sep="\t", quote=F, row.names=T, col.names=T)
