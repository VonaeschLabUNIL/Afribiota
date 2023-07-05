###   ITS analysis using Dada2 out of R full dataset Afribiota

  ## based on the pipeline described here: https://benjjneb.github.io/dada2/ITS_workflow.html
  ## adapted by Pascale Vonaesch in October 2019

####Environment Setup####
theme_set(theme_bw())
setwd("xxxx") # CHANGE to your home directory


#### load packages ####
          source("https://bioconductor.org/biocLite.R")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.9")
          packageVersion("dada2")
          library(ShortRead)
          packageVersion("ShortRead")
          library(Biostrings)
          packageVersion("Biostrings")
          library(ggplot2)
          library(phyloseq)
          library(tidyverse)
          library(ShortRead)
          library(dada2)

#### go to the files you need and look for primer read-through ####

          path <- "xxxxx" # CHANGE ME to your path
          list.files(path)
          
          
          # Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
          fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
          fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
          
          sample.names <- sapply(strsplit(basename(fnFs), "_L"), `[`, 1)
          
          # look at quality of your reads
          plotQualityProfile(fnFs[1:15]) # these are the control samples, water
          plotQualityProfile(fnFs[16:30]) # the real samples look much better
          plotQualityProfile(fnFs[1:2])
          plotQualityProfile(fnFs[3:4]) # values start  decreasing after 200 bp, some already before.  but remain largely above the common quality threshold of 30: hence suggested to trim after 250 bp for the forward reads
          
          
          plotQualityProfile(fnRs[15:30])
          plotQualityProfile(fnRs[1:2])
          plotQualityProfile(fnRs[3:4])
          plotQualityProfile(fnRs[5:6])# overall worse seuqencing than for forward primers. values start rapidly decreasing after ca. 150 bp, but remain largely above 30 until 200 bp? Trim at 200 bp suggested, ev. adapt later if too shaby results?
          
          # control that you did not seuqence by mistake the reverse primer
          
              #The BITS3 (forward) and B58S3 (reverse) primers were used to amplify this dataset. We record the DNA sequences, including ambiguous nucleotides, for those primers.
          
          FWD <- "GCATCGATGAAGAACGCAGC"  ## forward primer
          REV <- "TCCTCCGCTTATTGATATGC"  ## reverse primer 
          
          
          allOrients <- function(primer) {
            # Create all orientations of the input sequence
            require(Biostrings)
            dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
            orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                         RevComp = reverseComplement(dna))
            return(sapply(orients, toString))  # Convert back to character vector
          }
          FWD.orients <- allOrients(FWD)
          REV.orients <- allOrients(REV)
          FWD.orients
          
          # now filter for N's
          fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
          fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
          
          sample.namesF <- sapply(strsplit(basename(fnFs), "_L"), `[`, 1)
          sample.namesF #1160 samples
          sample.namesR <- sapply(strsplit(basename(fnRs), "_L"), `[`, 1)
          sample.namesR # 1157 samples
          FnotR <- setdiff(sample.namesF, sample.namesR) 
          FnotR #there are three files with now reverse that we have to delete manually in the original folder
          
          RnotF<-setdiff(sample.namesR, sample.namesF) 
          RnotF
          
          fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
          fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
          fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
          fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
          
          filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
          
          
          
          # now see how often the primer appears in the forward and reverse read
          primerHits <- function(primer, fn) {
            # Counts number of reads in which the primer is found
            nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
            return(sum(nhits > 0))
          }
          rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
                FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
                REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
                REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
              # comment: the primers have already been chopped-off by Microbiome Insights and the forward primer therefore does not appear! The reverse primer seems to read into the forward reads
          
          REV <- REV.orients[["RevComp"]] ## there was a problem here as it was in wrong order!, need to restart here!
          allOrients <- function(primer) {
            # Create all orientations of the input sequence
            require(Biostrings)
            dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
            orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                         RevComp = reverseComplement(dna))
            return(sapply(orients, toString))  # Convert back to character vector
          }
          FWD.orients <- allOrients(FWD)
          REV.orients <- allOrients(REV)
          FWD.orients
         
          primerHits <- function(primer, fn) {
            # Counts number of reads in which the primer is found
            nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
            return(sum(nhits > 0))
          }
          rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
                FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
                REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
                REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
          

          # RERUN install cutadapt to get rid of the primers. Do this with cutadapt dowloaded from  here using pip install and conda: https://cutadapt.readthedocs.io/en/stable/installation.html#quick-installation
          
          cutadapt <- "/home/wslnx/.local/bin/cutadapt" 
          system2(cutadapt, args = "--version") # Run shell commands from R
          
          # now run cutadapt
          path.cut <- file.path(path, "cutadapt")
          if(!dir.exists(path.cut)) dir.create(path.cut)
          fnFs.cut <- file.path(path.cut, basename(fnFs))
          fnRs.cut <- file.path(path.cut, basename(fnRs))
          
          FWD.RC <- dada2:::rc(FWD)
          REV.RC <- dada2:::rc(REV)
          # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
          R1.flags <- paste("-g", FWD, "-a", REV.RC) 
          # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
          R2.flags <- paste("-G", REV, "-A", FWD.RC) 
          # Run Cutadapt Carefull: any spaces in the path will make it crash here. If you change the folder name, you need to restart from the begining (it stores the path somewhere and does not work anymore if not)
          for(i in seq_along(fnFs)) {
            system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
            "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
            fnFs.filtN[i], fnRs.filtN[i]))} # input files
          
          
          # check if it worked out well 
          rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
                FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
                REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
                REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
          
#### Now prepare to feed into Dada2 pipeline and control quality before moving ahead ####
          
          # Forward and reverse fastq filenames have the format:
          cutFs <- sort(list.files(path.cut, pattern = "*R1_001.fastq", full.names = TRUE))
          cutRs <- sort(list.files(path.cut, pattern = "*R2_001.fastq", full.names = TRUE))
          
          # Extract sample names, assuming filenames have format:
          get.sample.name <- function(fname) strsplit(basename(fname), "_L")[[1]][1] ## careful, we do not get what is after the _; do not know how to do this!!
          sample.names <- unname(sapply(cutFs, get.sample.name))
          head(sample.names)
          
          # Inspect and read now the quality profiles 
          plotQualityProfile(cutFs[1:2])
          plotQualityProfile(cutRs[1:15])
          
          # now filter and trim
          filtFs <- file.path(path.cut, "filtered", basename(cutFs))
          filtRs <- file.path(path.cut, "filtered", basename(cutRs))
          
          out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(6, 8), 
                               truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
          retained <- as.data.frame(out)
          retained$percentage_retained <- retained$reads.out/retained$reads.in*100
          View(retained)
          head(out)
          dim(out)
          filtFs[5]
          out[5]
          
          samples_to_keep <- as.numeric(out[,"reads.out"]) > 0 # Or other cutoff #will need to be adjusted per-experiment
          #need to redefine objects from above to exclude samples you excluded in this section
          cutFs <- sort(list.files(path.cut, pattern="_R1_001.fastq", full.names = TRUE))[samples_to_keep]
          cutRs <- sort(list.files(path.cut, pattern="_R2_001.fastq", full.names = TRUE))[samples_to_keep]
        
          path.filtered <- file.path(path.cut, "filtered")
          filtFs <- sort(list.files(path.filtered, pattern="_R1_001.fastq", full.names = TRUE))[samples_to_keep]
          filtRs <- sort(list.files(path.filtered, pattern="_R2_001.fastq", full.names = TRUE))[samples_to_keep]
          out <- out[samples_to_keep,]
          
          
          # learn the error rates
          errF <- learnErrors(filtFs, multithread=TRUE)
          errR <- learnErrors(filtRs, multithread = TRUE)

          
          # Visualize the estimated error rates
          plotErrors(errF, nominalQ=TRUE)
          plotErrors(errR, nominalQ=TRUE)
          
          # Dereplicate your sequences
          derepFs <- derepFastq(filtFs, verbose=TRUE)
          derepRs <- derepFastq(filtRs, verbose=TRUE)
          
          # Name the derep-class objects by the sample names, we first need to re-adress the right sample name!
          get.sample.name <- function(fname) strsplit(basename(fname), "_L")[[1]][1] ## careful, we do not get what is after the _; do not know how to do this!!
          sample.names <- unname(sapply(filtFs, get.sample.name))
          names(derepFs) <- sample.names
          names(derepRs) <- sample.names
          
          # Apply now the core sample inference algorithm of Dada2
          dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
          dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
          
          dadaFs[[1]] ## 120 real sample sequences inferred
          dadaRs[[1]] ## 100 real sample sequences inferred
          
          # save the R objects as to lead them on a more powerfull computer
          saveRDS(derepFs, file = "derepFs.rds")
          saveRDS(derepRs, file = "derepRs.rds")
          saveRDS(dadaFs, file = "dadaFs.rds")
          saveRDS(dadaRs, file = "dadaR.rds")
          saveRDS(out, file = "out.rds")
          
          #samples_to_keep <- as.numeric(out[,"reads.out"]) > 100 #simple method
          getN <- function(x) sum(getUniques(x)) #keeping track of read retention, number of unique sequences after ASV inference
          track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
          samples_to_keep <- track[,4] > 50 #method that accounts for dereplication/ASVs left after inference
          length(which(track[,4] > 50)) #how many samples retained with this filter in place
          
         #create phyloseq otu_table
          mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=TRUE)
          
          # Inspect the merger data.frame from the first sample
          head(mergers[[1]])
          
          # construct sequence table 
          seqtab <- makeSequenceTable(mergers)
          dim(seqtab)
          
          ####View Sequence Length Distribution Post-Merging####
          #most useful with merged data. this plot will not show you much for forward reads only, which should have a nearly uniform length distribution.
          length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
          plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot
          
            otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)
          
          #some metrics from the sequence table
          otu_pres_abs <- otus
          otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
          otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV
          View(otu_pres_abs_rowsums)
          length(otu_pres_abs_rowsums) #how many ASVs 
          length(which(otu_pres_abs_rowsums == 1)) #how many ASVs only present in one sample 
          length(which(otu_pres_abs_rowsums < 5)) #how many ASVs present in less than five samples 
          length(which(otu_pres_abs_rowsums < 10)) #how many ASVs present in less than ten samples 
          
          #what are the counts of each ASV
          otu_rowsums <- rowSums(otus) #raw counts per ASV
          otu_singleton_rowsums <- as.data.frame(otu_rowsums[which(otu_pres_abs_rowsums == 1)]) #raw read counts in ASVs only presesnt in one sample
          hist(otu_singleton_rowsums[,1], breaks=10000, xlim = c(0,2000), xlab="# Reads in ASV", main="Read number of ASV in Afribiota ITS dataset")  #plot of above          		  length(which(otu_singleton_rowsums <= 50)) #how many are there with less than N reads #data exploration, not really going to inform ASV removal if we are 			  using relative abundances to filter
          length(otu_rowsums)
          
          #IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance     		  table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
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
          dim(otus_filt) #how many ASVs did you retain? 8243 OTUs were retained
          seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline
          
          ## now remove chimeras
          seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="consensus", multithread=TRUE, verbose=TRUE)
          table(nchar(getSequences(seqtab.nosingletons.nochim)))
          
          ##track reads through the pipeline 
          getN <- function(x) sum(getUniques(x))
          track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                                 getN), rowSums(seqtab.nosingletons.nochim))
          # If processing a single sample, remove the sapply calls: e.g. replace
          # sapply(dadaFs, getN) with getN(dadaFs)
          colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                               "nonchim")
          rownames(track) <- sample.names
          head(track)
          write.table(track, "trackingITS2.txt", sep="\t", quote=F, row.names=T, col.names=T)
          
          
          # now assign taxonomy
          unite.ref <- "~/Desktop/sh_general_release_dynamic_s_01.12.2017.fasta"  # CHANGE ME to location on your machine
          taxa <- assignTaxonomy(seqtab.nosingletons.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
          
          #NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
          unique(taxa[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
          NAs <- is.na(taxa[,1]) #test for NA
          NAs <- which(NAs == TRUE) #get indices of NA values
          taxa[NAs,1] <- "Unassigned" #apply new label to identified indices
          
          write.table(taxa, "taxonomy_table.afribiota_ITS_MI.txt", sep="\t", quote=F, row.names=T, col.names=T)
          write.table(seqtab.nosingletons.nochim, "sequence_table.afribiota_ITS_MI.txt", sep="\t", quote=F, row.names=T, col.names=T)
          
          # now try to re-assign the unassigned by pulling them out and blasting them on the server against NCBI, this has to be done on the server, see specific memo of how this was done
          taxa <- read.csv("taxonomy_table.afribiota_ITS_MI.csv", 
                              header = TRUE)
          
          Subs2<-subset(taxa, (is.na(taxa[,3])))
          View(Subs2)
          row.names(Subs2)<-Subs2$V1
          tobeblasted <-row.names(Subs2)
          tobeblasted
          nrow(Subs2) ## 4668 do not have any assignment at second level
          nrow(taxa)
          100/2625*963  # 36.7% are not assigned at Rank2, but nothing really usefull coming out of blast
          
          write.csv(tobeblasted, "tobeblastedITS.csv")
          
          # import the sequences that have been found to have a match on NCBI: note: there is almost no fungus that could be further identified          
          
          
          