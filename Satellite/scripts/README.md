This folder include all the scripts used on the HPC cluster to processes and analyze the data from the Satellite project. 
Project: Satellite - Afribiota
The Satellite project samples include oral, gastric aspirates, duodenal aspirates and stools from stunted and non-stunted children age 2-5 years, from Bangui, Central African Republic.
DNA extraction was done with either promega or qiagen kit. Extracted DNA was send for shotgun metagenomic sequencing 12Gb at Novogene (NVUK2024042951).
Reads were cleaned locally using:
    fastp for quality trimming and adapter removal (v.0.23.4)
Then reads were transferred on the Urblauna cluster for human host reads removal using:
    bowtie2 aligner using the human reference genome (GRCh38, downloaded on the 21st of Nov. 2022) for host reads removal. (v.2.5.1)
Samples:
    - Oral: 39
    - Gastric: 31
    - Duodenum: 15 (remove AGN043DU--> 14 samples)
    - Stool: 44

### Sequencing depth calculations ########
In bowtie2 log folder on the nas /nas/FAC/FBM/DMF/pvonaesc/default/D2c/Simon\ Yersin/Satellite/BioInfo/MetaG/Bowtie2_logs
Ran: 
    echo -e "SampleName\tTotal\tRetained\tRemoved\tRetained%" > bowtie_summary.tsv
    for file in *.log; do
    awk -v filename="$file" '
        BEGIN { sample=filename; sub(/_bowtie2\.log$/, "", sample) }
        /reads; of these:/ { total=$1 }
        /aligned concordantly 0 times/ && !/pairs/ { unaligned=$1 }
        END { printf "%s\t%d\t%d\t%d\t%.2f%%\n", sample, total, unaligned, total - unaligned, (unaligned/total)*100 }
    ' "$file" >> bowtie_summary.tsv
    done
Showing total number of reads, how many were retained after removal of the human reads, how many were removed due to alignment on human reads and the percentage retained.
This worked for the 110 samples that I sent to sequencing to Novogene 
Missing 19 stool samples previously sequenced and cleaned in March 2024. The bowtie2 reports are missing.
04.06.2025
Recleaned raw reads from the 19 remaining stool samples previously sequenced
Obtained bowtie2 logs and recomputed the summary for bowtie2: bowtie_summary.tsv
    Few thousands difference between the two cleaning runs
N.B. except AGN001FE all stool samples retained >99% reads after host reads removal

#############---- Read-based analysis ------###########################################################################################
Taxonomic annotations and abundance profiling of the clean reads is done with MetaPhlAn4 (v. 4.1.1) and ChocoPhlAn database SGB_202307 : metaphlan_profile.sh
The profiles were merged into two tables using: merge_profile.sh
    Table 1: merge relative abundance profile
    Table 2: merge read count profile
    Table 3: merge relative abundance profile with GTDB taxonomy
The resulting merged table was downloaded and analysis done in RStudio.

09.05.2025
Merging profiles with Afribiota profiles (see Afribiota directory) using merge_profile2.sh

The functional potential was profiled using HUMAnN v3.9. All metacyc, go terms, kegg orthologs, pfam, level 4 enzyme and eggnog tables were constructed from the original uniref90 table.
    All tables were normalized both in copy per millions and relative abundance.
    Tables from the different samples were then merged. 
    Output saved on the NAS recherche


StrainPhlAn: tracking individual strains accross a large set of samples
    (Inout: metagenomic samples, outpu: MSA file)
01. Taxonomic profiling with MetaPhlAn4 (01_sams.sh)
    Output: sampleX.sam.bz2 
        Reads mapped against the MetaPhlAn marker database
    Including isolates reads (cleaned reads for WGS short reads)
    MetaPhlAn4 does not suppport the long reads (not included)
02. Extracting sample markers (02_cons_markers.sh)
    Output: sampleX.json.bz2 
        Consensus marker file production using sample2markers.py
    Including short reads isolates consensus markers as well

2025-04-24 Redoing the next steps only for samples (without isolates) nad removing AGN043DU
Previous results saved on nas in /2024_all_samples_isolates

2025-05-06 Redoing the next steps for samples and isolates: 
    Removing isolates: SAT350 (contaminated), SAT622 (contaminated) and isoaltes not belonging to salivarius or parasanguinis

02.2. Run then in terminal in the outdir:
    activate metaphlan environment
    strainphlan -s ./*.json.bz2 --print_clades_only -o . --nproc 8 > ./clades.txt
        Extract the clade names and store them in a txt file
03. Extract marker from database and do MSA (03_msa.sh)
    (Array 1-480 as 480 markers) - extract_markers.py skipped when redone as same .fna file for CladeMarkers
    For markers: t__SGB4706, t__SGB14834, t__SGB19825, t__SGB5273, t__SGB15154, t__SGB15237
        [Error] Phylogeny can not be inferred. Less than 4 samples remained after filtering.
We decided to use the method presented in Feehily et al, 2022 which normalise the best tree for each clade, compute the distance betweem sample using the ape package and 
then using a 0.001 threshold, decide on a sharing event between samples
Therefore, we do not use the scripts 04 and 05 for strainphlan. We do not understand how strain_transmission.py assigned the threshold.

#############---- Upstream analysis - MAGS ------###########################################################################################
MAGs generation pipeline
pre-processing: reads were first copied from the nas storage to the scratch directory in /Satellite/data
    each sample was then moved into it's unique folder in /data/sampleX/
        sbatch standard_scripts/create_dir.sh

01. Metagenomic assembly:
The metagenomic assembly is done with metaSPAdes: 01_metaspades.sh
    2 samples stop during metaspades assembly: increased memory and threads --> add -t and -m in spades command

02. Filtering:
Filter the scaffolds and output the stats: 02_scaffolds_filter.sh

03.1. Indexing:
Index the scaffolds with bowtie2-build: 03_1_index_scaffolds.sh

03.2. Backmapping:
First we split the samples into groups as it is unecessary to have more than 50 backmapping on each samples.
We use the config file (config.txt) to create different batch of job:
    -: Oral samples
    - : Gastric samples + Duodenal samples
    - : Stool samples 
In metaspade dir: ls -d *FE > ~/Satellite/scripts/A_upstream_analysis/config_feces.txt
Downloaded file, opened in Excel to add colum with number and colnames
Save file on cluster and use: sed -e "s/\r//g" /scratch/syersin2/Satellite_scratch/tmp/config_feces.txt > /scratch/syersin2/Satellite_scratch/tmp/config_feces2.txt
To remove \r at the end of each line (otherwise return line in the config file)

Then backmapping is performed with samtools: 03_2_backmapping.sh (FE, SA, GA_DU)

04. Depth:
Depth calculation is done with MetaBAT2 utility script jgi_summarize_bam_contig_depths: 04_depth.sh

05. Binning:
Binning of scaffolds is done with MetaBAT2: 05_binning.sh

06. QC:
Quality control of the MAGs is done with CheckM: 06_checkM.sh

07. Taxonomy annotations:
Using GTDB-TK: 07_gtdktk_mags.sh

08. zipping files
09. copy contigs 

Merging tables and analysis in RStudio
Select HQ MAGs, In RStudio:
Combine checkM and GTDB results to select high quality MAGs and MAGs of interest based on the taxonomy
Copy the high quality MAGs in a new folder HQ_MAGs using a config file: hqmags_config.txt: 08_copy_HQ_mags.sh

We now have high quality MAGs ready for downstream analysis

10. CoverM
Quantification of the MAGs abundance in the sample. Using the dereplicated MAGs from inStrain (strain level).
Tried 2 methods: 
    - Trimmed mean: need further normalization per sequencing depth
    - TPM

#############---- Isolates assembly ------###########################################################################################
Strain were isolated from the different samples and sent for sequencing, either with PacBio Hifi or Illumina short reads.
The sequenced genomes were assembled either with:
    - SPAdes for short reads sequenced genomes
    - Flye for PacBio long reads genomes

For Illumina short reads, the steps include the following:
- QC using fastqc
- Read cleaning using fastp
- Reads count using grep
    cat /scratch/syersin2/Satellite_scratch/Isolates/read_count/*.txt > /scratch/syersin2/Satellite_scratch/Isolates/Merged_read_counts.txt
- Assembly using SPAdes
- Contigs and scaffolds filtering using scaffold_filter.py

For the PacBio, the sequenced were assembled by the sequencing facility. The data were then directly copied using: 
cp -r /nas/FAC/FBM/DMF/pvonaesc/default/D2c/Simon\ Yersin/Satellite/PacBio_Strept_Isolates/Fly_assemblies /scratch/syersin2/Satellite_scratch/Isolates/
Then renamed with the sample name and moved to a folder using: PacBio/00_copy_files

I combined all the assemblies files under one folder wgs_assemblies to run checkM and GTDB on all instead of separatly:
- Assembly QC using CheckM
- Taxonomic assignation using GTDB-TK

Saved files and assemblies in the NAS

Coverage calculation
Removed AGN043DU
Alignment of reads from isolates on assembled genomes using bowtie2 and samtools:
1. Index isolated genome: iso_indexing.sh
2. Align and depth calculation: iso_coverage.sh
Average coverage for each isolates with:
    find . -type f -name "*.depth" | while read file; do
        avg=$(awk '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }' "$file")
        echo "${file},${avg}" >> average_coverage.csv
    done
Overall general coverage:
    awk -F, 'NR > 1 { sum += $2 } END { print "Overall average coverage = ", sum / (NR - 1) }' average_coverage.csv
    Overall average coverage =  838
Saved .depth, .stats.txt and .log files on the nas and average_coverage.csv on the nas

Pangenome analysis
Ran prokka on S. salivarius genomes - 88 genomes from Satellite (removed contaminated strains, parasanguinis and afribiota strains)
Took .gff file produced by prokka (v.1.14.6)
Run roary (v 3.13.0) to analyse pangenome of S salivarius strains
    Use seqkit to translate pan_genome reference into amino acid:
        seqkit translate --frame 1 --trim pan_genome_reference.fa > pan_genome_reference_proteins.faa
Use scoary (v.1.6.16) to evaluate difference betwenn groups

Run egg-nog on file pan_genome_reference_proteins.faa
    Annotation saved in eggnogg folder
    emapper v.2.1.12, Eggnog db v.5.0.2, diammond v2.1.11, MMseqs2 v 17.b804f, novel families db v.1.0.1

#############---- Downstream analysis ------###########################################################################################
    Comparative genomic:
FastANI: I computed the fatsANI distance between the isolated genome and the Streptococcus MAGs (with completeness >90% and conta < 5%)
    Use fastANI script: fastANI.sh
    Input: genomes list (fastANI_list.txt)
    Output: SAT_fastANI.tsv file with the list of pairwise comparison and the ANI distance
I computed the fastANI distance between genomes of the Veillonella, Rothia and Haemophilus genus (with completeness >90% and conta < 5%)
    Use fastANI script: fastANI.sh
    Input: genomes list (fastANI_list_rvh.txt)
    Output: RVH_fastANI.tsv file with the list of pairwise comparison and the ANI distance
    
    01. Gene annotations:
    Annotation of the genes was done with prokka: 01_prokka.sh

    02. Phylogeny:
A phylogenetic tree with the genomes of interest can be produced using orthofinder and IQTree
For isolated genomes:
2025-05-06
Restart the phylogeny analysis using S. salivarius and S. parasanguinis genomes only (96 genomes)
Renamed file with rename function keep only SAT### and renaming SAT892 to SAT692.
Removing SAT350 and SAT622 as assembly are contaminated
    0. Gunzip the genomes
    1. Annotate genome using Prokka: 01_prokka.sh
        Isolates/prokka
        Copy .faa file in Isolates/prokka/faa_files using :
        for files in $(ls); do cp /scratch/syersin2/Satellite_scratch/Isolates/prokka/${files}/${files}.faa /scratch/syersin2/Satellite_scratch/Isolates/prokka/faa_files; done
        gzip genomes again
    2. Run orthofinder using: 02_orthofinder.sh
    3. Run IQ tree for phylogenetic tree reconstruction: 03_iq_tree.sh
        Timed out after 1 hour --> relaunched (4 hours time limit) and it starts at the previous time point
        Timed out after 4 hour --> relaunched (4 hours time limit) and it starts at the previous time point 
Visualize tree and evaluate if we should include MAGs?
2025-08-11
Included 13 MAGs in the phylogeny: MAGs had 75% completeness and <10% contamination 
    0. Run prodigal to obtain genome in .faa format
    1. Re-runing Orthofinder
    2. IQ-tree2
    --> Results in Aug11

    03. Transmission 
TRANSMISSION ANALYSIS - STRAIN-LEVEL DIVERSITY ANALYSIS
Once we obtained all the MAGs from our assembly, we used dRep to dereplicate the MAGs.
First all the MAGs from all samples were grouped together using: 00_group_mags.sh
Then, we moved all the assembly files to the folder with the MAGs to derep everything together

In R (MAGs_analysis.Rmd), loaded the checkM tables and merged them together. Saved a csv file for dRep with the genome name, completeness and contamination: checkm_drepV2.csv
In excel, changed .fasta to .fa names. Then for the isolates, changed the files from .fasta to .fa to have all the same extensions for MAGs and isolates, using: 
    for file in *.fasta; do mv "$file" "${file%.fasta}.fa"; done

1. Dereplication
For the dereplication, I provided the checkM results but simplified csv with 3 columns: genome, completeness, contamination (see line above).
I selected the first mash threshold at 95% (species) and the secondary ANI threshold at 98% (strains), using fastANI as algorithm.
I also added a coverage of 0.3 and completeness at 75 and contamination at 10%.
Each genomes in the dereplicated genome folder after dRep is representing a strains.
    01_1_drep_strains.sh
            At the strain level, we have 2180 dereplicated genomes
            Copying log, figure, data table and dereplicated genomes on the nas (zipping genomes)
    01_2_drep_species is the same with 90% and 95% threshold to have the representative genome at the species level.
            At the species level, we have 900 dereplicated genomes
            Copying log, figure, data table and dereplicated genomes on the nas (zipping genomes)

2. Preparing file for inStrain
Issue with contig name in dereplicated genome files, for isolates sequenced using PacBio, the assembly give the name of a contig within the files as contig_1, contig_2 etc 
In the MAGs the names are >MAGs_name_contig_1 and short reads assembly of isolates the name are Isolate_name_contig. 
When merging the dereplicated genomes into a single file allGenomes.fasta, several contigs have the same names contig1, contig2 etc but they are from different isolates 
It creates issues with the indexing and aligning using bowtie2
Using 02_00_rename_contigs.sh I renamed the contig name to include the isolate name in front. File rename_contigs_48711606.out show which file have been changed

Once the MAGs and isolates are dereplicated, they need to be combine into a single fasta file for inStrain:
    cat /scratch/syersin2/Satellite_scratch/transmission_dir/dereplication/dereplicated_genomes/*.fa > /scratch/syersin2/Satellite_scratch/transmission_dir/tables/genomes_db/allGenomes.fasta
Then, inStrain need a file telling from which genome each scaffold is coming from. For this we use dRep python script: parse_stb.py (use --help for explanations)
    Use --reverse option to generate the stb file from a list of genome.
    conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/dRep/drep
    parse_stb.py --reverse -f /scratch/syersin2/Satellite_scratch/transmission_dir/dereplication/dereplicated_genomes/*.fa -o /scratch/syersin2/Satellite_scratch/transmission_dir/tables/genomes.stb
Map reads to genomes (MAGs) using bowtie2
    Build index using bowtie2: note, use enough memory (8 was not enough, used 20 and it worked)
Map reads to the reference genomes index using bowtie2 and use sametools to go from a sam file to a bam file (mapping_output)
Run prodigal for gene-level profiling (produce .fna file)

3. Run inStrain profile (in database mode)
Profile needs: the bam file produced above, the .fna file from prodigal, the sample reads, the .stb file produced with dRep
Option: --data-base mode and --skip_plot_generation
    - Save profiles on Recherche nas

4. Run inStrain compare 
Compare needs: a list of profile directories to compare, the .stb file produced with dRep
Parameters: -ani 0.99999 for the genome identity and -cov 0.5 for coverage threshold
Output: instrain_compare/output/SAMPLE/instrain_compare_genomeWide_compare.tsv
    1. Do comparison for each individuals 
    2. Merging genome wide comparison table (in inStrain compare directory): 
        awk 'FNR==1 && NR!=1 {next} {print}' ./AGN*/output/*_genomeWide_compare.tsv > instrain_compare_genomeWide_compare_merged.tsv
    Copied both instrain profiles and compares on the nas (only compare for MAGs comparison)
