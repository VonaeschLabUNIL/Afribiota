#!usr/bin/env bash

cd /media/wslnx/18TERA/SIMON/Metagenomic/Satellite/Satellite_full_dataset24/BioInfo

mkdir TRIMMED_READS
mkdir REPORTS_TRIMMED_READS

# Install fastp
# conda install -c bioconda fastp

# Paired-end:
for prefix in $(ls *.fastq.gz | sed -E 's/_R[12]_001[.]fastq.gz//' | uniq)
	do
	echo "performing cleaning on files :" "${prefix}_R1_001.fastq.gz" "${prefix}_R2_001.fastq.gz"
	fastp -i "${prefix}_R1_001.fastq.gz" -I "${prefix}_R2_001.fastq.gz" \
	-o TRIMMED_READS/"${prefix}_FASTP_R1_001.fastq.gz"  -O TRIMMED_READS/"${prefix}_FASTP_R2_001.fastq.gz" \
	--report_title "${prefix}_fastp_report" \
	--thread 34 \
	-j REPORTS_TRIMMED_READS/"${prefix}_fastp".json -h REPORTS_TRIMMED_READS/"${prefix}_fastp".html

	# if single-ends :
	#echo "performing cleaning on files :" "${prefix}_R1_001.fastq.gz"
	#fastp -i "${prefix}_R1_001.fastq.gz" -o \
	#TRIMMED_READS/"${prefix}_FASTP_R1_001.fastq.gz" \
	#--report_title "${prefix}_fastp_report" \
	#--thread 38 \
	#-j REPORTS_TRIMMED_READS/"${prefix}_fastp".json \
	#-h REPORTS_TRIMMED_READS/"${prefix}_fastp".html
done

multiqc ./REPORTS_TRIMMED_READS/ --ignore-symlinks --outdir ./REPORTS_TRIMMED_READS \
--filename MULTIQC_ALL_SAMPLE_REPORT --fullnames --title MULTIQC_ALL_SAMPLE_REPORT

#cd ./TRIMMED_READS

