#!/bin/bash
#SBATCH --job-name=STARsolo.sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 5
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=100G
#SBATCH --mail-user=john.briseno@uconn.edu
#SBATCH -o /home/FCAM/jbriseno/JB_Nyholm_Lab/Scripts/Error.Output/STARsolo.sh_%j.out
#SBATCH -e /home/FCAM/jbriseno/JB_Nyholm_Lab/Scripts/Error.Output/STARsolo.sh_%j.err

###Loading modules for analysis pipeline###
module load star/2.7.11a


########################################################
###Generating and Setting Variables/Directories/Paths###
########################################################

#CPUs
CPUs=5

#genome/data paths/info
Es_Genome_V2=/labs/Nyholm/Es_V2_Genome/Lachesis_assembly.fasta
Es_V2_Annotation_GTF=/labs/Nyholm/Es_V2_Genome/eupsc_models_v2.2.tags.gtf
STAR_Genome_BRAKER_SCSN=/labs/Nyholm/Es_V2_Genome/STAR_Genome_BRAKER_SCSN/
scsnRNASeq=/labs/Nyholm/scsnRNASeq/
PRE_TRIM=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/QC/PRE_TRIM/


#STAR flags
CBLEN=16
UMILEN=12
STRAND=Forward
whitelist=/labs/Nyholm/scsnRNASeq/3M-february-2018.txt

#setting sample vairables
Male_Samples="squid_male_S2"
squid_sample="${Male_Samples}"

R1=R1_001
R2=R2_001

hostname

####################
###END USER INPUT###
####################


cd ${scsnRNASeq}

###Indexing genome for sc/sn STAR alignment###
STAR \
--runThreadN ${CPUs} \
--runMode genomeGenerate \
--genomeDir ${STAR_Genome_BRAKER_SCSN} \
--genomeFastaFiles ${Es_Genome_V2} \
--sjdbGTFfeatureExon ${FEATURE} \
--sjdbGTFfile ${Es_V2_Annotation_GTF}


###STARsolo Parameters###
#for sequences with v3.1 chemistry from 10X
#using braker gene annotations with good 3' UTR annotation
#standard 10X runs have cDNA as Read2 and barcode as Read1, the 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read
#all-read/write/execute (same as chmod 777)
#may need to run UMIstart at 17
#This was adapted from <- https://github.com/cellgeni/STARsolo
	#Will read in GeneFull output (recommended for snRNAseq with pre-mRNAs) into R Seuart
for sample in ${squid_sample}
	do
		echo "STARsolo run for ${sample}"
		
		/core/cbc/nreid/STAR-2.7.11a/source/STAR \
		--runThreadN ${CPUs} \
		--genomeDir ${STAR_Genome_BRAKER_SCSN} \
		--readFilesIn ${sample}_${R2}.fastq.gz ${sample}_${R1}.fastq.gz \
		--readFilesCommand gunzip -c \
		--runDirPerm All_RWX \
		--outFileNamePrefix ${sample}_ \
		--soloType CB_UMI_Simple --soloCBwhitelist ${whitelist} --soloBarcodeReadLength 0 \
		--soloCBlen ${CBLEN} --soloUMIstart $((CBLEN+1)) --soloUMIlen ${UMILEN} --soloStrand ${STRAND} \
		--soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
		--soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
		--soloFeatures Gene GeneFull Velocyto --soloOutFileNames ${sample}_SOLO.output/ features.tsv barcodes.tsv matrix.mtx \
		--outSAMtype BAM SortedByCoordinate \
		--soloMultiMappers EM --outReadsUnmapped Fastx \
		--soloCellReadStats Standard
	done