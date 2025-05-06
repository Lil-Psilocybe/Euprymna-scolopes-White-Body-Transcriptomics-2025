#!/bin/bash
#SBATCH --job-name=Initial_Trim_QC.sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mem=100G
#SBATCH --mail-user=john.briseno@uconn.edu
#SBATCH -o /home/FCAM/jbriseno/JB_Nyholm_Lab/Scripts/Error.Output/Initial_Trim_QC_%j.out
#SBATCH -e /home/FCAM/jbriseno/JB_Nyholm_Lab/Scripts/Error.Output/Initial_Trim_QC_%j.err


###Loading modules for analysis pipeline###
module load Trimmomatic/0.39 #trim
module load fastqc/0.11.7 #qc
module load MultiQC/1.12 #qc
module load star/2.7.11a #aligmment
module load subread/2.0.3 #quantification 


######################
##Variables & Paths##
#####################

###########
#TRIM & QC#
###########
#Raw reads trim and QC of demultiplexed data recieved in 2 batches from the Simakov lab
#where first batch (sequenced DATE) is 
FIRST_BATCH_RAW_READS=/labs/Nyholm/Oleg_project/HHH72DSX2_1_R11994_20210815/demultiplexed/
#and second batch (sequenced DATE) is 
SECOND_BATCH_RAW_READS=/labs/Nyholm/Oleg_project/HCVFYDSX3_1_R13227_20220410/demultiplexed/
#Raw data #Illumina read variables
FP=L001_R1_001.fastq.gz
RP=L001_R2_001.fastq.gz
#Other QC paths
PRE_TRIM=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/QC/PRE_TRIM/
POST_TRIM=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/QC/POST_TRIM/
MULTIQC=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/QC/MULTIQC/
LOGS=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/QC/LOGS/
#Output cleaned reads
CLEANED_READS=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/READS/CLEANED_READS/
Unpaired_Reads=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/READS/Unpaired_Reads/

###############
###ALIGNMENT###
###############
#E. scolopes V2 Gene Annotation for STAR Genome Aligner
Es_Genome_V2=/labs/Nyholm/Es_V2_Genome/Lachesis_assembly.fasta
Es_V2_Annotation_GTF=/labs/Nyholm/Es_V2_Genome/eupsc_models_v2.2.tags.gtf
STAR_Genome_BRAKER_BULK=/labs/Nyholm/Es_V2_Genome/STAR_Genome_BRAKER_BULK/
TAG=Parent
FEATURE=exon
CPUs=20
BAM=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/ALIGNMENT/BAM/

####################
###Quantification###
####################
#Quantifying reads to genes (exons) with featureCounts 
FEATURECOUNTS=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/WB/COUNTS/FEATURECOUNTS/


#Organ Sample Prefixes
MNTLs="
167717_S2
190853_S1
190910_S9
190922_S21
190931_S28
"
SKINs="
167720_S5
190904_S3
190912_S11
190924_S23
190933_S30
"
HECT0_ARMs="
167719_S4
190854_S2
190911_S10
190923_S22
190932_S29
"
CBs="
190909_S8
190914_S13
190921_S20
167744_S29
190939_S36
190954_S51
"
OLLs="
167722_S7
167728_S13
190906_S5
167736_S21
167745_S30
190935_S32
"
WBLs="
167730_S15
190913_S12
190916_S15
190918_S17
190907_S6
190925_S24
167737_S22
167747_S32
190951_S48
190936_S33
"
OVRYs="
167749_S34
167741_S26
190934_S31
190965_S53
"
TSTs="
167721_S6
167734_S19
190905_S4
190915_S14
"
LO_CCs="
190908_S7
190917_S16
190920_S19
167743_S28
167751_S36
190938_S35
190953_S50
"
#Combining mutiple variables into one
#https://stackoverflow.com/questions/8446146/combining-multiple-variables-into-another-variable-in-unix
ORGANS="${HECT0_ARMs} ${MNTLs} ${SKINs} ${WBLs} ${CBs} ${OLLs} ${OVRYs} ${TSTs} ${LO_CCs}"

####################
###END USER INPUT###
####################




###############
###TRIM & QC###
###############

#Doing fastqc reports of raw reads, 
#trimming reads with Trimmomatic, 
#then doing a multiqc report on cleaned reads

for organ in ${ORGANS} #This can be subbed for specific organs -> to use just use the organ! like "for wb in ${WB}"
	do
		#Pre-trim FastQC
		fastqc -o ${PRE_TRIM} ${FirstBatch_DEMULTIPEXED}"$organ"_${FP}
		fastqc -o ${PRE_TRIM} ${FirstBatch_DEMULTIPEXED}"$organ"_${RP}
		fastqc -o ${PRE_TRIM} ${SecondBatch_DEMULTIPEXED}"$organ"_${FP}
		fastqc -o ${PRE_TRIM} ${SecondBatch_DEMULTIPEXED}"$organ"_${RP}
		###This Wurks###
		
		#Need to clean of adapters
		echo "IlluminaClippin' raw $organ seqs of adapters from Batch1"
		java -jar $Trimmomatic PE \
		${FirstBatch_DEMULTIPEXED}"$organ"_${FP} \
		${FirstBatch_DEMULTIPEXED}"$organ"_${RP} \
		${CLEANED_READS}"$organ"_FP.fq.gz \
		${Unpaired_Reads}"$organ"_FuP.fq.gz \
		${CLEANED_READS}"$organ"_RP.fq.gz \
		${Unpaired_Reads}"$organ"_RuP.fq.gz \
		ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 
		MINLEN:150
		#4 files -> uP is unpaired good quality reads that got their pairs removed for low quality
		###This Works###
		
		#Need to clean of adapters
		echo "IlluminaClippin' raw $organ seqs of adapters from Batch2"
		java -jar $Trimmomatic PE \
		${SecondBatch_DEMULTIPEXED}"$organ"_${FP} \
		${SecondBatch_DEMULTIPEXED}"$organ"_${RP} \
		${CLEANED_READS}"$organ"_FP.fq.gz \
		${Unpaired_Reads}"$organ"_FuP.fq.gz \
		${CLEANED_READS}"$organ"_RP.fq.gz \
		${Unpaired_Reads}"$organ"_RuP.fq.gz \
		ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 
		MINLEN:150
		4 files -> uP is unpaired good quality reads that got their pairs removed for low quality
		###This Works###
	done
fastqc -o ${POST_TRIM}/5.2024.QC/ ${CLEANED_READS}*.fq.gz
multiqc -n multiqc_report_rnaseq -n ${MULTIQC}POST_TRIM_MULTIQC \ #Output in ${MULTIQC} as .html
${POST_TRIM}/5.2024.QC/*.zip \
${LOGS}*Log.final.out \
${COUNTS}WB_featureCounts.txt.summary



###############
###ALIGNMENT###
###############

#Indexing genome for Bulk STAR alignment
STAR \
--runThreadN ${CPUs} \
--runMode genomeGenerate \
--genomeDir ${STAR_Genome_BRAKER_BULK} \
--genomeFastaFiles ${Es_Genome_V2} \
--sjdbGTFfile ${Es_V2_Annotation_GTF} \
--sjdbGTFfeatureExon ${FEATURE} \
--sjdbOverhang 149 #150 bp read length so N - 1 = 149 for Overhang

#Aligning cleaned reads from TRIM & QC section to this newly indexed genome
cd ${BAM}
for organ in ${ORGANS} #This can be subbed for specific organs -> to use just use the organ! like "for wb in ${WB}"
	do	
	
		#Mappin to V2 of the genome using STAR 2pass alignment 
		echo "Basic two-pass mapping mode to ID splice junctions"
		STAR \
		--runThreadN ${CPUs} \
		--runMode alignReads \
		--runDirPerm All_RWX \
		--sjdbOverhang 149 \
		--twopassMode Basic \
		--sjdbGTFfeatureExon ${FEATURE} \
		\
		--genomeDir ${STAR_Genome_BRAKER_BULK} \
		--sjdbGTFfile ${Es_V2_Annotation_GTF} \
		\
		--readFilesIn ${CLEANED_READS}"$organ"_FP.fq.gz ${CLEANED_READS}"$organ"_RP.fq.gz \
		--readFilesCommand gunzip -c \
		\
		--outFileNamePrefix ${BAM}"$organ" \
		--outTmpKeep None \
		--outSAMtype BAM SortedByCoordinate \
		--outFilterMultimapNmax 1 \
	
	done



####################
###Quantification###
####################

##Quantification with FeatureCounts##
featureCounts -T 14 \
-a ${Es_V2_Annotation_GTF} \
-t exon -g gene_id \
-p --countReadPairs \
-o ${FEATURECOUNTS}WBL.CB.OLL.LO_CC.OVRY.TST.MNTL.SKIN.HECTO_ARM.BRAKER.gIDs.featureCounts.matrix \
--donotsort \































