#!/bin/bash
## -------------------------------------
## RUNNING ON Pre-GATK, GATK, and more
## 2015/05/10
## BHOOM SUKTITIPAT, MD, PhD
## based on GATK_CISCO_IPGG_20150507.sh
## -- create script to run all workflow at the end.
## -------------------------------------

## --------------------- INSTRUCTION -----------------------## 
# Please religiously follow the followings 
#-------------------------------------------------------------
# PUT THIS SCRIPT IN $PROJECT_FOLDER
#-------------------------------------------------------------

# 1) Install bwa, picard, and gatk in $HOME/bin
# 2) Download gatk_bundle files from GATK ftp sites in $GATK_RESOURCE_PATH
#   - this includes reference genome (hg19), known indels, dbsnp138
# 3) set up bwa index
# 4) put fastq files in $PROJECT_FOLDER/$FASTQ_DIR 
# 5) This script will create a subscript for each process in $SAMPLE/$SCRIPT
# 6) Change the number of threads to reflect your CPUs
# 7) Change the INPUT DATA below to match with your input
## =============================================================
## INPUT DATA
FASTQ_DIR=$1 # folder of all fastq files
FASTQ_PREFIX=$2 # Prefix in the format ${FASTQ_DIR}/${FASTQ_PREFIX}_R{1,2}.fastq.gz 
SAMPLE=$3
PROJECT=$4 ## use in read group creation
LB=$5 ## ?? use in read group creation e.g. WES/WGS
PLATFORM=$6 ## use in read group creation i.e. Illumina
NUMBER_THREADS=40
NT=10 ## paramter in GATK
NCT=10 ## parameter in GATK

## PROGRAM INFO PATH
BWA_PATH=$HOME/bin
PICARD_PATH=$HOME/bin/picard-tools-1.93 ## FOR NGS.CISCO SERVER
#* PICARD_PATH=~/bin/picard-tools-1.119 ## path on archimedes & MBP
GATK_DIR=$HOME/bin/gatk3.3-0
GATK_RESOURCE_PATH=/mnt/data/data/hg19/gatk_bundle


REF_GENOME=$GATK_RESOURCE_PATH/ucsc.hg19.fasta
RG_NAME="@RG\tID:${PROJECT}_${SAMPLE}\tSM:${SAMPLE}\tPL:${PLATFORM}\tLB:${LB}\tPU:unit1"
KNOWN_INDEL_1=$GATK_RESOURCE_PATH/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
KNOWN_INDEL_2=$GATK_RESOURCE_PATH/1000G_phase1.indels.hg19.sites.vcf
DBSNP=$GATK_RESOURCE_PATH/dbsnp_138.hg19.vcf
##BAIT_BED=

mkdir -p $SAMPLE/{BAM,SNV/QC/FILTERED_ON_BAIT,TEMP,LOG,GATK,GVCF,VCF/QC/FILTERED_ON_BAIT,SCRIPT}
## -------------------------------------
## SEQUENCE ASSEMBLY
## -------------------------------------

# 1) Align read with reference
cat <<EOL > $SAMPLE/SCRIPT/01_${PROJECT}_${SAMPLE}_bwa.sh
#!/bin/bash
$BWA_PATH/bwa mem -t $NUMBER_THREADS -M \
-R "$RG_NAME" \
$REF_GENOME \
${FASTQ_DIR}/${FASTQ_PREFIX}_R1.fastq.gz \
${FASTQ_DIR}/${FASTQ_PREFIX}_R2.fastq.gz > ${SAMPLE}/BAM/${FASTQ_PREFIX}_aligned.sam
    
EOL
## -------------------------------------

# 2) Sort Bam file
cat <<EOL > $SAMPLE/SCRIPT/02_${PROJECT}_${SAMPLE}_SortSam.sh
#!/bin/bash
java -Xmx${JAVA_MEM} -jar ${PICARD_PATH}/SortSam.jar \
INPUT=${SAMPLE}/BAM/${FASTQ_PREFIX}_aligned.sam \
OUTPUT=${SAMPLE}/BAM/${FASTQ_PREFIX}_sorted.bam \
SORT_ORDER=coordinate \
TMP_DIR=${SAMPLE}/TEMP
EOL
## -------------------------------------
# 3) Mark duplicates
cat <<EOL > $SAMPLE/SCRIPT/03_${PROJECT}_${SAMPLE}_MarkDuplicates.sh
#!/bin/bash
java -Xmx${JAVA_MEM} -jar ${PICARD_PATH}/MarkDuplicates.jar \
INPUT=${SAMPLE}/BAM/${FASTQ_PREFIX}_sorted.bam \
OUTPUT=${SAMPLE}/BAM/${FASTQ_PREFIX}_dedup.bam \
METRICS_FILE=${SAMPLE}/BAM/${FASTQ_PREFIX}_dedup.matrics \
CREATE_INDEX=TRUE \
TMP_DIR=${SAMPLE}/TEMP 
EOL
## -------------------------------------
# 4) Build bam index
cat <<EOL > $SAMPLE/SCRIPT/04_${PROJECT}_${SAMPLE}_BuildBamIndex.sh
#!/bin/bash
java -Xmx${JAVA_MEM} -jar ${PICARD_PATH}/BuildBamIndex.jar \
INPUT=${SAMPLE}/BAM/${FASTQ_PREFIX}_dedup.bam \
TMP_DIR=${SAMPLE}/TEMP 
if [[ -e ${SAMPLE}/BAM/${FASTQ_PREFIX}_dedup.bam ]]; then
  rm ${SAMPLE}/BAM/${FASTQ_PREFIX}_aligned.sam
fi
EOL

## -------------------------------------
# 5) realigner target creator
cat <<EOL > $SAMPLE/SCRIPT/05_${PROJECT}_${SAMPLE}_RealignerTargetCreator.sh
#!/bin/bash
java -Xmx48g -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
--disable_auto_index_creation_and_locking_when_reading_rods \
-known $KNOWN_INDEL_1 \
-known $KNOWN_INDEL_2 \
-R ${REF_GENOME} \
-I ${SAMPLE}/BAM/${FASTQ_PREFIX}_dedup.bam \
-dt NONE \
-nt 24 \
-o ${SAMPLE}/BAM/${FASTQ_PREFIX}_indel_target_intervals.list \
-log ${SAMPLE}/LOG/05_${PROJECT}_${SAMPLE}_indel_target_intervals.log  
EOL

## 6) Indel realigner
cat <<EOL > $SAMPLE/SCRIPT/06_${PROJECT}_${SAMPLE}_IndelRealigner.sh
#!/bin/bash
java -Xmx4g -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
-T IndelRealigner \
--disable_auto_index_creation_and_locking_when_reading_rods \
-I ${SAMPLE}/BAM/${FASTQ_PREFIX}_dedup.bam \
-R $REF_GENOME \
-known $KNOWN_INDEL_1 \
-known $KNOWN_INDEL_2 \
-targetIntervals ${SAMPLE}/BAM/${FASTQ_PREFIX}_indel_target_intervals.list \
-dt NONE \
-o ${SAMPLE}/BAM/${FASTQ_PREFIX}_realigned.bam \
-log ${SAMPLE}/LOG/06_${PROJECT}_${SAMPLE}_realigned_indel.log

EOL

## 7) BQSR
cat <<EOL > $SAMPLE/SCRIPT/07_${PROJECT}_${SAMPLE}_BaseRecalibrator.sh
#!/bin/bash
java -Xmx4g -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
-I ${SAMPLE}/BAM/${FASTQ_PREFIX}_realigned.bam \
-knownSites $KNOWN_INDEL_1 \
-knownSites $KNOWN_INDEL_2 \
-knownSites $DBSNP \
-nct 8 \
-o ${SAMPLE}/BAM/${FASTQ_PREFIX}_BQSR.bqsr \
-log ${SAMPLE}/LOG/07_${PROJECT}_${SAMPLE}_BQSR.log 
EOL

## 8) 

# ## 8) Final BAM
cat <<EOL > $SAMPLE/SCRIPT/08_${PROJECT}_${SAMPLE}_PrintReads.sh
#!/bin/bash
java -Xmx4g -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
-T PrintReads \
-R $REF_GENOME \
-I ${SAMPLE}/BAM/${FASTQ_PREFIX}_realigned.bam \
-BQSR ${SAMPLE}/BAM/${FASTQ_PREFIX}_BQSR.bqsr \
-dt NONE \
-EOQ \
-nct 8 \
-o ${SAMPLE}/GATK/${FASTQ_PREFIX}_GATK.bam \
-log ${SAMPLE}/LOG/08_${PROJECT}_${SAMPLE}_final_bam.log

if [[ -e ${SAMPLE}/BAM/${FASTQ_PREFIX}_GATK.bam ]]; then
  rm ${SAMPLE}/BAM/${FASTQ_PREFIX}_{dedup,sorted}.bam
fi 
EOL

## 9) Call variants
cat <<EOL > $SAMPLE/SCRIPT/09_${PROJECT}_${SAMPLE}_HaplotypeCaller.sh 
#!/bin/bash
java -Xmx32g -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $REF_GENOME \
--input_file ${SAMPLE}/GATK/${FASTQ_PREFIX}_GATK.bam \
--emitRefConfidence GVCF \
--variant_index_type LINEAR \
--variant_index_parameter 128000 \
-A DepthPerSampleHC \
-A ClippingRankSumTest \
-A MappingQualityRankSumTest \
-A ReadPosRankSumTest \
-A FisherStrand \
-A GCContent \
-A AlleleBalanceBySample \
-A AlleleBalance \
-A QualByDepth \
-pairHMM VECTOR_LOGLESS_CACHING \
-nct 3 \
-o ${SAMPLE}/GVCF/${FASTQ_PREFIX}_GATK.vcf \
-log ${SAMPLE}/LOG/09_${PROJECT}_${SAMPLE}_HaplotypeCaller.log
EOL

## -- CIDR used HaplotypeCaller on -L $BAIT_BED

## 10) Now creating a normal VCF file
cat <<EOL > $SAMPLE/SCRIPT/10_${PROJECT}_${SAMPLE}_GenotypeGVCF.sh
#!/bin/bash
java -Xmx${JAVA_MEM} -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $REF_GENOME \
--variant ${SAMPLE}/GVCF/${FASTQ_PREFIX}_GATK.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-nt 1 \
-o ${SAMPLE}/VCF/${FASTQ_PREFIX}_QC_RAW_OnBait.vcf \
-log ${SAMPLE}/LOG/10_${PROJECT}_${SAMPLE}_GenotypeGVCF.log
EOL

#11 Extract SNPs from the above file
cat <<EOL > $SAMPLE/SCRIPT/11_${PROJECT}_${SAMPLE}_SelectSNV.sh
#!/bin/bash
java -Xmx${JAVA_MEM} -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
--variant ${SAMPLE}/VCF/${FASTQ_PREFIX}_QC_RAW_OnBait.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-selectType SNP \
--excludeFiltered \
-nt 1 \
-o ${SAMPLE}/SNV/QC/${FASTQ_PREFIX}_QC_RAW_OnBait_SNV.vcf \
-log ${SAMPLE}/LOG/11_${PROJECT}_${SAMPLE}_SelectSNV.log
EOL

#12) Annotate the above file with dbSNP variant type
cat <<EOL > $SAMPLE/SCRIPT/12_${PROJECT}_${SAMPLE}_SNV_Annotation.sh
#!/bin/bash
java -Xmx${JAVA_MEM} -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R $REF_GENOME \
--variant ${SAMPLE}/SNV/QC/${FASTQ_PREFIX}_QC_RAW_OnBait_SNV.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--dbsnp $DBSNP \
-L ${SAMPLE}/SNV/QC/${FASTQ_PREFIX}_QC_RAW_OnBait_SNV.vcf \
-A GCContent \
-A VariantType \
-dt NONE \
-nt 1 \
-o ${SAMPLE}/SNV/QC/FILTERED_ON_BAIT/${FASTQ_PREFIX}_QC_RAW_OnBait_SNV_ANNOTATED.vcf \
-log ${SAMPLE}/LOG/12_${PROJECT}_${SAMPLE}_SNV_Annotation.log
## -nt 1 finished in 20 minutes while -nt 4 takes 4 hours and crashed
EOL
## -- replace -L with bait_bed ??

#13 Filtering or flagging the variants in the above file.
cat <<EOL > $SAMPLE/SCRIPT/13_${PROJECT}_${SAMPLE}_FilterSNV.sh
#!/bin/bash
java -Xmx${JAVA_MEM} -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $REF_GENOME \
--variant ${SAMPLE}/SNV/QC/FILTERED_ON_BAIT/${FASTQ_PREFIX}_QC_RAW_OnBait_SNV_ANNOTATED.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--filterExpression 'QD < 2.0' \
--filterName 'QD' \
--filterExpression 'MQ < 30.0' \
--filterName 'MQ' \
--filterExpression 'FS > 40.0' \
--filterName 'FS' \
--filterExpression 'MQRankSum < -12.5' \
--filterName 'MQRankSum' \
--filterExpression 'ReadPosRankSum < -8.0' \
--filterName 'ReadPosRankSum' \
--filterExpression 'DP < 8.0' \
--filterName 'DP' \
--logging_level ERROR \
-o ${SAMPLE}/SNV/QC/FILTERED_ON_BAIT/${FASTQ_PREFIX}_QC_FILTERED_OnBait_SNV.vcf \
-log ${SAMPLE}/LOG/13_${PROJECT}_${SAMPLE}_QC_FilterSNV.log
EOL

#14 Now removing any SNVs that were flagged in the above file
cat <<EOL > $SAMPLE/SCRIPT/14_${PROJECT}_${SAMPLE}_cleanSNV.sh
#!/bin/bash
java -Xmx${JAVA_MEM} -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
--variant ${SAMPLE}/SNV/QC/FILTERED_ON_BAIT/${FASTQ_PREFIX}_QC_FILTERED_OnBait_SNV.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--excludeFiltered \
-nt 1 \
-o ${SAMPLE}/SNV/QC/FILTERED_ON_BAIT/${FASTQ_PREFIX}_Clean_OnBait.vcf \
-log ${SAMPLE}/LOG/14_${FASTQ_PREFIX}_QC_cleanSNV.log
EOL

#15 Doing the same as above, except now also filtering on the target bed file

##!/bin/bash
# java -jar $GATK_DIR/GenomeAnalysisTK.jar \
#-T SelectVariants \
#-R $REF_GENOME \
#--variant ${SAMPLE}/SNV/QC/FILTERED_ON_BAIT/${FASTQ_PREFIX}_QC_OnBait.vcf \
#-L $TARGET_BED \
#-nt ${NT} \
#-o ${SAMPLE}/SNV/QC/FILTERED_ON_TARGET/${FASTQ_PREFIX}_QC_OnTarget.vcf \
#-log ${SAMPLE}/LOG/${FASTQ_PREFIX}_15QC3_filterSNV.log

# Master Script to run all the subscripts
cat <<EOF > GATK_${PROJECT}_${SAMPLE}.sh
#bash $SAMPLE/SCRIPT/01_${PROJECT}_${SAMPLE}_bwa.sh 2>&1 | tee $SAMPLE/LOG/01_${PROJECT}_${SAMPLE}_bwa.log
#bash $SAMPLE/SCRIPT/02_${PROJECT}_${SAMPLE}_SortSam.sh 2>&1 | tee $SAMPLE/LOG/02_${PROJECT}_${SAMPLE}_SortSam.log
bash $SAMPLE/SCRIPT/03_${PROJECT}_${SAMPLE}_MarkDuplicates.sh 2>&1 |tee ${SAMPLE}/LOG/03_${PROJECT}_${SAMPLE}_MarkDuplicates.log
bash $SAMPLE/SCRIPT/04_${PROJECT}_${SAMPLE}_BuildBamIndex.sh 2>&1 | tee ${SAMPLE}/LOG/04_${PROJECT}_${SAMPLE}_BuildBamIndex.log
bash $SAMPLE/SCRIPT/05_${PROJECT}_${SAMPLE}_RealignerTargetCreator.sh ## log files of GATK in GATK command -- in LOG folder
bash $SAMPLE/SCRIPT/06_${PROJECT}_${SAMPLE}_IndelRealigner.sh ## log files of GATK in GATK command -- in LOG folder
bash $SAMPLE/SCRIPT/07_${PROJECT}_${SAMPLE}_BaseRecalibrator.sh ## log files of GATK in GATK command -- in LOG folder
bash $SAMPLE/SCRIPT/08_${PROJECT}_${SAMPLE}_PrintReads.sh ## log files of GATK in GATK command -- in LOG folder
bash $SAMPLE/SCRIPT/09_${PROJECT}_${SAMPLE}_HaplotypeCaller.sh  ## log files of GATK in GATK command -- in LOG folder
bash $SAMPLE/SCRIPT/10_${PROJECT}_${SAMPLE}_GenotypeGVCF.sh ## log files of GATK in GATK command -- in LOG folder
bash $SAMPLE/SCRIPT/11_${PROJECT}_${SAMPLE}_SelectSNV.sh ## log files of GATK in GATK command -- in LOG folder
bash $SAMPLE/SCRIPT/12_${PROJECT}_${SAMPLE}_SNV_Annotation.sh ## log files of GATK in GATK command -- in LOG folder
bash $SAMPLE/SCRIPT/13_${PROJECT}_${SAMPLE}_FilterSNV.sh ## log files of GATK in GATK command -- in LOG folder
bash $SAMPLE/SCRIPT/14_${PROJECT}_${SAMPLE}_cleanSNV.sh ## log files of GATK in GATK command -- in LOG folder
EOF
