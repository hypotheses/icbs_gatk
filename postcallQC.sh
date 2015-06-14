#!/bin/bash
## -------------------------------------
## RUNNING ON Pre-GATK, GATK, and more
## 2015/04/29
## BHOOM SUKTITIPAT, MD, PhD
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
NUMBER_THREADS=12
NT=12 ## paramter in GATK
NCT=12 ## parameter in GATK
JAVA_MEM=30g ## for java -Xmx${JAVA_MEM} -jar 

## PROGRAM INFO PATH
BWA_PATH=$HOME/bin
PICARD_PATH=$HOME/bin/picard-tools-1.119 ## FOR NGS.CISCO & SIBC SERVER
GATK_DIR=$HOME/bin
GATK_RESOURCE_PATH=/Volumes/Raid1/ResearchData/hg19/gatk_bundle
#* PICARD_PATH=~/bin/picard-tools-1.119 ## path on archimedes & MBP
SAMTOOLS_DIR=$HOME/bin

# BED FILES
TI_TV_BED
BAIT_BED
TARGET_BED
EXON_BED
VERIFY_VCF

REF_GENOME=$GATK_RESOURCE_PATH/ucsc.hg19.fasta
RG_NAME="@RG\tID:${PROJECT}_${SAMPLE}\tSM:${SAMPLE}\tPL:${PLATFORM}\tLB:${LB}\tPU:unit1"
KNOWN_INDEL_1=$GATK_RESOURCE_PATH/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
KNOWN_INDEL_2=$GATK_RESOURCE_PATH/1000G_phase1.indels.hg19.sites.vcf
DBSNP=$GATK_RESOURCE_PATH/dbsnp_138.hg19.vcf


#------------------------------------# 
# POST HAPLOTYPE CALLING QC 
#____________________________________#

# from I.01_TITV_ALL.sh
# 13) Calculate overal ti/tv ratio

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
--excludeFiltered \
-L $TI_TV_BED \
--variant $SAMPLE/SNV/QC/FILTERED_ON_BAIT/$FASEQ_PREFIX"_QC_OnBait.vcf" \
-o $SAMPLE/TEMP/$FASEQ_PREFIX"_QC_FILTERED_TiTv_SNV.vcf"

# 14) Ti/Tv of known variants
# from I.02_TITV_KNOWN.sh 
## Below is writing out a vcf file of "known" variants (as defined by dbsnp_138.b37.excluding_sites_after_129.vcf).

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
--variant $SAMPLE/SNV/QC/FILTERED_ON_BAIT/$FASEQ_PREFIX"_QC_OnBait.vcf" \
--excludeFiltered \
-L $TI_TV_BED \
--concordance /isilon/sequencing/GATK_resource_bundle/2.8/b37/dbsnp_138.b37.excluding_sites_after_129.vcf \
-o $SAMPLE/TEMP/$FASEQ_PREFIX"_QC_FILTERED_TiTv_SNV_KNOWN.vcf"

# Now TiTv on those

$SAMTOOLS_DIR/bcftools/vcfutils.pl qstats $SAMPLE/TEMP/$FASEQ_PREFIX"_QC_FILTERED_TiTv_SNV_KNOWN.vcf" \
>| $SAMPLE/REPORTS/TI_TV/$FASEQ_PREFIX"_Known_.titv.txt"

# 15) Ti/Tv of novel variants
# from I.03_TITV_NOVEL.sh
## Below is writing out a vcf file of "novel" variants (as defined by dbsnp_138.b37.excluding_sites_after_129.vcf).

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
--variant $SAMPLE/SNV/QC/FILTERED_ON_BAIT/$FASEQ_PREFIX"_QC_OnBait.vcf" \
--excludeFiltered \
-L $TI_TV_BED \
--discordance /isilon/sequencing/GATK_resource_bundle/2.8/b37/dbsnp_138.b37.excluding_sites_after_129.vcf \
-o $SAMPLE/TEMP/$FASEQ_PREFIX"_QC_FILTERED_TiTv_SNV_NOVEL.vcf"

# Now TiTv on those

$SAMTOOLS_DIR/bcftools/vcfutils.pl qstats $SAMPLE/TEMP/$FASEQ_PREFIX"_QC_FILTERED_TiTv_SNV_NOVEL.vcf" \
 >| $SAMPLE/REPORTS/TI_TV/$FASEQ_PREFIX"_Novel_.titv.txt"

# 16) BQSR Summary Table and Plots
# from H.02_BQSR_PLOT.sh
## --Generate post BQSR table--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-I $SAMPLE/BAM/AGGREGATE/$FASEQ_PREFIX".bam" \
-R $REF_GENOME \
-knownSites $KNOWN_INDEL_1 \
-knownSites $KNOWN_INDEL_2 \
-knownSites $DBSNP \
-BQSR $SAMPLE/REPORTS/COUNT_COVARIATES/GATK_REPORT/$FASEQ_PREFIX"_PERFORM_BQSR.bqsr" \
-o $SAMPLE/REPORTS/COUNT_COVARIATES/GATK_REPORT/$FASEQ_PREFIX"_AFTER_BQSR.bqsr" \
-nct 8

## --Generate BQSR plots--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R $REF_GENOME \
-before $SAMPLE/REPORTS/COUNT_COVARIATES/GATK_REPORT/$FASEQ_PREFIX"_PERFORM_BQSR.bqsr" \
-after $SAMPLE/REPORTS/COUNT_COVARIATES/GATK_REPORT/$FASEQ_PREFIX"_AFTER_BQSR.bqsr" \
-plots $SAMPLE/REPORTS/COUNT_COVARIATES/PDF/$FASEQ_PREFIX".BQSR.pdf"

# 17) H.03_DOC_EXON_BED.sh
## --Depth of Coverage ALL UCSC EXONS--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $REF_GENOME \
-I $SAMPLE/BAM/AGGREGATE/$FASEQ_PREFIX".bam" \
-geneList:REFSEQ $GENE_LIST \
-L $EXON_BED \
-mmq 20 \
-mbq 10 \
--outputFormat csv \
-omitBaseOutput \
-o $SAMPLE/REPORTS/GENES_COVERAGE/UCSC_EXONS/$FASEQ_PREFIX".ALL.UCSC.EXONS" \
-ct 8 \
-ct 15 \
-ct 20

# 18) H.05_DOC_SUPERSET_BED.sh
### --Depth of Coverage On Bait--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $REF_GENOME \
-I $SAMPLE/BAM/AGGREGATE/$FASEQ_PREFIX".bam" \
-geneList:REFSEQ $GENE_LIST \
-L $BAIT_BED \
-mmq 20 \
-mbq 10 \
--outputFormat csv \
-omitBaseOutput \
-o $SAMPLE/REPORTS/GENES_COVERAGE/BED_SUPERSET/$FASEQ_PREFIX".SUPERSET.BED" \
-ct 8 \
-ct 15 \
-ct 20

# 19) H.06_ALIGNMENT_SUMMARY_METRIC.sh
## --Alignment Summary Metrics--

java -jar $PICARD_DIR/CollectAlignmentSummaryMetrics.jar \
INPUT=$SAMPLE/BAM/AGGREGATE/$FASEQ_PREFIX".bam" \
OUTPUT=$SAMPLE/REPORTS/ALIGNMENT_SUMMARY/$FASEQ_PREFIX"_alignment_summary_metrics.txt" \
R=$REF_GENOME \
VALIDATION_STRINGENCY=SILENT

# 20) H.07_BASECAL_Q_SCORE_DISTRIBUTION.sh
## --Base Call Quality Score Distribution--

java -jar $PICARD_DIR/QualityScoreDistribution.jar \
INPUT=$SAMPLE/BAM/AGGREGATE/$FASEQ_PREFIX".bam" \
OUTPUT=$SAMPLE/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/METRICS/$FASEQ_PREFIX"_quality_score_distribution.txt" \
CHART=$SAMPLE/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/PDF/$FASEQ_PREFIX"_quality_score_distribution_chart.pdf" \
R=$REF_GENOME \
VALIDATION_STRINGENCY=SILENT \
INCLUDE_NO_CALLS=true

# 21) H.08_GC_BIAS.sh
## --GC Bias Metrics--

java -jar $PICARD_DIR/CollectGcBiasMetrics.jar \
INPUT=$SAMPLE/BAM/AGGREGATE/$FASEQ_PREFIX".bam" \
OUTPUT=$SAMPLE/REPORTS/GC_BIAS/METRICS/$FASEQ_PREFIX"_gc_bias_metrics.txt" \
CHART_OUTPUT=$SAMPLE/REPORTS/GC_BIAS/PDF/$FASEQ_PREFIX"_gc_bias_metrics.pdf" \
SUMMARY_OUTPUT=$SAMPLE/REPORTS/GC_BIAS/SUMMARY/$FASEQ_PREFIX"_gc_bias_summary.txt" \
R=$REF_GENOME \
VALIDATION_STRINGENCY=SILENT

# 22 Insert size H.09_INSERT_SIZE.sh
## --Insert Size--

java -jar $PICARD_DIR/CollectInsertSizeMetrics.jar \
INPUT=$SAMPLE/BAM/AGGREGATE/$FASEQ_PREFIX".bam" \
OUTPUT=$SAMPLE/REPORTS/INSERT_SIZE/METRICS/$FASEQ_PREFIX"__insert_size_metrics.txt" \
H=$SAMPLE/REPORTS/INSERT_SIZE/PDF/$FASEQ_PREFIX"_insert_size_metrics_histogram.pdf" \
R=$REF_GENOME \
VALIDATION_STRINGENCY=SILENT

# 23 H.10_MEAN_QUALITY_BY_CYCLE
## --Mean Quality By Cycle--

java -jar $PICARD_DIR/MeanQualityByCycle.jar \
INPUT=$SAMPLE/BAM/AGGREGATE/$FASEQ_PREFIX".bam" \
OUTPUT=$SAMPLE/REPORTS/MEAN_QUALITY_BY_CYCLE/METRICS/$FASEQ_PREFIX"_mean_quality_by_cycle.txt" \
CHART=$SAMPLE/REPORTS/MEAN_QUALITY_BY_CYCLE/PDF/$FASEQ_PREFIX"_mean_quality_by_cycle_chart.pdf" \
R=$REF_GENOME \
VALIDATION_STRINGENCY=SILENT

# 24) H.11_CALCULATE_HS_METRICS.sh
# Calculate HS metrics bed files

($SAMTOOLS_DIR/samtools view -H $SAMPLE/BAM/AGGREGATE/$FASEQ_PREFIX".bam" \
| grep "@SQ" ; sed 's/\r//g' $BAIT_BED | awk '{print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' | sed 's/ /\t/g') \
>| $SAMPLE/TEMP/$FASEQ_PREFIX".OnBait.picard.bed"

($SAMTOOLS_DIR/samtools view -H $SAMPLE/BAM/AGGREGATE/$FASEQ_PREFIX".bam" \
| grep "@SQ" ; sed 's/\r//g' $TARGET_BED | awk '{print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' | sed 's/ /\t/g') \
>| $SAMPLE/TEMP/$FASEQ_PREFIX".OnTarget.picard.bed"

## --Calculate HS metrics and also calculate a per target coverage report--

java -jar /isilon/sequencing/BDC/Programs/PicardModifications/CalculateHSMetrics_CIDR_Mod_Jars_Deploy/1.57/CalculateHsMetrics_1.57_build1027_CIDR_130129.jar \
INPUT=$SAMPLE/BAM/AGGREGATE/$FASEQ_PREFIX".bam" \
OUTPUT=$SAMPLE/REPORTS/HYB_SELECTION/$FASEQ_PREFIX"__hybridization_selection_metrics.txt" \
PER_TARGET_COVERAGE=$SAMPLE/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE/$FASEQ_PREFIX"_per_target_coverage.txt" \
R=$REF_GENOME \
BI=$SAMPLE/TEMP/$FASEQ_PREFIX".OnBait.picard.bed" \
TI=$SAMPLE/TEMP/$FASEQ_PREFIX".OnTarget.picard.bed" \
VALIDATION_STRINGENCY=SILENT

#25) H.12_VerifyBamID.sh
# --Creating an on the fly VCF file to be used as the reference for verifyBamID--

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
--variant $VERIFY_VCF \
-L $TARGET_BED \
-o $SAMPLE/TEMP/$FASEQ_PREFIX".VerifyBamID.vcf"

## --Running verifyBamID--

$VERIFY_DIR/verifyBamID \
--bam $SAMPLE/BAM/AGGREGATE/$FASEQ_PREFIX".bam" \
--vcf $SAMPLE/TEMP/$FASEQ_PREFIX".VerifyBamID.vcf" \
--out $SAMPLE/REPORTS/verifyBamID/$FASEQ_PREFIX \
--precise \
--verbose \
--maxDepth 1000

# 26 I.01_TITV_ALL.sh
 filter the on bait SNV file to TiTv

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
--excludeFiltered \
-L $TI_TV_BED \
--variant $SAMPLE/SNV/QC/FILTERED_ON_BAIT/$FASEQ_PREFIX"_QC_OnBait.vcf" \
-o $SAMPLE/TEMP/$FASEQ_PREFIX"_QC_FILTERED_TiTv_SNV.vcf"

## so now you write out the TiTv ratio for all of those jokes;

$SAMTOOLS_DIR/bcftools/vcfutils.pl qstats $SAMPLE/TEMP/$FASEQ_PREFIX"_QC_FILTERED_TiTv_SNV.vcf" \
>| $SAMPLE/REPORTS/TI_TV/$FASEQ_PREFIX"_All_.titv.txt"

# 27 I.02_TITV_KNOWN.sh
## Below is writing out a vcf file of "known" variants (as defined by dbsnp_138.b37.excluding_sites_after_129.vcf).

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
--variant $SAMPLE/SNV/QC/FILTERED_ON_BAIT/$FASEQ_PREFIX"_QC_OnBait.vcf" \
--excludeFiltered \
-L $TI_TV_BED \
--concordance /isilon/sequencing/GATK_resource_bundle/2.8/b37/dbsnp_138.b37.excluding_sites_after_129.vcf \
-o $SAMPLE/TEMP/$FASEQ_PREFIX"_QC_FILTERED_TiTv_SNV_KNOWN.vcf"

# Now TiTv on those

$SAMTOOLS_DIR/bcftools/vcfutils.pl qstats $SAMPLE/TEMP/$FASEQ_PREFIX"_QC_FILTERED_TiTv_SNV_KNOWN.vcf" \
>| $SAMPLE/REPORTS/TI_TV/$FASEQ_PREFIX"_Known_.titv.txt"


#28 I.03_TITV_NOVEL.sh
# Below is writing out a vcf file of "novel" variants (as defined by dbsnp_138.b37.excluding_sites_after_129.vcf).

java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
--variant $SAMPLE/SNV/QC/FILTERED_ON_BAIT/$FASEQ_PREFIX"_QC_OnBait.vcf" \
--excludeFiltered \
-L $TI_TV_BED \
--discordance /isilon/sequencing/GATK_resource_bundle/2.8/b37/dbsnp_138.b37.excluding_sites_after_129.vcf \
-o $SAMPLE/TEMP/$FASEQ_PREFIX"_QC_FILTERED_TiTv_SNV_NOVEL.vcf"

# Now TiTv on those

$SAMTOOLS_DIR/bcftools/vcfutils.pl qstats $SAMPLE/TEMP/$FASEQ_PREFIX"_QC_FILTERED_TiTv_SNV_NOVEL.vcf" \
 >| $SAMPLE/REPORTS/TI_TV/$FASEQ_PREFIX"_Novel_.titv.txt"
