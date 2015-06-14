# finna-be-octo-bear
## GATK Workflow

- custom gatk workflow adopted from CIDR
- Modify to suite the setting at local servers

1) Install bwa, picard, and gatk in $HOME/bin
2) Download gatk_bundle files from GATK ftp sites in $GATK_RESOURCE_PATH
   - this includes reference genome (hg19), known indels, dbsnp138
3) set up bwa index
4) put fastq files in $PROJECT_FOLDER/$FASTQ_DIR 
5) This script will create a subscript for each process in $SAMPLE/$SCRIPT
6) Change the number of threads to reflect your CPUs
7) Change the INPUT DATA below to match with your input
