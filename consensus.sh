#!/bin/bash
# This scripts builds consensus sequence out of separated minor and major bam files
set -e
echo "##### consensus script started #####"
id=$1
trimQ=$2

mainBam=${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam
mkdir -p $SAMPLE_PATH/output/${id}_${trimQ}/typing

# get output reads file
reads=( `ls $SAMPLE_PATH/output/${id}_${trimQ}/*.reads` )
echo ${reads[0]}
echo ${reads[1]}

# create major.bam and minor.bam
python rmreads.py ${reads[0]} $mainBam ${SAMPLE_PATH}/output/${id}_${trimQ}/major.bam
python rmreads.py ${reads[1]} $mainBam ${SAMPLE_PATH}/output/${id}_${trimQ}/minor.bam
# index bams
samtools index ${SAMPLE_PATH}/output/${id}_${trimQ}/major.bam
samtools index ${SAMPLE_PATH}/output/${id}_${trimQ}/minor.bam

# sort bam by name (needed for bedtools)
samtools sort -n -o ${SAMPLE_PATH}/output/${id}_${trimQ}/major.namesorted.bam ${SAMPLE_PATH}/output/${id}_${trimQ}/major.bam
samtools sort -n -o ${SAMPLE_PATH}/output/${id}_${trimQ}/minor.namesorted.bam ${SAMPLE_PATH}/output/${id}_${trimQ}/minor.bam

# bam to paired end fastq
bedtools bamtofastq -i ${SAMPLE_PATH}/output/${id}_${trimQ}/major.namesorted.bam \
    -fq ${SAMPLE_PATH}/output/${id}_${trimQ}/major_1.fq -fq2 ${SAMPLE_PATH}/output/${id}_${trimQ}/major_2.fq

bedtools bamtofastq -i ${SAMPLE_PATH}/output/${id}_${trimQ}/minor.namesorted.bam \
    -fq ${SAMPLE_PATH}/output/${id}_${trimQ}/minor_1.fq -fq2 ${SAMPLE_PATH}/output/${id}_${trimQ}/minor_2.fq



# do calls on isolates
bcftools mpileup -Ou -f $REF_PATH ${SAMPLE_PATH}/output/${id}_${trimQ}/major.bam | bcftools call -mv -Oz -o ${SAMPLE_PATH}/output/${id}_${trimQ}/major-calls.vcf.gz
bcftools mpileup -Ou -f $REF_PATH ${SAMPLE_PATH}/output/${id}_${trimQ}/minor.bam | bcftools call -mv -Oz -o ${SAMPLE_PATH}/output/${id}_${trimQ}/minor-calls.vcf.gz
bcftools index ${SAMPLE_PATH}/output/${id}_${trimQ}/major-calls.vcf.gz
bcftools index ${SAMPLE_PATH}/output/${id}_${trimQ}/minor-calls.vcf.gz
# filter out indels
bcftools view --exclude-types indels ${SAMPLE_PATH}/output/${id}_${trimQ}/major-calls.vcf.gz -Ob -o ${SAMPLE_PATH}/output/${id}_${trimQ}/major-calls.noindel.bcf
bcftools view --exclude-types indels ${SAMPLE_PATH}/output/${id}_${trimQ}/minor-calls.vcf.gz -Ob -o ${SAMPLE_PATH}/output/${id}_${trimQ}/minor-calls.noindel.bcf

bcftools index ${SAMPLE_PATH}/output/${id}_${trimQ}/major-calls.noindel.bcf
bcftools index ${SAMPLE_PATH}/output/${id}_${trimQ}/minor-calls.noindel.bcf
# create consensus seq
cat $REF_PATH | bcftools consensus ${SAMPLE_PATH}/output/${id}_${trimQ}/major-calls.noindel.bcf > ${SAMPLE_PATH}/output/${id}_${trimQ}/major-consensus.fa
cat $REF_PATH | bcftools consensus ${SAMPLE_PATH}/output/${id}_${trimQ}/minor-calls.noindel.bcf > ${SAMPLE_PATH}/output/${id}_${trimQ}/minor-consensus.fa
