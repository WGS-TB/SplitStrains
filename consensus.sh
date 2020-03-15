#!/bin/bash
# This scripts builds consensus sequence out of separated minor and major bam files
set -e
id=$1
trimQ=$2
dir=data/mixed_synth_samples
original_ref=refs/tuberculosis.fna
mainBam=${dir}/aligned/trimmed-${id}.sorted.bam

# get output reads file
reads=( `ls $dir/output/${id}_${trimQ}/*.reads` )
echo ${reads[0]}
echo ${reads[1]}

# create major.bam and minor.bam
python rmreads.py ${reads[0]} $mainBam ${dir}/output/${id}_${trimQ}/major.bam
python rmreads.py ${reads[1]} $mainBam ${dir}/output/${id}_${trimQ}/minor.bam
# index bams
samtools index ${dir}/output/${id}_${trimQ}/major.bam
samtools index ${dir}/output/${id}_${trimQ}/minor.bam
# do calls on isolates
bcftools mpileup -Ou -f $original_ref ${dir}/output/${id}_${trimQ}/major.bam | bcftools call -mv -Oz -o ${dir}/output/${id}_${trimQ}/major-calls.vcf.gz
bcftools mpileup -Ou -f $original_ref ${dir}/output/${id}_${trimQ}/minor.bam | bcftools call -mv -Oz -o ${dir}/output/${id}_${trimQ}/minor-calls.vcf.gz
bcftools index ${dir}/output/${id}_${trimQ}/major-calls.vcf.gz
bcftools index ${dir}/output/${id}_${trimQ}/minor-calls.vcf.gz
# filter out indels
bcftools view --exclude-types indels ${dir}/output/${id}_${trimQ}/major-calls.vcf.gz -Ob -o ${dir}/output/${id}_${trimQ}/major-calls.noindel.bcf
bcftools view --exclude-types indels ${dir}/output/${id}_${trimQ}/minor-calls.vcf.gz -Ob -o ${dir}/output/${id}_${trimQ}/minor-calls.noindel.bcf

bcftools index ${dir}/output/${id}_${trimQ}/major-calls.noindel.bcf
bcftools index ${dir}/output/${id}_${trimQ}/minor-calls.noindel.bcf
# create consensus seq
cat $original_ref | bcftools consensus ${dir}/output/${id}_${trimQ}/major-calls.noindel.bcf > ${dir}/output/${id}_${trimQ}/major-consensus.fa
cat $original_ref | bcftools consensus ${dir}/output/${id}_${trimQ}/minor-calls.noindel.bcf > ${dir}/output/${id}_${trimQ}/minor-consensus.fa
