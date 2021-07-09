#!/bin/bash
# gsc.sh stands for Generate Separate and Consensus
# This is a master script that does the following:
#   1) Generate alternative fasta genomes
#   2) generate samples of different proportions or use real data samples
#   3) run splitStrains on each sample
#   4) create consensus seq and finally check separation results
set -e

####################################
#### Paths to refs and software ####
####################################
export TRIMMOMATIC_PATH=~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar

export refDir=refs                                    # directory for all references references
# export REF_PATH=${refDir}/tuberculosis-mixInfect-ref.fasta
export REF_PATH=${refDir}/tuberculosis.fna
export GFF_PATH=${refDir}/tuberculosis.filtered-intervals.gff
export included_intervals=${refDir}/included_intervals.txt

####################################
#### Read generation parameters ####
####################################
export numSNP=99         # number of snps to alter in the reference genome
export depth=60           # fold coverage
export start=4000         # where splitstrains will start on the ref
export end=4400000         # where splitstrains ends on the ref


####################################
#### Paths to data directory    ####
####################################
export SAMPLE_PATH=data/mixed_synth_samples_${numSNP}snps             # path where to output all synth fastQ files
# export SAMPLE_PATH=data/mixed_data                                    # path to real samples.
# export SAMPLE_PATH=data/mixInfect_data/original_fastq                 # path to mixInfect samples.
# export SAMPLE_PATH=data/mentalist_data_1                              # path to mentalist samples.
# export SAMPLE_PATH=data/mentalist_data_2                              # path to mentalist samples.
# export SAMPLE_PATH=data/mentalist_data_3                              # path to mentalist samples.
# export SAMPLE_PATH=data/mentalist_data_pure_1                           # path to mentalist pure samples.
# export SAMPLE_PATH=data/mentalist_data_pure_2                           # path to mentalist pure samples.
# export SAMPLE_PATH=data/mentalist_data_pure_3                           # path to mentalist pure samples.


####################################
#### Strains naming vars        ####
####################################
export strainNameA=strainA_${numSNP}snp      # name of the altered strain (major)
export strainNameB=strainB_${numSNP}snp      # name of the altered strain (minor)

export ref_strainA=${SAMPLE_PATH}/refs/${strainNameA}.fasta   # path to a reference of a major strain
export ref_strainB=${SAMPLE_PATH}/refs/${strainNameB}.fasta   # path to a reference of a minor strain
export artOutput=${SAMPLE_PATH}/strain_                          # prefix of the output from art read generator, these files are tmp
mkdir -p ${SAMPLE_PATH}/refs


#######################################################
#### Trim, alignment and splitStrains parameters ######
#######################################################

export READLENGTH=40    # trimmomatic minimum read length

alterRef=1             # create altered references. Skip if 0
genReads=1             # generate reads. Skip if 0
doAlignment=1          # align generated reads. Skip if 0
trimQ=16               # parameter for Trimmomatic
reuse=0                # splitStrains will reuse the csv from prev run. Set to 1 after the first run!
entropyFilter=0        # default 0.7
minDepth=0.7            # default 0.75
entropyStep=60         # good value 60
split=0                # attempt to split strains
components=2           # number of strains
bowtie=0               # use bowtie2 as an aligner (not recommended)
model='gmm'            # model for clustering. Two options: bmm and gmm (recommended)
alpha=0.05              # default 0.05

# check if need to generate new fasta references
if [ $alterRef -eq 1 ]; then
    echo "altering reference"
    python alter_ref.py $strainNameA $start $end $numSNP $REF_PATH $included_intervals ${SAMPLE_PATH}/refs
    python alter_ref.py $strainNameB $start $end $numSNP $REF_PATH $included_intervals ${SAMPLE_PATH}/refs
    echo "bwa index refs"
    bwa index $ref_strainA
    bwa index $ref_strainB
fi


mkdir -p $SAMPLE_PATH
mkdir -p $SAMPLE_PATH/aligned   # directory for bam files
mkdir -p $SAMPLE_PATH/trimmed   # directory for trimmed files
mkdir -p $SAMPLE_PATH/output    # directory for splitStrains.py output
mkdir -p $SAMPLE_PATH/calls     # directory for bcftools calls

# mix=( 10 20 30 40 50 60 70 80 90 100 )     # Pedro pure data different coverage
# mix=( 50 55 60 65 70 75 80 85 90 95 )      # Pedro mixed synth data
# mix=( 95 90 )            # my synth data
mix=( 95 90 85 80 70 60 55 50 )            # my synth data
# mix=( `seq 60` )                           # Inaki data
# mix=( `seq 48` )                           # mixInfect data
# mix=( 6 7 8 9 10 12 16 17 19 20 )            # mixInfect samples for separation

for id in ${mix[@]}
do
    # output log file from splitStrains stdout
    resultFile="results-id${id}_trim${trimQ}.txt"

    # The project masks fastq files as sample{id}_1.fastq.gz and sample{id}_2.fastq.gz
    sampleR1="art-sample_${id}_R1.fq"         # synth sample
    sampleR2="art-sample_${id}_R2.fq"         # synth sample
    # sampleR1="sample${id}_1.fastq.gz"       # real data sample, mixInfect and mentalist data
    # sampleR2="sample${id}_2.fastq.gz"       # real data sample, mixInfect and mentalist data

    # generate reads using art_illumina
    if [ $genReads -eq 1 ]; then
        ./gen_reads.sh $id
    fi

    # run the script that does trimming, alignment and splitStrains.py
    ./runSplitStrains.sh $trimQ $id $reuse $entropyFilter $doAlignment $resultFile $minDepth $start $end $entropyStep $split $sampleR1 $sampleR2 $components $bowtie $model $alpha
    # ./consensus.sh $id $trimQ
done
