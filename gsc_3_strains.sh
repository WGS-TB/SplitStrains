#!/bin/bash
# gsc.sh stands for Generate Separate and Consensus
# This is a master script to generate samples of different proportions, run splitStrains on each sample,
# create consensus seq and finally check separation results
set -e

####################################
#### Read generation parameters ####
####################################


export depth=300        # fold coverage
export numSNP=30        # number of snps to alter in the reference genome, changing this requires running alter_ref.py
export start=100000     # where to start on the ref
export end=3500000

export strainNameA=3_strainA_${numSNP}      # name of the altered strain (major)
export strainNameB=3_strainB_${numSNP}      # name of the altered strain (minor)
export strainNameC=3_strainC_${numSNP}

export dataDir=data/mixed_synth_3_samples_${numSNP}snps             # path where to output all fastQ files
export refDir=refs      # directory for all references references
export ref_strainA=${dataDir}/${strainNameA}.fasta   # path to a reference of a major strain
export ref_strainB=${dataDir}/${strainNameB}.fasta   # path to a reference of a minor strain
export ref_strainC=${dataDir}/${strainNameC}.fasta

export artOutput=${dataDir}/strain_                                 # prefix of the output from art read generator, these files are tmp


#######################################################
#### Trim, alignment and splitStrains parameters ######
#######################################################

export TRIMMOMATIC_PATH=~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar
export REF_PATH=${refDir}/tuberculosis.fna
export GFF_PATH=${refDir}/tuberculosis.gff
export REF_BOWTIE_PATH={refDir}/bowtie-index-tuberculosis
export SAMPLE_PATH=$dataDir
# export SAMPLE_PATH=./data/mixed_data

alterRef=0       # create altered references
genReads=0       # generate reads
trimQ=16
reuse=1
entropyFilter=1
doAlignment=0
depthScale=0.75
entropyStep=60
split=0
components=3
bowtie=0
model='gmm'

# check if need to generate new references
if [ $alterRef -eq 1 ]; then
    echo "altering reference"
    python alter_ref.py $strainNameA $start $end $numSNP $REF_PATH $dataDir
    python alter_ref.py $strainNameB $start $end $numSNP $REF_PATH $dataDir
    python alter_ref.py $strainNameC $start $end $numSNP $REF_PATH $dataDir

    echo "bwa index refs"
    bwa index $ref_strainA
    bwa index $ref_strainB
    bwa index $ref_strainC
fi


mix=( 35-50-15 ) # distance between the means increased
# mix=( 10-25-65 15-30-55 20-35-45 25-40-35 30-45-25 35-50-15 ) # distance between the means increased

for id in ${mix[@]}
do
    resultFile="results-id${id}_trim${trimQ}.txt"
    sampleR1="art-sample_${id}_R1.fq"         # synth sample
    sampleR2="art-sample_${id}_R2.fq"         # synth sample

    if [ $genReads -eq 1 ]; then
        ./gen_reads_3_strain.sh $id
    fi

    ./runSplitStrains.sh $trimQ $id $reuse $entropyFilter $doAlignment $resultFile $depthScale $start $end $entropyStep $split $sampleR1 $sampleR2 $components $bowtie $model
    # ./consensus.sh $id $trimQ
done
