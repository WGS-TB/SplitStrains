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
export REF_PATH=${refDir}/tuberculosis.fna
export GFF_PATH=${refDir}/tuberculosis.filtered.gff
export REF_BOWTIE_PATH={refDir}/bowtie-index-tuberculosis
export SAMPLE_PATH=$dataDir

####################################
#### Read generation parameters ####
####################################

alterAndGenReads=0       # create altered references and generate new reads

export depth=80          # fold coverage
export numSNP=100        # number of snps to alter in the reference genome
export start=50000       # where to start on the ref
export end=4400000

export strainNameA=strainA_${numSNP}snp      # name of the altered strain (major)
export strainNameB=strainB_${numSNP}snp      # name of the altered strain (minor)

export refDir=refs                                  # directory for all references references
export ref_strainA=${refDir}/${strainNameA}.fasta   # path to a reference of a major strain
export ref_strainB=${refDir}/${strainNameB}.fasta   # path to a reference of a minor strain

export dataDir=data/mixed_synth_samples_${numSNP}snps        # path where to output all synth fastQ files
# export dataDir=data/mixed_data                             # path to real samples. Comment this to generate synth samples
export artOutput=${dataDir}/strain_                          # prefix of the output from art read generator, these files are tmp


#######################################################
#### Trim, alignment and splitStrains parameters ######
#######################################################

export READLENGTH=40    # trimmomatic minimum read length

alterRef=0             # create altered references. Skip if 0
genReads=1             # generate reads. Skip if 0
doAlignment=1          # align generated reads
trimQ=16               # parameter for Trimmomatic
reuse=0                # splitStrains will reuse the csv from prev run. Set to 1 after the first run!
entropyFilter=0.70     # default 0.7
depthScale=0.75        # default 0.75
entropyStep=100        # good value 60
split=0                # attempt to split strains
components=2           # number of strains
bowtie=0               # use bowtie2 as an aligner.
model='gmm'            # model for clustering. Two options: bmm and gmm (recommended)

# check if need to generate new fasta references
if [ $alterRef -eq 1 ]; then
    echo "altering reference"
    python alter_ref.py $strainNameA $start $end $numSNP $REF_PATH $refDir
    python alter_ref.py $strainNameB $start $end $numSNP $REF_PATH $refDir
    echo "bwa index refs"
    bwa index $ref_strainA
    bwa index $ref_strainB
fi


# The project masks fastq files as sample{id}_1.fastq.gz and sample{id}_2.fastq.gz
# mix is an array of ids for fastq files

# mix=( 95 90 80 70 60 55 50 )
# mix=( `seq 51 60` )
mix=( `seq 50` )

for id in ${mix[@]}
do
    resultFile="results-id${id}_trim${trimQ}.txt"

    sampleR1="art-sample_${id}_R1.fq"         # synth sample
    sampleR2="art-sample_${id}_R2.fq"         # synth sample
    # sampleR1="sample${id}_1.fastq.gz"       # real data sample
    # sampleR2="sample${id}_2.fastq.gz"       # real data sample

    if [ $genReads -eq 1 ]; then
        ./gen_reads.sh $id
    fi

    ./runSplitStrains.sh $trimQ $id $reuse $entropyFilter $doAlignment $resultFile $depthScale $start $end $entropyStep $split $sampleR1 $sampleR2 $components $bowtie $model
    # ./consensus.sh $id $trimQ
done
