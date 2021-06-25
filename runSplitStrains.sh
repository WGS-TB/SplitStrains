#!/bin/bash
set -e
# This scripts defines splitStrains.py pipeline.
# The script is not meant to be ran directly; instead batchData.sh should be used to call this script.
trimThreshold=$1
id=$2
reuse=$3
entropyFilter=$4
doAlignment=$5
resultFile=$6
mindDepth=$7
startRegion=$8
endRegion=$9
entropyStep=${10}
split=${11}
sampleR1=${12}
sampleR2=${13}
components=${14}
bowtie=${15}
model=${16}
alpha=${17}

if [ "$#" != "17" ]
then
    echo "Illegal number of parameters, provided: $# needed: 16"
    exit
fi

sampleR1_trimmed=${SAMPLE_PATH}/trimmed/${sampleR1}.fastq.gz
sampleR2_trimmed=${SAMPLE_PATH}/trimmed/${sampleR2}.fastq.gz

sampleR1_BWA_out=${SAMPLE_PATH}/aligned/trimmed-1-${id}.sai
sampleR2_BWA_out=${SAMPLE_PATH}/aligned/trimmed-2-${id}.sai
sampleBAM_out=${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam

if [ $doAlignment == 1 ]
then
    adapter=/home/user1/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa
    echo "########################### trimming start: ###############################"
    java -jar $TRIMMOMATIC_PATH PE -phred33 $SAMPLE_PATH/$sampleR1 $SAMPLE_PATH/$sampleR2 \
        $sampleR1_trimmed $sampleR1_trimmed.se \
        $sampleR2_trimmed $sampleR2_trimmed.se \
        ILLUMINACLIP:$adapter:3:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:$trimThreshold MINLEN:$READLENGTH


    echo "########################### alignment start: ###############################"
    if [ $bowtie == 1 ]
    then
        bowtie2 -x $REF_BOWTIE_PATH --very-sensitive -1 $sampleR1_trimmed -2 $sampleR2_trimmed -S ${SAMPLE_PATH}/aligned/bowtie-sample-${id}.sam
        samtools view -S -b ${SAMPLE_PATH}/aligned/bowtie-sample-${id}.sam | samtools sort -o ${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam -
        samtools index ${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam
    else
        bwa aln -q 15 -t 6 -R '@RG\tID:mixed\tSM:mixed\tLB:None\tPL:Illumina' $REF_PATH $sampleR1_trimmed > $sampleR1_BWA_out

        bwa aln -q 15 -t 6 -R '@RG\tID:mixed\tSM:mixed\tLB:None\tPL:Illumina' $REF_PATH $sampleR2_trimmed > $sampleR2_BWA_out

        bwa sampe $REF_PATH $sampleR1_BWA_out $sampleR2_BWA_out $sampleR1_trimmed $sampleR2_trimmed | \
        samtools view -bhS - | \
        samtools sort -n -o ${SAMPLE_PATH}/aligned/trimmed-${id}.name.sorted.bam -
        samtools fixmate -m -r ${SAMPLE_PATH}/aligned/trimmed-${id}.name.sorted.bam ${SAMPLE_PATH}/aligned/trimmed-${id}.fixm.sorted.bam
        samtools sort -o ${SAMPLE_PATH}/aligned/trimmed-${id}.pos.sorted.bam ${SAMPLE_PATH}/aligned/trimmed-${id}.fixm.sorted.bam
        samtools markdup -r ${SAMPLE_PATH}/aligned/trimmed-${id}.pos.sorted.bam ${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam

        rm $sampleR1_BWA_out \
            $sampleR2_BWA_out \
            ${SAMPLE_PATH}/aligned/trimmed-${id}.fixm.sorted.bam \
            ${SAMPLE_PATH}/aligned/trimmed-${id}.pos.sorted.bam \
            ${SAMPLE_PATH}/aligned/trimmed-${id}.name.sorted.bam

        # bwa mem -w 0 -t 8 -R '@RG\tID:mixed\tSM:mixed\tLB:None\tPL:Illumina' $REF_PATH \
        #     $sampleR1_trimmed $sampleR2_trimmed \
        #     | samtools view -b -S - | samtools sort -o ${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam -

        samtools index ${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam
    fi
fi

echo "########################### splitStrain section: ###############################"

outputDir=${id}_${trimThreshold}

if [ ! -d $SAMPLE_PATH/output/$outputDir ]
then
    mkdir $SAMPLE_PATH/output/$outputDir
fi

echo start = $startRegion, end = $endRegion
depthFile=bam_depth.txt

# if reuse is 1 then reuse a csv file from the prev run
if [ $reuse == 1 ]
then
    # Compute avg depth of the bam file
    samtools depth $SAMPLE_PATH/aligned/trimmed-${id}.sorted.bam > $SAMPLE_PATH/output/$outputDir/$depthFile
    avgDepth=`awk '{ total += $3; count++ } END { print int(total/count) }' $SAMPLE_PATH/output/$outputDir/$depthFile`
    echo 'bam avg depth: ' $avgDepth | tee $SAMPLE_PATH/output/$outputDir/$depthFile

    # look up the depth from the txt file
    avgDepth=`awk 'BEGIN {FS=":  "}{print($2)}' $SAMPLE_PATH/output/$outputDir/$depthFile`
    depth=`awk "BEGIN {print(int($avgDepth*$mindDepth))}"`

    if [ $split == 1 ]
    then
        python splitStrains.py -a $alpha -b $GFF_PATH -c -z -mo $model -fe $entropyFilter -fes $entropyStep -g $components -f sample-${id} -s $startRegion -e $endRegion -r $REF_PATH -o $SAMPLE_PATH/output/$outputDir -fd $depth ${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam | tee $SAMPLE_PATH/output/$resultFile
    else
        python splitStrains.py -a $alpha -b $GFF_PATH -z -mo $model -fe $entropyFilter -fes $entropyStep -g $components -f sample-${id} -s $startRegion -e $endRegion -r $REF_PATH -b $GFF_PATH -o $SAMPLE_PATH/output/$outputDir -fd $depth ${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam | tee $SAMPLE_PATH/output/$resultFile
    fi
else
    # Compute avg depth of the bam file
    samtools depth $SAMPLE_PATH/aligned/trimmed-${id}.sorted.bam > $SAMPLE_PATH/output/$outputDir/$depthFile
    avgDepth=`awk '{ total += $3; count++ } END { print int(total/count) }' $SAMPLE_PATH/output/$outputDir/$depthFile`
    echo 'bam avg depth: ' $avgDepth | tee $SAMPLE_PATH/output/$outputDir/$depthFile

    # if [ $avgDepth -lt 5 ]
    # then
    #     echo > $SAMPLE_PATH/output/$outputDir/skipped
    #     exit
    # fi

    depth=`awk "BEGIN {print(int($avgDepth*$mindDepth))}"`

    if [ $split == 1 ]
    then
        # no GFF
        # python splitStrains.py -a $alpha -c -mo $model -fe $entropyFilter -fes $entropyStep -g $components -f sample-${id} -s $startRegion -e $endRegion -r $REF_PATH -o $SAMPLE_PATH/output/$outputDir/ -fd $depth ${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam | tee $SAMPLE_PATH/output/$resultFile
        python splitStrains.py -a $alpha -b $GFF_PATH -c -mo $model -fe $entropyFilter -fes $entropyStep -g $components -f sample-${id} -s $startRegion -e $endRegion -r $REF_PATH -o $SAMPLE_PATH/output/$outputDir/ -fd $depth ${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam | tee $SAMPLE_PATH/output/$resultFile
    else
        # no GFF
        # python splitStrains.py -a $alpha -mo $model -fe $entropyFilter -fes $entropyStep -g $components -f sample-${id} -s $startRegion -e $endRegion -r $REF_PATH -b $GFF_PATH -o $SAMPLE_PATH/output/$outputDir -fd $depth ${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam | tee $SAMPLE_PATH/output/$resultFile
        python splitStrains.py -a $alpha -b $GFF_PATH -mo $model -fe $entropyFilter -fes $entropyStep -g $components -f sample-${id} -s $startRegion -e $endRegion -r $REF_PATH -b $GFF_PATH -o $SAMPLE_PATH/output/$outputDir -fd $depth ${SAMPLE_PATH}/aligned/trimmed-${id}.sorted.bam | tee $SAMPLE_PATH/output/$resultFile
    fi
fi
