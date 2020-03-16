#!/bin/bash
set -e
id=$1

echo "#### Working on proportion $id"
proportionA=`python -c "print(int($depth*0.01*$id))"`
proportionB=`python -c "print(int($depth*(1-0.01*$id)))"`
# proportionC=`python -c "print(int($depth*0.5))"`

echo "generate reads"
art_illumina -ss HSXn -i $ref_strainA -o ${artOutput}${strainNameA}_R -l 150 -qL 21 -qU 31 -qs -9 -qs2 -8 -na -f $proportionA -p -s 10 -m 300
art_illumina -ss HSXn -i $ref_strainB -o ${artOutput}${strainNameB}_R -l 150 -qL 21 -qU 31 -qs -9 -qs2 -8 -na -f $proportionB -p -s 10 -m 300

# concat
cat ${artOutput}${strainNameA}_R1.fq ${artOutput}${strainNameB}_R1.fq > ${SAMPLE_PATH}/art-sample_${id}_R1.fq
cat ${artOutput}${strainNameA}_R2.fq ${artOutput}${strainNameB}_R2.fq > ${SAMPLE_PATH}/art-sample_${id}_R2.fq


# get rid of original reads
rm ${artOutput}${strainNameA}_R1.fq ${artOutput}${strainNameA}_R2.fq
rm ${artOutput}${strainNameB}_R1.fq ${artOutput}${strainNameB}_R2.fq
