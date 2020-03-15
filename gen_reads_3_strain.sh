#!/bin/bash
set -e
id=$1

p1=${id:0:2}
p2=${id:3:2}
p3=${id:6:2}

proportionA=`python -c "print(int($depth*0.01*$p1))"`
proportionB=`python -c "print(int($depth*0.01*$p2))"`
proportionC=`python -c "print(int($depth*0.01*$p3))"`

echo "depth for 3 strains: " $proportionA $proportionB $proportionC

echo "generate reads"
art_illumina -ss HSXn -i $ref_strainA -o ${artOutput}${strainNameA}_R -l 150 -qL 21 -qU 31 -qs -9 -qs2 -8 -na -f $proportionA -p -s 10 -m 300
art_illumina -ss HSXn -i $ref_strainB -o ${artOutput}${strainNameB}_R -l 150 -qL 21 -qU 31 -qs -9 -qs2 -8 -na -f $proportionB -p -s 10 -m 300
art_illumina -ss HSXn -i $ref_strainC -o ${artOutput}${strainNameC}_R -l 150 -qL 21 -qU 31 -qs -9 -qs2 -8 -na -f $proportionC -p -s 10 -m 300

# concat
cat ${artOutput}${strainNameA}_R1.fq ${artOutput}${strainNameB}_R1.fq ${artOutput}${strainNameC}_R1.fq > ${dataDir}/art-sample_${p1}-${p2}-${p3}_R1.fq
cat ${artOutput}${strainNameA}_R2.fq ${artOutput}${strainNameB}_R2.fq ${artOutput}${strainNameC}_R2.fq > ${dataDir}/art-sample_${p1}-${p2}-${p3}_R2.fq

# get rid of original reads
rm ${artOutput}${strainNameA}_R1.fq ${artOutput}${strainNameA}_R2.fq
rm ${artOutput}${strainNameB}_R1.fq ${artOutput}${strainNameB}_R2.fq
rm ${artOutput}${strainNameC}_R1.fq ${artOutput}${strainNameC}_R2.fq
