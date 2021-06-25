

export SAMPLE_PATH=data/mixInfect_data
export REF_PATH=$SAMPLE_PATH/tuber-ref-mixInf.fasta
export GFF_PATH=${SAMPLE_PATH}/tuberculosis.filtered.gff

export startRegion=100000      # where to start on the ref
export endRegion=4400000

id=ERR221663
# id=ERR245801
sample="${id}_gi_41353971_emb_AL123456.2_.bam"
outputDir=$id
model='gmm'
resultFile="result.txt"
reuse=0
entropyFilter=0.7  # default 0.7
depthScale=0.75     # default 0.75
entropyStep=100      # good value 60
components=2
model='gmm'

depthFile=${id}_bam_depth.txt

# #Compute avg depth of the bam file
samtools depth $SAMPLE_PATH/$sample > $SAMPLE_PATH/output/$outputDir/$depthFile
avgDepth=`awk '{ total += $3; count++ } END { print int(total/count) }' $SAMPLE_PATH/output/${outputDir}/$depthFile`
# echo 'bam avg depth: ' $avgDepth | tee $SAMPLE_PATH/output/${outputDir}/$depthFile
depth=`awk "BEGIN {print(int($avgDepth*$depthScale))}"`


if [ $reuse == 1 ]
then
    python splitStrains.py -z -b $GFF_PATH -mo $model -fe $entropyFilter -fes $entropyStep -g $components -f sample-${id} -s $startRegion -e $endRegion -r $REF_PATH -b $GFF_PATH -o $SAMPLE_PATH/output/$outputDir -fd $depth $SAMPLE_PATH/$sample | tee -a $SAMPLE_PATH/output/$resultFile
else
    python splitStrains.py -b $GFF_PATH -mo $model -fe $entropyFilter -fes $entropyStep -g $components -f sample-${id} -s $startRegion -e $endRegion -r $REF_PATH -b $GFF_PATH -o $SAMPLE_PATH/output/$outputDir -fd $depth $SAMPLE_PATH/$sample | tee -a $SAMPLE_PATH/output/$resultFile
fi
