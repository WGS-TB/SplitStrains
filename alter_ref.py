import os
import numpy as np
import random
import sys

if len(sys.argv) < 7:
    print('please provide the name of the alternative strain')
    print('Usage: python alter_ref.py new_strainName start end numSNPs path_to_ref path_to_intervals ouput_dir_path')
    exit()

strainName = sys.argv[1]
strain = f'>gi|{strainName}|'

start_genome = int(sys.argv[2])
end_genome = int(sys.argv[3])
numSNPs = int(sys.argv[4])

refPath = sys.argv[5]

intervalPath = sys.argv[6]

outputDirPath = sys.argv[7]

alteredRefPath = f'{outputDirPath}/{strainName}.fasta'
alteredRefIndexPath = f'{outputDirPath}/{strainName}.subs.index'



# open reference file
with open(refPath, 'r') as refFile:
    next(refFile)
    refSeq = refFile.read().replace('\n', '')


intervals = []
# open intervals file and populate inervals list with included intervals
with open(intervalPath, 'r') as intervalFile:
    lines = intervalFile.readlines()
    for line in lines:
        line_split = line.split(',')
        start = int(line_split[0])
        end = int(line_split[1].rstrip())
        intervals.append([start, end])


# make sequence a list and prepare indexList
refSeq = list(refSeq)
indexList = []

# do substitutions
i = 0
while i < numSNPs:
    is_inside_legal_interval = False
    ind = int(random.uniform(start_genome, end_genome))
    allele = random.choice(['A', 'C', 'T', 'G'])

    # for every interval check if the ind is inside
    print('checking if is inside legal interval')
    for interval in intervals:
        if interval[0] < ind and ind < interval[1]:
            print('ind is inside legal interval. ok.')
            is_inside_legal_interval = True

    # do the check if the ind falls into the ignored gff region
    if not is_inside_legal_interval:
        continue

    # check if char from ref fasta is letter if not skip
    if not refSeq[ind].isalpha() or refSeq[ind] == allele:
        continue
    # if allele is the same letter as reference skip
    elif refSeq[ind] == allele:
        continue

    else:
        refSeq[ind] = allele
        indexList.append(ind)
        i += 1

refSeq = ''.join(refSeq)
indexList = sorted(indexList)
indexList = ''.join(str(item)+'\n' for item in indexList)

lineLength = 81

with open(alteredRefPath, 'w') as f:
    # write genome name
    f.write(strain + '\n')
    # write lines of genome
    for i in range(0, len(refSeq), lineLength):
        f.write(refSeq[i:i+lineLength]+'\n')
    # write the remainding line
    f.write(refSeq[i+lineLength:])

with open(alteredRefIndexPath, 'w') as f:
    f.write(indexList)
