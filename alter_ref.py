import os
import numpy as np
import random
import sys

if len(sys.argv) < 7:
    print('please provide the name of the alternative strain')
    print('Usage: python alter_ref.py new_strainName start end numSNPs path_to_ref ouput_dir_path')
    exit()

strainName = sys.argv[1]
strain = f'>gi|{strainName}|'

start = int(sys.argv[2])
end = int(sys.argv[3])
numSNPs = int(sys.argv[4])

refPath = sys.argv[5]
outputDirPath = sys.argv[6]
alteredRefPath = f'{outputDirPath}/{strainName}.fasta'
alteredRefIndexPath = f'{outputDirPath}/{strainName}.subs.index'



# open reference file
with open(refPath, 'r') as refFile:
    next(refFile)
    refSeq = refFile.read().replace('\n', '')

# make sequence a list and prepare indexList
refSeq = list(refSeq)
indexList = []

# do substitutions
i = 0
while i < numSNPs:

    ind = int(random.uniform(start, end))
    allele = random.choice(['A', 'C', 'T', 'G'])

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
