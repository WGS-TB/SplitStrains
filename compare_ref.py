import os
import sys


def compare_ref(id, originalRefPath, altRefPath, consensusPath):
    """ This function compares consensus seq from splitStains and bcftools with the target genome """


    with open(originalRefPath, 'r') as f:
        next(f)
        originalRefData = f.read().replace('\n', '')

    with open(altRefPath, 'r') as f:
        next(f)
        altRefData = f.read().replace('\n', '')

    altRefIndexPath = altRefPath.split('.')[0] + '.subs.index'
    with open(altRefIndexPath, 'r') as f:
        altRefSubsIndex = f.read()

    with open(consensusPath, 'r') as f:
        next(f)
        consensusData = f.read().replace('\n', '')

    altRefSubsIndex = altRefSubsIndex.split()
    altRefSubsIndex = [int(x) for x in altRefSubsIndex]

    refDifList = []
    i = 0
    if len(originalRefData) != len(altRefData):
        print('original ref length and alt reference length do not match')
        exit()

    """ Compare altered and original genome references. The number of differences should be exactly
        the number of introduced SNPs """

    for x,y in zip(originalRefData, altRefData):
        if x != y:
            refDifList.append(i)
        i += 1

    # print('# of differences between alt ref and original ref:', len(refDifList))

    """ Compare a consensus sequence obtained from splitStrain.py and bcftools with the altered reference.
        Idealy, this should match. However, since splitStrains filters out bad positions based on different samflags and
        quality settings it is very likely that there would be some some differences. """

    consenDifList = []
    i = 0

    for x,y in zip(altRefData, consensusData):
        if x != y:
            consenDifList.append(i)
        i += 1

    # print('# of differences between consensus and alt ref:', len(consenDifList))
    # print(consenDifList)

    """ Check if consensus non matching positions are of those which were introduced by aler_ref.py. Here we find if the SNPs introduced by alter_ref.property are resolved. This is the main metric of success. The error will increase as proportions approach 50:50 ir 10:90 splits """
    count = 0
    unresolvedList = []
    for pos in consenDifList:
        # if the pos belongs to the true SNP by alterRef.py then count it as unresolved.
        if pos in altRefSubsIndex:
            unresolvedList.append(pos)
            count += 1

    # print('# of introduced SNPs that were not resolved by splitStrains.py:', count)
    # print(unresolvedList)
    return unresolvedList

if __name__ == '__main__':
    id = sys.argv[1]
    originalRefPath = sys.argv[2]
    altRefPath = sys.argv[3]
    consensusPath = sys.argv[4]
    compare_ref(id, originalRefPath, altRefPath, consensusPath)
