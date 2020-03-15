# this script filters a gff file by removing entries based on locus_tag field
# if locus_tag field is in excluded list, the entry is ignored

fileEx = '/home/user1/Documents/lab/DrugResistance/splitStrains/refs/excluded2.txt'
fileIn = '/home/user1/Documents/lab/DrugResistance/splitStrains/refs/tuberculosis.gff'
fileOut = '/home/user1/Documents/lab/DrugResistance/splitStrains/refs/tuberculosis.filtered.gff'

excluded = []
processedGenes = []

# open the list of genes to exclude
with open(fileEx, 'r') as e:
    for line in e:
        splitLine = line.split()
        excluded.append(splitLine[0])

print(excluded, ' ,size: ', len(excluded))

count = 0

# open files for reading and writing
with open(fileIn, 'r') as f:
    with open(fileOut, 'w') as o:

        for line in f:

            splitLine = line.split()
            # write headers into output
            if splitLine[0][0] == '#':
                o.write(line)
                continue

            # if not the header do other checks
            else:
                regionStart = splitLine[3]
                regionEnd = splitLine[4]
                fields = splitLine[8].split(';')    # get variable like fields
                locus_tag = fields[-1].split('=')[-1]   # get locus tag

                if locus_tag not in excluded:
                    o.write(line)
                else:
                    processedGenes.append(locus_tag)
                    count += 1

print('excluded count: ', count)

notProcessed = [set(excluded) - set(processedGenes)]
print(notProcessed)
