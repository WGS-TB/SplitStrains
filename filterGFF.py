# this script filters a gff file by removing entries based on locus_tag field
# if locus_tag field is in excluded list, the entry is ignored


def isCloseInterval(start, end, excludedInterval, startDist=50, endDist=50):
    """ This function takes an interval and compares it with intervals in the excludedIntervals list.
    If the end points are close enough to some interval in the list return True """

    for interval in excludedIntervals:
        if abs(interval[0]-start) < startDist and abs(interval[1]-end) < endDist:
            # print('close to the excluded interval')
            return True

    return False



fileEx = '/home/user1/Documents/lab/SplitStrains/refs/excluded2.txt'
fileExIntervals = '/home/user1/Documents/lab/SplitStrains/refs/excluded-intervals.txt'

fileIn = '/home/user1/Documents/lab/SplitStrains/refs/tuberculosis.gff'
fileOut = '/home/user1/Documents/lab/SplitStrains/refs/tuberculosis.filtered-intervals.gff'

# this txt file will only contain included region coordinates
file_only_included_intervals = '/home/user1/Documents/lab/SplitStrains/refs/included_intervals.txt'

excluded = []
excludedIntervals = []


processedGenes = []

# open the list of genes to exclude
with open(fileEx, 'r') as e:
    for line in e:
        splitLine = line.split()
        excluded.append(splitLine[0])

# open the list of intervals to exclude
with open(fileExIntervals, 'r') as e:
    for line in e:
        splitLine = line.rstrip().split('\t')
        excludedIntervals.append([int(splitLine[0]), int(splitLine[1])])

count = 0

# open files for reading and writing
with open(fileIn, 'r') as f:
    with open(fileOut, 'w') as o:
            with open(file_only_included_intervals, 'w') as o_intervals:

                for line in f:

                    splitLine = line.split()
                    # write headers into output
                    if splitLine[0][0] == '#':
                        o.write(line)
                        continue

                    # if not the header do other checks
                    else:
                        interval_type = splitLine[2]
                        regionStart = int(splitLine[3])
                        regionEnd = int(splitLine[4])
                        fields = splitLine[8].split(';')    # get variable like fields
                        locus_tag = fields[-1].split('=')[-1]   # get locus tag

                        if locus_tag not in excluded and not isCloseInterval(regionStart, regionEnd, excludedIntervals):
                                o.write(line)

                                # write all included intervals into a txt file
                                if interval_type=='gene':
                                    print('writing gene:', interval_type)
                                    o_intervals.write(f'{regionStart}, {regionEnd}\n')

                        else:
                            processedGenes.append(locus_tag)
                            count += 1

print('excluded count: ', count)

notProcessed = [set(excluded) - set(processedGenes)]
# print(notProcessed)
