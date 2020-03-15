import pysam
import numpy as np
import argparse

description = """ Given the output of splitStrains.py, remove reads from bam file.
                  Usage: python rmreads.py <reads_file> <input_bam> <output_bam> """


parser = argparse.ArgumentParser(description=description)
parser.add_argument(dest='filenames', metavar='files', nargs=3, help='reads_file input_bam output_bam')

args = parser.parse_args()

readsfile = args.filenames[0]
inputsam = args.filenames[1]
outputsam = args.filenames[2]


# open reads list as a numpy matrix
try:
    readsData = np.genfromtxt(readsfile, dtype='str', delimiter=':')
except Exception as err:
    raise SystemExit('ERROR:{}'.format(str(err)))

# create a set for fast access
readsSet = set(readsData[:,0])    # splitStrains.py output reads must be organized as columns
# open sam/bam file to do reading writing
try:
    samfile = pysam.AlignmentFile(inputsam, 'rb')
    filteredFile = pysam.AlignmentFile(outputsam, "wb", template=samfile)
except Exception as err:
    raise SystemExit('ERROR: Check bam file name or path.')

print('read removal is started')

for read in samfile.fetch():

    if read.query_name not in readsSet:
        filteredFile.write(read)

samfile.close()
filteredFile.close()

print('read removal is complete')
