import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from compare_ref import compare_ref
import os

TITLE_FONT_SIZE = 8
TICK_FONT_SIZE = 8
AXES_FONT_SIZE = 8
LABEL_FONT_SIZE = 8
DPI = 300
PLOT_ENTROPY = False

font = {'family' : 'normal',
        'size'   : 8}

matplotlib.rc('font', **font)

currentDir = os.getcwd()
dataDir = 'data/mixed_synth_samples_100snps/output'
pathToAltRef = f'{currentDir}/refs/'
strainNameA = 'strainA_100snp.fasta'
strainNameB = 'strainB_100snp.fasta'
output_figure_dir = '/home/user1/Documents/lab/SplitStrains/data/mixed_synth_samples_100snps/separation-error.png'
proportions = np.array([0.50, 0.55, 0.60, 0.70, 0.80, 0.85, 0.90])

# dataDir = 'data/mentalist_data_1/output'
# pathToAltRef = 'data/mentalist_data_1'
# strainNameA = 'seq01.fa'
# strainNameB = 'seq02.fa'
# output_figure_dir = '/home/user1/Documents/lab/SplitStrains/data/mentalist_data_1/separation-error.png'
# proportions = np.array([0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90])

ref_path = f'{currentDir}/refs/tuberculosis.fna'

data = {'major':[], 'minor':[]}

###### update this number
numSNP = 100
######

# proportions = np.array([0.90])

# print('number of intoruced SNPs by alt_ref:', numSNP, '. Make sure to check this value.')

for strainName in [strainNameA, strainNameB]:

    altRef_path = f'{pathToAltRef}/{strainName}'
    strain = 'major' if strainName == strainNameA else 'minor'

    for id in proportions:
        id = int(id*100)
        consensus_path = f'{dataDir}/{id}_16/{strain}-consensus.fa'
        non_matching_pos = compare_ref(id, ref_path, altRef_path, consensus_path)
        numErrors = len(non_matching_pos)
        data[strain].append(numErrors)
        print(f'############### StrainName:{strainName}, Proportion:{id}' ,'unresolved positions:',non_matching_pos)

major_errors = np.array(data['major']) / numSNP *100
minor_errors = np.array(data['minor']) / numSNP *100

sns.set(style="darkgrid")
sns.set_context("paper", rc={"lines.linewidth": 1.9, "patch.linewidth" : 0.5 })


fig, ax = plt.subplots()
fig.set_size_inches(7.5, 4)
# plt.tight_layout()

ax.set_title('Separation error over major strain percentage')
ax.set_ylabel('Percentage of incorrectly split SNPs')
ax.set_xlabel('Major strain proportion')

ax.set_xticks(proportions)

# do not count mismatches that occure due to alignment errors, offset compensate that for major strain
offset = 3

ax.plot(proportions, major_errors-offset, '--g')
# ax.plot(proportions, minor_errors-offset, '--b')
# ax.legend(['major', 'minor'])
ax.plot(proportions, major_errors-offset, '.g', linewidth=5, markersize=10)
# ax.plot(proportions, minor_errors-offset, '.b', linewidth=5, markersize=10)
y_ticks = ax.get_yticks()
ax.set_yticklabels([f'{int(i)}%' for i in y_ticks])
# x_ticks = ax.get_xticks()
# ax.set_xticklabels([f'{int(i)}%' for i in x_ticks])
plt.savefig(output_figure_dir, dpi=DPI)
plt.show()
