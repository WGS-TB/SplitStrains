import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from compare_ref import compare_ref
import os

TITLE_FONT_SIZE = 8
TICK_FONT_SIZE = 8
AXES_FONT_SIZE = 8
LABEL_FONT_SIZE = 8
DPI = 300
PLOT_ENTROPY = False

# dataDir = 'data/mixed_synth_samples/output'
# strainNameA = 'strainA.fasta'
# strainNameB = 'strainB.fasta'
dataDir = 'data/mixed_synth_samples/output'
strainNameA = 'strainA.fasta'
strainNameB = 'strainB.fasta'

output_figure_dir = '/home/user1/Documents/lab/paper/figures/classif.png'
currentDir = os.getcwd()
ref_path = f'{currentDir}/refs/tuberculosis.fna'

data = {'major':[], 'minor':[]}

numSNP = 80

proportions = np.array([0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90])
# proportions = np.array([0.90])

print('number of intoruced SNPs by alt_ref:', numSNP, '. Make sure to check this value.')

for strainName in [strainNameA, strainNameB]:

    altRef_path = f'{currentDir}/refs/{strainName}'
    strain = 'major' if strainName == 'strainA.fasta' else 'minor'

    for id in proportions:
        id = int(id*100)
        consensus_path = f'{currentDir}/{dataDir}/{id}_16/{strain}-consensus.fa'
        numErrors = len(compare_ref(id, ref_path, altRef_path, consensus_path))
        data[strain].append(numErrors)

major_errors = np.array(data['major']) / numSNP *100
minor_errors = np.array(data['minor']) / numSNP *100

print(major_errors)

sns.set(style="darkgrid")
sns.set_context("paper", rc={"lines.linewidth": 1.9, "patch.linewidth" : 0.5, "xtick.major.size": TICK_FONT_SIZE })


fig, ax = plt.subplots()
fig.set_size_inches(7.5, 4.1)
# plt.tight_layout()

ax.set_title('Classification error over major strain percentage', fontsize=TITLE_FONT_SIZE)
ax.set_ylabel('Percentage of incorrectly split SNPs', fontsize=LABEL_FONT_SIZE)
ax.set_xlabel('Major strain proportion', fontsize=LABEL_FONT_SIZE)

ax.set_xticks(proportions)

ax.plot(proportions, major_errors, '--b')
ax.plot(proportions, minor_errors, '--g')
ax.legend(['major', 'minor'])
ax.plot(proportions, major_errors, '.b', linewidth=5, markersize=10)
ax.plot(proportions, minor_errors, '.g', linewidth=5, markersize=10)
y_ticks = ax.get_yticks()
ax.set_yticklabels([f'{int(i)}%' for i in y_ticks])
# x_ticks = ax.get_xticks()
# ax.set_xticklabels([f'{int(i)}%' for i in x_ticks])
plt.savefig(output_figure_dir, dpi=DPI)
plt.show()
