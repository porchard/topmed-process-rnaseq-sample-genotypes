#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# don't use kinship -- use InfType
# since kinship estimate isn't as accurate overall...

KIN0, FAM = sys.argv[1:]

#kin0 = '/net/topmed10/working/porchard/rnaseq/work/unrelated-samples/freeze-alpha-truncated/work/b0/69de6221190c5d5ba68f2470236599/test.kin0'
#fam = '/net/topmed10/working/porchard/rnaseq/work/unrelated-samples/freeze-alpha-truncated/work/b0/69de6221190c5d5ba68f2470236599/plink.merged.fam'
kin0 = pd.read_csv(KIN0, sep='\t')

fig, ax = plt.subplots()
sns.histplot(x='Kinship', hue='InfType', data=kin0, ax=ax)
fig.tight_layout()
fig.savefig('kinship.png', dpi=300)
fig.clf()

kin0 = kin0[kin0.InfType != 'UN']

# for each individual:
# get number of non-relatives...
# rank, and shuffle within each rank
all_individuals = pd.read_csv(FAM, sep=' ', header=None, names=['fid', 'iid', 'x1', 'x2', 'x3', 'x4'])
all_individuals = all_individuals.iid.to_list()
relatives = {i: set() for i in all_individuals}

for id1, id2 in zip(kin0.ID1, kin0.ID2):
    relatives[id1].add(id2)
    relatives[id2].add(id1)

relative_counts = pd.DataFrame([[i, (len(all_individuals) - 1) - len(relatives[i])] for i in relatives.keys()], columns=['iid', 'nonrelated_count'])
relative_counts = relative_counts.sample(frac=1.0, random_state=2022).reset_index(drop=True)
relative_counts = relative_counts.sort_values('nonrelated_count', ascending=False)


unrelated = set()
drop = set()
for i in relative_counts.iid.to_list():
    if i not in drop:
        unrelated.add(i)
        # drop all relatives
        for r in relatives[i]:
            drop.add(r)
    else:
        continue


kin0_keep = kin0[(kin0.ID1.isin(unrelated)) & (kin0.ID2.isin(unrelated))]
assert(len(kin0_keep) == 0)

print('\n'.join(list(unrelated)))