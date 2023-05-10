#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

#EIG_PCA_FILE = '/net/topmed10/working/porchard/rnaseq/work/genotype-pca/unrelated/results/pca/filtered.pca.evec'
#EIG_EIGENVALUES_FILE = '/net/topmed10/working/porchard/rnaseq/work/genotype-pca/unrelated/results/pca/filtered.eval'
#ANCESTRIES = '/net/topmed10/working/porchard/rnaseq/data/ancestry/NWD_ancestry.txt'
EIG_PCA_FILE, EIG_EIGENVALUES_FILE, ANCESTRIES = sys.argv[1:]


def plot_pca(pc_scores, explained_variance=None, labels={}, kwargs={}, plot_top=None, panel_width=4, panel_height=4):
    """
    pc_scores: pd.DataFrame (columns = PC1, PC2, etc)
    explained_variance: pandas series or dict, PC1 --> fraction variance explained, etc
    kwargs will be passed through to sns.scatterplot
    """
    plot_top = plot_top if plot_top is not None else len(pc_scores.columns)
    number_plots = plot_top / 2
    nrows = int(round(np.sqrt(number_plots)))
    ncols = int(np.ceil(number_plots/nrows))
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(panel_width*ncols, panel_height*nrows), squeeze=False)
    all_texts = []
    for count, ax in enumerate(axs.flatten(), 1):
        if (2*(count-1)+1) > plot_top:
            ax.remove()
            continue
        first_pc = 'PC{}'.format(2*(count-1)+1)
        second_pc = 'PC{}'.format(2*(count-1)+2)
        xlab = '{} ({}%)'.format(first_pc, round(explained_variance[first_pc]*100, 2)) if explained_variance is not None else first_pc
        ylab = '{} ({}%)'.format(second_pc, round(explained_variance[second_pc]*100, 2)) if explained_variance is not None else second_pc
        ax = sns.scatterplot(x=first_pc, y=second_pc, data=pc_scores, ax=ax, **kwargs)
        if len(labels) > 0:
            labelme = pc_scores[pc_scores.index.isin(labels.keys())]
            texts = [ax.text(x=r[first_pc], y=r[second_pc], s=labels[i]) for i, r in labelme.iterrows()]
            all_texts.append(texts)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
    for texts, ax in zip(all_texts, axs.flatten()):
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='red'))
    #handles, labels = ax.get_legend_handles_labels()
    #fig.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left')
    #for ax in axs.flatten():
    #    ax.legend().remove()
    fig.tight_layout()
    return (fig, axs)


def load_eigensoft_pca(evec_file, eval_file):
    pca = None
    with open(evec_file, 'r') as fh:
        pca = [i.strip() for i in fh.readlines()]
    eigenvalues = [float(i) for i in pca[0].split()[1:]]
    pc_scores = pd.DataFrame([i.split() for i in pca[1:]])
    pc_scores.columns = ['id'] + [f'PC{i}' for i in range(1, len(pc_scores.columns) - 1)] + ['pop']
    pc_scores = pc_scores.drop(columns='pop').set_index('id').astype(float).reset_index()

    eigenvalues = pd.read_csv(eval_file, header=None, names=['eigenvalue'])
    eigenvalues['PC'] = [f'PC{i}' for i in range(1, len(eigenvalues)+1)]
    eigenvalues['fraction_variance_explained'] = eigenvalues.eigenvalue / eigenvalues.eigenvalue.sum()
    eigenvalues = eigenvalues[['PC', 'fraction_variance_explained']]

    return (pc_scores, eigenvalues)



pc_scores, variance_explained = load_eigensoft_pca(EIG_PCA_FILE, EIG_EIGENVALUES_FILE)


ancestry = pd.read_csv(ANCESTRIES, sep='\t').set_index('NWD_ID').rename_axis(index='id')
all_ancestries = ancestry.columns.to_list()
ancestry['ancestry'] = ancestry.idxmax(axis=1)


pc_scores.id = pc_scores.id.str.split(':', expand=True)[1]
pc_scores.to_csv('PC-scores.txt', sep='\t', index=False)
variance_explained.to_csv('PC-variance-explained.txt', sep='\t', index=False)

# plot by ancestries
pc_scores = pc_scores.merge(ancestry.reset_index(), how='left')
pc_scores.ancestry= pc_scores.ancestry.fillna('???')


fig, axs = plot_pca(pc_scores, variance_explained.set_index('PC').fraction_variance_explained, kwargs={'alpha': 0.1, 'hue': 'ancestry'}, plot_top=14)
fig.savefig('genotype-PCA.by-ancestry.png', dpi=300)
fig.clf()

for ancestry in all_ancestries:
    fig, axs = plot_pca(pc_scores, variance_explained.set_index('PC').fraction_variance_explained, kwargs={'alpha': 0.1, 'hue': ancestry}, plot_top=14)
    fig.savefig(f'genotype-PCA.by-{ancestry}.png', dpi=300)
    fig.clf()