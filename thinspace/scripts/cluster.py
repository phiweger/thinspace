#!/usr/bin/env python


'''
python cluster.py --models ../../models/ --outdir ../../results --force --genomes /Volumes/container/genomic_files/

# --4vis vis.tsv
# better just include in output

TODO: do we need metadata? can we not just dereplicate and done? for plotting
for the publication we can code stuff by hand, like looking up clusters etc.:

> The clustering and dereplication based on `nanotext` is entirely
_taxonomy-free_ -- i.e. we need know nothing about these genomes and can
nevertheless dereplicate them confidently.
'''


import argparse
from collections import defaultdict
from glob import glob
import json
import os
from pathlib import Path

# Filter the multitude of warnings thrown by the HDBSCAN library
# stackoverflow.com/questions/879173
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

import numpy as np
import pandas as pd
# from sklearn.manifold import TSNE
from tqdm import tqdm
import umap
import hdbscan
# from openTSNE import TSNE, TSNEEmbedding, affinity, initialization

from nanotext.classes import GenomeModel
from nanotext.io import dbopen
from nanotext.utils import get_taxa_from_names, get_names_from_taxon, cosine
from nanotext.utils import strip_name


parser = argparse.ArgumentParser(description='Process some integers.')

parser.add_argument(
    '--models', type=str, required=True, 
    help='Path to models (directory)')

# parser.add_argument(
#     '--metadata', type=str, required=True, 
#     help='Path to metadata SQLite database')

parser.add_argument(
    '--config', type=str, default='config.json',
    help='Path to config (file)')

parser.add_argument(
    '--outdir', type=str, required=True,
    help='Path to results (directory)')

parser.add_argument(
    '--genomes', type=str, required=True,
    help='Path to genomes (directory)')

# default=False is implied by action='store_true'
# https://stackoverflow.com/questions/8259001
parser.add_argument(
    '--reduce', action='store_true',
    help='Reduce dimensions before clustering')

parser.add_argument(
    '--force', action='store_true',
    help='Overwrite results if exist')

# TODO: rm? bc/ this is part of some bash script
# parser.add_argument(
#     '--min_identity', dest='ani', default=0.95, type=float,
#     help='ANI above which clustered genomes are dereplicated')

args = parser.parse_args()


in_ = Path(args.genomes).resolve(strict=True)
# TODO: this has to be the absolute path, not relative!
# https://stackoverflow.com/questions/42513056
# p.resolve(strict=True)
# returns a Path object
# Note: strict=True needs Python >= v3.6
# print(in_)

out = Path(args.outdir).resolve()  # does not exist yet, so no strict=True
if args.force:
    out.mkdir(exist_ok=True)
else:
    out.mkdir()


with open(args.config, 'r') as file:
    config = json.load(file)


model = GenomeModel(args.models, mode='core', norm='l2')
genomes = list(in_.glob('[!^\.]*.fasta'))
# Don't include hidden files (e.g. the OS index on a Mac)

print('Collecting embedding for each genome ...')
names, m = [], []
for g in genomes:
    name = g.name.replace('.fasta', '')
    try:
        v = model.embedding[name]
        names.append(name)
        m.append(v)
    except KeyError:
        continue
# names, m = zip(*[(k, v) for k, v in model.embedding.items()])
percent = int(100*len(m)/len(genomes))
print(f'Found {len(m)} embeddings for {len(genomes)} genomes ({percent}%)')

# For testing:
# m = m[:1000]
'''

'''

# on min_cluster_size and min_samples
# (3, 3) -- 3260 clusters, 29548 outliers
# (2, 2) -- 8744 clusters, 38053 outliers
# (3, 2) -- 6537 clusters, 35306 outliers

# (3, 3, 100 dim) -- 3800 clusters, 18k outliers
# (3, 1, 100 dim) -- 4344 clusters, 16109 outliers


'''
Calling clusters --------------------------------------------------------------

On reducing the dimensionality of the embedding before clustering:

> UMAP can be used as an effective preprocessing step to boost the performance
of density based clustering. This is somewhat controversial, and should be
attempted with care. For a good discussion of some of the issues involved in
this please see the various answers in this stackoverflow thread (1) on
clustering the results of t-SNE. 
-- https://umap-learn.readthedocs.io/en/latest/clustering.html

(1) https://stats.stackexchange.com/questions/263539/clustering-on-the-output-of-t-sne

We used "enhanced clustering".
-- https://umap-learn.readthedocs.io/en/latest/clustering.html#umap-enhanced-clustering

>  [B]ut euclidean distance on l2 normalized vectors is equivalent to angular
distance, so you can normalize the vectors and then use euclidean distance.
-- https://github.com/scikit-learn-contrib/hdbscan/issues/69

https://hdbscan.readthedocs.io/en/latest/parameter_selection.html
https://hdbscan.readthedocs.io/en/latest/outlier_detection.html

Param selection:

"umap": {
    "random_state": 42,
    "n_components": 10,
    "metric": "cosine",
    "n_neighbors": 30,
    "min_dist": 0.5
}
Clusters: 5304
Outliers: 15266

- more mindist less clusters/ outliers

"hdbscan": {
    "min_cluster_size": 3,
    "min_samples": 1,
    "metric": "euclidean",
    "cluster_selection_method": "eom"
}
Clusters: 4344
Outliers: 16109
'''

clusterer = hdbscan.HDBSCAN(**config['hdbscan'])

if args.reduce:
    print('Reducing dimensions ...')
    redux = umap.UMAP(**config['umap']).fit_transform(m)
    # There is no random_state for the current HDBSCAN implementation
    print(f'Clustering {len(m)} points ...')
    cluster_labels = clusterer.fit_predict(redux)
else:
    print('No dimension reduction applied')
    print(f'Clustering {len(m)} points ...')
    cluster_labels = clusterer.fit_predict(m)

print(f'Clusters: {len(set(cluster_labels))}')
print(f'Outliers: {sum(cluster_labels == -1)}')


'''
Vis ---------------------------------------------------------------------------
'''

# print(m)

# TODO: optimize
print('Projecting points for visualization ...')
# projection = TSNE(n_components=2, metric='cosine', n_jobs=8).fit(np.array(m))
from sklearn.manifold import TSNE
projection = TSNE(n_components=2).fit_transform(np.array(m))


# TODO: Add vis? Simply uncomment.
l = []
outdirs = defaultdict(list)
for i, label in tqdm(enumerate(cluster_labels)):
    name = names[i]
    
    c1, c2 = projection[i]
    # c1, c2 = 0, 0  # hack until proper vis flag

    # Collect coords for visualization
    l.append([name, label, c1, c2])
    # We will dereplicate the clusters in parallel, so each cluster needs
    # to get a separate directory in the results.
    outdirs[label].append(name)

coords = pd.DataFrame.from_records(l)
coords.columns = 'acc label c1 c2'.split()
coords.to_csv(out/'4vis.csv', header=True, index=False)


'''
Create symlinks

https://kite.com/python/examples/1484/os-create-a-symbolic-link
'''

outc = out/'clusters'
outc.mkdir(exist_ok=True)  # will clobber if exists

print('Adding softlinks ...')
for label, names in outdirs.items():
    # Create a folder for each cluster (outliers are labelled -1).
    # When using this script in nextflow, it wants to feed a dir "-1" to
    # the dereplicator.py script, which is wrongly interpreted as a command
    # line argument -- rename.
    if label == -1:
        label = "outliers"

    fp = outc/str(label)
    fp.mkdir(exist_ok=True)
    # Add symbolic links to input genomes to each cluster directory
    for name in names:
        genome = f'{name}.fasta'
        target = fp/genome
        target.symlink_to(in_/genome)
'''
We will make use of the scripts from Ryan Wick to create the necessary files
to compute a Centrifuge index. For this, we need the GTDB taxonomy and the
genomes meant to go into the index. The genome IDs in the GTDB (taxonomy) are
prefixed w/ "GB_GCA..." (Genbank) and "RS_GCF" (Refseq). However, the index
scripts already prune those names, and so the filenames of the genomes need to
be formatted like "GCA_xxx.fasta".
'''




