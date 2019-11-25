'''
python workflow/scripts/annotate.py --taxonomy db/taxonomy.tsv --clusters foo/4vis.csv --rank species

Metrics first seen here:
 https://umap-learn.readthedocs.io/en/latest/clustering.html


# TODO?
# print(f'Loading metadata ...')
# db = args.metadata
# df = get_taxa_from_names(db, names)
# names_w_metadata = list(df['name'])

# selection = list(df[df.order == 'Clostridiales']['name'].unique())
# selection = list(df[df.genus == 'Campylobacter']['name'].unique())
# selection = list(df['name'].unique())
'''
import argparse

import pandas as pd
from nanotext.utils import strip_name

from sklearn.metrics import adjusted_rand_score as ARS
from sklearn.metrics import adjusted_mutual_info_score as AMIS
# https://umap-learn.readthedocs.io/en/latest/clustering.html


parser = argparse.ArgumentParser(description='How good are these taxonomic clusters, really?')

parser.add_argument('--rank', type=str, default='genus',
    help='At which taxonomic rank to assess cluster homogeneity?')

parser.add_argument('--taxonomy', type=str, required=True,
    help='GTDB-style formatted taxonomy: GCF_000979745.1\td__Archaea;p__Halobacterota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina mazei\n')

parser.add_argument('--clusters', type=str, required=True,
    help='Cluster labels and coordinates: GCF_000979745.1,262,-31.92,9.89\n')

args = parser.parse_args()


ranks = {
    'domain': 'd',
    'phylum': 'p',
    'class': 'c',
    'order': 'o',
    'family': 'f',
    'genus': 'g',
    'species': 's',
}


qr = ranks[args.rank]


labels = {}

with open(args.taxonomy, 'r') as file:
    for line in file:
        name, taxa = line.strip().split('\t')
        d = dict([rank.split('__') for rank in taxa.split(';')])
        
        labels[strip_name(name)] = d[qr]
        

x, y = [], []
with open(args.clusters, 'r') as file:
    _ = next(file)  # skip header
    for line in file:
        name, cluster, c1, c2 = line.strip().split(',')
        x.append(cluster)
        y.append(labels[name])


print('{0:.4f}'.format(ARS(x, y)), 'Adjusted Rand score')
print('{0:.4f}'.format(AMIS(x, y, average_method='arithmetic')),
    'Adjusted mutual information score')
# FutureWarning: The behavior of AMI will change in version 0.22. To match the behavior of 'v_measure_score', AMI will use average_method='arithmetic' by default.



