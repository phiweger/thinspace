#!/bin/bash


# time ./thinspace.sh <models> <genomes> <taxonomy> <tmpdir> <outdir>
# cd /Users/phi/tmp/thinspace_testrun
# time ./thinspace.sh models/ /Volumes/container/genomes_gtdb_r89/ taxonomy.tsv tmp results

# The results have to be written to the local disk, otherwise dereplicator.py
# complains:
# ERROR: could not open "/var/folders/[...]/tmp5bw0tggl/mash.msh" for reading.

printf "\nStarting to thin genome space ..."
mkdir -p ${4}/clustering
mkdir ${4}/clustering2
mkdir $5
touch ${5}/log
THRESHOLD=0.05
# 0.1 is aggressive, 0.05 is standard procedure


# 1. Cluster the nanotext space using HDBSCAN
# There is a --reduce option if dimension reduction is required.
python scripts/cluster.py \
    --config config.json \
    --models $1 \
    --outdir ${4}/clustering \
    --force \
    --genomes $2 \
    2>&1 | tee -a ${5}/log
cp ${4}/clustering/4vis.csv $5


echo "Cluster evaluation:"
python scripts/annotate.py \
    --taxonomy $3 \
    --clusters ${4}/clustering/4vis.csv \
    --rank species \
    2>&1 | tee -a ${5}/log


# 2. Dereplicate the individual clusters at <threshold>
CLUSTERS=${4}/clustering/clusters
DEREPLICA=${4}/dereplica
mkdir $DEREPLICA

for i in $(ls -1 $CLUSTERS); do
    if [[ $i -ne "outliers" ]]
    then
        echo "Cluster ${i} ..."
        scripts/dereplicator/dereplicator.py \
            --threads 8 \
            --threshold $THRESHOLD \
            ${CLUSTERS}/${i} ${DEREPLICA}/ \
            2>&1 | tee -a ${5}/log
    fi
done

# cp ${CLUSTERS}/outliers/* ${DEREPLICA}
# find /Volumes/box/tmp/results -name "*.fasta" | wc -l
find ${CLUSTERS}/outliers -name *.fasta -exec cp {} ${DEREPLICA} \;


# Second round
printf "\nSecond round of dereplication ..."
python scripts/cluster.py \
    --config config.json \
    --models $1 \
    --outdir ${4}/clustering2 \
    --force \
    --genomes ${DEREPLICA} \
    2>&1 | tee -a ${5}/log
cp ${4}/clustering2/4vis.csv $5/4vis2.csv
# Found 28456 embeddings for 28456 genomes (100%)
# No dimension reduction applied
# Clustering 145817 points ...
# Clusters: 1467
# Outliers: 8614
# TODO: change cluster.py so we don't have to write $PWD/${DEREPLICA}


# 3. Assess clusters
echo "Cluster evaluation:"
python scripts/annotate.py \
    --taxonomy $3 \
    --clusters ${4}/clustering2/4vis.csv \
    --rank species \
    2>&1 | tee -a ${5}/log
# species
# round 1
# 0.6511 Adjusted Rand score
# 0.7595 Adjusted mutual information score
# round 2
# 0.0113 Adjusted Rand score
# 0.1558 Adjusted mutual information score
# genus
# round 1
# 0.7139 Adjusted Rand score
# 0.7964 Adjusted mutual information score

CLUSTERS2=${4}/clustering2/clusters
DEREPLICA2=${4}/dereplica2
mkdir $DEREPLICA2

for i in $(ls -1 $CLUSTERS2); do
    if [[ $i -ne "outliers" ]]
    then
        echo "Cluster ${i} ..."
        scripts/dereplicator/dereplicator.py \
            --threads 8 \
            --threshold $THRESHOLD \
            ${CLUSTERS2}/${i} ${DEREPLICA2}/ \
            2>&1 | tee -a ${5}/log
    fi
done
# 13184 genomes after dereplication of clusters + 8614 outliers = 21k genomes
# Before derep round 2, there were 28456 genomes, so about 25% reduction still.


# TODO: The following is not possible bc/ dereplica2 outliers are soft linked
# to dereplica1 -- this copying is silly, find some solution.
# rm -r $DEREPLICA
# frees like 100 GB
find ${CLUSTERS2}/outliers -name *.fasta -exec cp {} ${DEREPLICA2} \;


# 4. Process files so we can build an index for Centrifuge
# scripts/Metagenomics-Index-Correction/tax_from_gtdb.py \
#     --gtdb $3 \
#     --assemblies $DEREPLICA2 \
#     --nodes ${5}/ex.tree \
#     --names ${5}/ex.name \
#     --conversion ${5}/ex.conv \
#     --cat_fasta ${5}/ex.fa \
#     2>&1 | tee -a ${5}/log


# 5. Remove tmpdir
# mkdir $5/dereplica
mv ${DEREPLICA2} ${5}/dereplica
rm -r $4


# 6. Pack
# tar cf - results | pigz -p 8 > 4centrifuge.tar.gz
echo "Done."