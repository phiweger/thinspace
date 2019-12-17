## thinspace

Many recent microbial genome collections curate hundreds of thousands of genomes. This volume complicates many genomic analyses such as taxon assignment because the associated computational burden is substantial. However, the number of representatives of each species is highly skewed towards human pathogens and model organisms. Thus many genomes contain little additional information and could be removed. We created a frugal dereplication method that can reduce massive genome collections based on genome sequence alone, without the need for manual curation nor taxonomic information.


### Installation

There are a couple of dependencies, we apologize -- on the other hand, `thinspace` will embed genomes, cluster them, assess cluster quality, dereplicate them using k-mer hashes and finally create a `Centrifuge` index for direct use.


```bash
# Create a new conda environment to experiment in
conda create -y -n nanotext && conda activate nanotext

# Then we need a couple of useful programs that are used internally by
# nanotext and thinspace, as well as the Open Science Foundation (OSF) 
# client to fetch test data
conda install -y -c pytorch faiss-cpu
conda install -y -c bioconda biopython pybedtools pysam mash
conda install -y -c conda-forge umap-learn hdbscan osfclient 

# Install the nanotext library
git clone https://github.com/phiweger/nanotext
cd nanotext
pip install -e .
pytest  # run tests
cd ..
```


### Run

We are now ready dereplicate some genomes. Put them in a directory `genomes`. The filename needs to match the first column in the `taxonomy.csv` metadata. We provide some _Pseudomonas_ sample genomes in this repo for a quick trial (source: Lopes _et al._, 2018. “Genome Variations between Rhizosphere and Bulk Soil Ecotypes of a Pseudomonas Koreensis Population.” Environmental Microbiology 20 (12): 4401–14.).


```bash
# Get thinspace
git clone https://github.com/phiweger/thinspace && cd thinspace/thinspace

# Get nanotext model(s)
osf -p pjf7m fetch models.zip && unzip models.zip
osf -p pjf7m fetch genomes.zip && unzip genomes.zip

# Get taxonomy metadata from GTDB
wget -O taxonomy.tsv https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/ar122_taxonomy_r89.tsv

# And go
./thinspace.sh models genomes taxonomy.tsv tmp results
```


###  Run via Docker

#### Using example data

```bash
# Get taxonomy metadata from GTDB
wget -O taxonomy.tsv https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/ar122_taxonomy_r89.tsv

# Start docker
docker run --rm -it -v $PWD:/input nanozoo/thinspace:latest /bin/bash

# Move into the git dir
cd /thinspace/thinspace

# Use the preloaded models and genomes
thinspace.sh /thinspace/thinspace/models /thinspace/thinspace/genomes /input/taxonomy.tsv tmp /input/results
```

* **done**, close docker via `ctrl + D`
    * results are located in `results/`

#### Use your own data

```bash
# prepare models and genomes dirs (these are example data sets)
osf -p pjf7m fetch models.zip && unzip models.zip && rm models.zip
osf -p pjf7m fetch genomes.zip && unzip genomes.zip && rm genomes.zip

# Get taxonomy metadata from GTDB
wget -O taxonomy.tsv https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/ar122_taxonomy_r89.tsv
```

* Mount your current dir containing `models/`, `genomes/` and `taxonomy.tsv` into the docker thinspace container via:

```bash
docker run --rm -it -v $PWD:/input nanozoo/thinspace:latest /bin/bash
```

* inside the docker container do the following for your own data:

```bash
# move into the git dir
cd /thinspace/thinspace

# execute thinspace
thinspace.sh /input/models /input/genomes /input/taxonomy.tsv tmp /input/results
```

* **done**, close docker via `ctrl + D`
    * results are located in `results/`



