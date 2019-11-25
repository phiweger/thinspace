'''
conda activate thinspace
'''
import shutil
from pathlib import Path
# https://stackoverflow.com/questions/33625931

from tqdm import tqdm

from nanotext.utils import strip_name


# GTDB
genomes = Path('/Volumes/container/genomic_files').glob('*.tsv')
out = Path('/Volumes/container/genomes_gtdb_r89')
out.mkdir(exist_ok=True)

for i in tqdm(genomes):
    target = out/strip_name(i.name.replace('_genomic.tsv', '.fasta'))
    shutil.copy(i, target)


# GTDB redux
genomes = Path('/Volumes/container/genomes_gtdb_r89_redux_rename').glob('[!^\.]*.fna')
out = Path('/Volumes/container/genomes_gtdb_r89_redux_forindex')
out.mkdir(exist_ok=True)

for i in tqdm(genomes):
    target = out/strip_name(i.name.replace('_genomic.fna', '.fasta'))
    shutil.copy(i, target)