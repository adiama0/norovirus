from Bio import SeqIO
import sys
import pandas as pd
from tqdm import tqdm

nextstrain_fasta = sys.argv[1]
nextstrain_metadata = pd.read_csv(sys.argv[2], sep="\t")
output_file = sys.argv[3]

downsampled_records = []

for record in tqdm(SeqIO.parse(nextstrain_fasta, "fasta")):
    if record.id in nextstrain_metadata["name"].values:
        downsampled_records.append(record)

SeqIO.write(downsampled_records, output_file, "fasta")
