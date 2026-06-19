import sys
from Bio import SeqIO

input_gb = sys.argv[1]
output_fasta = sys.argv[2]

record = SeqIO.read(input_gb, "genbank")
record.id = record.id.split(".")[0]
record.name = record.id
record.description = ""

SeqIO.write(record, output_fasta, "fasta")
