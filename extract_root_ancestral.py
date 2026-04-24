import sys
from Bio import Phylo, SeqIO
 
tree_file = sys.argv[1]
anc_fasta = sys.argv[2]
out_fasta = sys.argv[3]
 
tree = Phylo.read(tree_file, "newick")
root_name = tree.root.name
 
records = SeqIO.to_dict(SeqIO.parse(anc_fasta, "fasta"))
root_record = records[root_name]
 
root_record.id = "norovirus_all_VP1_ancestral"
root_record.name = root_record.id
root_record.description = ""
 
SeqIO.write(root_record, out_fasta, "fasta")
print(f"Root node: {root_name}")
print(f"Wrote ancestral reference to: {out_fasta}")