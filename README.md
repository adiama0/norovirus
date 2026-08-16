To create the mock samples, bygul was used with the following command inside the sequences/ folder: 
```
bygul simulate-proportions \
  source_fastas/GII2_fasta,source_fastas/GII17_fasta \
  --proportions 0.50,0.50 \
  --outdir sample_43 \
  --simulation_mode metagenomics \
  --readcnt 100000 \
  --error_rate 0.001
```
This is an example for the generation of sample_44. Individual sample proportions can be found in 'data/sample_proportions.xlsx'

run the following command to individuate medoid sequences:
```
run_medoid_reference.sh
```

run the following command to generate sample proportions:
```
bash run_completitive_mapping.sh
```