# Dissertation - CURRENTLY UNDER CONSTRUCTION /\
All scripts required to run the analysis for my final individual project

conda activate pipeline_pckgs

To get all nucleotide and protein fasta files run this commmand
python3 ~/scripts/get_nucleotide_and_protein_fasta.py -i ~/all_gammaproteobacteria_data/assemblies -d1 ~/all_gammaproteobacteria_data/assemblies/nucleotide_fasta_files -d2 ~/all_gammaproteobacteria_data/assemblies/protein_fasta_files

conda deactivate 



All the bacteria in our study are gammaproteobateria and so we need to use that specific lineage dataset when running busco. To download the lineage dataset follow the steps below

Activate the conda environemnt with busco 
conda activate pipeline_pckgs

Make a directory to store the lineage dataset 
```bash
mkdir -p ~/all_gammaproteobacteria_data/busco_lineage_dataset/
```

Download the gammaproteobacteria lineage dataset with busco and direct the output to the datasets directory created earlier
```bash
busco --download_path ~/all_gammaproteobacteria_data/busco_lineage_dataset/ --download gammaproteobacteria_odb10
```
