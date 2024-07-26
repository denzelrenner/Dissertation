# Dissertation - LIFE4137 CURRENTLY UNDER CONSTRUCTION /\

All scripts required to run the analysis for my final individual project can be found here

Note: Unless on the cloud or stated otherwise, all command line code was ran on a Mac Intel i5

conda activate pipeline_pckgs

## Project Title: An investigation into the association between host genetic elements and antibiotic resistant phenotype

## What are the aims of our analysis?
The aim of this study was to examine the relationship between antibiotic resistance, which is implied through the presence of an antibiotic resistance gene (ARG) family, and other unique genes within host genomes of diverse Gammaproteobacteria, focusing on how the presence of ARG families correlates with the presence or absence of various other genes in the genome. The use of ARG families and not a single ARG to determine our resistant phenotypes ensures we are not labelling resistance based on a single gene. We also investigated if there were COG categories or GO terms that were enriched in genes positively and negatively correlated with the presence of ARG families relative to the entire pan-genome.

## What are the expected outcomes?
We hypothesize that there will be several genes, belonging to different compartments of the cell, that are correlated with the presence of ARG families, and these correlated genes will be enriched for different GO terms and COG categories relative to the entire pangenome that ultimatley affect the resistant phenotype

# Prerequisites

## Tool versions and links
These are all the tools that were used in our analysis with versions and links provided where applicable. Dependencies for certain packages, and their versions, are placed in parentheses. Some references were chosen based on what was recommended on the tool's online help page/documentation.

| Tool | Version | Reference(Harvard Style) |
|------|---------|-----------|
|[ncbi blast](https://blast.ncbi.nlm.nih.gov/)|NA| Johnson, M., Zaretskaya, I., Raytselis, Y., Merezhuk, Y., McGinnis, S. and Madden, T.L., 2008. NCBI BLAST: a better web interface. Nucleic acids research, 36(suppl_2), pp.W5-W9.|
|[TAIR BLAST](https://www.arabidopsis.org/Blast/index.jsp)|version 2.9.0+ (uses NCBI BLAST 2.9.0+)| Garcia-Hernandez, M., Berardini, T., Chen, G., Crist, D., Doyle, A., Huala, E., Knee, E., Lambrecht, M., Miller, N., Mueller, L.A. and Mundodi, S., 2002. TAIR: a resource for integrated Arabidopsis data. Functional & integrative genomics, 2, pp.239-253.|
|[orf finder](https://www.ncbi.nlm.nih.gov/orffinder/)|NA| NA|
|[uniprot](https://www.uniprot.org/align)|clustalO version 1.2.4 (clustal was ran through the align function on the uniprot website)|UniProt Consortium, 2019. UniProt: a worldwide hub of protein knowledge. Nucleic acids research, 47(D1), pp.D506-D515.|
|[alphafold collab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)|(ColabFold v1.5.5: AlphaFold2 using MMseqs2|Mirdita, M., Schütze, K., Moriwaki, Y., Heo, L., Ovchinnikov, S. and Steinegger, M., 2022. ColabFold: making protein folding accessible to all. Nature methods, 19(6), pp.679-682.|
|[InterPro](https://www.ebi.ac.uk/interpro/)|version 98.0|Hunter, S., Apweiler, R., Attwood, T.K., Bairoch, A., Bateman, A., Binns, D., Bork, P., Das, U., Daugherty, L., Duquenne, L. and Finn, R.D., 2009. InterPro: the integrative protein signature database. Nucleic acids research, 37(suppl_1), pp.D211-D215.|
|[ITOL:Interactive Tree Of Life](https://itol.embl.de/)|version 6.8.2|Letunic, I. and Bork, P., 2021. Interactive Tree Of Life (iTOL) v5: an online tool for phylogenetic tree display and annotation. Nucleic acids research, 49(W1), pp.W293-W296.|
|[MEGA](https://www.megasoftware.net/)|version 11|Tamura, K., Dudley, J., Nei, M. and Kumar, S., 2007. MEGA4: molecular evolutionary genetics analysis (MEGA) software version 4.0. Molecular biology and evolution, 24(8), pp.1596-1599.|
|[Jalview](https://www.jalview.org/)|version 2.11.3.2|Waterhouse, A.M., Procter, J.B., Martin, D.M., Clamp, M. and Barton, G.J., 2009. Jalview Version 2—a multiple sequence alignment editor and analysis workbench. Bioinformatics, 25(9), pp.1189-1191.|
|[PyMOL](https://pymol.org/)|version 2.5.8|The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC.|
|[GATK](https://github.com/broadinstitute/gatk) (HTSJDK,Picard)|version 4.2.2.0 (version 2.24.

## Tool intallation 

To install the environments with the different tools necessary to reproduce this analysis you could install the different yaml files we used for the analysis following the command below.

Firstly copy the yaml file to the HPC using rsync or scp and then type in the following command.

```bash
conda env create -f environment.yaml
```

Alternatively, the guidance below outlines the scripts,steps or commands that have to be ran to install some of the tools needed to reproduce our results.


# THE ANALYSIS

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
