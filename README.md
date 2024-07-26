# Dissertation - LIFE4137 CURRENTLY UNDER CONSTRUCTION /\

All scripts required to run the analysis for my final individual project can be found here

Note: Unless on the cloud or stated otherwise, all command line code was ran on a Mac Intel i5

The order of the analysis carried out here is mostly the same as in the report so if there is a specific script or piece of code you may be interested in you can look for the header name for that analysis in the report, and find a similar heading name here.

Dont forget input data

## Project Title: An investigation into the association between host genetic elements and antibiotic resistant phenotype

## Why carry out this analysis?

Antibiotic resistance is a global crisis that continues to worsen as the rate of discovery for novel antimicrobial therapies reduces and antibiotics continue to be misused (misuse of antibiotics continues to rise in healthcare and agriculture), allowing for bacteria to develop resistance mechanisms to multiple different antibiotics. Currently, a single antibiotic resistance gene is considered responsible for the resistance phenotype of bacteria, and there is ongoing interest in using novel tools to discover and predict new antibiotic resistant loci thought to contribute to a resistant phenotype. However, it is more likely that there are many genes that are not ARGs but have a role in confering a resistance phenotype.

## What are the aims of our analysis?
The aim of this study was to examine the relationship between antibiotic resistance, which is implied through the presence of an antibiotic resistance gene (ARG) family, and other unique genes within host genomes of diverse Gammaproteobacteria, focusing on how the presence of ARG families correlates with the presence or absence of various other genes in the genome. The use of ARG families and not a single ARG to determine our resistant phenotypes ensures we are not labelling resistance based on a single gene. We also investigated if there were COG categories or GO terms that were enriched in genes positively and negatively correlated with the presence of ARG families relative to the entire pan-genome.

## What are the expected outcomes?
We hypothesize that there will be several genes, belonging to different compartments of the cell, that are correlated with the presence of ARG families, and these correlated genes will be enriched for different GO terms and COG categories relative to the entire pangenome that ultimatley affect the resistant phenotype

# Prerequisites

## Tool versions and links
These are all the tools that were used in our analysis with versions and links provided where applicable. Dependencies for certain packages, and their versions, are placed in parentheses. Some references were chosen based on what was recommended on the tool's online help page/documentation.

| Tool | Version | Reference(Harvard Style) |
|------|---------|-----------|
|[PanTA](https://github.com/amromics/panta)|NA| Le, D.Q., Nguyen, T.A., Nguyen, S.H., Nguyen, T.T., Nguyen, C.H., Phung, H.T., Ho, T.H., Vo, N.S., Nguyen, T., Nguyen, H.A. and Cao, M.D., 2023. Efficient inference of large pangenomes with PanTA. bioRxiv, pp.2023-07.|
|[RGI](https://github.com/arpcard/rgi)|version 2.9.0+ (uses NCBI BLAST 2.9.0+)| Alcock, B.P., Huynh, W., Chalil, R., Smith, K.W., Raphenya, A.R., Wlodarski, M.A., Edalatmand, A., Petkau, A., Syed, S.A., Tsang, K.K. and Baker, S.J., 2023. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic acids research, 51(D1), pp.D690-D699.|
|[Scoary](https://github.com/AdmiralenOla/Scoary)|NA| Brynildsrud, O., Bohlin, J., Scheffer, L. and Eldholm, V., 2016. Rapid scoring of genes in microbial pan-genome-wide association studies with Scoary. Genome biology, 17, pp.1-9.|
|[Eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)|NA|Cantalapiedra, C.P., Hernández-Plaza, A., Letunic, I., Bork, P. and Huerta-Cepas, J., 2021. eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. Molecular biology and evolution, 38(12), pp.5825-5829., Huerta-Cepas, J., Szklarczyk, D., Heller, D., Hernández-Plaza, A., Forslund, S.K., Cook, H., Mende, D.R., Letunic, I., Rattei, T., Jensen, L.J. and von Mering, C., 2019. eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Nucleic acids research, 47(D1), pp.D309-D314.|
|[GOATOOLS](https://github.com/tanghaibao/goatools)|(ColabFold v1.5.5: AlphaFold2 using MMseqs2|Klopfenstein, D.V., Zhang, L., Pedersen, B.S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C.J., Yunes, J.M., Botvinnik, O., Weigel, M. and Dampier, W., 2018. GOATOOLS: A Python library for Gene Ontology analyses. Scientific reports, 8(1), pp.1-17.|
|[BUSCO](https://gitlab.com/ezlab/busco)|version 98.0|Simão, F.A., Waterhouse, R.M., Ioannidis, P., Kriventseva, E.V. and Zdobnov, E.M., 2015. BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics, 31(19), pp.3210-3212., Manni, M., Berkeley, M.R., Seppey, M., Simão, F.A. and Zdobnov, E.M., 2021. BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. Molecular biology and evolution, 38(10), pp.4647-4654.|
|[PyANI](https://github.com/widdowquinn/pyani)|version 6.8.2|Pritchard, L., Glover, R.H., Humphris, S., Elphinstone, J.G. and Toth, I.K., 2016. Genomics and taxonomy in diagnostics for food security: soft-rotting enterobacterial plant pathogens. Analytical methods, 8(1), pp.12-24.|
|[NCBI-Datasets](https://github.com/ncbi/datasets)|version 11|O’Leary, N.A., Cox, E., Holmes, J.B., Anderson, W.R., Falk, R., Hem, V., Tsuchiya, M.T., Schuler, G.D., Zhang, X., Torcivia, J. and Ketter, A., 2024. Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. Scientific Data, 11(1), p.732.|
|[Prokka](https://github.com/tseemann/prokka)|version 11|Seemann, T., 2014. Prokka: rapid prokaryotic genome annotation. Bioinformatics, 30(14), pp.2068-2069.|

## Tool intallation 

To install the environments with the different tools necessary to reproduce this analysis you could install the different yaml files we used for the analysis following the command below.

Firstly copy the yaml file to the HPC using rsync or scp and then type in the following command.

```bash
conda env create -f environment.yaml
```

Alternatively, the guidance below outlines the scripts,steps or commands that have to be ran to install some of the tools needed to reproduce our results.

Add eggnog stuff to your bash profile


# THE ANALYSIS

You will notice a log_files directory, an all_gammaproteobacteria directory and a scripts directory.

## Part1 - Data Acquisition
Before conducting any sort of analysis we need to get the nucleotide and protein fastas for different Gammaproteobacteria, as well as their associated metadata. Running the `data_acquisition_pipeline.sh` script will run the whole data acquisition pipeline for the data acquisition stage of the analysis. It first gets all the different complete and refseq (GCF) assemblies for different gammaproteobcteria from ncbi. This includes The accession code , organism name, assembly release data, completeness, and a whole battery of other metrics are for each species is written to an output file. Then the resulting file is filtered to remove duplicate assemblies, and between different strains of the same species, take the most recent assembly. All the different assemblies for our final list of gammaproteobacteria are then downloaded as zip files, and we put metrics like assembly length into a txt file. The R script `distribution_plotting.r` is also ran to produce a plot showing the distribution of the GC content and assembly lengths which would become useful when doing phylogenetic trees because sometimes considering GC content is important.

1. Activate conda environment with the tools the script uses
```bash
conda activate pipeline_pckgs
```

2. Run the script
```bash
sbatch ~/scripts/data_acquisition_pipeline.sh
```

3. Run another script to unzip the zip files and put all nucleotide and protein fasta files into separate directories 
```bash
python3 ~/scripts/get_nucleotide_and_protein_fasta.py \
-i ~/all_gammaproteobacteria_data/assemblies -d1 ~/all_gammaproteobacteria_data/assemblies/nucleotide_fasta_files -d2 ~/all_gammaproteobacteria_data/assemblies/protein_fasta_files
```

This should produce three different directories. One directory is called `assemblies` and contains zip files for each genome. Another directory is called `output_data` and this contains the list of gammaproteobacteria we accesed from NCBI. The final directory is called `metric_files` and this contains files with different metrics about every genome in the data liek N50 and assembly length.

## Part2 - Quality Filtering

Now that we have successfully downloaded our genomes, we want to filter out genomes that have bad quality scores, and also remove any two genomes that appear genetically identical. To do this we first deduplicated the dataset using average nucelotide identity (ANI) and then BUSCO.

### ANI Filtering

### BUSCO Filtering

All the bacteria in our study are gammaproteobateria and so we need to use that specific lineage dataset when running busco. To download the lineage dataset follow the steps below

Activate the conda environemnt with busco 
```bash
conda activate pipeline_pckgs
```

1. Activate the right conda env
```bash
conda activate pipeline_pckgs
```

2.Make a directory to store the lineage dataset 
```bash
mkdir -p ~/all_gammaproteobacteria_data/busco_lineage_dataset/
```

3.Download the gammaproteobacteria lineage dataset with busco and direct the output to the datasets directory created earlier
```bash
busco --download_path ~/all_gammaproteobacteria_data/busco_lineage_dataset/ --download gammaproteobacteria_odb10
```

4.Run busco and filter the genomes with a busco score of 90%. This script contains two other scripts which are important for plotting
```bash
sbatch run_busco.sh
```

The important files form this is a directory called `busco_plotting_output_and_summary_files` which contains a pdf file showing the different busco scors for completeness, duplicated etc.


## Part3 - Finding ARG families

Now that we have  our final list of genomes to use in the analysis, we need to find the ARG families across all the 1137 genomes and create a list of all the unique ARG families found.

Before we actually run rgi, we need to put the protein fasta files into their own unique directory. That is accomplished following the steps below

1. Activate the conda env
```bash
conda activate pipeline_pckgs
```

2. We initally stored the protein fasta files for all 1348 files in the [Data Acquisition](#part1---data-acquisition) portion of this analysis. We now want to create a new directory which only has protein fasta files for the final dataset of 1173 genomes. These files will act as input for RGI to then identify ARG families. To do this enter the command below into the command line.
```bash
python3 filtering_data_directories.py -id ~/all_gammaproteobacteria_data/assemblies/protein_fasta_files --input_file ~/all_gammaproteobacteria_data/busco_plotting_output_and_summary_files/busco_filtered_gammaproteobacteria_correction_names.txt -od ~/all_gammaproteobacteria_data/rgi_input_protein_fasta_1173_genomes
```

3. Now we can actually run RGI to get the different ARG families. The script we will use is called `rgi_main_pipeline_1173_genomes.sh` and it uses some other python scripts to create a final list of ARG families. To accomplish this, follow the commands below.
```bash
sbatch ~/scripts/rgi_scripts/rgi_main_pipeline_1173_genomes.sh
```

This command will produce a directory called `~/all_gammaproteobacteria_data/rgi_output_1173_genomes/gene_family_files`, and within this directory there will be a file called `unique_gene_families.txt` and it contains a list of all ARG families used in this study

## Part4 - Building Pan-genome

## Part5 - Running Scoary
Now that we have the unique ARG families in our dataset, we essentially have the traits that will be used by scoary. For the purpose of our analysis we are interested in correlation and coincidence, rather than causative stuff. To do this you can simply enter the two command into the command line.

```bash
conda activate pipeline_pckgs

sbatch ~/scripts/scoary_scripts/main_scoary_pipeline_groups_1173genomes_NOPAIRWISE_with_groups
```

This script does a lot of different things and uses a lot of intermediate scripts within it. Each script is commented but here is what they do

## Part5A - Correlation

## Part5B - Cytoscape

## Part6 - Annotating the Pan-genome

We have successfully found genes that were positively or negatively correlated with the presence of the ARG families in this study. The next piece of analysis involves looking for enrichemnt of GO terms, and over-representation of COG categories in our positively and negatively correlated genes compared to the entire pan-genome. To do this, we first need to annotate the protein sequences for each gene in the pangenome with COG categories and GO terms. Then we can use GOATOOLS, and the SciPY python package to look for enrichment and over-representation.

# Part 6A - EggNOG-mapper
To annotate the pan-genome with GO terms and COG categories run the command below.

```bash
sbatch ~/scripts/eggnog_scripts/emapper_pangenome_input.sh
```
This will produce a file called something.annotations, and it contains a tsv of genes and their GO terms, and COG annotations.

# Part6B - Goatools

To accomplish this you can simply run the command below.

```bash
sbatch ~/scripts/goatools_scripts/run_goatools_1173genomes_UPDATED.sh
```

Again, to breakdown the script a bit and the output you get from running it

# Part6C - COG Categories
For the COG categories we want to . This is accomplished by running the commands below

```bash
conda activate pipeline_pckgs

python3 ~/scripts/cog_scripts/cog_calculations.py -od ~/all_gammaproteobacteria_data/COG_calculation_1173_genomes --annotation_file ~/all_gammaproteobacteria_data/eggnog_test_WITH_GROUPS/whole_dataset.emapper.annotations --positively_correlated_genes ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE_WITH_GROUPS/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families/master_positive_correlated_genes.tsv --negatively_correlated_genes ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE_WITH_GROUPS/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families/master_negative_correlated_genes.tsv --Rtab_file ~/all_gammaproteobacteria_data/panta_output/panta_030pid_e7_LD07_split_1173_genomes_with_sflag_sotruenosplit/gene_presence_absence.Rtab
```

Again, this script does a lot of different things. It takes the annotations file produced from EggNOG-mapper in [the EggNOF-mapper step of the analysis](#part6A---EggNOG-mapper)


# Conclusion
