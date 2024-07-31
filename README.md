# Dissertation - LIFE4137 CURRENTLY UNDER CONSTRUCTION /\

All scripts required to run the analysis for my final individual project can be found here

Note: Unless on the cloud or stated otherwise, all command line code was ran on a Mac Intel i5. In all sbatch scripts you would need to change the #SBATCH --output and
#SBATCH --error lines at the beginning of the script such that output and errors are written to your own personal directory.

The order of the analysis carried out here is mostly the same as in the report so if there is a specific script or piece of code you may be interested in you can look for the header name for that analysis in the report, and find a similar heading name here.

Dont forget input data

## Project Title: An investigation into the association between host genetic elements and antibiotic resistant phenotype

## Why carry out this analysis?

Antibiotic resistance is a global crisis that continues to worsen as the rate of discovery for novel antimicrobial therapies reduces and antibiotics continue to be misused (misuse of antibiotics continues to rise in healthcare and agriculture), allowing for bacteria to develop resistance mechanisms to multiple different antibiotics. Currently, a single antibiotic resistance gene is considered responsible for the resistance phenotype of bacteria, and there is ongoing interest in using novel tools to discover and predict new antibiotic resistant loci thought to contribute to a resistant phenotype. However, it is more likely that there are many genes that are not ARGs but have a role in confering a resistance phenotype.

## What are the aims of our analysis?
The aim of this study was to examine the relationship between antibiotic resistance, which is implied through the presence of an antibiotic resistance gene (ARG) family, and other unique genes within host genomes of diverse Gammaproteobacteria, focusing on how the presence of ARG families correlates with the presence or absence of various other genes in the genome. The use of ARG families and not a single ARG to determine our resistant phenotypes ensures we are not labelling resistance based on a single gene. We also investigated if there were COG categories or GO terms that were enriched in genes positively and negatively correlated with the presence of ARG families relative to the entire pan-genome.

## What are the expected outcomes?
We hypothesize that there will be several genes, belonging to different compartments of the cell, that are correlated with the presence of ARG families, and these correlated genes will be enriched for different GO terms and COG categories relative to the entire pangenome that ultimatley affect the resistant phenotype

# Table of Contents

[Part1 - Data Acquisition](#part1---data-acquisition)

[Part2 - Quality Filtering](#part2---quality-filtering)

[Part3 - Finding ARG families](#part3---finding-arg-families)


# Prerequisites

## Tool versions and links
These are all the tools that were used in our analysis with versions and links provided where applicable. Dependencies for certain packages, and their versions, are placed in parentheses. Some references were chosen based on what was recommended on the tool's online help page/documentation.


| Tool | Version | Reference(Harvard Style) |
|------|---------|-----------|
|[PanTA](https://github.com/amromics/panta)|Version 1.0.0| Le, D.Q., Nguyen, T.A., Nguyen, S.H., Nguyen, T.T., Nguyen, C.H., Phung, H.T., Ho, T.H., Vo, N.S., Nguyen, T., Nguyen, H.A. and Cao, M.D., 2023. Efficient inference of large pangenomes with PanTA. bioRxiv, pp.2023-07.|
|[RGI](https://github.com/arpcard/rgi)|Version 6.0.3| Alcock, B.P., Huynh, W., Chalil, R., Smith, K.W., Raphenya, A.R., Wlodarski, M.A., Edalatmand, A., Petkau, A., Syed, S.A., Tsang, K.K. and Baker, S.J., 2023. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic acids research, 51(D1), pp.D690-D699.|
|[Scoary](https://github.com/AdmiralenOla/Scoary)|Version 1.6.16| Brynildsrud, O., Bohlin, J., Scheffer, L. and Eldholm, V., 2016. Rapid scoring of genes in microbial pan-genome-wide association studies with Scoary. Genome biology, 17, pp.1-9.|
|[Eggnog-mapper (DIAMOND)](https://github.com/eggnogdb/eggnog-mapper)|Version 5.0.2 (Version 2.0.11) |Cantalapiedra, C.P., Hernández-Plaza, A., Letunic, I., Bork, P. and Huerta-Cepas, J., 2021. eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. Molecular biology and evolution, 38(12), pp.5825-5829., Huerta-Cepas, J., Szklarczyk, D., Heller, D., Hernández-Plaza, A., Forslund, S.K., Cook, H., Mende, D.R., Letunic, I., Rattei, T., Jensen, L.J. and von Mering, C., 2019. eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Nucleic acids research, 47(D1), pp.D309-D314.|
|[GOATOOLS](https://github.com/tanghaibao/goatools)|Version 1.4.12|Klopfenstein, D.V., Zhang, L., Pedersen, B.S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C.J., Yunes, J.M., Botvinnik, O., Weigel, M. and Dampier, W., 2018. GOATOOLS: A Python library for Gene Ontology analyses. Scientific reports, 8(1), pp.1-17.|
|[BUSCO](https://gitlab.com/ezlab/busco)|Version 5.4.3|Simão, F.A., Waterhouse, R.M., Ioannidis, P., Kriventseva, E.V. and Zdobnov, E.M., 2015. BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics, 31(19), pp.3210-3212., Manni, M., Berkeley, M.R., Seppey, M., Simão, F.A. and Zdobnov, E.M., 2021. BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. Molecular biology and evolution, 38(10), pp.4647-4654.|
|[PyANI](https://github.com/widdowquinn/pyani)|Version 0.2.12|Pritchard, L., Glover, R.H., Humphris, S., Elphinstone, J.G. and Toth, I.K., 2016. Genomics and taxonomy in diagnostics for food security: soft-rotting enterobacterial plant pathogens. Analytical methods, 8(1), pp.12-24.|
|[NCBI-Datasets](https://github.com/ncbi/datasets)|version 16.14.0|O’Leary, N.A., Cox, E., Holmes, J.B., Anderson, W.R., Falk, R., Hem, V., Tsuchiya, M.T., Schuler, G.D., Zhang, X., Torcivia, J. and Ketter, A., 2024. Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. Scientific Data, 11(1), p.732.|
|[Prokka](https://github.com/tseemann/prokka)|Version 1.14.6|Seemann, T., 2014. Prokka: rapid prokaryotic genome annotation. Bioinformatics, 30(14), pp.2068-2069.|
|[Cytoscape](https://github.com/tseemann/prokka)|Version 3.10.2|Shannon, P., Markiel, A., Ozier, O., Baliga, N.S., Wang, J.T., Ramage, D., Amin, N., Schwikowski, B. and Ideker, T., 2003. Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome research, 13(11), pp.2498-2504.|

## Tool intallation 

To install the environments with the different tools necessary to reproduce this analysis you could install the different yaml files we used for the analysis following the command below.

The yaml files themselves were downloaded from my environments by following the format below:
```bash
conda activate env
conda env export > env_name.yaml
```

The yaml files are located in this github repository in a folder called `yaml_file`. Firstly copy the yaml files to the HPC using rsync or scp and then type in the following commands. This will give you the same conda environment I have used in my analysis without having to install any packages yourself.

```bash
conda env create -f goatoolsenv.yaml
conda env create -f eggnogenv.yaml
conda env create -f pantaenv.yaml
conda env create -f rgienv.yaml
conda env create -f pyanienv.yaml
conda env create -f pipeline_pckgs.yaml
conda env create -f Renv.yaml
```

Alternatively, you may want to install all the tools manually and the guidance below outlines the scripts,steps or commands that have to be ran to install the tools needed to reproduce our results.

To install prokka,busco,ncbi datasets, and scoary run this command

```bash
sbatch pipeline_tools_install.sh
```

To install pyani run the command below

```bash
sbatch installing_pyani.sh
```

To install rgi run the command below

```bash
sbatch installing_RGI.sh 
```


To install panta follow the steps below

1.Download the [github repository](https://github.com/amromics/panta) for panta onto your local machine. The directory should be called panta-main

2.Unzip the file and scp it to the scripts folder on the cloud

3.Then move into the `panta-main` directory and run these commands in the command line:

```bash
conda create -y -c conda-forge -c defaults --name panta python=3.10 mamba
conda activate panta
mamba install -y -c conda-forge -c bioconda -c anaconda -c defaults  --file requirements.txt
pip install .
```


To install goatools we followed the series of steps outlined below:

1. Install goatools using this script.

```bash
sbatch installing_goatools.sh 
```

2. Activate the goatools env and install this module which caused an error when trying to run one of the goatools scripts

```bash
pip install docopt
```

3. Update goatools with the command below

```bash
pip install -U goatools
```


To set up an environment we use when plotting figures using R on the cloud run the command below:
```bash
sbatch installing_Renv.sh 
```


To install eggnog and all the associated files needed for the tool follow the steps below:

1. Go on the [github page for eggnog](https://github.com/eggnogdb/eggnog-mapper) and clone the repository, then send to ADA with scp or Rsync

2. Unzip the cloned repository on ADA

3. create a conda env for eggnog by running the command below

```bash
conda create --name eggnog python=3.8 -y
```

4. Go into the unzipped github directory and run

```bash
conda install --file requirements.txt
```

5. Add these two lines to your bash profile on ADA `~/.bash_profile`.

```bash
export PATH=~/scripts/eggnog_scripts/eggnog-mapper-master:~/scripts/eggnog_scripts/eggnog-mapper-master/eggnogmapper/bin:"$PATH"
export EGGNOG_DATA_DIR=~/all_gammaproteobacteria_data/eggnog-mapper-data
```

6.Make the directory that the database gets put into when running emapper

```bash
mkdir -p ~/all_gammaproteobacteria_data/eggnog-mapper-data
```


To install cytoscape follow the guidance on the [cytoscape web page](https://cytoscape.org/)



# THE ANALYSIS

You will notice a log_files directory, an all_gammaproteobacteria directory and a scripts directory.

## Part1 - Data Acquisition
Before conducting any sort of analysis we need to get the nucleotide and protein fastas for different Gammaproteobacteria, as well as their associated metadata. Running the `data_acquisition_pipeline.sh` script will run the whole data acquisition pipeline for the data acquisition stage of the analysis. It first gets all the different complete and refseq (GCF) assemblies for different gammaproteobcteria from ncbi. This includes The accession code , organism name, assembly release data, completeness, and a whole battery of other metrics are for each species is written to an output file. Then the resulting file is filtered to remove duplicate assemblies, and between different strains of the same species, take the most recent assembly. All the different assemblies for our final list of gammaproteobacteria are then downloaded as zip files, and we put metrics like assembly length into a txt file. The R script `distribution_plotting.r` is also ran to produce a plot showing the distribution of the GC content and assembly lengths which would become useful when doing phylogenetic trees because sometimes considering GC content is important.

1. Activate conda environment with the tools the script uses
```bash
conda activate pipeline_pckgs
```

2. This script contains a lot of intermediate python scripts within it so i will mention what they all do. The bash script itself  will obtain the names, GCF accession codes, and other metrics for all gammaproteobacteria with a refseq annotation (there were 17009). The `filtering_gammaprotobacteria.py` script then filters this dataset to remove multiple strains of the same species, this reduces the dataset down to 1348 species to work with. The `assembly_and_metrics_grab.py` script will then run ncbi-datasets again but this time we are downloading the actual nucloeitde and protein fasta files as a zip file, rather than just the species name and metrics. The most important part of this script below is that it produces a directory with the path `~/all_gammaproteobacteria_data/assemblies` and this is where the zip files are stored.
 
```bash
sbatch ~/scripts/data_acquisition_pipeline.sh
```

3. Run another script to unzip the zip files and put all nucleotide and protein fasta files into separate directories 
```bash
python3 ~/scripts/get_nucleotide_and_protein_fasta.py \
-i ~/all_gammaproteobacteria_data/assemblies -d1 ~/all_gammaproteobacteria_data/assemblies/nucleotide_fasta_files -d2 ~/all_gammaproteobacteria_data/assemblies/protein_fasta_files
```


Running this part of the analysis also produces a directory called `output_data`. Here you can find some plots showing the metrics of the inital 1348 species before any filtering. There is also names and NCBI accession codes for the initial 1348 species in the file called `filtered_gammaproteobacteria.txt `. The file called `all_gammaprotobacteria.txt` contains the GCF codes and other metrics for the 17006 Gammaproteobacteria with refseq annotations and this was filtered by the `filtering_gammaproteobacteria.py` script.

## Part2 - Quality Filtering

Now that we have successfully downloaded our genomes, we want to filter out genomes that have bad quality scores, and also remove any two genomes that appear genetically identical. To do this we first deduplicated the dataset using average nucelotide identity (ANI) and then remove low quality genomes using BUSCO.

### Part2A - Average Nucleotide Identity (ANI) Filtering
So now we want to calculate the ANI for all of the species we have in our dataset. This will then allow us to deduplicate the dataset by removing any 2 species that are found to be too similar. The similarity cutoff we use is 95% so where two species are found to have an ANI of 95%, one of them is chosen by the script to be removed from the dataset. Values of 100% that are a species being compared against itself are ignored by the script. To carry out this analysis follow the steps below

1. Run the main script that actually calculates ANI between species by running pyANI. The output will be sent to a directory called `~/all_gammaproteobacteria_data/ANIm_output `. This should take 1 to 3 days to complete depending on the resources provided
```bash
sbatch ~/scripts/calculating_ANI_mummer.sh
```

2. Now that all the ANI similarity matrices have been created you can run the script below to deduplicate the dataset and give you your list of ANI-filtered gammaproteobacteria. Enter the command below into the command line. This script has a python script with the path `~/scripts/deduplicate_dataset.py` inside of it. This script takes the ANI matrix produced by pyANI in the command above, and it uses some logic and code i wrote to remove one of any two species that are more than 95% similar in terms of ANI values. This then produces a list of ANI-filter gammaproteobacteria. This step reduced the initial genome set from 1348 to 1237. The file with ANI-filtered species has the path `~/all_gammaproteobacteria_data/output_data/final_deduplicated_gammaproteobacteria.txt`
```bash
sbatch ~/scripts/deduplicate_dataset_command.sh
```


### Part2B - BUSCO Filtering

All the bacteria in our study are gammaproteobateria and so we need to use that specific lineage dataset when running busco. To download the lineage dataset follow the steps below

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

This `run_busco.sh` file contains two other scripts called `busco_plotting.py` and `busco_quality_check.py`. The `busco_quality_check.py` is what actually runs busco by taking the list of nucloetide fasta files downloade as input, and using the gammaproteobacteria lineage dataset. This will produce an output directory called `~/all_gammaproteobacteria_data/busco_output` with the output files from busco. Next, the `busco_quality_check.py` script will sort through these output files form busco and using the list of ANI filtered genomes as a starting point, will start removing all of the ANI-fitered genomes that also have a busco cut off score less than what is specified by the user. This will produce a new output directory as well called `~/all_gammaproteobacteria_data/busco_plotting_output_and_summary_files`, and there will be a file called `busco_filtered_gammaproteobacteria.txt` in this directory. This wiill serve as input when we run prokka later on.


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

To begin with we will actually be building the pan-genome. We tried using the roary tool but it just didnt work with our data so we tried using panta and the run was able to finish. Keep in mind this takes 5 days to complete with max resources. To build the pan-genome panta needs the gff3 files produced from prokka so the first thing we need to do is run prokka to get gff files for our 1173 genomes

1.Running the command below will run prokka for all the 1173 genomes in the study. The gff files will be located in the directory with this path `~/all_gammaproteobacteria_data/prokka_output_files_busco_filtered_input/all_species_gff`.

```bash
sbatch ~/scripts/run_prokka.sh
```

2. Now that is done you can run Panta by running the command below. It will use the gff files just produced to build the pan-genome

```bash
sbatch ~/scripts/panta_scripts/panta_030pid_e7_LD07_split_1173genomes.sh
```

3. To create the plot for the pan-genome YOU MUST BE ON YOUR LOCAL MACHINE. You must also have matplot lib version 3.8.3 installed and numpy 1.25.0 installed as well. Once you have those installed go to the command line and run the command below. This will open a new window with the pan-genome plot present and you can take a screenshot or save it to your machine.

```bash
python3 pangenome_plot.py
```

## Part5 - Running Scoary
Now that we have the unique ARG families in our dataset, we essentially have the traits that will be used by scoary. For the purpose of our analysis we are interested in correlation and coincidence, rather than causative interpretations. To do this you can simply enter the two command into the command line.

## Part5A - Correlation

```bash
conda activate pipeline_pckgs

sbatch ~/scripts/scoary_scripts/main_scoary_pipeline_groups_1173genomes_NOPAIRWISE_with_groups.sh
```

This script does a lot of different things and uses a lot of intermediate scripts within it. Each script is commented but here is what they do. The two most important directories are the `~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/parsed_output_and_plotting_files_bonferroni/across_all_gene_families` and `all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families`. These contain files called `master_positive_correlated_genes.tsv` and `master_negative_correlated_genes.tsv`. Each of these files contain genes that were positively or negatively correlated with different the different ARG families we used, and they represent our pool of genes we use in downstream analysis. 


## Part5B - Cytoscape
Because the script above already produced the input needed for cytoscape there is no actual code for this section. Only a series of steps on how to reproduce the figures used in the paper.


### Figure X
To make figure X you first load the tsv file into cytoscape

when in cytoscape there is a search bar near the top-right of the window. Select this and search 'positive' , this will highlight all the genes with a positive correlation to the ARG. Then drag these to one side of the screen. Go the search bar again and search 'negative', this will highlight all the genes with a negative correlation to the ARG. 

Go to the styles menu bar

Tick Lock node width and height

Set size to 400 by selecting the Def. column, for Map. column select correlation type as the column, discrete mapping type, set size to 40 for positive and negative

For shape size select 'ellipse' for Def. column

For fill colour you can select the colour you want, but to get different pallettes for positively and negatively correlated genes select the Map. column and choose odds raito, select mapping type as continuous and then set the limits such that value from 0 to 1 are coloured based and value from 0 to 1 get lighter . and values from 1 and above get darker.

Then go to menu bar on mac and run . go to layout , bundle edges, and then All nodes and edges.

THen save the image and take the figure into powerpoint to add labels


Search CORRELATEDARGSW:N and collapse all the genes with similar coincidence patterns
Select all nodes in the star

## Part6 - Annotating the Pan-genome

We have successfully found genes that were positively or negatively correlated with the presence of the ARG families in this study. The next piece of analysis involves looking for enrichemnt of GO terms, and over-representation of COG categories in our positively and negatively correlated genes compared to the entire pan-genome. To do this, we first need to annotate the protein sequences for each gene in the pangenome with COG categories and GO terms. Then we can use GOATOOLS, and the SciPY python package to look for enrichment and over-representation.

# Part6A - EggNOG-mapper
To annotate the pan-genome with GO terms and COG categories run the command below. Verified by me

```bash
sbatch ~/scripts/eggnog_scripts/emapper_pangenome_input.sh
```
This will produce a file called something.annotations, and it contains a tsv of genes and their GO terms, and COG annotations.

# Part6B - Goatools

To accomplish this you can simply run the command below. Verified by me

```bash
sbatch ~/scripts/goatools_scripts/run_goatools_1173genomes_UPDATED.sh
```

Again, to breakdown the script a bit and the output you get from running it.

# Part6C - COG Categories
For the COG categories we want to perform a fishers exact test and then an odds ratio test to check for overrepresentation. This is accomplished by running the commands below. It contains the script called `cog_calculations_and_plots.py` which takes the genes in the pangenome, the significant positively correlated genes after applying Benjamini Hochberg correction method, and the negatively correlated genes after applying Benjamini Hochberg correction method. The script then looks for a significant difference in COG categories between the pangenome and both correlated gene sets. The most important files from running this script are `~/all_gammaproteobacteria_data/COG_calculation_1173_genomes/stats_output.txt` and `~/all_gammaproteobacteria_data/COG_calculation_1173_genomes/COG_plotting.R`. The first file mentioned contains pvalues and odds ratio values for all the statistics carried out in python and is the Supplementary Table S3 in our supplementary figures. The R script produces the plots we use in the main paper showing the percentage of COGs in the pangenome and differentially correlated gene sets. The plots produced by the R script are called `COG_category_plot_positive.pdf` and `COG_category_plot_negative.pdf`

```bash
sbatch run_COG_analysis.sh
```

Again, this script does a lot of different things. It takes the annotations file produced from EggNOG-mapper in [the EggNOG-mapper step of the analysis](#part6a---eggnog-mapper)

# Part7 - Supplementary Figures
There were a few supplmentary figures we produced for the report. Below is code and steps on how to get those figures.

### Plot for Pipeline
This bit of analysis needs to be carried out on your local machine. The first thing we need for this is a file called `project_pipeline.tsv`. Each line of this file details the 'flow' of the pipeline where one step (the source) is mapped to the next (the target). The layer of the heirachy we want the step to be on is also included on the line. We will also need to use Cytoscape for this and a python script that converts the `project_pipeline.tsv` file to a format acceptable py Cytoscape. To accomplish all this follow the steps below.

1. Create input for Cytoscape. Move into a directory with the project_pipeline.tsv file and run the command below.
```bash
python3 create_cytoscape_pipeline.py
```
This creates a new file called `cytoscape_project_pipeline.tsv` and will act as the input file for cytoscape

2. Open Cytoscape and load the `cytoscape_project_pipeline.tsv` file by selecting `Import Network From File system` then navigating through your file structure till you find the right file.

3. A pop up will appear showing which the first two columns columns in the tsv file are assigned as the target and source. The final column which details the position in the heirachy we want the different steps to be has the header `Layer` and we clicked on the black arrow next to the header name. From the drop-down menu that appears we select the `Source Node Attribute` option and this makes it so we can use one of Cytoscape layouts to organise the steps in our pipeline based on this column. 

4. When the network has loaded in go to `Layout` in the menu bar on your local PC and select yfiles Hierarchic Layout.

5. On the left of the cytoscape window there should be a menu bar with option like `Network` and `Style`. Select the `Style` option and change the style option from default to Minimal. 
   
6. After completing that, we should still be in the Style options for Nodes. At the bottom of the properties page uncheck the box saying `Lock node width and height`. We should be able to see the different properties for the nodes and these specific property values should be changed by selecting the `Def.` column and entering the values mentioned below.
   
Set Border Width to 1

Set Label Font Size 10

Set Width to 235

Set Height to 80
 
7. To adjust the colour of each node to your liking, select the `map.` column for the `Fill Color` property. There should be two rows called `Column` and `Mapping Type`. For `Column` select `name` from the dropdown menu, and for `Mapping Type` select `Discrete Mapping`. This will produce more rows for each step in the analysis and you can change the colours of their boxes here.

### Plot for busco stats and assembly length
This script will also have to be ran from your LOCAL MACHINE. Ensure you have R downloaded, and these version for these packages:
"ggplot2" -> ‘3.4.4’
tidyr -> ‘1.3.0’
gridExtra -> 2.3
grid -> ‘4.3.1’

Once those are all downloaded, you want to use scp or rsync to copy the file with all the busco metrics data for our 1173 genomes from ADA. The file you need has the path `~/all_gammaproteobacteria_data/busco_plotting_output_and_summary_files/sorted_by_complete_all_species_busco_metrics.tsv`. Once you have done that you can  run this command from your local machine. Make sure the R script is in your current directory. The file produced from this script will be called ``. The first argument of the command is the working directory you want to set, and the second is the path to the `sorted_by_complete_all_species_busco_metrics.tsv` file on your machine.

```bash
Rscript plot_busco_completeness.R ./ sorted_by_complete_all_species_busco_metrics.tsv
```

### Plot for ANI values
To create the figure for the ANI values between species you need to copy the file with this path to your LOCAL MACHINE. The file path is `all_gammaproteobacteria_data/ANIm_output/ANIm_percentage_identity.tab`. Once on your local machine, you need to run the commands below and you will produce a plot called `ANI_distribution.pdf`.

1. First run this python command which takes the ANIm identity matrix produced by ANI.

```bash
python3 ANImatrix_to_list.py --output_directory ./ -i ANIm_percentage_identity.tab 
```

2. Using the `ANI_values.txt` file we can plot the distribution of ANI values in R by running the command below. Make sure the R script and the `ANI_values.txt` file are in the same directory before running the command. This will produce the plot called `ANI_distribution.pdf`

```bash
Rscript ani_plots.R
```


# Conclusion
