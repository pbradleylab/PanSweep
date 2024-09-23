# PanSweep
PanSweep is a tool to evaluate gene-level metagenomic analysis for differentially present genes by species, evaluate to determine if the gene are contamination from different species and provide functional information these genes. PanSweep is designed to work with MIDAS2 output and parquet databases from figshare *HERE LINK*. The package is designed to run two parts separately, an analysis and user interface. The analysis should be run on a HPC cluster while the user interface is can be run on personal computers. 

The analysis and results evaluation is fully described in *PAPER HERE*

# Installing PanSweep Package

PanSweep can be installed by using dev.tools to install from github. The following packages need to be installed for each part of the analysis.

~~~~
devtools::install_github()
~~~~

PanSweep Analysis:


PanSweep Shiny:


# Running PanSweep Analysis

The PanSweep Analysis is designed to be run on a HPC cluster. It will need at least 15GB of memory to run. 

First, download the parquet databases from figshare *HERE* and place them in a directory the R analysis can read.

Next, prepare the JSON config file with the below paths. An empty config file is provided on this repository. If any paths are missing the analysis will be unable to run. For an example please reference the example.json file on this repository.

|Path Name                         |Description                                                            |
|----------------------------------|-----------------------------------------------------------------------|
|Species_Set                       |An array of target species using MIDAS IDs                             |
|path_uhgp_50_cluster              |Path to the folder that holds the uhgp_50_cluster parquet file         |
|path_uhgp_90_cluster              |Path to the folder that holds the uhgp_90_cluster parquet file         |
|path_uhgp_50_eggNOG               |Path to the folder that holds the uhgp_50_eggNOG parquet file          |
|path_uhgp_90_eggNOG               |Path to the folder that holds the uhgp_90_eggNOG parquet file          |
|path_for_genomes_all_metadata     |metadata.tsv file for genomes from UHGG                                |
|path_to_sample_metadata           |Path to metadata tsv for samples                                       |
|path_to_species_abundance         |species_marker_read_counts.tsv from MIDAS                              |
|path_to_presabs                   |path to folder that holds the species files for presabs and gene_counts|
|save_folder_location              |path to user designated save folder                                    |  

After the config file is saved, the analysis is ready to be run. Use the PanSweep_Analysis function outlined below:

~~~~    
PanSweep_Analysis(Json_Config_Path = "Your/JSON/Path/here.JSON")
~~~~

**Note:** The Corr_lower_limit variable is set to 3 automatically. This sets the number of significant genes needed per species to run the co-oocurance analysis, which includes the heatmap and ordination plots. This limit can be altered by changing the varriable but is is not recommended to go lower than 3 for the analysis.

The results of the analysis will be saved in a date stamped folder called "PanSweep_Analysis_Output_YYYY-MM-DD" as the file called "PanSweep_Analysis_Output.rds"

# Running PanSweep Shiny

The PanSweep Shiny GUI is designed to use the "PanSweep_Analysis_Output.rds" from the PanSweep_Analysis function to allow for local evaluation of the results. The GUI is broken into "Analysis Report", "eggNOG & Correlation Report", "Ordination & Heatmap", and "NMDS". The PanSweep Shiny GUI is run using the below function:

~~~~
PanSweep_Shiny(loadData_Path = "Path/to/file/PanSweep_Analysis_Output_YYYY-MM-DD/PanSweep_Analysis_Output.rds")
~~~~
### Analysis Report

The analysis report provides information on the number of genes found to be significant, the number of species that have significant genes, the number of genes per species, and if there are any repeated UHGP-90 and UHGP-50 cluster ids. 

### eggNOG & Correlation Report

This report provides infromation on the individual genes from eggNOG from the UHGP-90 cluster ids. The linage test repots if the highest correlated species matches the genes originating pangenomes to the family level. The FDR for the Fisher's exact test multiple testing correction are provided. 

### Ordination & Heatmap

Ordination plots for UMAP, NMDS, and PCoA are made from the Jaccard co-occurrence matrix of genes by sample. The numerical sliders allow the adjustment of the n_neighbors and min_dist values for the UMAP ordination plot. The Jaccard similarity heatmap is based off of co-occurance of genes by sample.  If the heatmap or ordination plots are clicked the gene on the heatmap will be highlighted. 

### NMDS
The NMDS provides the NMDS ordination plot and stress plot by species for evaluation of NMDS results.
