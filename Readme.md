# PanSweep

PanSweep is a tool to evaluate gene-level metagenomic analysis for differentially present genes by species, to determine whether these genes are likely to be pangenome contaminants, and to provide functional information about these genes.

PanSweep is designed to work with the output of MIDAS2 (https://github.com/czbiohub-sf/MIDAS) run on the UHGG database, plus [Parquet databases](https://zenodo.org/records/14852853) from Zenodo. For an example workflow, see https://github.com/pbradleylab/cirrhosis-pansweep/. Example test data is also available at Zenodo: https://zenodo.org/uploads/13891285

The package is designed to run two parts separately, 1. an analysis of the MIDAS2 output, and 2. the Shiny visualization of the results. The analysis should be run on a HPC cluster or a workstation with sufficient RAM, while the visualization can be run on personal computers. 

The analysis and results evaluation are fully described in [the PanSweep manuscript](https://www.biorxiv.org/content/10.1101/2024.10.11.617902v1).

# Installing PanSweep Package

PanSweep can be installed directly from Github using devtools.
~~~~
devtools::install_github("https://github.com/pbradleylab/pansweep")
~~~~

PanSweep depends on the following additional packages:

    * tidyverse (>= 2.0.0)
    * progress (>= 1.2.3)
    * readr (>= 2.1.4)
    * viridisLite (>= 0.4.2)
    * arrow (>= 13.0.0)
    * pillar (>= 1.9.0)
    * dbplyr (>= 1.1.2)
    * parallelDist (>= 0.2.6)
    * umap (>= 0.2.10.0)
    * vegan (>= 2.6-4)
    * jsonlite (>= 1.8.8)
    * shiny (>= 1.7.5.1)
    * plotly (>= 4.10.3)
    * DT (>= 0.33)
    * knitr (>= 1.43)
    * kableExtra (>= 1.3.4)
    * DiscreteFDR

# Running PanSweep Analysis

The PanSweep Analysis is designed to be run on a HPC cluster or on a workstation with at least 15GB of memory.

First, download the [Parquet databases](https://zenodo.org/uploads/13891285) from Zenodo and unzip them. They can be unzipped into any directory, as long as the correct path is added to the config file.

Next, prepare the JSON config file with the below paths. An empty config file is provided on this repository. If any paths are missing, the analysis will be unable to run. For an example please reference the `Config.json` file on this repository, or the `cirrhosis.json` file in the cirrhosis-pansweep repository: https://github.com/pbradleylab/cirrhosis-pansweep/tree/main/pansweep


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

Note: the metadata file should be a tab-delimited file with three columns named "sample", "subject", and "env" (for environment). The "subject" column allows metagenomic samples from the same subject to be grouped/averaged together. The "env" column gives the status of the subject, e.g., "case" or "control." Currently, only two conditions can be compared.

After the config file is saved, the analysis is ready to be run. Use the PanSweep_Analysis function outlined below:

~~~~    
PanSweep_Analysis(Json_Config_Path = "Your/JSON/Path/here.JSON")
~~~~

**Note:** The Corr_lower_limit parameter is set to 3 automatically. This sets the number of significant genes needed per species to run the co-occurrence analysis, which includes the heatmap and ordination plots. This limit can be altered by changing the appropriate function argument, but it is not recommended to go lower than three for the analysis as the results will not be very meaningful.

**Note:** The analysis may take some time to finish. You can obtain more detailed progress information by adding the flag `verbose=TRUE` to the function call.

The results of the analysis will be saved in a date stamped folder called "PanSweep_Analysis_Output_YYYY-MM-DD" as the file called "PanSweep_Analysis_Output.rds"

# Running PanSweep Shiny

The PanSweep Shiny GUI is designed to use the "PanSweep_Analysis_Output.rds" from the PanSweep_Analysis function to allow for local evaluation of the results. The GUI is broken into "Analysis Report", "eggNOG & Correlation Report", "Ordination & Heatmap", and "NMDS". The PanSweep Shiny GUI is run using the below function:

~~~~
PanSweep_Shiny(loadData_Path = "Path/to/file/PanSweep_Analysis_Output_YYYY-MM-DD/PanSweep_Analysis_Output.rds")
~~~~
### Analysis Report

The analysis report provides information on the number of genes found to be significant, the number of species that have significant genes, the number of genes per species, and if there are any repeated UHGP-90 and UHGP-50 cluster ids. 

### eggNOG & Correlation Report

This report provides infromation on the individual genes from eggNOG from the UHGP-90 cluster ids. The lineage test reports if the most-correlated species matches the pangenome's annotated species at the family level or lower. The FDR-corrected p-values for the Fisher's exact test are provided. 

### Ordination & Heatmap

Ordination plots for UMAP, NMDS, and PCoA are made from the Jaccard co-occurrence matrix of genes by sample. The numerical sliders allow the adjustment of the n_neighbors and min_dist values for the UMAP ordination plot. The Jaccard similarity heatmap is based off of co-occurrence of genes by sample.  If the heatmap or ordination plots are clicked the gene on the heatmap will be highlighted. 

### NMDS
The NMDS provides the NMDS ordination plot and stress plot by species for evaluation of NMDS results.
