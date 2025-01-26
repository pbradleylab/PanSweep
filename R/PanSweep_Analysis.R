#'@import dplyr
#'@import progress
#'@import readr
#'@import purrr
#'@import viridisLite
#'@import arrow
#'@import pillar
#'@import dbplyr
#'@import parallelDist
#'@import umap
#'@import vegan
#'@import jsonlite
#'@import tibble
#'@import DiscreteFDR
#'@import stringr
#'@import tidyr
#'@import rlang
NULL

#' PanSweep Analysis:
#'
#' This function run the analysis of MIDAS2 data and outputs a save file to be
#' used with the PanSweep_Shiny() function. The inputs are Json_Config_Path and
#' Co_occurrence_lower_limit, which is set to 3. Refer to the PanSweep github for
#' the Json file and explanations.
#'
#' To run the function:
#' PanSweep_Analysis(Json_Config_Path = "Path/To/Json/file.json")
#'
#' To change the lower bound for number of genes needed for co-occurence
#'  analysis:
#' PanSweep_Analysis(Json_Config_Path = "Path/To/Json/file.json",
#'                    Co_occurrence_lower_limit = 5)
#'
#' @param Json_Config_Path Path to the completed json config file
#' @param Co_occurrence_lower_limit  lower bound for number of genes needed for
#' co-occurrence analysis. Set to 3.
#' @param Pres_abs_lower_limit  lower bound for number of times a gene must be
#' present and absent in order for it to be counted. Set to 0 (no filtering) by
#' default since DiscreteFDR automatically handles these cases well.
#' @param Max_FDR  False discovery rate threshold. Set to 0.05 by default.
#' @param merge_fn_binary  Function to use for merging binary data from different samples but the same subject. Default is `base::max`.
#' @param merge_fn_counts  Function to use for merging count data from different samples but the same subject. Default is `base::sum`.
#' @param signif_test_function  Function to use to take a given table of gene counts and return p-values, corrected p-values, and sample sizes; this allows users to override the built-in statistical test. Default is `analyze_tbl`.
#' @param correlation_function  Function to use to correlate gene and species matrices; this allows users to override the built-in correlation test. Default is `spearman_cor_wrapper`.
#' @param save_folder_location  String. Path to save output results. Default is NULL; will override JSON value only if not NULL.
#' @param return_not_save  Boolean. If FALSE, save PanSweep output as .rds files; if TRUE, simply return the data structure. Useful for running PanSweep within a larger pipeline. Default is FALSE.
#' @param verbose  Boolean. Print messages; default is FALSE.
#' @return Returned a date stamped folder called PanSweep_Analysis_Output_YYYY-MM-DD
#' containing the file "PanSweep_Analysis_Output.rds".
#'
#'@export
PanSweep_Analysis <- function(Json_Config_Path,
                              Co_occurrence_lower_limit = 3,
                              Pres_abs_lower_limit = 0,
                              Max_FDR = 0.05,
                              merge_fn_binary = base::max,
                              merge_fn_counts = base::sum,
                              signif_test_function = analyze_tbl,
                              correlation_function = spearman_cor_wrapper,
                              save_folder_location = NULL,
                              return_not_save = FALSE,
                              verbose = FALSE
                              ) {
  ####################################################################################################################
  #Load in Paths and variables via JSON#
  Paths_and_Variables <-fromJSON(Json_Config_Path)
  #Variables:
  Species_Set <- Paths_and_Variables$Variables$Species_Set
  #Databases:
  path_uhgp_50_cluster <- normalizePath(Paths_and_Variables$Databases$path_uhgp_50_cluster)
  path_uhgp_90_cluster <- normalizePath(Paths_and_Variables$Databases$path_uhgp_90_cluster)
  path_uhgp_50_eggNOG <- normalizePath(Paths_and_Variables$Databases$path_uhgp_50_eggNOG)
  path_uhgp_90_eggNOG <- normalizePath(Paths_and_Variables$Databases$path_uhgp_90_eggNOG)
  #Metadata:
  path_for_genomes_all_metadata <- normalizePath(Paths_and_Variables$Metadata$path_for_genomes_all_metadata)
  path_phylo_md2 <- normalizePath(Paths_and_Variables$Metadata$path_to_sample_metadata)
  path_to_presabs <- normalizePath(Paths_and_Variables$Metadata$path_to_presabs)
  path_to_species_abundance <- normalizePath(Paths_and_Variables$Metadata$path_to_species_abundance)
  path_to_read_counts <- normalizePath(Paths_and_Variables$Metadata$path_to_read_counts)
  is_compressed <- Paths_and_Variables$Metadata$gene_data_is_compressed
  #Output:
  # If the user did not provide anything, we use the value from the JSON file, otherwise we use what they provided
  if (is.null(save_folder_location)) {
    save_folder_location <- normalizePath(Paths_and_Variables$Output$save_folder_location)
  }
  if (is.null(save_folder_location)) { error("Neither user nor JSON file provided a save location")}
  if (!dir.exists(save_folder_location)) { dir.create(save_folder_location) }
  #Change variable name:
  Corr_lower_limit <- Co_occurrence_lower_limit

  arg_to_string <- function(f) { rlang::as_label(rlang::enexpr(f)) }

  ####################################################################################################################

  if (length(Species_Set) == 0) { stop("No species provided") }

  if (verbose) message("Reading in presence-absence data...")
  phylo_md2 <- read_tsv(path_phylo_md2, show_col_types = FALSE)
  #Check for files:

  gene_paths <- purrr::map_chr(Species_Set, ~ {
    gpath <- file.path(path_to_read_counts, .x, paste0(.x, ".genes_presabs.tsv"))
    if (is_compressed) {
      gpath <- paste0(gpath, ".lz4")
    }
  }) %>% setNames(Species_Set)
  gene_paths_exist <- purrr::map_lgl(gene_paths, file.exists)
  if (any(!gene_paths_exist)) {
    iwalk(gene_paths_exist, \(x, idx) if (!x) {
      warning(paste0("Missing gene variant file: ", gene_paths[idx]))
    })
    Species_Set <- names(which(gene_paths_exist))
    gene_paths <- gene_paths[Species_Set]
  }
  if (length(Species_Set) == 0) { stop("No species had gene variant files") }

  test_tbls <- purrr::map(gene_paths, \(gpath) {
    if (is_compressed) {
      tsv <- arrow::read_tsv_arrow(gpath)
    } else {
      tsv <- readr::read_tsv(gpath, col_types=cols())
    }
    merge_columns_tbl(tsv, md=phylo_md2, fn=merge_fn_binary)
  })

  if (verbose) message("Calculating tests...")
  # run Fisher test analysis
  pb <- progress::progress_bar$new(total = length(test_tbls))
  test_results <- purrr::map(test_tbls, ~ {
    pb$tick()
    analyze_tbl(.x, md=phylo_md2)
  })
  names(test_results) <- Species_Set

  # ...next:
  # extract for each fdrs in each member of test_results, names of genes with fdr <= Max_FDR
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Extracting genes with fdrs<=Max_FDR#
  test_results_fdrs <- purrr::map(test_results,
                                  function(x) tibble::enframe(x$fdrs,
                                                              name = "Gene_id",
                                                              value = "Fdrs"))
  #enframe is used in named vector
  #creates a list under species_id and removes other unneeded lists (ie pvals)
  test_results_lt_tbl <- tibble::enframe(test_results_fdrs, name = "Species_id")
  #creates a list of tibbles under species
  test_results_tbl <- tidyr::unnest(test_results_lt_tbl, cols = c(value))
  #makes a tbl of the test results ALL fdrs are included

  Gene_extract_tbl <- dplyr::filter(test_results_tbl, Fdrs<=Max_FDR)
  #Extract all genes by (with species) that are Fdrs<=Max_FDR
  Genes_of_intr <- unique(Gene_extract_tbl$Gene_id)


  Genes_intr_extr <- dplyr::filter(test_results_tbl, Gene_id %in% Genes_of_intr)
  # Extract all genes from the test results tbl (ALL fdrs) that are in the "Genes_of_interest" vector.

  #start list for report:
  Number_of_Significant_Genes <- nrow(Gene_extract_tbl)
  if (is.null(Number_of_Significant_Genes)) {
    stop("Error: No significant genes identified")
  }
  if (Number_of_Significant_Genes == 0) {
    stop("Error: No significant genes identified")
  }



  if (verbose) message("Getting extra info from Parquet databases...")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Adding cluster ids#
  Gut_Trans <- gsub("UHGG","", Gene_extract_tbl$Gene_id) %>%
    gsub("_", "", .) %>%
    as.numeric()

  prqt_uhgp_90 <- arrow::open_dataset(path_uhgp_90_cluster, format = "parquet")
  prqt_uhgp_50 <- arrow::open_dataset(path_uhgp_50_cluster, format = "parquet")

  Gene_intr_uhgp_90 <- prqt_uhgp_90 %>%
    filter(gene_id %in% Gut_Trans) %>%
    dplyr::collect() %>%
    dplyr::select(-c(group)) %>%
    rename(cluster_id_n = cluster_id)

  Gene_intr_uhgp_90 <- Gene_intr_uhgp_90 %>%
    map_dfc(function(x) as.character(x)) %>%
    map_dfc(function(x) str_pad(x, 11, "left", pad = "0"))

  Gene_intr_uhgp_50 <- prqt_uhgp_50 %>%
    filter(gene_id %in% Gut_Trans) %>%
    dplyr::collect() %>%
    dplyr::select(-c(group)) %>%
    rename(cluster_id_n = cluster_id)

  Gene_intr_uhgp_50 <- Gene_intr_uhgp_50 %>%
    map_dfc(function(x) as.character(x)) %>%
    map_dfc(function(x) str_pad(x, 11, "left", pad = "0"))

  Genes_intr_extr <- mutate(Genes_intr_extr,
                            "gene_id_n"=gsub("UHGG",
                                             "",
                                             Gene_extract_tbl$Gene_id) %>%
                                gsub("_", "", .))
  UHGP_50_genes_of_int <- left_join(Genes_intr_extr,
                                    Gene_intr_uhgp_50,
                                    by=c("gene_id_n" = "gene_id"))
  UHGP_90_genes_of_int <- left_join(Genes_intr_extr,
                                    Gene_intr_uhgp_90,
                                    by=c("gene_id_n" = "gene_id"))
  #______________________________________________________________________________#
  #Return Gut_Genome Syntax#

  UHGP_90_genes_of_int <- mutate(UHGP_90_genes_of_int, "cluster_id" =
                                   gsub('^(.{6})(.*)$','GUT_GENOME\\1_\\2',
                                        UHGP_90_genes_of_int$cluster_id_n))

  UHGP_50_genes_of_int <- mutate(UHGP_50_genes_of_int, "cluster_id" =
                                   gsub('^(.{6})(.*)$','GUT_GENOME\\1_\\2',
                                        UHGP_50_genes_of_int$cluster_id_n))
  #______________________________________________________________________________#
  #Look for cluster repeats#
  UHGP_90_cluster_id_summ <-UHGP_90_genes_of_int %>%
    group_by(cluster_id) %>%
    filter(n()>1) %>%
    count()

  Num_rep_UHGP_90_clus_id <- sum(UHGP_90_cluster_id_summ$n)

  UHGP_50_cluster_id_summ <-UHGP_50_genes_of_int %>%
    group_by(cluster_id) %>%
    filter(n()>1) %>%
    count()

  Num_rep_UHGP_50_clus_id <- sum(UHGP_50_cluster_id_summ$n)
  #______________________________________________________________________________#
  if (verbose) message("Adding EggNOG information...")
  #Assign eggNOG Taxonomy and gene info#
  prqt_uhgp_50_eggNOG <- arrow::open_dataset(path_uhgp_50_eggNOG, format = "parquet")
  prqt_uhgp_90_eggNOG <- arrow::open_dataset(path_uhgp_90_eggNOG, format = "parquet")

  uhgp_90_eggNOG <- prqt_uhgp_90_eggNOG %>%
    filter(query_name %in% UHGP_90_genes_of_int$cluster_id) %>%
    select(query_name, Predicted_taxonomic_group, Predicted_protein_name, eggNOG_free_text_description) %>%
    dplyr::collect() %>%
    right_join(UHGP_90_genes_of_int, by = c("query_name" = "cluster_id"), keep = TRUE)

  uhgp_50_eggNOG <- prqt_uhgp_50_eggNOG %>%
    filter(query_name %in% UHGP_50_genes_of_int$cluster_id) %>%
    select(query_name, Predicted_taxonomic_group, Predicted_protein_name, eggNOG_free_text_description) %>%
    dplyr::collect() %>%
    right_join(UHGP_90_genes_of_int, by = c("query_name" = "cluster_id"), keep = TRUE)
  #______________________________________________________________________________#
  #Add in Taxonomy info#
  genome_metadata <- read_tsv(file= path_for_genomes_all_metadata, show_col_types=FALSE)



  meta_genome_sep_taxa <- genome_metadata %>% separate_taxonomy_with_s("Lineage") %>%
    mutate(species_id = as.character(species_id))

  uhgp_90_eggNOG <- uhgp_90_eggNOG %>%
    left_join(meta_genome_sep_taxa, by = c("Species_id" = "species_id"), keep = TRUE)

  # uhgp_50_eggNOG <- uhgp_50_eggNOG %>%
  #   mutate("Genome" = gsub('^(.{6})(.*)$','GUT_GENOME\\1', uhgp_50_eggNOG$Species_id)) %>%
  #   left_join(meta_genome_sep_taxa, by = "Genome", keep = TRUE)
  #Clean up extra columns.
  #______________________________________________________________________________#
  #Number of genes per species#
  Num_Sig_Genes_per_sp <-Gene_extract_tbl %>%
    group_by(Species_id) %>%
    count()

  Test <- sum(Num_Sig_Genes_per_sp$n)
  #______________________________________________________________________________#
  #Separate by species#
  test_tbls_sp <- test_tbls
  names(test_tbls_sp) <- Species_Set
  Runs_by_Genes_intr<- lapply(test_tbls_sp, function(x)  subset(x,x$gene_id %in% Genes_of_intr))
  Sp_to_corr_Runs <- keep(Runs_by_Genes_intr, function(x) nrow(x)>Corr_lower_limit) %>%
    lapply(., function(x) column_to_rownames(x, var = 'gene_id'))

  ################################################################################
  #Shiny Prep#
  #______________________________________________________________________________#
  if (verbose) message("Computing distance matrices...")
  #Correlate by runs#
  Sp_corr <- lapply(Sp_to_corr_Runs, function(x) parDist(as.matrix(x), method = 'binary'))
  M.Sp_corr <- lapply(Sp_corr, function(x) as.matrix(x))
  if (verbose) message("Computing ordinations...")
  #______________________________________________________________________________#
  #UMAP_IT#
  Species <- Num_Sig_Genes_per_sp$Species_id

  U.Sp_corr <- lapply(M.Sp_corr, function(x) {
    max_n_neighbs <- ceiling(nrow(x) / 3)
    if ((max_n_neighbs - 2) > 10) {
      umap_step_size = ceiling((max_n_neighbs-2)/10)
    } else {
      umap_step_size = 1
    }
    n_n <- seq.int(from = 2, to = max_n_neighbs, by = umap_step_size)
    result_nn <- lapply(n_n, function(y) {
      min_dist_seq  <- seq(from = 0.1, to = 0.9, by = 0.1)
      umap_seq_result <- lapply(min_dist_seq, function(m){
        umap(x, n_neighbors = y,min_dist = m, input = 'dist')
      })
      names(umap_seq_result) <- min_dist_seq
      return(umap_seq_result)
    })
    names(result_nn) <- n_n
    return(result_nn)
  })
  #______________________________________________________________________________#
  #NMDS#
  N.Sp_corr <- lapply(M.Sp_corr, function(x) {
    metaMDS(x,
            distance = jaccard,
            trace = (as.numeric(verbose) - 1))
  })
  N.Stress <- lapply(N.Sp_corr, function(x) round(x$stress, 3))
  N.StressPlot <- lapply(N.Sp_corr, function(x) stressplot(x))
  #______________________________________________________________________________#
  #PCoA#
  P.Sp_corr <- lapply(M.Sp_corr, function(x) cmdscale(x, k = 2, eig = TRUE))
  #______________________________________________________________________________#
  #Significant genes to species correlation#
  #!!!Read in of docs needs to be fixed!!#
  #Now merges by default#
  if (verbose) message("Performing gene-to-species correlation lineage test...")
  #Determine significant species:
  species_set_corr <- unique(Genes_intr_extr$Species_id)
  #Read in and merge gene counts:
  Gene_reads <- purrr::map(species_set_corr, ~ {
    gpath <- file.path(path_to_read_counts, .x, paste0(.x, ".genes_reads.tsv"))
    if (is_compressed) {
      gpath <- paste0(gpath, ".lz4")
      tsv <- read_tsv_arrow(gpath)
    } else {
      tsv <- read_tsv(gpath, col_types=cols())
    }
    merge_columns_tbl(tsv, md=phylo_md2, fn=merge_fn_counts)
  }) %>% setNames(species_set_corr)

  Sig_Gene_reads <- purrr::map(species_set_corr, ~ {
    dplyr::filter(Gene_reads[[.x]], gene_id %in% Genes_intr_extr$Gene_id) %>%
      column_to_rownames("gene_id") %>%
      as.matrix()
  }) %>% setNames(species_set_corr)
 #Prepare species for correlation
  Species_Abd <-read_tsv(path_to_species_abundance, show_col_types=FALSE) %>%
    merge_columns_tbl(md=phylo_md2, fn=merge_fn_counts) %>%
    column_to_rownames("species_id") %>%
    as.matrix()

  Corr_Results_All <- species_to_gene_correlations(Sig_Gene_reads, Species_Abd, meta_genome_sep_taxa)

  #______________________________________________________________________________#
  #Add lineage and max species correlation to eggNOG table#
  uhgp_90_eggNOG <- left_join(uhgp_90_eggNOG,
                              Corr_Results_All$Cor_Sp_Ln_max_sp_DF %>%
                                select("Gene",
                                       "Lineage_Shared",
                                       "cor_max_species"),
                              by = c("Gene_id" = "Gene"))
  uhgp_90_eggNOG <- left_join(uhgp_90_eggNOG,
                              Corr_Results_All$Species_Cor_DF %>%
                                filter(mark == "max_f") %>%
                                select(Gene, rank),
                              by = c("Gene_id" = "Gene")) %>%
    rename(Family_max_rank = rank)
  uhgp_90_eggNOG <- left_join(uhgp_90_eggNOG,
                              Corr_Results_All$Species_Cor_DF %>%
                                filter(mark == "target") %>%
                                select(Gene, rank),
                              by = c("Gene_id" = "Gene")) %>%
    rename(Sp_rank = rank)
  #______________________________________________________________________________#
  #Create Report#
  if (verbose) message("Creating report and outputting results...")
  Analysis_report <- tibble::as_tibble(cbind(c("Number of Significant Genes", "Number of Species with Significant Genes", "Number of repeated UHGP-90 ids", "Number of repeated UHGP-50 ids"),
                                     c(Number_of_Significant_Genes, nrow(Num_Sig_Genes_per_sp), Num_rep_UHGP_90_clus_id, Num_rep_UHGP_50_clus_id)))
  Analysis_output_names <- c('Analysis_report','Number_of_Significant_Genes', 'UHGP_90_cluster_id_summ', 'UHGP_50_cluster_id_summ', 'Num_Sig_Genes_per_sp', 'uhgp_90_eggNOG')
  Analysis_output <- setNames(mget(Analysis_output_names), Analysis_output_names)
  # Allow us to track how the analysis was run
  PanSweep_Analysis_Parameters <- rlang::call_match(defaults=TRUE)
  ################################################################################
  #Create save files#
  All_RDS_to_Save <- c("M.Sp_corr", "U.Sp_corr", "N.Sp_corr", "N.Stress", "P.Sp_corr", "Analysis_output", "PanSweep_Analysis_Parameters")
  save_folder_name <- paste0("PanSweep_Analysis_Output_", format(Sys.Date(), "%Y-%m-%d"))
  save_folder_location_full <- file.path(save_folder_location, save_folder_name)
  if (!return_not_save) {
    dir.create(save_folder_location_full)
    MIDAS_Analysis_Output <- setNames(mget(All_RDS_to_Save), All_RDS_to_Save)
    saveRDS(MIDAS_Analysis_Output, file.path(save_folder_location_full, paste0("PanSweep_Analysis_Output", ".rds")))
    saveRDS(Analysis_output, file.path(save_folder_location_full, paste0("PanSweep_Analysis_TablesOnly", ".rds")))
  } else {
    return(MIDAS_Analysis_Output)
  }
  ################################################################################
}

## Helper functions

#' Merge all columns in a matrix that correspond to the same subject.
#'
#' This function will "merge" all columns that correspond to the same subject in a metadata table. Merge behavior can be controlled by passing a function. The default merge behavior is to take the maximum (`max()`).
#'
#' @param mtx Matrix to merge.
#' @param md Metadata tbl/df. Must have the columns "sample", "env" (environment, e.g., "case" and "control"), and "subject". All "sample" columns for the same "subject" will be merged.
#' @param fn Function to use for merging. Default is `base::max()`.
#' @export
merge_columns <- function(mtx, md, fn=base::max) {
  unique_subj <- unique(dplyr::filter(md, sample %in% colnames(mtx))$subject)
  output_mtx <- matrix(nrow=nrow(mtx), ncol=length(unique_subj),
                       dimnames = list(rownames(mtx), unique_subj))
  for (s in unique_subj) {
    samples <- md$sample[md$subject==s]
    which_cols <- intersect(samples, colnames(mtx))
    output_mtx[, s] <- apply(mtx[, which_cols, drop=FALSE], 1, fn)
  }
  output_mtx
}

#' Alternative way to merge -- tidyverse friendly but slower
merge_by_pivoting <- function(tbl, md, fn=base::max) {
  long_tbl <- pivot_longer(tbl, -1, names_to="sample") %>%
    left_join(., md, by=c("sample"))
  merged <- long_tbl %>%
    group_by(gene_id, subject) %>%
    summarize(V=fn(value), .groups="keep")
  pivot_wider(merged, names_from=subject, values_from=V)
}


#' Merge all columns in a tbl that correspond to the same subject.
#'
#' Wrapper around merge_columns() that allows seamless use with tbls. This function will "merge" all columns that correspond to the same subject in a metadata table. Merge behavior can be controlled by passing a function. The default merge behavior is to take the maximum (`max()`).
#'
#' @param tbl Tbl to merge.
#' @param md Metadata tbl/df. Must have the columns "sample", "env" (environment, e.g., "case" and "control"), and "subject". All "sample" columns for the same "subject" will be merged.
#' @param fn Function to use for merging. Default is `base::max()`.
#' @param rowcol Which column contains the row names? Default is 1.
#' @export
merge_columns_tbl <- function(tbl, md, fn=base::max, rowcol=1) {
  mtx <- as.matrix(tbl[, -rowcol])
  rownames(mtx) <- tbl[[rowcol]]
  output_mtx <- merge_columns(mtx, md, fn)
  output_tbl <- tibble::as_tibble(output_mtx, rownames=colnames(tbl)[rowcol])
  output_tbl
}

#'Analyze gene tables with Fisher's exact test and report discrete FDR.
#'
#'This function analyzes a gene table with an associated metadata table by performing a Fisher's exact test per gene, then correcting the resulting p-values using DiscreteFDR. Genes may also optionally be filtered by the minimum number of presences and absences.
#'
#' @param tbl Tbl or data.frame of the gene presence/absence table. Assumes that the first column is the gene ID and that all other columns are 1/0.
#' @param md  Metadata tbl. Must have the columns "sample", "env" (environment, e.g., "case" and "control"), and "subject". All "sample" columns for the same "subject" will be merged by taking the max value.
#' @param min_obs  Integer. Minimum number of times a gene must be
#' present and absent in order for it to be counted. Set to 0 (no filtering) by
#' default since DiscreteFDR automatically handles these cases well.
#' @param merge Boolean. If TRUE, merge samples from the same subject. If FALSE, assume samples are already merged. Default is FALSE.
#' @param merge_fn Function. How to merge samples from the same subject? Default is `base::max`.
#' @param verbose Boolean or integer. If not FALSE/0, print messages to indicate where we are in the process.
#' @return Returns a list with p-values (`pvals`), adjusted p-values (`fdrs`), a list of subjects per condition (`subjects_per_condition`), and the filtered and merged data (`clean_mtx`).
#'
#'@export
analyze_tbl <- function(tbl, md, min_obs = 0, merge=FALSE, merge_fn = base::max, verbose=FALSE) {
  # convert tibble to matrix with rownames
  mtx <- as.matrix(tbl[,-1]) #-1 = do not grab first column into matrix
  rownames(mtx) <- tbl$gene_id #name the first row based off of the tibble
  if (merge) {
    if (verbose) message("Merging...")
    merge_mtx <- merge_columns(mtx, md, max)
  } else {
    merge_mtx <- mtx
  }
  if (verbose) message("Filtering...")
  # filter out genes that don't have at least min_obs observations (and min_obs non-observations)
  which_rows <- apply(merge_mtx, 1, function(x) (sum(x) >= min_obs) & (sum(!x) >= min_obs))
  clean_mtx <- merge_mtx[which_rows, ]
  conditions <- unique(md$env)
  if (verbose) message(paste0(length(conditions), " different conditions detected: ", paste0(conditions, collapse=", ")))
  if (length(conditions) != 2) { stop("Currently PanSweep only works when there are two conditions") }
  subjects_per_condition <- lapply(conditions, \(this_cond) {
    intersect(colnames(clean_mtx), md$subject[md$env==this_cond])
  })
  for (spc in subjects_per_condition) {
    if (length(spc)==0) { stop("Some conditions have no subjects represented in the data matrix")}
  }
  if (verbose) message("Generating contingencies...")
  contingency_rows <- t(apply(clean_mtx, 1, \(x) {
    # "sum" here gives "how many TRUEs" ## apply analysis across gene for species
    #     contingency_tbl <- matrix(nr=2, byrow=TRUE, data=c(sum(x[samples_per_condition[[1]]]),
    c(sum(x[subjects_per_condition[[1]]]),
      sum(x[subjects_per_condition[[2]]]),
      sum(!x[subjects_per_condition[[1]]]),
      sum(!x[subjects_per_condition[[2]]]))
  }))
  rownames(contingency_rows) <- rownames(clean_mtx)

  if (verbose) message("Performing Discrete FDR...")
  fisher_results <- DiscreteFDR::direct.discrete.BH(
    data.frame(t(contingency_rows)),
    "fisher",
    direction="sd"
  ) # we use the step-down method to report accurate adjusted p-values

  if (verbose) message("Returning results...")
  pvals <- fisher_results$Data$Raw.pvalues
  names(pvals) <- rownames(clean_mtx)
  fdrs <- fisher_results$Adjusted
  names(fdrs) <- rownames(clean_mtx)

  #fdrs <- p.adjust(pvals, 'BH')
  return(list(pvals=pvals,
              fdrs=fdrs,
              subjects_per_condition=subjects_per_condition,
              clean_mtx=clean_mtx))
}

#' Perform and analyze gene-to-species correlation lineage tests.
#'
#' This function correlates genes in a given species' pangenome to all species abundances, then tests whether the highest-correlating taxa are in from the same taxonomic family as the species whose pangenome was tested.
#'
#' @param Sig_Gene_reads A named list of matrices with one matrix per pangenome tested. Names should correspond to a species ID in meta_genome_sep_taxa. Rows of the matrices should be gene names, columns should be sample IDs, and values should be gene abundance.
#' @param Species_Abd A matrix with species IDs as rows and samples as columns; values are the abundance of that species in that column.
#' @param meta_genome_sep_taxa  Tibble containing taxonomic information for each species ID. Must contain the columns "species_id" and "Domain" through "Species". See `PanSweep::separate_taxonomy_with_s()` if you have lineage strings instead.
#' @param cor_fxn Function to actually correlate the two matrices. Defaults to `spearman_cor_wrapper()`.
#' @return Returns a table with a report on the results.
#'
#' @export
species_to_gene_correlations <- function(Sig_Gene_reads, Species_Abd, meta_genome_sep_taxa, cor_fxn = spearman_cor_wrapper) {

  meta_genome_sep_taxa <- mutate(meta_genome_sep_taxa,
                                 species_id = as.character(species_id))

  # Correlate with cor (much faster, returns matrix):
  pb <- progress::progress_bar$new(total = length(Sig_Gene_reads))
  Cor_Results <- purrr::map(Sig_Gene_reads, function(sgr) {
    pb$tick()
    n <- intersect(colnames(sgr), colnames(Species_Abd))
    spearman_cor_wrapper(sgr[, n, drop=FALSE],
                         Species_Abd[, n, drop=FALSE])
  })

  # Find max and the species of the pangenome from the results (added multiple maxes), and turn this into "long" data
  Cor_Results_max_targ <-  lapply(names(Cor_Results), function(s) {
    lapply(rownames(Cor_Results[[s]]), function(g){
      tbl <- enframe(Cor_Results[[s]][g, ], name="Species_Cor", value="Rho")
      target_family <- meta_genome_sep_taxa %>% dplyr::filter(species_id == s) %>% pull(Family)
      family_df <- left_join(
        tbl,
        meta_genome_sep_taxa %>% rename(Species_Cor_Name = Species),
        by = c("Species_Cor" = "species_id"))
      #Cannot do which.max here!!
      withCallingHandlers({
        Corr_Order <- family_df %>%
          mutate(Rho_n = as.numeric(Rho)) %>%
          filter(!is.na(Rho_n)) %>%
          arrange(desc(Rho_n)) %>%
          mutate(rank = row_number()) %>%
          mutate(mark = "") %>%
          mutate(f_match = (Family == target_family),
                 s_match = (Species_Cor == s))
        max_any <- Corr_Order %>%
          slice_max(Rho_n, n=1, with_ties = TRUE) %>%
          mutate(mark = "max") %>%
          arrange(rank)
        max_any_fam <- Corr_Order %>%
          filter(f_match) %>%
          slice_max(Rho_n, n=1, with_ties = FALSE) %>%
          mutate(mark = "max_f") %>%
          arrange(rank)
        s_target <- Corr_Order %>%
          filter(s_match) %>%
          mutate(mark = "target")
      }, warning = function(w){
        if(grepl("NAs introduced by coercion", conditionMessage(w))){
          invokeRestart("muffleWarning")
        }
      })
      df <- rbind(s_target, max_any, max_any_fam) %>%
        select(Species_Cor,
               Species_Cor_Name,
               Rho, mark, rank, f_match, s_match) %>%
        mutate(Gene = g, Species_id = s)
      return(df)
    }) %>% setNames(rownames(Cor_Results[[s]]))
  }) %>% setNames(names(Cor_Results))

  # Completely flatten (easier because we've already added species, gene)
  Species_Cor_DF <- Reduce(rbind, lapply(Cor_Results_max_targ,
                                         \(.) Reduce(rbind, .)))

  # Report key values per species and gene
  report <- Species_Cor_DF %>%
    dplyr::group_by(Species_id, Gene) %>%
    dplyr::summarize(
      Lineage_Shared = any(f_match[mark=="max"]),
      cor_max_species = paste0(Species_Cor_Name[mark=="max"], collapse=";"),
      cor_max_value = unique(Rho[mark=="max"])[1],
      cor_max_fam_value = max(as.numeric(Rho[mark=="max_f"])),
      cor_max_sp_value = as.numeric(Rho[mark=="target"]),
      .groups="keep") %>%
    dplyr::mutate(
      cor_ratio = cor_max_value / cor_max_fam_value
    )

  return(list(Cor_Sp_Ln_max_sp_DF=report,
              Species_Cor_DF=Species_Cor_DF))

}

#' Separate a column containing GTDB taxonomy strings to a standardized set of separate columns.
#'
#' @param inpt A table with a column that has GTDB taxonomy strings.
#' @param taxa_col The column to separate.
#' @return The table in `inpt` with new columns named "Domain", "Phylum", ..., "Species"
#'
#' @export
separate_taxonomy_with_s <- function(inpt, taxa_col){
  inpt <- inpt %>%
    tidyr::separate_wider_delim({{ taxa_col }},
                                names = c("d", "p", "c", "o", "f", "g", "s"),
                                delim = ";") %>%
    mutate(across(c("d", "p", "c", "o", "f", "g", "s"), ~ gsub("[dpcofgs]__","", .))) %>%
    rename_with(~ case_when(
      . == "d" ~ "Domain",
      . == "p" ~ "Phylum",
      . == "c" ~ "Class",
      . == "o" ~ "Order",
      . == "f" ~ "Family",
      . == "g" ~ "Genus",
      . == "s" ~ "Species",
      TRUE ~ .
    ))
  return(inpt)
}

#' Wrapper function for performing Spearman correlation between the genes and species matrices.
#'
#' This function is used by `species_to_gene_correlations()` and should not need to be called directly. It is made available mainly so that users can test that their own replacements accept and return data in the same format.
#'
#' @param genes_mtx Matrix of genes (rows) by samples (columns); values indicate gene abundance.
#' @param species_mtx Matrix of genes (rows) by samples (columns); values indicate species abundance.
#' @return Returns a matrix of Spearman's correlation values (possibly with some NAs) between genes (rows) and species (columns).
#' @export
spearman_cor_wrapper <- function(genes_mtx, species_mtx) {
  if (colnames(genes_mtx) != colnames(species_mtx)) {
    stop("Error: gene and species matrices need to have the same columns")
  }
  cor(t(genes),
      t(species),
      method = 'spearman')
}
