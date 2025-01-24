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
#'@import progress
#'@import DiscreteFDR
#'@import stringr
#'@import tidyr
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
                              verbose=FALSE
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
  save_folder_location <- normalizePath(Paths_and_Variables$Output$save_folder_location)
  if (!dir.exists(save_folder_location)) { dir.create(save_folder_location) }
  #Change variable name:
  Corr_lower_limit <- Co_occurrence_lower_limit

  ####################################################################################################################
  if (verbose) message("Reading in presence-absence data...")
  phylo_md2 <- read_tsv(path_phylo_md2)
  test_tbls <- purrr::map(Species_Set, ~ {
    gpath <- file.path(path_to_read_counts, .x, paste0(.x, ".genes_presabs.tsv"))

    if (is_compressed) {
      gpath <- paste0(gpath, ".lz4")
      tsv <- read_tsv_arrow(gpath)
    } else {
      tsv <- read_tsv(gpath, col_types=cols())
    }
    merge_columns_tbl(tsv, md=phylo_md2, fn=merge_fn_binary)
  })
  if (verbose) message("Calculating tests...")
  # run Fisher test analysis
  pb <- progress_bar$new(total = length(test_tbls))
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
  genome_metadata <- read_tsv(file= path_for_genomes_all_metadata)

  separate_taxonomy_with_s <- function(inpt, taxa_col){
    inpt <- inpt %>%
      tidyr::separate(taxa_col, into = c("d", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
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
    n_n <- seq.int(from = 2, to = ceiling(nrow(x) / 3))
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
  N.Sp_corr <- lapply(M.Sp_corr, function(x) metaMDS(x, distance = jaccard))
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
  Species_Abd <-read_tsv(path_to_species_abundance) %>%
    merge_columns_tbl(md=phylo_md2, fn=merge_fn_counts) %>%
    column_to_rownames("species_id") %>%
    as.matrix()

  #Correlate with try catch & cor.test:
  pb <- progress_bar$new(total = Reduce(sum, lapply(Sig_Gene_reads, nrow)) * nrow(Species_Abd))
  Cor_Results <-
    lapply(names(Sig_Gene_reads), function(sp){
      lapply(rownames(Sig_Gene_reads[[sp]]), function(g){
        lapply(rownames(Species_Abd), function(s){
          pb$tick()
          tryCatch({
            g_x <- Sig_Gene_reads[[sp]][g,]
            sp_y <- Species_Abd[s,colnames(Sig_Gene_reads[[sp]])]
            cor.test(g_x, sp_y[names(g_x)], method = "spearman", exact = FALSE) %>% .$estimate
          },
          warning = function(w){
            if (grepl("NaNs produced", w$message)){
              return("NaNs")
            } else if (grepl("the standard deviation is zero", w$message)){
              return("SD0")
            }
          }
          )
        }) %>% setNames(rownames(Species_Abd))
      }) %>% setNames(rownames(Sig_Gene_reads[[sp]]))
    }) %>% setNames(names(Sig_Gene_reads))



  # clean up the cor results
  Species_family_only <- meta_genome_sep_taxa %>% select("species_id", "Family")

  Cor_Results_df <-lapply(names(Cor_Results), function(s) {
    lapply(names(Cor_Results[[s]]), function(g){
      return(data.frame(Rho = Cor_Results[[s]][[g]] %>% unlist(use.names = FALSE), Species_Cor = names(Cor_Results[[s]][[g]]), stringsAsFactors = FALSE))
    }) %>% setNames(names(Cor_Results[[s]]))
  }) %>% setNames(names(Cor_Results))
  # Find max and the species of the pangenome from the results (added multiple maxes)
  Cor_Results_max_targ <-  lapply(names(Cor_Results_df), function(s) {
    lapply(names(Cor_Results_df[[s]]), function(g){
      target_family <- meta_genome_sep_taxa %>% dplyr::filter(species_id == s) %>% pull(Family)
      family_df <- left_join(Cor_Results_df[[s]][[g]],  Species_family_only, by = c("Species_Cor" = "species_id"))
                          #Cannot do which.max here!!
      withCallingHandlers({

        max_n <- max(as.numeric(Cor_Results_df[[s]][[g]][["Rho"]]), na.rm = TRUE)
        max_i <- which(Cor_Results_df[[s]][[g]][["Rho"]] == max_n)
        max <- Cor_Results_df[[s]][[g]][max_i,]
        max_c <- length(max_i)

        Corr_Order <- family_df %>%
          mutate(Rho_n = as.numeric(Rho)) %>%
          filter(!is.na(Rho_n)) %>%
          arrange(desc(Rho_n)) %>%
          mutate(rank = row_number())
        # print(colnames(Corr_Order))
        max_fam <- Corr_Order %>%
          filter(Family == target_family) %>%
          mutate(Rho_n = as.numeric(Rho)) %>%
          filter(!is.na(Rho_n)) %>%
          slice_max(Rho_n, n=1, with_ties = FALSE) %>%
          select(Species_Cor, rank, Rho)

        max_f <- max_fam %>%
          pull(Species_Cor)
        max_f_r <- max_fam %>%
          pull(rank)
      }, warning = function(w){
        if(grepl("NAs introduced by coercion", conditionMessage(w))){
          invokeRestart("muffleWarning")
        }
      })
      target <- Cor_Results_df[[s]][[g]][Cor_Results_df[[s]][[g]]$Species_Cor == s, ]
      target_r <- Corr_Order %>%
        filter(Species_Cor == s) %>%
        pull(rank)
      # print(max_f_r)
      df1 <- rbind(target, max_f, max)
      mark <- c("target", "max_f", rep("max", max_c))
      rank <- c(target_r, max_f_r, rep("max", max_c))
      df <- cbind(df1, mark, rank)
      return(df)
    }) %>% setNames(names(Cor_Results_df[[s]]))
  }) %>% setNames(names(Cor_Results_df))
  #Completely flatten
  Layer_1 <- Cor_Results_max_targ %>% tibble::enframe(name =  "Species", value = "Data")
  Layer_2 <- Layer_1 %>%
    mutate(Data = purrr::map(Data, ~ tibble::enframe(.x, name = "Gene", value = "Data2"))) %>%
    tidyr::unnest(Data)
  Species_Cor_DF <- Layer_2 %>% tidyr::unnest(cols = Data2)
  #______________________________________________________________________________#
  #Determine if lineage to Family is shared#
  #Extract species:
  Cor_Sp_label <- lapply(names(Cor_Results_max_targ),function(s){
    lapply(names(Cor_Results_max_targ[[s]]), function(g){
      m <- Cor_Results_max_targ[[s]][[g]] %>% dplyr::filter(mark == "max") %>% pull(Species_Cor)
      m_sp <- as.data.frame(meta_genome_sep_taxa %>% dplyr::filter(species_id == m) %>% select("Domain", "Phylum", "Class", "Order", "Family"))
      rownames(m_sp) <- paste(rep("max", nrow(m_sp)), 1:nrow(m_sp))
      t <-  Cor_Results_max_targ[[s]][[g]] %>% dplyr::filter(mark == "target") %>% pull(Species_Cor)
      t_sp <- as.data.frame(meta_genome_sep_taxa %>% dplyr::filter(species_id == t) %>% select("Domain", "Phylum", "Class", "Order", "Family"))
      rownames(t_sp) <- "target"
      df <- bind_rows(t_sp, m_sp)
      return(df)
    }) %>% setNames(names(Cor_Results_max_targ[[s]]))
  }) %>% setNames(names(Cor_Results_max_targ))
  #Determine if lineage is shared (TRUE) or not (FALSE):
  #All genes with multiple correlation maximums (or no max) return "Problem with correlation".
  Cor_Sp_agree <- lapply(names(Cor_Sp_label), function(s){
    lapply(names(Cor_Sp_label[[s]]), function(g){
      if (nrow(Cor_Sp_label[[s]][[g]]) == 2){
        max <- Cor_Sp_label[[s]][[g]] %>% t() %>% as.data.frame() %>% .$max
        target <- Cor_Sp_label[[s]][[g]] %>% t() %>% as.data.frame() %>% .$target
        all(max == target)
      } else{
        return("Problem with correlation")
      }
    }) %>% setNames(names(Cor_Sp_label[[s]]))
  }) %>% setNames(names(Cor_Sp_label))
  #Flatten to dataframe:
  Layer_1 <- Cor_Sp_agree %>% tibble::enframe(name =  "Species", value = "Data")
  Cor_Sp_agree_DF <- Layer_1 %>% mutate(Data = purrr::map(Data, ~ tibble::enframe(.x, name = "Gene", value = "Lineage_Shared"))) %>% tidyr::unnest(Data)
  Cor_Sp_agree_DF$Lineage_Shared <- Cor_Sp_agree_DF$Lineage_Shared %>% unlist()
  #______________________________________________________________________________#
  #Report species with max spearman correlation#
  sp_lable <- lapply(names(Cor_Results_max_targ),function(s){
    lapply(names(Cor_Results_max_targ[[s]]), function(g){
      m <- Cor_Results_max_targ[[s]][[g]] %>% dplyr::filter(mark == "max") %>% pull(Species_Cor)
      m_sp <- as.data.frame(meta_genome_sep_taxa %>% dplyr::filter(species_id == m) %>% select("Species"))
      if(m_sp == ''){
        m_sp <- as.data.frame(meta_genome_sep_taxa %>% dplyr::filter(species_id == m) %>% select("Genus"))
      }
      return(as.data.frame(m_sp))
    }) %>% setNames(names(Cor_Results_max_targ[[s]]))
  }) %>% setNames(names(Cor_Results_max_targ))
  #Need to keep information when unnesting:
  Layer_1 <- sp_lable  %>% tibble::enframe(name =  "Species_OG", value = "Data")
  Layer_2 <- Layer_1 %>% mutate(Data = purrr::map(Data, ~ tibble::enframe(.x, name = "Gene", value = "Data2"))) %>% tidyr::unnest(Data)
  sp_report_DF <- Layer_2 %>% tidyr::unnest(cols = Data2)
  #Transfer post DF being made:
  sp_report_DF <- sp_report_DF %>%  mutate(Species = ifelse(is.na(Species), Genus, Species)) %>%
    select("Species_OG", "Gene", "Species") %>% rename(cor_max_species = Species)
  # Left_join Lineage agree:
  Cor_Sp_Ln_max_sp_DF <- Cor_Sp_agree_DF %>% left_join(sp_report_DF, by = "Gene") %>%
    select("Species", "Gene", "Lineage_Shared", "cor_max_species")
  #______________________________________________________________________________#
  #Add lineage and max species correlation to eggNOG table#
  uhgp_90_eggNOG <- left_join(uhgp_90_eggNOG, Cor_Sp_Ln_max_sp_DF %>% select("Gene", "Lineage_Shared", "cor_max_species"), by = c("Gene_id" = "Gene"))
  uhgp_90_eggNOG <- left_join(uhgp_90_eggNOG, Species_Cor_DF %>% filter(mark == "max_f") %>%
                                select(Gene, rank), by = c("Gene_id" = "Gene")) %>% rename(Family_max_rank = rank)
  uhgp_90_eggNOG <- left_join(uhgp_90_eggNOG, Species_Cor_DF %>% filter(mark == "target") %>%
                                select(Gene, rank), by = c("Gene_id" = "Gene")) %>% rename(Sp_rank = rank)
  #______________________________________________________________________________#
  #Create Report#
  if (verbose) message("Creating report and outputting results...")
  Analysis_report <- tibble::as_tibble(cbind(c("Number of Significant Genes", "Number of Species with Significant Genes", "Number of repeated UHGP-90 ids", "Number of repeated UHGP-50 ids"),
                                     c(Number_of_Significant_Genes, nrow(Num_Sig_Genes_per_sp), Num_rep_UHGP_90_clus_id, Num_rep_UHGP_50_clus_id)))
  Analysis_output_names <- c('Analysis_report','Number_of_Significant_Genes', 'UHGP_90_cluster_id_summ', 'UHGP_50_cluster_id_summ', 'Num_Sig_Genes_per_sp', 'uhgp_90_eggNOG')
  Analysis_output <- setNames(mget(Analysis_output_names), Analysis_output_names)

  ################################################################################
  #Create save files#
  All_RDS_to_Save <- c("M.Sp_corr", "U.Sp_corr", "N.Sp_corr", "N.Stress", "P.Sp_corr", "Analysis_output")
  save_folder_name <- paste0("PanSweep_Analysis_Output_", format(Sys.Date(), "%Y-%m-%d"))
  save_folder_location_full <- file.path(save_folder_location, save_folder_name)
  dir.create(save_folder_location_full)
  MIDAS_Analysis_Output <- setNames(mget(All_RDS_to_Save), All_RDS_to_Save)
  saveRDS(MIDAS_Analysis_Output, file.path(save_folder_location_full, paste0("PanSweep_Analysis_Output", ".rds")))
  saveRDS(Analysis_output, file.path(save_folder_location_full, paste0("PanSweep_Analysis_TablesOnly", ".rds")))
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
#' @param verbose Boolean. If TRUE, print messages to indicate where we are in the process.
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


