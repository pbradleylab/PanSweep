#'PanSweep Analysis:
#'
#'This function run the analysis of MIDAS2 data and outputs a save file to be
#'used with the PanSweep_Shiny() function. The inputs are Json_Config_Path and
#'Co_occurrence_lower_limit, which is set to 3. Refer to the PanSweep github for 
#'the Json file and explanations.
#'
#' To run the function:
#' PanSweep_Analysis(Json_Config_Path = "Path/To/Json/file.json")
#'
#'To change the lower bound for number of genes needed for co-occurence 
#'  analysis: 
#' PanSweep_Analysis(Json_Config_Path = "Path/To/Json/file.json", 
#'                    Co_occurrence_lower_limit = 5)
#'
#' @param Json_Config_Path Path to the completed json config file
#' @param Co_occurrence_lower_limit  lower bound for number of genes needed for 
#' co-occurrence analysis. Set to 3.
#' @return Returned a date stamped folder called PanSweep_Analysis_Output_YYYY-MM-DD
#' containing the file "PanSweep_Analysis_Output.rds".            
#'
#'@import tidyverse
#'@import progress
#'@import readr
#'@import viridisLite
#'@import arrow
#'@import pillar
#'@import dbplyr
#'@import parallelDist
#'@import umap
#'@import vegan
#'@import jsonlite
#'@export
PanSweep_Analysis <- function(Json_Config_Path, Co_occurrence_lower_limit = 3){
  ####################################################################################################################
  #Load in Paths and variables via JSON#
  Paths_and_Varriables <-fromJSON(Json_Config_Path)
  #Variables:
  Species_Set <- Paths_and_Varriables$Varriables$Species_Set
  #Databases:
  path_uhgp_50_cluster <- Paths_and_Varriables$Databases$path_uhgp_50_cluster
  path_uhgp_90_cluster <- Paths_and_Varriables$Databases$path_uhgp_90_cluster
  path_uhgp_50_eggNOG <- Paths_and_Varriables$Databases$path_uhgp_50_eggNOG
  path_uhgp_90_eggNOG <- Paths_and_Varriables$Databases$path_uhgp_90_eggNOG
  #Metadata:
  path_for_genomes_all_metadata <- Paths_and_Varriables$Metadata$path_for_genomes_all_metadata
  path_phylo_md2 <- Paths_and_Varriables$Metadata$path_to_sample_metadata
  path_to_presabs <- Paths_and_Varriables$Metadata$path_to_presabs
  path_to_species_abundance <- Paths_and_Varriables$Metadata$path_to_species_abundance
  path_to_read_counts <- Paths_and_Varriables$Metadata$path_to_read_counts
  #Output:
  save_folder_location <- Paths_and_Varriables$Output$save_folder_location
  
  #Change variable name:
  Corr_lower_limit <- Co_occurrence_lower_limit
  
  ####################################################################################################################
  phylo_md2 <- read_tsv(path_phylo_md2)
  test_tbls <- map(Species_Set, ~ read_tsv(paste0(path_to_presabs,
                                                  .x,
                                                  "/",
                                                  .x,
                                                  ".genes_presabs.tsv"),
                                           col_types=cols()))
  
  analyze_tbl <- function(tbl, md=phylo_md2, min_obs = 3) {
    # convert tibble to matrix with rownames
    mtx <- as.matrix(tbl[,-1]) #-1 = do not grab first column into matrix
    rownames(mtx) <- tbl$gene_id #name the first row based off of the tibble
    # filter out genes that don't have at least 3 observations (and 3 non-observations) => be there three times and not there three times
    which_rows <- apply(mtx, 1, function(x) (sum(x) >= min_obs) & (sum(!x) >= min_obs))
    clean_mtx <- mtx[which_rows, ]
    ctrl <- intersect(colnames(clean_mtx), md$sample[md$env=="control"]) #just the samples where the env is control
    case <- intersect(colnames(clean_mtx), md$sample[md$env!="control"]) #just the samples where the env is not control
    pvals <- apply(clean_mtx, 1, function(x) {
      # "sum" here gives "how many TRUEs" ## apply analysis across gene for species
      contingency_tbl <- matrix(nr=2, byrow=TRUE, data=c(sum(x[ctrl]),
                                                         sum(x[case]),
                                                         sum(!x[ctrl]),
                                                         sum(!x[case])))
      # looking at distribution could be defined by random chance
      fisher.test(contingency_tbl)$p.value
    })
    fdrs <- p.adjust(pvals, 'BH')
    return(list(pvals=pvals, fdrs=fdrs, ctrl=ctrl, case=case, clean_mtx=clean_mtx))
  }
  
  # run analysis
  pb <- progress_estimated(length(test_tbls))
  test_results <- purrr::map(test_tbls, ~ { pb$tick()$print(); analyze_tbl(.x) })
  names(test_results) <- Species_Set
  
  # ...next:
  # extract for each fdrs in each member of test_results, names of genes with fdr <= 0.05
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Extracting genes with fdrs<=0.05#
  test_results_fdrs <- map(test_results, function(x) enframe(x$fdrs, name = "Gene_id", value = "Fdrs"))
  #enframe is used in named vector
  #creates a list under species_id and removes other unneeded lists (ie pvals)  
  test_results_lt_tbl <-enframe(test_results_fdrs, name = "Species_id")
  #creates a list of tibbles under species
  test_results_tbl <- unnest(test_results_lt_tbl, cols = c(value))
  #makes a tbl of the test results ALL fdrs are included
  
  Gene_extract_tbl <- filter(test_results_tbl, Fdrs<=0.05)
  #Extract all genes by (with species) that are Fdrs<=0.05
  Genes_of_intrest <- unique(Gene_extract_tbl$Gene_id) 
  
  
  Genes_intrest_extr <- filter(test_results_tbl, Gene_id %in% Genes_of_intrest)
  # Extract all genes from the test results tbl (ALL fdrs) that are in the "Genes_of_interest" vector.
  
  #start list for report:
  Number_of_Significant_Genes <- nrow(Gene_extract_tbl)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Adding cluster ids#
  Gut_Trans <- gsub("UHGG","", Gene_extract_tbl$Gene_id) %>%
    gsub("_", "", .)
  
  prqt_uhgp_90 <- arrow::open_dataset(path_uhgp_90_cluster, format = "parquet")
  prqt_uhgp_50 <- arrow::open_dataset(path_uhgp_50_cluster, format = "parquet")
  
  Gene_intr_uhgp_90 <- prqt_uhgp_90 %>%
    filter(gene_id %in% Gut_Trans) %>%
    collect() %>%
    dplyr::select(-c(group)) %>%
    rename(cluster_id_n = cluster_id)
  
  Gene_intr_uhgp_90 <- Gene_intr_uhgp_90 %>%
    map_dfc(function(x) as.character(x)) %>%
    map_dfc(function(x) str_pad(x, 11, "left", pad = "0"))
  
  Gene_intr_uhgp_50 <- prqt_uhgp_50 %>%
    filter(gene_id %in% Gut_Trans) %>%
    collect() %>%
    dplyr::select(-c(group)) %>%
    rename(cluster_id_n = cluster_id)
  
  Gene_intr_uhgp_50 <- Gene_intr_uhgp_50 %>%
    map_dfc(function(x) as.character(x)) %>%
    map_dfc(function(x) str_pad(x, 11, "left", pad = "0"))
  
  Genes_intrest_extr<- mutate(Genes_intrest_extr, "gene_id_n"=gsub("UHGG","", Gene_extract_tbl$Gene_id) %>%
                                gsub("_", "", .))
  UHGP_50_genes_of_int <- left_join(Genes_intrest_extr, Gene_intr_uhgp_50, by=c("gene_id_n" = "gene_id"))
  UHGP_90_genes_of_int <- left_join(Genes_intrest_extr, Gene_intr_uhgp_90, by=c("gene_id_n" = "gene_id"))
  #______________________________________________________________________________#
  #Return Gut_Genome Syntax#
  
  UHGP_90_genes_of_int <- mutate(UHGP_90_genes_of_int, "cluster_id" = 
                                   gsub('^(.{6})(.*)$','GUT_GENOME\\1_\\2', UHGP_90_genes_of_int$cluster_id_n))
  
  UHGP_50_genes_of_int <- mutate(UHGP_50_genes_of_int, "cluster_id" = 
                                   gsub('^(.{6})(.*)$','GUT_GENOME\\1_\\2', UHGP_50_genes_of_int$cluster_id_n))
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
  #Assign eggNOG Taxonomy and gene info#
  prqt_uhgp_50_eggNOG <- arrow::open_dataset(path_uhgp_50_eggNOG,format = "parquet")
  prqt_uhgp_90_eggNOG <- arrow::open_dataset(path_uhgp_90_eggNOG, format = "parquet")
  
  uhgp_90_eggNOG <- prqt_uhgp_90_eggNOG %>%
    filter(query_name %in% UHGP_90_genes_of_int$cluster_id) %>%
    select(query_name, Predicted_taxonomic_group, Predicted_protein_name, eggNOG_free_text_description) %>%
    collect() %>%
    right_join(UHGP_90_genes_of_int, by = c("query_name" = "cluster_id"), keep = TRUE)
  
  uhgp_50_eggNOG <- prqt_uhgp_50_eggNOG %>%
    filter(query_name %in% UHGP_50_genes_of_int$cluster_id) %>%
    select(query_name, Predicted_taxonomic_group, Predicted_protein_name, eggNOG_free_text_description) %>%
    collect() %>%
    right_join(UHGP_90_genes_of_int, by = c("query_name" = "cluster_id"), keep = TRUE)
  #______________________________________________________________________________#
  #Add in Taxonomy info#
  genome_metadata <- read_tsv(file= path_for_genomes_all_metadata)
  
  separate_taxonomy_with_s <- function(inpt, taxa_col){
    inpt <- inpt %>%
      separate(taxa_col, into = c("d", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
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
  Runs_by_Genes_intrest<- lapply(test_tbls_sp, function(x)  subset(x,x$gene_id %in% Genes_of_intrest))
  Sp_to_corr_Runs <- keep(Runs_by_Genes_intrest, function(x) nrow(x)>Corr_lower_limit) %>%
    lapply(., function(x) column_to_rownames(x, var = 'gene_id'))
  
  ################################################################################
  #Shiny Prep#
  #______________________________________________________________________________#
  #Correlate by runs#
  Sp_corr <- lapply(Sp_to_corr_Runs, function(x) parDist(as.matrix(x), method = 'binary'))
  M.Sp_corr <- lapply(Sp_corr, function(x) as.matrix(x))
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
  Species_Abd <-read_tsv(path_to_species_abundance) %>% 
    pivot_longer(!species_id, names_to = "run", values_to = "species_count") %>% 
    pivot_wider(names_from = "species_id", values_from = "species_count") %>% column_to_rownames("run")
  #Determine significant species:
  species_set_corr <- unique(Genes_intrest_extr$Species_id)
  #Read in gene counts:
  Gene_reads <- map(species_set_corr, ~ read_tsv(paste0(path_to_read_counts,
                                                        .x,
                                                        "/",
                                                        .x,
                                                        ".genes_reads.tsv"),
                                                 col_types=cols())) %>% setNames(species_set_corr)
  #Extract reads for significant genes:
  Sig_Gene_copyNum <- map(species_set_corr, ~dplyr::filter(Gene_reads[[.x]], gene_id %in% Genes_intrest_extr$Gene_id) %>% pivot_longer(!gene_id, names_to = "run", values_to = "Gene_count") %>% pivot_wider(names_from = "gene_id", values_from = "Gene_count") %>% column_to_rownames("run")) %>% setNames(species_set_corr)
  #Correlate with try catch & cor.test:
  Cor_Results <- 
    lapply(names(Sig_Gene_copyNum), function(sp){
      lapply(names(Sig_Gene_copyNum[[sp]]), function(g){
        lapply(colnames(Species_Abd), function(s){
          tryCatch({
            g_x <- as.numeric(Sig_Gene_copyNum[[sp]][[g]])
            sp_y <- as.numeric(Species_Abd %>% filter(rownames(.) %in% rownames(Sig_Gene_copyNum[[sp]])) %>% pull(s))
            cor.test(g_x, sp_y, method = "spearman", exact = FALSE) %>% .$estimate
          },
          warning = function(w){
            if (grepl("NaNs produced", w$message)){
              return("NaNs")
            } else if (grepl("the standard deviation is zero", w$message)){
              return("SD0")
            }
          }
          )
        }) %>% setNames(colnames(Species_Abd))
      }) %>% setNames(names(Sig_Gene_copyNum[[sp]]))
    }) %>% setNames(names(Sig_Gene_copyNum))
  # clean up the cor results
  Cor_Results_df <-lapply(names(Cor_Results), function(s) {
    lapply(names(Cor_Results[[s]]), function(g){
      return(data.frame(Rho = Cor_Results[[s]][[g]] %>% unlist(use.names = FALSE), Species_Cor = names(Cor_Results[[s]][[g]]), stringsAsFactors = FALSE))
    }) %>% setNames(names(Cor_Results[[s]]))
  }) %>% setNames(names(Cor_Results))
  # Find max and the species of the pangenome from the results (added multiple maxes)
  Cor_Results_max_targ <-  lapply(names(Cor_Results_df), function(s) {
    lapply(names(Cor_Results_df[[s]]), function(g){
      withCallingHandlers({
        max_n <- max(as.numeric(Cor_Results_df[[s]][[g]][["Rho"]]), na.rm = TRUE)
        max_i <- which(Cor_Results_df[[s]][[g]][["Rho"]] == max_n)
        max <- Cor_Results_df[[s]][[g]][max_i,]
        max_c <- length(max_i)
      }, warning = function(w){
        if(grepl("NAs introduced by coercion", conditionMessage(w))){
          invokeRestart("muffleWarning")
        }
      })
      
      target <- Cor_Results_df[[s]][[g]][Cor_Results_df[[s]][[g]]$Species_Cor == s, ]
      df1 <- rbind(target, max)
      mark <- c("target", rep("max", max_c))
      df <- cbind(df1, mark)
      return(df)
    }) %>% setNames(names(Cor_Results_df[[s]]))
  }) %>% setNames(names(Cor_Results_df))
  #Completely flatten
  Layer_1 <- Cor_Results_max_targ %>% enframe(name =  "Species", value = "Data")
  Layer_2 <- Layer_1 %>% mutate(Data = map(Data, ~enframe(.x, name = "Gene", value = "Data2"))) %>% unnest(Data)
  Species_Cor_DF <- Layer_2 %>% unnest(cols = Data2)
  #______________________________________________________________________________#
  #Determine if lineage to Family is shared#
  #Extract species:
  Cor_Sp_label <- lapply(names(Cor_Results_max_targ),function(s){
    lapply(names(Cor_Results_max_targ[[s]]), function(g){
      m <- Cor_Results_max_targ[[s]][[g]] %>% filter(mark == "max") %>% pull(Species_Cor)
      m_sp <- as.data.frame(meta_genome_sep_taxa %>% filter(species_id == m) %>% select("Domain", "Phylum", "Class", "Order", "Family"))
      rownames(m_sp) <- paste(rep("max", nrow(m_sp)), 1:nrow(m_sp))
      t <-  Cor_Results_max_targ[[s]][[g]] %>% filter(mark == "target") %>% pull(Species_Cor)
      t_sp <- as.data.frame(meta_genome_sep_taxa %>% filter(species_id == t) %>% select("Domain", "Phylum", "Class", "Order", "Family"))
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
  Layer_1 <- Cor_Sp_agree %>% enframe(name =  "Species", value = "Data")
  Cor_Sp_agree_DF <- Layer_1 %>% mutate(Data = map(Data, ~enframe(.x, name = "Gene", value = "Lineage_Shared"))) %>% unnest(Data)
  Cor_Sp_agree_DF$Lineage_Shared <- Cor_Sp_agree_DF$Lineage_Shared %>% unlist()
  #______________________________________________________________________________#
  #Report species with max spearman correlation#
  sp_lable <- lapply(names(Cor_Results_max_targ),function(s){
    lapply(names(Cor_Results_max_targ[[s]]), function(g){
      m <- Cor_Results_max_targ[[s]][[g]] %>% filter(mark == "max") %>% pull(Species_Cor)
      m_sp <- as.data.frame(meta_genome_sep_taxa %>% filter(species_id == m) %>% select("Species"))
      if(m_sp == ''){
        m_sp <- as.data.frame(meta_genome_sep_taxa %>% filter(species_id == m) %>% select("Genus"))
      }
      return(as.data.frame(m_sp))
    }) %>% setNames(names(Cor_Results_max_targ[[s]]))
  }) %>% setNames(names(Cor_Results_max_targ))
  #Need to keep information when unnesting:
  Layer_1 <- sp_lable  %>% enframe(name =  "Species_OG", value = "Data")
  Layer_2 <- Layer_1 %>% mutate(Data = map(Data, ~enframe(.x, name = "Gene", value = "Data2"))) %>% unnest(Data)
  sp_report_DF <- Layer_2 %>% unnest(cols = Data2)
  #Transfer post DF being made:
  sp_report_DF <- sp_report_DF %>%  mutate(Species = ifelse(is.na(Species), Genus, Species)) %>% 
    select("Species_OG", "Gene", "Species") %>% rename(cor_max_species = Species)
  # Left_join Lineage agree:
  Cor_Sp_Ln_max_sp_DF <- Cor_Sp_agree_DF %>% left_join(sp_report_DF, by = "Gene") %>% 
    select("Species", "Gene", "Lineage_Shared", "cor_max_species")
  #______________________________________________________________________________#
  #Add lineage and max species correlation to eggNOG table#
  uhgp_90_eggNOG <- left_join(uhgp_90_eggNOG, Cor_Sp_Ln_max_sp_DF %>% select("Gene", "Lineage_Shared", "cor_max_species"), by = c("Gene_id" = "Gene"))
  #______________________________________________________________________________#
  #Create Report#
  Analysis_report <- as_tibble(cbind(c("Number of Significant Genes", "Number of Species with Significant Genes", "Number of repeated UHGP-90 ids", "Number of repeated UHGP-50 ids"), 
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
  ################################################################################
}