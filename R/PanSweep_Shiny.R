#'PanSweep Shiny Application:
#'
#'This is the function to run the shiny analysis of the output from PanSweep_Analysis()
#'function. The file to be imputed is the "PanSweep_Analysis_Output.rds" located
#'in the PanSweep_Analysis_Output_YYYY-MM-DD folder created by the PanSweep_Analysis()
#'function.
#'
#'To run the application use the below function:
#'
#'PanSweep_Shiny(loadData_Path = "Path/to/PanSweep_Analysis_Output_YYYY-MM-DD/PanSweep_Analysis_Output.rds")
#'
#'Note: The PanSweep_Analysis_Output.rds is not required to be in the folder
#'called PanSweep_Analysis_Output_YYYY-MM-DD for the application to be run
#'
#'Once the function is called a shiny application will run and the R console cannot
#'be used until the application is closed.
#'
#'@param loadData_Path Path to the PanSweep_Analysis_Output.rds file from the
#'PanSweep_Analysis function.
#'@return This will not return any files but run a shiny UI for visual analysis
#'of results.
#'
#'@import shiny
#'@import umap
#'@import vegan
#'@import dplyr
#'@import purrr
#'@import readr
#'@import tidyr
#'@import stringr
#'@import plotly
#'@import DT
#'@import knitr
#'@import kableExtra
#'
#'@export
PanSweep_Shiny <- function(loadData_Path){
  loadData <- read_rds(loadData_Path)
  if(interactive()){

    ui <- navbarPage( "PanSweep",

                      tabPanel("Analysis Report",
                               sidebarPanel(selectInput(inputId = "report",
                                                        label = "Choose Report:",
                                                        choices = c("Overall Report", "Species Report", "UHGP-90 Repeat ids", "UHGP-50 Repeat ids"))
                               ),
                               mainPanel(DTOutput("AnaR")),
                      ),

                      tabPanel("eggNOG & Correlation Report",
                               mainPanel(DTOutput("eNR"))
                      ),
                      tabPanel("Ordination & Heatmap",
                               fluidRow(column(2,
                                               selectInput(inputId = "ord_plt",
                                                           label = "Choose analysis:",
                                                           choices = c("UMAP","NMDS","PCoA")),
                                               sliderInput("n_n",
                                                           label = "Number of n_neighbors:",
                                                           min = 2,
                                                           max = 10,
                                                           value = 2),
                                               sliderInput("min_dist",
                                                           label = "min_dist:",
                                                           min = 0.1,
                                                           max = 0.9,
                                                           value = 0.1,
                                                           step = 0.1),
                                               selectInput(inputId = "species_c",
                                                           label = "Choose species:",
                                                           choices = names(loadData$N.Sp_corr)),
                                               textOutput("Ord_Species"),
                                               actionButton("reset", "Reset")
                                               

                               ),
                               column(10, plotlyOutput("ordination_plot"),
                                      plotlyOutput("corrMax2"),
                                      tableOutput("mtData")
                               ),


                               )
                      ),
                      tabPanel("NMDS",
                               sidebarPanel(
                                 selectInput(inputId = "species_c2",
                                             label = "Choose species:",
                                             choices = names(loadData$N.Sp_corr))
                               ),
                               mainPanel(plotlyOutput("NMDS"),
                                         plotOutput("stressPlot"))
                      )
    )

    server <- function(input, output, session) {

      clickData <- reactiveValues(x = vector(), y = vector(), cd = vector())
      #Builds up click data:
      observeEvent(event_data("plotly_click"), {
        new_clickData <- event_data("plotly_click")
        if (!is.null(new_clickData)){
          clickData$x <- c(clickData$x, new_clickData$x)
          clickData$y <- c(clickData$y, new_clickData$y)
          clickData$cd <- c(clickData$cd, new_clickData$customdata)
        }
      })
      #Resets click data:
      observeEvent(input$reset,{
        clickData$x <- vector()
        clickData$y <- vector()
        clickData$cd <- vector()
      })


      observeEvent(input$species_c, {
        req(loadData)
        possible_n <- names(loadData$U.Sp_corr[[paste0(input$species_c, sep='')]])
        n_stepsize <- 1
        if (!is.null(possible_n)) {
          n_num <- as.numeric(possible_n)
          if (length(n_num) > 1) { n_stepsize <- n_num[2] - n_num[1] }
        } else {
          n_num <- 2
        }
        updateSliderInput(session = session,
                          "n_n",
                          value = 2,
                          min = min(n_num),
                          max = max(n_num),
                          step = n_stepsize
                          )
      }
      )

      output$ordination_plot <-  renderPlotly({
        if(input$ord_plt == "UMAP"){
          n_n = input$n_n
          n_num <- as.numeric(names(loadData$U.Sp_corr[[paste0(input$species_c, sep='')]]))
          names(n_num) <- as.character(n_num)
          n_closest <- names(sort(abs(n_num-n_n)))[1]
          p <- loadData$U.Sp_corr %>%
            .[[paste0(input$species_c, sep='')]] %>%
            .[[n_closest]] %>%
            .[[paste0(input$min_dist, sep='')]] %>%
            .$"layout" %>%
            as.data.frame()%>%
            rownames_to_column() %>%
            rename("Gene_id" = "rowname") %>%
            left_join(., loadData$Analysis_output$uhgp_90_eggNOG, "Gene_id") %>%
            plot_ly(., x = ~V1, y = ~V2, type = 'scatter', mode = 'markers',
                    color= ~replace(.$Predicted_taxonomic_group, is.na(.$Predicted_taxonomic_group), "NA"),
                    text = paste(.$Gene_id),
                    customdata = ~paste(.$Gene_id)) %>%
            layout(showlegend = TRUE,
                   plot_bgcolor = "#e5ecf6",
                   xaxis = list(
                     title = "0"),
                   yaxis = list(
                     title = "1"),
                   annotations = list(text = ~paste("n_neighbors:", input$n_n, "min_dist", input$min_dist), showarrow=FALSE ), ##MAKE PRERDY##
                   ggplot2::theme(plot.title.position = ggplot2::element_text(vjust = 0.5))
            )
          event_register(p, 'plotly_click')


        }
        else if(input$ord_plt == "NMDS"){
          loadData$N.Sp_corr %>%
            .[[paste0(input$species_c, sep='')]] %>%
            scores() %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            rename("Gene_id" = "rowname") %>%
            left_join(., loadData$Analysis_output$uhgp_90_eggNOG, "Gene_id") %>%
            plot_ly(x = ~NMDS1, y = ~NMDS2, type = "scatter", mode = "markers",  color= ~replace(.$Predicted_taxonomic_group, is.na(.$Predicted_taxonomic_group), "NA"),
                    text = paste(.$Gene_id), customdata = ~paste(.$Gene_id), showlegend = TRUE, legendgroup = "markers")%>%
            #add_text(textfont = list(color = "black"), textposition = "top right", showlegend = FALSE) %>%
            layout(plot_bgcolor = "#e5ecf6",
                   annotations = list(text = ~paste("stress", loadData$N.Stress[[paste0(input$species_c2, sep='')]]), showarrow=FALSE))
        }
        else if(input$ord_plt == "PCoA"){
          loadData$P.Sp_corr %>%
            .[[paste0(input$species_c, sep='')]] %>%
            .$points %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            rename("Gene_id" = "rowname") %>%
            left_join(., loadData$Analysis_output$uhgp_90_eggNOG, "Gene_id") %>%
            plot_ly(x = ~V1, y = ~V2, type = 'scatter', mode = 'markers',
                    color= ~replace(.$Predicted_taxonomic_group, is.na(.$Predicted_taxonomic_group), "NA"),
                    text = paste(.$Gene_id),
                    customdata = ~paste(.$Gene_id)) %>%
            layout(showlegend = TRUE, plot_bgcolor = "#e5ecf6")
        }
      })

      output$eNR <- renderDT({
        loadData$Analysis_output$uhgp_90_eggNOG %>%
          select("Gene_id", "Species_id", "Species","Lineage_Shared", "cor_max_species", "Sp_rank", "Family_max_rank", "Fdrs", "cluster_id", "Predicted_taxonomic_group",
                 "Predicted_protein_name", "eggNOG_free_text_description") %>%
          mutate(Fdrs = format(signif(Fdrs, 3), scientific = TRUE))%>%
          arrange(Species_id, desc(Lineage_Shared), desc(Fdrs)) %>%
          rename("Gene ID" = "Gene_id", "Species ID" = "Species_id", "Lineage Shared" = "Lineage_Shared",
                 "Species with Max Correlation Value" = "cor_max_species", "Species Correlation Rank" = "Sp_rank", "Family Correlation Rank" = "Family_max_rank",
                 "Cluster ID" = "cluster_id","EggNOG Predicted taxonomic group" = "Predicted_taxonomic_group",
                 "Predicted protein name" = "Predicted_protein_name",
                 "eggNOG free text description" = "eggNOG_free_text_description")
      })

      output$AnaR <- renderDT({
        if (input$report == "Overall Report"){
          loadData$Analysis_output$Analysis_report %>% datatable(colnames = NULL)
        } else if (input$report == "UHGP-90 Repeat ids"){
          loadData$Analysis_output$UHGP_90_cluster_id_summ
        } else if (input$report == "UHGP-50 Repeat ids"){
          loadData$Analysis_output$UHGP_50_cluster_id_summ
        } else if (input$report == "Species Report"){
          loadData$Analysis_output$Num_Sig_Genes_per_sp %>%
            select("Species_id", "Species", "n") %>%
            rename("Species Id" = "Species_id") %>%
            rename("Number of Significant Genes in Species" = "n")
        }
      })
      
      output$Ord_Species <- renderText({
        species <- loadData$Analysis_output$Num_Sig_Genes_per_sp %>% filter(Species_id == input$species_c2) %>% pull(Species)
      })


      output$NMDS <- renderPlotly({
        loadData$N.Sp_corr %>%
          .[[paste0(input$species_c2, sep='')]] %>%
          scores() %>%
          as.data.frame() %>%
          rownames_to_column() %>%
          rename("Gene_id" = "rowname") %>%
          left_join(., loadData$Analysis_output$uhgp_90_eggNOG, "Gene_id") %>%
          plot_ly(x = ~NMDS1, y = ~NMDS2, type = "scatter", mode = "markers",  color= ~replace(.$Predicted_taxonomic_group, is.na(.$Predicted_taxonomic_group), "NA"),
                  text = paste(.$Gene_id))%>%
          layout(showlegend = TRUE,
                 annotations = list(text = ~paste("stress", loadData$N.Stress[[paste0(input$species_c2, sep='')]]), showarrow=FALSE))
      })

      output$stressPlot <- renderPlot({
        loadData$N.Sp_corr %>%
          .[[paste0(input$species_c2, sep='')]] %>%
          stressplot()
      })

      output$corrMax2 <- renderPlotly({
        corMax <- loadData$M.Sp_corr%>%
          .[[paste0(input$species_c, sep='')]] %>%
          as.data.frame() %>%
          mutate(across(everything(), ~  1 - .)) %>%
          as.matrix()

        if (length(clickData$cd > 0)){

          highlight_row <- isolate(clickData$cd)

          heatmap_plot <- plot_ly(x = rownames(corMax), y = colnames(corMax),
                                  z = corMax, zmin = 0, zmax = 1,
                                  type = "heatmap",
                                  colorscale = "Greys",
                                  colorbar = list(title = "Greys", y = 0.45, len = 0.45)
          )

          HiLite_mtx <-matrix(NA, nrow = nrow(corMax), ncol = ncol(corMax))
          rownames(HiLite_mtx) <- rownames(corMax)
          colnames(HiLite_mtx) <- colnames(corMax)
          HiLite_mtx[,paste0(highlight_row, sep='')] <- corMax[,paste0(highlight_row, sep='')]
          HiLite_mtx[paste0(highlight_row, sep=''),] <- corMax[paste0(highlight_row, sep=''),]

          heatmap_plot <- heatmap_plot %>%
            add_trace(
              z = HiLite_mtx,
              type = "heatmap",
              colorscale = "Viridis",
              colorbar = list(title = "Viridis", y = 1, len = 0.45)
            ) %>%
            layout(
              title = "Jaccard Similarity"
            )
        }
        else{
          heatmap_plot <-
            plot_ly(x = rownames(corMax), y = colnames(corMax),
                    z = corMax, zmin = 0, zmax = 1,
                    type = "heatmap") %>%
            layout(
              title = "Jaccard Similarity"
            )
        }
      })
      output$mtData <- function()({
        if(length(clickData$cd) > 0){

          Test_c <- isolate(clickData$cd)

          loadData$Analysis_output$uhgp_90_eggNOG %>%
            select("Gene_id", "Species_id", "Species","Lineage_Shared", "cor_max_species", "Fdrs", "cluster_id", "Predicted_taxonomic_group",
                   "Predicted_protein_name", "eggNOG_free_text_description") %>%
            mutate(Fdrs = format(signif(Fdrs, 3), scientific = TRUE))%>%
            filter(Gene_id %in% Test_c) %>%
            rename("Gene ID" = "Gene_id", "Species ID" = "Species_id", "Lineage Shared" = "Lineage_Shared",
                   "Species with Max Correlation Value" = "cor_max_species", "Cluster ID" = "cluster_id","EggNOG Predicted taxonomic group" = "Predicted_taxonomic_group",
                   "Predicted protein name" = "Predicted_protein_name",
                   "eggNOG free text description" = "eggNOG_free_text_description") %>%
            knitr::kable("html") %>%
            kable_styling("striped", full_width = F)

        }
        else if (length(clickData$x) >0){ #x is the print out from the histogram and then print nothing if nothing else

          Test_c <- c(isolate(clickData$x), isolate(clickData$y))

          loadData$Analysis_output$uhgp_90_eggNOG %>%
            select("Gene_id", "Species_id", "Species","Lineage_Shared", "cor_max_species", "Fdrs", "cluster_id", "Predicted_taxonomic_group",
                   "Predicted_protein_name", "eggNOG_free_text_description") %>%
            mutate(Fdrs = format(signif(Fdrs, 3), scientific = TRUE))%>%
            filter(Gene_id %in% Test_c) %>%
            rename("Gene ID" = "Gene_id", "Species ID" = "Species_id", "Lineage Shared" = "Lineage_Shared",
                   "Species with Max Correlation Value" = "cor_max_species", "Cluster ID" = "cluster_id","EggNOG Predicted taxonomic group" = "Predicted_taxonomic_group",
                   "Predicted protein name" = "Predicted_protein_name",
                   "eggNOG free text description" = "eggNOG_free_text_description") %>%
            knitr::kable("html") %>%
            kable_styling("striped", full_width = F)
        }
      })
    }

    shinyApp(ui = ui, server = server)}
}
