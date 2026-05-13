library(shiny)
library(shinyjs)
library(shinydashboard)
library(rintrojs)

library(iSEE)
library(iSEEu)

library(BasicPlots) # Funciones propias del departamento VHIOinformatics
library(BasicFunctions) # Funciones propias del departamento VHIOinformatics

library(SummarizedExperiment)
library(SingleCellExperiment)
library(RColorBrewer)
library(ggplot2)
library(edgeR) 

options(shiny.maxRequestSize = 50 * 1024^2)

ui <- dashboardPage(
  dashboardHeader(title = "VHIO's VISUAL OMICS (VVO)", titleWidth = 240),
  
  dashboardSidebar(
    width = 240,
    sidebarMenu(
      id = "tabs",
      menuItem("Global Analysis", tabName = "basic_tab", icon = icon("eye")),
      menuItem("Genomic Exploration", tabName = "df_tab", icon = icon("dna")),
      menuItem("Help", tabName = "help_tab", icon = icon("question-circle"))
    ),
    
    hr(style = "border-top: 2px solid white; margin-top:4px; margin-bottom:4px;"),
    
    # Per a RanggedSummarizedExperiment
    fileInput("rse_file", "Upload .rds (RSE object):", accept = ".rds"),
    
    uiOutput("contrast"),
    uiOutput("cluster_var"), 
    uiOutput("num_genes"), 
    uiOutput("mostres"),
    uiOutput("gens")
  ),
  
  dashboardBody(
    includeCSS(system.file(package="iSEE", "www", "iSEE.css")),
    useShinyjs(),
    introjsUI(),
    
    shiny::tags$head(
      shiny::tags$style(HTML("
        iframe.shiny-frame {
          height: 1200px !important;
        }
      "))
    ),
    
    tabItems(
      tabItem(tabName = "basic_tab", uiOutput("isee_ui")),
      tabItem(tabName = "df_tab", uiOutput("isee_ui2")),
      tabItem(tabName = "help_tab",      
              HTML('
          <h3 style="color:#2c3e50;">User Guide - RSE Explorer</h3>
          <hr>
          <h4 style="color:#2c3e50;">1. Data Upload</h4>
          <div style="background-color: #f9f9f9; padding:10px;">
            <ul>
              <li>Upload a valid <code>.RDS</code> file containing a <b>RangedSummarizedExperiment</b> object.</li>
              <li>The object should contain an assay named <code>counts</code>.</li>
            </ul>
          </div>
          <h4 style="color:#2c3e50;">2. DEA Visualization</h4>
          <div style="background-color: #f9f9f9; padding:10px;">
            <ul>
              <li>Select the contrast from the sidebar to update the Volcano Plots and tables.</li>
              <li>Ensure your RSE object has DEA results in the <code>rowData</code> or <code>metadata</code>.</li>
            </ul>
          </div>
        ')
      )
    ) 
  ) 
)

server <- function(input, output, session) {
  datos_reset <- reactiveVal(NULL)
  #1. CARREGA I PROCESSAMENT_______________________________________________________________________________________
  
  # Primer reactive() --> Lectura pura del fitxer
  rse_raw <- reactive({
    req(input$rse_file) # Espera a que hi hagi un fitxer
    readRDS(input$rse_file$datapath)
  })
  # Segon reactive() --> filtrar/normalitzar
  rse_processed <- reactive({
    req(rse_raw())
    #  Llegim i convertim
    obj <- as(rse_raw(), "SingleCellExperiment")
    
    # Validació del format
    if (!inherits(obj, "RangedSummarizedExperiment")) {
      showNotification("Error: This file is not a RangedSummarizedExperiment type object", closeButton = TRUE, type = "error")
      return(NULL)
    }
    
    # Validació: Existeix la columna "Cond"? 
    if (!"Cond" %in% colnames(colData(obj))) {
      showNotification("Error: Column 'Cond' not found", duration = 30, type = "error")
      return(NULL)
    }
    
    # FILTRATGE
    counts <- assay(obj, "counts")
    keep <- filterByExpr(counts, group = colData(obj)[["Cond"]]) # Mira quins gens tenen prou expressió per ser estadísticament útils segons la teva columna de condicions ("Cond")
    obj <- obj[keep, ]# ens quedem amb les files (gens) que han passat el filtre i eliminem la resta
    countsF <- assay(obj, "counts")
    
    # NORMALITZACIÓ
    #countsTMM <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(countsF)), log = TRUE), ens assegurem que el paquet edgeR esta carregat
    countsTMM <- normTMM(countsF, log = TRUE) # normalització
    assay(obj, "countsTMM") <- countsTMM
    
    # CÀLCUL DE LA PCA
    pca <- prcomp(t(countsTMM), scale. = TRUE)
    n <- ncol(countsTMM)
    pca_scores_makePCA <- sweep(pca$x, 2,pca$sdev * sqrt(n), FUN = "/")
    
    # guardem pca a obj perque isee la trobi
    reducedDims(obj)$PCA <- pca_scores_makePCA[, 1:2]
    
    return(obj)
  })
  
  #2. OUTPUTS DE LA UI________________________________________________________________________________________________
  
  # Variable per agrupar (metadata) --> cluster_var
  output$cluster_var <- renderUI({
    req(rse_processed())
    req(input$tabs)
    if (input$tabs == "basic_tab") {
      rse <- rse_processed()
      cols_fil <- colnames(colData(rse))
      selectizeInput(
        "cluster_var",
        "Color / Group by:",
        choices = cols_fil,
        selected = cols_fil[1]
      )
    } else {
      return(NULL)
    }
  })
  
  # Input numeric genes a graficar --> num_genes
  
  output$num_genes <- renderUI({
    req(rse_processed(), input$tabs)
    if (input$tabs == "basic_tab") {
      numericInput(
        "num_genes",
        "Number of genes to display (min 2 - max 2000):",
        50,
        min =2, max = 2000
      )
    } else {
      return(NULL)
    }
  })
  
  # Selecció de mostres --> mostres
  output$mostres <- renderUI({
    req(rse_processed())
    req(input$tabs)
    
    if (input$tabs == "basic_tab") {
      coldata <- as.data.frame(colData(rse_processed()))
      choices_mostres <- rownames(coldata)
      
      selectizeInput(
        inputId = "mostres",
        label = "Select samples to visualize:",
        choices = choices_mostres,
        selected = choices_mostres,
        multiple = TRUE,
        options = list(
          placeholder = 'Select sample...',
          plugins = list('remove_button')
        )
      )
    } else {
      return(NULL)
    }
  })
  
  # Selecció de gens a graficar --> gens
  output$gens <- renderUI({
    req(rse_processed())
    req(input$tabs)
    req(input$num_genes)
    
    if (input$tabs == "basic_tab") {
      choices_genes <- rownames(rse_processed())[1:input$num_genes]
      
      selectizeInput(
        inputId = "gens",
        label = "Select genes to visualize:",
        choices = choices_genes,
        multiple = TRUE,
        options = list(
          placeholder = 'Select genes...',
          plugins = list('remove_button')
        )
      )
    } else {
      return(NULL)
    }
  })
  
  # GENOMIC EXPLORATION
  
  # Selector de cromosoma
  output$chrom_selector <- renderUI({
    req(rse_processed())
    req(input$tabs)
    
    if (input$tabs == "df_tab") {
      chrs <- unique(as.character(seqnames(rowRanges(rse_processed()))))
      selectInput("chrom", "Select chromosome:", choices = chrs)
    } else {
      return(NULL)
    }
  })
  
  # Selector de regió genòmica
  output$region_selector <- renderUI({
    req(input$tabs)
    
    if (input$tabs == "df_tab") {
      numericRangeInput(
        "region",
        "Genomic region (bp):",
        value = c(1e6, 2e6)
      )
    } else {
      return(NULL)
    }
  })
  
  # Selector de tipus de tracks
  output$track_selector <- renderUI({
    req(input$tabs)
    
    if (input$tabs == "df_tab") {
      checkboxGroupInput(
        "tracks",
        "Tracks to display:",
        choices = c("Karyoplot", "Lolliplot", "Coverage", "Annotation", "Genes"),
        selected = c("Coverage", "Annotation", "Genes")
      )
    } else {
      return(NULL)
    }
  })
  
  #3. REACTIUS ISEE_UI________________________________________________________________________________________________________
  
  # Colors per levels de cluster_var
  color_palette <- reactive({
    req(rse_processed(), input$cluster_var)
    levs <- unique(as.character(colData(rse_processed())[[input$cluster_var]]))
    pal <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))
    setNames(pal[seq_along(levs)], levs)
  })
  
  # Selecció de mostres --> quines mostres (cols) conservem després de la selecció de l'usuari
  selected_samples <- reactive({
    req(rse_processed())
    
    samples <- colnames(rse_processed())
    sel <- input$mostres
    
    if (is.null(sel) || length(sel) == 0) {
      return(samples)
    }
    intersect(sel, samples)
  })
  
  # Selecció de gens
  selected_genes <- reactive({
    req(rse_processed())
    
    genes <- rownames(rse_processed())
    sel <- input$gens
    
    if (is.null(sel) || length(sel) == 0) {
      return(genes)
    }
    intersect(sel, genes)
  })
  
  # Objecte filtrat final -> rse_filtrat
  rse_filtrat <- reactive({
    req(rse_processed(), selected_samples(), selected_genes())
    
    rse <- rse_processed()
    samples <- selected_samples()
    genes <- selected_genes()
    
    rse[genes, samples, drop = FALSE]
  })
  
  # Funció per seleccionar columnes 
  get_selected_rse <- function(rse, columns = NULL) {
    sel_global <- colnames(rse)
    if (!is.null(columns) && length(columns) > 0) {
      sel_local <- unique(unlist(columns))
      sel_local <- sel_local[sel_local %in% sel_global]
      sel <- sel_local
    } else {
      sel <- sel_global
    }
    rse[, sel, drop = FALSE]
  }
  
  # Funció QC per panell iSEE
  QC_fun <- function(rse, rows = NULL, columns = NULL) {
    rse <- get_selected_rse(rse, columns)
    counts_mat <- assay(rse, "counts")
    sample_totals <- colSums(counts_mat)
    sample_order <- colnames(rse)
    
    cluster_raw <- as.character(colData(rse)[[input$cluster_var]])
    cluster_levels <- unique(cluster_raw)
    
    if (all(!is.na(suppressWarnings(as.numeric(cluster_levels))))) {
      cluster_levels <- cluster_levels[order(as.numeric(cluster_levels))]
    } else {
      cluster_levels <- sort(cluster_levels)
    }
    
    df <- data.frame(
      sample  = factor(sample_order, levels = sample_order),
      total   = sample_totals / 1e6,
      cluster = factor(cluster_raw, levels = cluster_levels)
    )
    
    cond_levels <- levels(df$cluster)
    gg_default_palette <- function(n) {
      hues <- seq(15, 375, length = (n + 1))
      hcl(h = hues, l = 65, c = 100)[seq_len(n)]
    }
    palette <- gg_default_palette(length(cond_levels))
    color_map <- setNames(palette, cond_levels)
    label_colors <- color_map[df$cluster]
    
    remove_grid <- ncol(rse) > 50
    
    p <- ggplot(df, aes(x = sample, y = total)) +
      geom_bar(aes(fill = total), stat = "identity") +
      geom_point(aes(colour = cluster), y = -Inf, alpha = 0) +
      scale_colour_manual(values = color_map, name = input$cluster_var) +
      guides(
        colour = guide_legend(
          override.aes = list(shape = 15, size = 5, alpha = 1)
        )
      ) +
      theme_bw() +
      theme(
        panel.grid.major.x = if (remove_grid) element_blank() else element_line(),
        panel.grid.minor.x = if (remove_grid) element_blank() else element_line(),
        axis.text.x = element_text(
          angle = 90, vjust = 0.5, hjust = 1,
          colour = label_colors, size = 8
        )
      ) +
      ylab("Million reads") +
      xlab(NULL) +
      labs(fill = NULL)
    
    return(p)
  }
  
  #4. LLANÇAMENT DE L'ISEE + PANELLS DE GLOBAL ANALYSIS TAB________________________________________________________________________________
  
  output$isee_ui <- renderUI({
    req(rse_filtrat())
    #Isolate per evitar que l'iSEE s'auto-reiniciï cada cop que canvies un slider
    rse_display <- isolate(rse_filtrat())
    current_cluster_var <- isolate(input$cluster_var)
    
    initial_panels <- list()
    
    #Configuració del panell del QC plot
    QC_plot <- createCustomPlot(
      QC_fun,
      restrict = NULL,
      className = "QCPlot1", #nova classe
      fullName = "Library Size")( PanelHeight = 400L, PanelWidth = 6L, ColumnSelectionDynamicSource = TRUE,
                                  ColumnSelectionSource = "ReducedDimensionPlot1", #default class
                                  RowSelectionRestrict = FALSE,
                                  ColumnSelectionRestrict = TRUE,
                                  SelectionHistory = list())
    
    #PANELLS1_____________________________
    
    initial_panels[["QCPlot1"]] <- QC_plot
    
    #PCA
    initial_panels[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot",
                                                     Type = "PCA",
                                                     XAxis = 1L,
                                                     YAxis = 2L,
                                                     ColorByColumnData = current_cluster_var,
                                                     ColorBy           = "Column data",
                                                     PanelHeight       = 400L,
                                                     PanelWidth        = 6L,
                                                     ColumnSelectionDynamicSource = FALSE,
                                                     RowSelectionDynamicSource    = FALSE,
                                                     SelectionHistory             = list())
    # aqui van mes panells
    
    # Output final
    iSEE(
      rse_display,
      initial  = initial_panels,
      appTitle = "Global Analysis of RSE"
    )
  })
  
  #5. REACTIUS ISEE_UI2 _____________________________________________________________________________________________________________
  
}

shinyApp(ui, server)