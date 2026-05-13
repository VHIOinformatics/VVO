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
library(GenomicRanges)
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
    fileInput("sce_file", "Upload .rds (SCE object):", accept = ".rds"),
    
    uiOutput("contrast"),
    uiOutput("cluster_var"), 
    uiOutput("num_genes"), 
    uiOutput("mostres"),
    uiOutput("gens"),
    uiOutput("chrom_basic")
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
                <h3 style="color:#2c3e50;">User Guide – VVO Explorer</h3>
                <hr>
                
                <h4 style="color:#2c3e50;">1. Data Upload</h4>
                <div style="background-color:#f9f9f9; padding:10px;">
                  <ul>
                    <li>Upload a valid <code>.RDS</code> file containing a <b>SingleCellExperiment (SCE)</b> object. This object must contain rowRanges (RangedSummarizedExperiment converted to SingleCellExperiment).</li>
                    <li>The object must contain an assay named <code>counts</code> and genomic coordinates in <code>rowRanges</code>.</li>
                    <li>The column <code>Cond</code> must exist in <code>colData</code>, as it is used for filtering and normalization.</li>
                    <li>After upload, the app automatically filters low‑expression genes, normalizes counts (TMM), and computes PCA.</li>
                  </ul>
                </div>
                
                <h4 style="color:#2c3e50;">2. Global Analysis</h4>
                <div style="background-color:#f9f9f9; padding:10px;">
                  <p>This tab provides an overview of sample‑level and gene‑level properties.</p>
                  <ul>
                    <li><b>QC Plot:</b> Visualizes library sizes per sample.</li>
                    <li><b>PCA:</b> Displays the first two principal components using normalized expression.</li>
                    <li><b>Row Data Table:</b> Inspect genomic features.</li>
                    <li><b>Sample and Gene Selectors:</b> Filter which samples and genes are displayed in the plots.</li>
                  </ul>
                </div>
                
                <h4 style="color:#2c3e50;">3. Genomic Exploration</h4>
                <div style="background-color:#f9f9f9; padding:10px;">
                  <p>This tab focuses on genomic visualization using the coordinates stored in <code>rowRanges</code>.</p>
                  <ul>
                    <li><b>Chromosome Selector:</b> Choose the chromosome to visualize.</li>
                    <li><b>Row Data Table:</b> Inspect genomic features and select them for visualization.</li>
                    <li><b>Karyoplot:</b> Displays genomic ranges as rectangles along the selected chromosome.</li>
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
  sce_raw <- reactive({
    req(input$sce_file) # Espera a que hi hagi un fitxer
    readRDS(input$sce_file$datapath)
  })
  # Segon reactive() --> filtrar/normalitzar
  sce_processed <- reactive({
    req(sce_raw())
    obj <- sce_raw()
    
    # Validació del format
    if (!inherits(obj, "SingleCellExperiment")) {
      showNotification("Error: This file is not a SingleCellExperiment type object", closeButton = TRUE, type = "error")
      return(NULL)
    }
    
    # Validació: Existeix la columna "Cond"? 
    if (!"Cond" %in% colnames(colData(obj))) {
      showNotification("Error: Column 'Cond' not found", duration = 30, type = "error")
      return(NULL)
    }
    
    # FILTRATGE
    counts <- assay(obj, "logcounts")
    #keep <- filterByExpr(counts, group = colData(obj)[["Cond"]]) # Mira quins gens tenen prou expressió per ser estadísticament útils segons la teva columna de condicions ("Cond")
    #obj <- obj[keep, ]# ens quedem amb les files (gens) que han passat el filtre i eliminem la resta
    countsF <- assay(obj, "logcounts")
    
    # NORMALITZACIÓ
    #countsTMM <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(countsF)), log = TRUE), ens assegurem que el paquet edgeR esta carregat
    #countsTMM <- normTMM(countsF, log = TRUE) # normalització
    #assay(obj, "countsTMM") <- countsTMM
    countsTMM <- countsF
    
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
    req(sce_processed())
    req(input$tabs)
    if (input$tabs == "basic_tab") {
      sce <- sce_processed()
      cols_fil <- colnames(colData(sce))
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
    req(sce_processed(), input$tabs)
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
    req(sce_processed())
    req(input$tabs)
    
    if (input$tabs == "basic_tab") {
      coldata <- as.data.frame(colData(sce_processed()))
      choices_mostres <- rownames(coldata)

      selectizeInput(
        inputId = "mostres",
        label = "Select samples to visualize:",
        choices = choices_mostres,
        selected = choices_mostres,
        multiple = TRUE,
        options = list(
          placeholder = "Select sample...",
          plugins = list("remove_button")
        )
      )
    } else {
      return(NULL)
    }
  })
  
  # Selecció de gens a graficar --> gens
  output$gens <- renderUI({
    req(sce_processed())
    req(input$tabs)
    req(input$num_genes)
    
    if (input$tabs == "basic_tab") {
      choices_genes <- rownames(sce_processed())[1:input$num_genes]
      
      selectizeInput(
        inputId = "gens",
        label = "Select genes to visualize:",
        choices = choices_genes,
        multiple = TRUE,
        options = list(
          placeholder = "Select genes...",
          plugins = list("remove_button")
        )
      )
    } else {
      return(NULL)
    }
  })
  
  # Genomic exploration
  
  output$chrom_basic <- renderUI({
    req(sce_processed())
    req(input$tabs == "df_tab")
    
    chrs <- unique(as.character(seqnames(rowRanges(sce_processed()))))
    selectInput("chrom_basic", "Chromosome:", choices = chrs)
  })
  

  #3. REACTIUS ISEE_UI________________________________________________________________________________________________________
  
  # Colors per levels de cluster_var
  color_palette <- reactive({
    req(sce_processed(), input$cluster_var)
    levs <- unique(as.character(colData(sce_processed())[[input$cluster_var]]))
    pal <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))
    setNames(pal[seq_along(levs)], levs)
  })
  
  # Selecció de mostres --> quines mostres (cols) conservem després de la selecció de l'usuari
  selected_samples <- reactive({
    req(sce_processed())
    
    samples <- colnames(sce_processed())
    sel <- input$mostres
    
    if (is.null(sel) || length(sel) == 0) {
      return(samples)
    }
    intersect(sel, samples)
  })
  
  # Selecció de gens
  selected_genes <- reactive({
    req(sce_processed())
    
    genes <- rownames(sce_processed())
    sel <- input$gens
    
    if (is.null(sel) || length(sel) == 0) {
      return(genes)
    }
    intersect(sel, genes)
  })
  
  # Objecte filtrat final -> sce_filtrat
  sce_filtrat <- reactive({
    req(sce_processed(), selected_samples(), selected_genes())
    
    input$chrom_basic
    sce <- sce_processed()
    samples <- selected_samples()
    genes <- selected_genes()
    
    sce[genes, samples, drop = FALSE]
  })
  
  # Funció per seleccionar columnes 
  get_selected_sce <- function(sce, columns = NULL) {
    sel_global <- colnames(sce)
    if (!is.null(columns) && length(columns) > 0) {
      sel_local <- unique(unlist(columns))
      sel_local <- sel_local[sel_local %in% sel_global]
      sel <- sel_local
    } else {
      sel <- sel_global
    }
    sce[, sel, drop = FALSE]
  }
  
  # Funció QC per panell iSEE
  QC_fun <- function(sce, rows = NULL, columns = NULL) {
    sce <- get_selected_sce(sce, columns)
    counts_mat <- assay(sce, "logcounts")
    sample_totals <- colSums(counts_mat)
    sample_order <- colnames(sce)
    
    cluster_raw <- as.character(colData(sce)[[input$cluster_var]])
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
    
    remove_grid <- ncol(sce) > 50
    
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
    req(sce_filtrat())
    # Isolate per evitar que l'iSEE s'auto-reiniciï cada cop que canvies un slider
    sce_display <- isolate(sce_filtrat())
    current_cluster_var <- isolate(input$cluster_var)
    
    initial_panels <- list()
    
    # Configuració del panell del QC plot
    QC_plot <- createCustomPlot(
      QC_fun,
      restrict = NULL,
      className = "QCPlot1", #nova classe
      fullName = "Library Size")( PanelHeight = 400L, PanelWidth = 6L, ColumnSelectionDynamicSource = TRUE,
                                  ColumnSelectionSource = "ReducedDimensionPlot1", #default class
                                  RowSelectionRestrict = FALSE,
                                  ColumnSelectionRestrict = TRUE,
                                  SelectionHistory = list())
    

  
    # PANELLS1_____________________________
    
    initial_panels[["QCPlot1"]] <- QC_plot
    
    # PCA
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
    
    # Output final
    iSEE(
      sce_display,
      initial  = initial_panels,
      appTitle = "Global Analysis"
    )
  })
    
  #5. REACTIUS ISEE_UI2 _____________________________________________________________________________________________________________
  
  karyo_fun <- function(sce, rows = NULL, columns = NULL) {
    
    #  Selecció de columnes
    if (!is.null(columns) && length(columns) > 0) {
      cols <- unique(unlist(columns))
      sce <- sce[, cols, drop = FALSE]
    }
    
    gr <- rowRanges(sce)
    
    # Cromosoma seleccionat
    chr_sel <- input$chrom_basic
    if (is.null(chr_sel)) {
      chr_sel <- unique(as.character(seqnames(gr)))[1]
    }
    
    # GRanges del cromosoma seleccionat
    idx_chr <- which(seqnames(gr) == chr_sel)
    gr_chr <- gr[idx_chr]
    
    # Crear karyoplot
    kp <- karyoploteR::plotKaryotype(
      genome = "hg38",
      chromosomes = chr_sel,
      plot.type = 1
    )
    
    # Pintar tots els rangs en blau
    karyoploteR::kpRect(
      kp,
      chr = chr_sel,
      x0 = start(gr_chr),
      x1 = end(gr_chr),
      y0 = 0,
      y1 = 0.5,
      col = "#2c7bb6",
      border = "#2c7bb6"
    )
    
    if (!is.null(rows)) {
      gens_seleccionats <- unlist(rows)
      
      # Busquem quins d'aquests gens estan al cromosoma que estem veient ara
      sel_gr <- gr_chr[rownames(gr_chr) %in% gens_seleccionats]
      
      if (length(sel_gr) > 0) {
        # Dibuixem una línia vertical vermella
        karyoploteR::kpAbline(kp, chr = chr_sel, v = start(sel_gr), col = "red", lwd = 2)
        
        # Dibuixem el rectangle vermell una mica més alt que el blau
        karyoploteR::kpRect(
          kp, chr = chr_sel,
          x0 = start(sel_gr), x1 = end(sel_gr),
          y0 = 0, y1 = 0.8, 
          col = "red", border = "red"
        )
      }
    }
  }
  
  
  #6. LLANÇAMENT DE L'ISEE + PANELLS DE GENOMIC EXPLORATION________________________________________________________________________________

  output$isee_ui2 <- renderUI({
    req(sce_processed(), input$chrom_basic)
    
    sce_display <- sce_processed()

    initial_panels <- list()
    
    ## PANELLS2______________________________
  
    # GRanges table
    initial_panels[["RowDataTable1"]] <- new("RowDataTable",
                                             PanelHeight = 400L,
                                             PanelWidth  = 6L,
                                             SelectionHistory = list())
    
    # Això obliga a iSEE a redibuixar el panell quan canvies el selector
    panel_id <- as.integer(as.factor(input$chrom_basic))
    
    initial_panels[["KaryoPlot1"]] <- createCustomPlot(karyo_fun,
                                                      className = "KaryoPlot1",
                                                      fullName = "Karyoplot"
                                                    )(
                                                      PanelId = panel_id,
                                                      PanelHeight = 800L,
                                                      PanelWidth = 6L,
                                                      ColumnSelectionDynamicSource = TRUE,
                                                      RowSelectionDynamicSource = TRUE,
                                                      RowSelectionSource = "RowDataTable1")
    
  
    iSEE(
      sce_display,
      initial = initial_panels,
      appTitle = "Genomic Exploration")
    })
  
  }

shinyApp(ui, server)
    