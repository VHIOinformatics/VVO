library(shiny)
library(shinyjs)
library(shinydashboard)
library(rintrojs)

library(iSEE)
library(iSEEu)

library(BasicPlots) # Funcions propies del departament VHIOinformatics
library(BasicFunctions) # Funcions propies del departament VHIOinformatics

library(SummarizedExperiment)
library(SingleCellExperiment)
library(GenomicRanges)
library(GenomeInfoDb)
library(IRanges)
library(RColorBrewer)
library(ggplot2)
library(edgeR)
library(karyoploteR)

options(shiny.maxRequestSize = 50 * 1024^2)


# DEFINICIÓ DE LA CLASSE KaryoPlot

# Classe nova que hereta directament de Panel (classe base d'iSEE)
# slot RowSelectionSource: nom del panel del qual llegirà la selecció
setClass("KaryoPlot",
         contains = "Panel",
         slots = c(RowSelectionSource = "character"))

# Constructor amb valor per defecte apuntant a la taula de files
KaryoPlot <- function(RowSelectionSource = "RowDataTable1", ...) {
  new("KaryoPlot", RowSelectionSource = RowSelectionSource, ...)
}

# Nom que mostrarà iSEE a la interfície
setMethod(".fullName", "KaryoPlot", function(x) "Genomic Karyoplot")

# Color del panel a la interfície
setMethod(".panelColor", "KaryoPlot", function(x) "#A84DA2")

# Indica que aquest panel treballa amb seleccions de files (gens)
setMethod(".multiSelectionDimension", "KaryoPlot", function(x) "row")

# Genera el contingut del panel
# contents = NULL perquè aquest panel no transmet seleccions a altres panels
setMethod(".generateOutput", "KaryoPlot",
          function(x, se, all_memory, all_contents) {
            list(
              contents = NULL,
              commands = list("# karyoploteR plot"),
              varname  = NULL
            )
          })

# Element UI del panel: un plotOutput de Shiny amb ID = "KaryoPlot1"
setMethod(".defineOutput", "KaryoPlot", function(x) {
  plotOutput(.getEncodedName(x), height = "600px")
})

# Renderitza el panel: tota la lògica de dibuix aquí
setMethod(".renderOutput", "KaryoPlot",
          function(x, se, ..., output, pObjects, rObjects) {
            
            # Nom codificat d'aquest panel
            panel_name <- .getEncodedName(x)
            # Nom del panel font de la selecció
            panel_src  <- x@RowSelectionSource
            
            output[[panel_name]] <- renderPlot({
              
              # Sense això el plot no es torna a dibuixar quan canvia la selecció
              force(rObjects[[paste0(panel_src, "_INTERNAL_single_select")]])
              
              # Llegim el gen seleccionat des de la memòria del panel font
              # mem@Selected conté el nom de la fila (gen) que s'ha clicat
              sel_rows <- character(0)
              mem <- pObjects$memory[[panel_src]]
              if (!is.null(mem@Selected) && nchar(mem@Selected) > 0) {
                sel_rows <- mem@Selected
              }
              
              # Construïm un GRanges net sense objectes S4
              # karyoploteR només necessita seqnames, start i end
              gr_raw <- SummarizedExperiment::rowRanges(se)
              gr <- GenomicRanges::GRanges(
                seqnames = as.character(GenomeInfoDb::seqnames(gr_raw)),
                ranges   = IRanges::IRanges(
                  start = as.integer(BiocGenerics::start(gr_raw)),
                  end   = as.integer(BiocGenerics::end(gr_raw))
                )
              )
              # Assignem noms dels gens per poder indexar per nom
              names(gr) <- rownames(se)
              
              # Dibuixem el cariotip base del genoma hg38
              kp <- karyoploteR::plotKaryotype(genome = "hg38", plot.type = 1)
              
              # Tots els gens en blau 
              karyoploteR::kpRect(kp,
                                  data   = gr,
                                  y0     = 0, y1 = 0.3,
                                  col    = "#3498db55",
                                  border = NA)
              
              # Gen seleccionat en vermell 
              if (length(sel_rows) > 0 && sel_rows %in% names(gr)) {
                karyoploteR::kpRect(kp,
                                    data   = gr[sel_rows],
                                    y0     = 0, y1 = 0.7,
                                    col    = "red",
                                    border = "red")
              }
            })
          })

# Sense paràmetres configurables per ara
setMethod(".defineInterface", "KaryoPlot", function(x, se, select_info) {
  list()
})


# --------- UI -------------

ui <- dashboardPage(
  dashboardHeader(title = "VHIO's VISUAL OMICS (VVO)", titleWidth = 240),
  
  dashboardSidebar(
    width = 240,
    sidebarMenu(
      id = "tabs",
      menuItem("Global Analysis",      tabName = "basic_tab", icon = icon("eye")),
      menuItem("Genomic Exploration",  tabName = "df_tab",    icon = icon("dna")),
      menuItem("Help",                 tabName = "help_tab",  icon = icon("question-circle"))
    ),
    
    hr(style = "border-top: 2px solid white; margin-top:4px; margin-bottom:4px;"),
    
    fileInput("sce_file", "Upload .rds (SCE object):", accept = ".rds"),
    
    uiOutput("contrast"),
    uiOutput("cluster_var"),
    uiOutput("num_genes"),
    uiOutput("mostres"),
    uiOutput("gens"),
    uiOutput("chrom_basic")
  ),
  
  dashboardBody(
    includeCSS(system.file(package = "iSEE", "www", "iSEE.css")),
    useShinyjs(),
    introjsUI(),
    
    tags$head(
      tags$style(HTML("
        iframe.shiny-frame { height: 1200px !important; }
      "))
    ),
    
    tabItems(
      tabItem(tabName = "basic_tab", uiOutput("isee_ui")),
      tabItem(tabName = "df_tab",    uiOutput("isee_ui2")),
      tabItem(tabName = "help_tab",
              HTML('
          <h3 style="color:#2c3e50;">User Guide – VVO Explorer</h3>
          <hr>
          <h4 style="color:#2c3e50;">1. Data Upload</h4>
          <div style="background-color:#f9f9f9; padding:10px;">
            <ul>
              <li>Upload a valid <code>.RDS</code> file containing a <b>SingleCellExperiment (SCE)</b> object.</li>
              <li>The object must contain an assay named <code>counts</code> and genomic coordinates in <code>rowRanges</code>.</li>
              <li>The column <code>Cond</code> must exist in <code>colData</code>.</li>
            </ul>
          </div>
          <h4 style="color:#2c3e50;">2. Global Analysis</h4>
          <div style="background-color:#f9f9f9; padding:10px;">
            <ul>
              <li><b>QC Plot:</b> Library sizes per sample.</li>
              <li><b>PCA:</b> First two principal components.</li>
            </ul>
          </div>
          <h4 style="color:#2c3e50;">3. Genomic Exploration</h4>
          <div style="background-color:#f9f9f9; padding:10px;">
            <ul>
              <li><b>Row Data Table:</b> Selecciona un gen fent click.</li>
              <li><b>Karyoplot:</b> Mostra el gen seleccionat en vermell.</li>
            </ul>
          </div>
        ')
      )
    )
  )
)


# ------------ SERVER ----------------

server <- function(input, output, session) {
  
  # 1. CARREGA I PROCESSAMENT_______________________________________________________________________________________
  
  # Primer reactive() --> Lectura pura del fitxer
  sce_raw <- reactive({
    req(input$sce_file) # Espera a que hi hagi un fitxer
    readRDS(input$sce_file$datapath)
  })
  
  # Segon reactive() --> filtrar/normalitzar
  sce_processed <- reactive({
    req(sce_raw())
    obj <- sce_raw()
    
    # Validació de format
    if (!inherits(obj, "SingleCellExperiment")) {
      showNotification("Error: not a SingleCellExperiment object",
                       closeButton = TRUE, type = "error")
      return(NULL)
    }
    
    #  Validació: Existeix la columna "Cond"? 
    if (!"Cond" %in% colnames(colData(obj))) {
      showNotification("Error: column 'Cond' not found",
                       duration = 30, type = "error")
      return(NULL)
    }
    
    countsTMM <- assay(obj, "logcounts")
    
    # Càlcul de la PCA sobre les dades normalitzades
    pca <- prcomp(t(countsTMM), scale. = TRUE)
    n   <- ncol(countsTMM)
    pca_scores <- sweep(pca$x, 2, pca$sdev * sqrt(n), FUN = "/")
    reducedDims(obj)$PCA <- pca_scores[, 1:2]
    
    return(obj)
  })
  
  # 2. OUTPUTS DE LA UI_______________________________________________________________________
  
  # Variable per agrupar (metadata) --> cluster_var
  output$cluster_var <- renderUI({
    req(sce_processed(), input$tabs == "basic_tab")
    cols <- colnames(colData(sce_processed()))
    selectizeInput("cluster_var", "Color / Group by:", choices = cols, selected = cols[1])
  })
  
  # Input numeric genes a graficar --> num_genes
  output$num_genes <- renderUI({
    req(sce_processed(), input$tabs == "basic_tab")
    numericInput("num_genes", "Number of genes (min 2 - max 2000):",
                 50, min = 2, max = 2000)
  })
  
  # Selecció de mostres --> mostres
  output$mostres <- renderUI({
    req(sce_processed(), input$tabs == "basic_tab")
    choices <- rownames(as.data.frame(colData(sce_processed())))
    selectizeInput("mostres", "Select samples:",
                   choices = choices, selected = choices, multiple = TRUE,
                   options = list(plugins = list("remove_button")))
  })
  
  #  Selecció de gens a graficar --> gens
  output$gens <- renderUI({
    req(sce_processed(), input$tabs == "basic_tab", input$num_genes)
    choices <- rownames(sce_processed())[seq_len(input$num_genes)]
    selectizeInput("gens", "Select genes:",
                   choices = choices, multiple = TRUE,
                   options = list(plugins = list("remove_button")))
  })
  
  # GENOMIC EXPLORATION
  output$chrom_basic <- renderUI({
    req(sce_processed(), input$tabs == "df_tab")
    chrs <- unique(as.character(seqnames(rowRanges(sce_processed()))))
    selectInput("chrom_basic", "Chromosome:", choices = chrs)
  })
  
  #3. REACTIUS ISEE_UI________________________________________________________________________________________________________
  
  # Paleta de colors per cluster_var 
  color_palette <- reactive({
    req(sce_processed(), input$cluster_var)
    levs <- unique(as.character(colData(sce_processed())[[input$cluster_var]]))
    pal  <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))
    setNames(pal[seq_along(levs)], levs)
  })
  
  # Selecció de mostres --> quines mostres (cols) conservem després de la selecció de l'usuari
  selected_samples <- reactive({
    req(sce_processed())
    sel <- input$mostres
    if (is.null(sel) || length(sel) == 0) return(colnames(sce_processed()))
    intersect(sel, colnames(sce_processed()))
  })
  
  # Selecció de gens
  selected_genes <- reactive({
    req(sce_processed())
    sel <- input$gens
    if (is.null(sel) || length(sel) == 0) return(rownames(sce_processed()))
    intersect(sel, rownames(sce_processed()))
  })
  
  # Objecte filtrat final -> sce_filtrat
  sce_filtrat <- reactive({
    req(sce_processed(), selected_samples(), selected_genes())
    sce <- sce_processed()
    sce[selected_genes(), selected_samples(), drop = FALSE]
  })
  
  # Funció per seleccionar columnes 
  get_selected_sce <- function(sce, columns = NULL) {
    sel <- colnames(sce)
    if (!is.null(columns) && length(columns) > 0) {
      sel_local <- unique(unlist(columns))
      sel <- sel_local[sel_local %in% sel]
    }
    sce[, sel, drop = FALSE]
  }
  
  # Funció QC per al panell d'iSEE
  QC_fun <- function(sce, rows = NULL, columns = NULL) {
    sce         <- get_selected_sce(sce, columns)
    counts_mat  <- assay(sce, "logcounts")
    sample_totals <- colSums(counts_mat)
    cluster_raw <- as.character(colData(sce)[[input$cluster_var]])
    cluster_levels <- sort(unique(cluster_raw))
    
    df <- data.frame(
      sample  = factor(colnames(sce), levels = colnames(sce)),
      total   = sample_totals / 1e6,
      cluster = factor(cluster_raw, levels = cluster_levels)
    )
    
    gg_palette <- function(n) {
      hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[seq_len(n)]
    }
    color_map    <- setNames(gg_palette(length(cluster_levels)), cluster_levels)
    label_colors <- color_map[df$cluster]
    remove_grid  <- ncol(sce) > 50
    
    ggplot(df, aes(x = sample, y = total)) +
      geom_bar(aes(fill = total), stat = "identity") +
      geom_point(aes(colour = cluster), y = -Inf, alpha = 0) +
      scale_colour_manual(values = color_map, name = input$cluster_var) +
      guides(colour = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
      theme_bw() +
      theme(
        panel.grid.major.x = if (remove_grid) element_blank() else element_line(),
        panel.grid.minor.x = if (remove_grid) element_blank() else element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   colour = label_colors, size = 8)
      ) +
      ylab("Million reads") + xlab(NULL) + labs(fill = NULL)
  }
  
  # 4. TAB GLOBAL ANALYSIS + PANELLS______________________________________________________________________________
  
  output$isee_ui <- renderUI({
    req(sce_filtrat())
    # Isolate per evitar que l'iSEE s'auto-reiniciï cada cop que canvies un slider
    sce_display         <- isolate(sce_filtrat())
    current_cluster_var <- isolate(input$cluster_var)
    
    initial_panels <- list()
    
    # Configuració del panell del QC plot
    QC_plot <- createCustomPlot(
      QC_fun,
      restrict  = NULL,
      className = "QCPlot1",
      fullName  = "Library Size")(
        PanelHeight                  = 400L,
        PanelWidth                   = 6L,
        ColumnSelectionDynamicSource = TRUE,
        ColumnSelectionSource        = "ReducedDimensionPlot1",
        RowSelectionRestrict         = FALSE,
        ColumnSelectionRestrict      = TRUE,
        SelectionHistory             = list()
      )
    
    # -------------- PANELLS ------------------
    
    # QC plot
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
  
  
  # 5. TAB GENOMIC EXPLORATION + PANELLS______________________________________________________________________
  
  output$isee_ui2 <- renderUI({
    req(sce_processed(), input$chrom_basic)
    sce_display <- sce_processed()
    
    initial_panels <- list()
    
    # ------------- PANELLS ---------------
    
    # GRanges table
    initial_panels[["RowDataTable1"]] <- new("RowDataTable",
                                             PanelHeight = 600L,
                                             PanelWidth  = 6L,
                                             SelectionHistory = list())
    
    # Karyoplot
    initial_panels[["KaryoPlot"]] <- new("KaryoPlot",
                                         PanelHeight = 600L,
                                         PanelWidth  = 6L,
                                         RowSelectionSource = "RowDataTable1")
    
    
    iSEE(
      sce_display,
      initial = initial_panels,
      appTitle = "Genomic Exploration")
  })
  
}

shinyApp(ui, server)