library(shiny)
library(shinyjs)
library(shinydashboard)
library(rintrojs)

library(iSEE)
library(iSEEu)
library(iSEEde) #per a les subclasses (VolcanoPlot)

library(S4Vectors)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(DeeDeeExperiment)

library(BasicPlots) # Funciones propias del departamento VHIOinformatics
library(BasicFunctions) # Funciones propias del departamento VHIOinformatics

library(RColorBrewer)
library(ggplot2)
library(edgeR)
library(ggvenn)
library(ggrepel)
library(karyoploteR)
library(circlize)


options(shiny.maxRequestSize = 50 * 1024^2)

# ------------ DEFINICIÓ DE LA CLASSE KaryoPlot -------------------

# Classe nova que hereta directament de Panel (classe base d'iSEE)
# slot RowSelectionSource: nom del panel del qual llegirà la selecció
setClass("KaryoPlot",
         contains = "Panel",
         slots = c(RowSelectionSource = "character",
                   ZoomChr = "logical"))

# Constructor amb valor per defecte apuntant a la taula de files
KaryoPlot <- function(RowSelectionSource = "RowDataTable1",ZoomChr = FALSE, ...) {
  new("KaryoPlot", RowSelectionSource = RowSelectionSource, ZoomChr = ZoomChr,...)
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

# Element UI del panel: un plotOutput de Shiny estàdard"
setMethod(".defineOutput", "KaryoPlot", function(x) {
  plotOutput(.getEncodedName(x), height = "600px")
})

# Sense paràmetres configurables per ara 
setMethod(".defineInterface", "KaryoPlot", function(x, se, select_info) {
  list()
})

# ----------- DEFINICIÓ DE LA CLASSE CircusPlot ----------------

# Definició de la Classe
setClass("CircosPlot",
         contains = "Panel",
         slots = c(RowSelectionSource = "character",
                   Contrast1 = "character",
                   Contrast2 = "character"))

# Constructor amb valors per defecte
CircosPlot <- function(RowSelectionSource = "RowDataTable1", 
                       Contrast1 = "FFXvsWT", 
                       Contrast2 = "GEMvsWT", ...) {
  new("CircosPlot", 
      RowSelectionSource = RowSelectionSource, 
      Contrast1 = Contrast1, 
      Contrast2 = Contrast2, ...)
}

# MÈTODES
# Nom que mostra en la interfície
setMethod(".fullName", "CircosPlot", function(x) "Circos Plot")

# Color del panell
setMethod(".panelColor", "CircosPlot", function(x) "#7EC7DE")

# Indica que aquest panel treballa amb seleccions de files (gens)
setMethod(".multiSelectionDimension", "CircosPlot", function(x) "row")

# .generateOutput retorna NULL perquè aquest panell no transmet dades a altres panells
setMethod(".generateOutput", "CircosPlot", function(x, se, all_memory, all_contents) {
  list(contents = NULL, commands = list("# circlize plot"), varname = NULL)
})

# Element UI del panel: un plotOutput de Shiny estàdard
setMethod(".defineOutput", "CircosPlot", function(x) {
  plotOutput(.getEncodedName(x), height = "600px")
})

setMethod(".defineDataInterface", "CircosPlot", function(x, se, select_info) {
  
  rd_cols <- colnames(rowData(se))
  contrastos <- gsub("_log2FoldChange$", "",
                     grep("_log2FoldChange$", rd_cols, value = TRUE))
  
  panel_name <- .getEncodedName(x)
  
  # Seguim el patró PANEL_SLOT que demana la documentació
  list(
    selectInput(
      inputId  = paste0(panel_name, "_Contrast1"),  # PANEL_SLOT
      label    = "Contrast 1:",
      choices  = contrastos,
      selected = x@Contrast1
    ),
    selectInput(
      inputId  = paste0(panel_name, "_Contrast2"),  # PANEL_SLOT
      label    = "Contrast 2:",
      choices  = contrastos,
      selected = x@Contrast2
    )
  )
})

setMethod(".createObservers", "CircosPlot",
          function(x, se, input, session, pObjects, rObjects) {
            
            callNextMethod()
            panel_name <- .getEncodedName(x)
            
            # Quan canvia Contrast1, actualitzem el slot de l'objecte
            observeEvent(input[[paste0(panel_name, "_Contrast1")]], {
              pObjects$memory[[panel_name]]@Contrast1 <- input[[paste0(panel_name, "_Contrast1")]]
              .requestUpdate(panel_name, rObjects)
            }, ignoreInit = TRUE)
            
            observeEvent(input[[paste0(panel_name, "_Contrast2")]], {
              pObjects$memory[[panel_name]]@Contrast2 <- input[[paste0(panel_name, "_Contrast2")]]
              .requestUpdate(panel_name, rObjects)
            }, ignoreInit = TRUE)
          })

setMethod(".hideInterface", "CircosPlot", function(x, field) {
  # Amaguem tots els camps relacionats amb seleccions
  hidden_fields <- "SelectionBoxOpen" 
  if (field %in% hidden_fields) TRUE else callNextMethod()
})

# ---------- UI ---------------

ui <- dashboardPage(
  
  dashboardHeader(title = "VHIO's VISUAL OMICS (VVO)", titleWidth = 240),
  
  dashboardSidebar(
    width = 240,
    sidebarMenu(
      id = "tabs",
      menuItem("Basic Exploration", tabName = "basic_tab", icon = icon("eye")),
      menuItem("DEA and FEA results", tabName = "df_tab", icon = icon("chart-line")),
      menuItem("Genomic Exploration",  tabName = "g_tab",    icon = icon("dna")),
      menuItem("Help", tabName = "help_tab", icon = icon("question-circle"))
    ),
    tags$hr(style = "border-top: 2px solid white; margin-top:4px; margin-bottom:4px;"),
    
    #uiOutput("Dades_input"), #select data source --> no cal perque només hi ha un
    fileInput("dde_file", "Upload .rds data:(DeeDeeExperiment)", accept = ".rds"),
    
    #uiOutput("Entrada_dades"), #carregar fitxer
    uiOutput("contrast"),
    uiOutput("cluster_var"), #cluster by
    uiOutput("num_genes"), #quantitat de gens
    uiOutput("mostres"),
    uiOutput("gens"),
    uiOutput("padj"),
    uiOutput("logFC"),
    uiOutput("chr_selector")
  ),
  
  dashboardBody(
    includeCSS(system.file(package="iSEE", "www", "iSEE.css")),
    useShinyjs(),
    introjsUI(),
    
    tags$head(
      tags$style(HTML("
        iframe.shiny-frame {
          height: 1200px !important;
        }
      "))
    ),
    
    tabItems(
      tabItem(tabName = "basic_tab", uiOutput("isee_ui")),
      tabItem(tabName = "df_tab", uiOutput("isee_ui2")),
      tabItem(tabName = "g_tab",    uiOutput("isee_ui3")),
      tabItem(tabName = "help_tab",    
              HTML('
          <h3 style="color:#2c3e50;">User Guide</h3>
          <hr>

          <h4 style="color:#2c3e50;">1. Data Upload</h4>
          <div style="background-color: #f9f9f9; padding:10px;">
            <ul>
              <li>Select the data type: <b>Bulk RNA-seq (.RDS)</b> or <b>Single-cell RNA-seq (.RDS)</b>.</li>
              <li>Upload a valid <code>.RDS</code> file containing a <code>SummarizedExperiment</code> or <code>SingleCellExperiment</code> object.</li>
              <li>Other formats will cause errors.</li>
            </ul>
          </div>

          <h4 style="color:#2c3e50;">2. Clustering Options</h4>
          <div style="background-color: #f9f9f9; padding:10px;">
            <ul>
              <li><b>Clustering by:</b> Select a column from <code>colData</code> to group samples or cells.</li>
              <li>Color schemes are automatically assigned to unique levels of the selected variable.</li>
              <li>Dendrograms are generated using correlation or Euclidean distance (configurable internally).</li>
            </ul>
          </div>

          <h4 style="color:#2c3e50;">3. Parameter Configuration</h4>
          <div style="background-color: #f9f9f9; padding:10px;">
            <ul>
              <li><b>Number of genes:</b> Choose between 2 and 2000 genes to display in heatmaps.</li>
              <li><b>Samples / Cells:</b> Filter which samples (bulk) or clusters (single-cell) to include in the analysis.</li>
              <li><b>Genes to visualize:</b> Optionally select specific genes to focus the plots.</li>
            </ul>
          </div>

          <h4 style="color:#2c3e50;">4. iSEE Panels Overview</h4>
          <div style="background-color: #f9f9f9; padding:10px;">
            <ul>
              <li><b>QCPlot:</b> Displays library size (million reads) per sample or cluster. Colors indicate cluster membership.</li>
              <li><b>ReducedDimensionPlot:</b> PCA of samples/cells. Color points by the selected cluster variable.</li>
              <li><b>ComplexHeatmapPlot:</b> Heatmap of top variable genes or selected genes. Rows clustered and scaled by default. Column selection linked to PCA plot.</li>
              <li><b>SampleAssayPlot:</b> Allows visualization of counts for individual genes or samples.</li>
              <li>Selections in one panel propagate to others if dynamic selection is enabled (ColumnSelectionDynamicSource = TRUE).</li>
            </ul>
          </div>

          <h4 style="color:#2c3e50;">5. Advanced Options & Customization</h4>
          <div style="background-color: #f9f9f9; padding:10px;">
            <ul>
              <li>Adjust visual parameters via panel boxes: color, size, text labels, contour overlays, violin plots, etc.</li>
              <li>Downsampling is available for large datasets.</li>
              <li>Panel dimensions and aspect ratios can be modified to improve visualization.</li>
            </ul>
          </div>

          <h4 style="color:#2c3e50;">6. Exporting & Saving Plots</h4>
          <div style="background-color: #f9f9f9; padding:10px;">
            <ul>
              <li>Plots can be downloaded directly from the iSEE interface using the download buttons.</li>
              <li>Changing analysis type or cluster variable will reset iSEE panels to reflect the new selection.</li>
            </ul>
          </div>

          <h4 style="color:#2c3e50;">7. Additional information</h4>
          <div style="background-color: #d5d2d4; padding:10px;">
            <ul>
              <li>If plots do not appear, verify that the RDS file contains the expected object type.</li>
              <li>For bulk RNA-seq, PCA will be computed automatically if missing.</li>
              <li>For large single-cell datasets, selecting too many genes or cells may slow down the app.</li>
              <li>Always select at least 2 samples/cells for clustering and heatmaps.</li>
              <li>If you have any questions or issues, please contact the Bioinformatics Unit directly.</li>
            </ul>
          </div>
        ')
      )
    )
  )
)

# --------------- SERVER ---------------------

server <- function(input, output, session) {
  
  # ---------- Karyoplot -----------
  
  # Valor reactiu per guardar el cromosoma seleccionat.
  chr_sel_reactive <- reactiveVal("chr1") 
  
  # Observador que escolta canvis a input$chr_sel.
  observe({
    req(input$chr_sel)
    chr_sel_reactive(input$chr_sel)
  })
  
  # Renderitza el panel: tota la lògica de dibuix aquí
  setMethod(".renderOutput", "KaryoPlot",
            function(x, se, ..., output, pObjects, rObjects) {
              
              # Nom codificat d'aquest panel
              panel_name <- .getEncodedName(x)
              # Nom del panel font de la selecció
              panel_src  <- x@RowSelectionSource
              is_zoom    <- x@ZoomChr
              
              output[[panel_name]] <- renderPlot({
                
                # Sense això el plot no es torna a dibuixar quan canvia la selecció
                force(rObjects[[paste0(panel_src, "_INTERNAL_single_select")]])
                
                par(mar = c(1, 1, 2, 1))   # redueix marges del plot
                
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
                
                if (is_zoom) {
                  # cromosoma concret del selector
                  chr_sel <- chr_sel_reactive()
                  # dimensions del plot 
                  pp <- karyoploteR::getDefaultPlotParams(plot.type = 1)
                  pp$topmargin      <- 0   
                  pp$bottommargin   <- 150   
                  pp$ideogramheight <- 25    
                  pp$leftmargin     <- 0.15
                  
                  # Dibuixem el cariotip base del genoma hg38
                  # dibuixem el cromosoma selecciont
                  kp <- karyoploteR::plotKaryotype(
                    genome      = "hg38",
                    chromosomes = chr_sel,
                    plot.type   = 1,
                    plot.params = pp
                  )
                } else {
                  # sino dibuixem tots els cromosomes
                  kp <- karyoploteR::plotKaryotype(genome = "hg38", plot.type = 1)
                }
                
                # tots els gens en blau
                karyoploteR::kpRect(kp, data = gr,
                                    y0 = 0, y1 = 0.3,
                                    col = "#3498db55", border = NA)
                
                # gen seleccionat a la taula en vermell
                if (length(sel_rows) > 0 && sel_rows %in% names(gr)) {
                  karyoploteR::kpRect(kp, data = gr[sel_rows],
                                      y0 = 0, y1 = 0.7,
                                      col = "red", border = "red", lwd = 3)
                }
              })  
            })
  
  # --------------- CircusPlot --------------
  
  setMethod(".renderOutput", "CircosPlot",
            function(x, se, ..., output, pObjects, rObjects) {
              
              # Agafem el nom únic d'aquest panell
              # iSEE el genera automàticament a partir del nom de la classe + ID
              panel_name <- .getEncodedName(x)
              
              # Creem el renderPlot associat a aquest panell
              output[[panel_name]] <- renderPlot({
                
                # Forcem reactivitat
                force(rObjects[[paste0(panel_name, "_INTERNAL_output_update")]])
                
                # Llegim dels slots, que ja estan actualitzats per .createObservers
                c1 <- pObjects$memory[[panel_name]]@Contrast1
                c2 <- pObjects$memory[[panel_name]]@Contrast2
                
                lfc1_nm <- paste0(c1, "_log2FoldChange")
                lfc2_nm <- paste0(c2, "_log2FoldChange")
                padj_nm <- paste0(c1, "_padj")
                padj_nm2 <- paste0(c2, "_padj")
                
                # Convertim rowData a dataframe
                rd <- as.data.frame(rowData(se))
                
                gr <- rowRanges(se)
                
                # Construïm el data.frame en format BED a partir de rowRanges
                bed_data <- data.frame(
                  chr   = as.character(seqnames(gr)),
                  start = start(gr),
                  end   = end(gr),
                  val1  = as.numeric(rowData(se)[[lfc1_nm]]),
                  val2  = as.numeric(rowData(se)[[lfc2_nm]]),
                  p1    = as.numeric(rowData(se)[[padj_nm]]),   # padj contrast 1
                  p2    = as.numeric(rowData(se)[[padj_nm2]])    # padj contrast 2
                )
                
                # Eliminem NAs
                bed_data <- bed_data[complete.cases(bed_data), ]
                
                # Cromosomes estàndards
                bed_data <- bed_data[bed_data$chr %in% paste0("chr", c(1:22, "X", "Y")), ]
                
                validate(need(nrow(bed_data) > 0, "No valid data to plot"))
                
                # Filtres independents per cada contrast
                bed_c1 <- bed_data[bed_data$p1 < 0.05 & abs(bed_data$val1) > 1, ]
                bed_c2 <- bed_data[bed_data$p2 < 0.05 & abs(bed_data$val2) > 1, ]
                bed_sig <- bed_data[bed_data$p1 < 0.05 | bed_data$p2 < 0.05, c("chr","start","end")]
                
                validate(need(nrow(bed_c1) > 0 | nrow(bed_c2) > 0, "No significant genes found"))
                
                # --- PLOT ---
                # Netegem qualsevol plot de circlize anterior
                circlize::circos.clear()
                
                # Paràmetres globals del circos
                circlize::circos.par(
                  start.degree = 90, # el cromosoma 1 comença a dalt
                  gap.degree = 2, # espai entre cromosomes
                  track.height = 0.2 # alçada de cada anell
                )
                
                # Inicialitzem hg38
                circlize::circos.initializeWithIdeogram(
                  species = "hg38",
                  plotType = c("ideogram", "labels")
                )
                
                # ------- TRACK 1: Densitat ----------
                # Seleccionem tots els gens amb padj < 0.05 (sense filtre de logFC)
                bed_sig <- bed_data[
                  (bed_data$p1 < 0.05 | bed_data$p2 < 0.05),
                  c("chr", "start", "end")
                ]
                
                # histograma de densitat, més gens significatius = pic més alt
                if (nrow(bed_sig) > 0) {
                  circlize::circos.genomicDensity(bed_sig, col = "#A84DA288", count_by = "number")
                }
                
                # ------- TRACK 2: logFC del contrast 1 ---------
                
                # logFC contrast 1 (només gens significatius del contrast 1)
                if (nrow(bed_c1) > 0) {
                  circlize::circos.genomicTrackPlotRegion(
                    
                    # Li passem chr, start, end i val1 (logFC del contrast 1)
                    bed_c1[, c("chr", "start", "end", "val1")],
                    
                    # L'eix Y va des del mínim fins al màxim de val1 
                    ylim = range(c(bed_c1$val1, 0), na.rm = TRUE),
                    
                    # panel.fun s'executa una vegada per cada cromosoma
                    panel.fun = function(region, value, ...) {
                      
                      v_raw <- as.numeric(value[[1]])
                      
                      # Vermell si el gen puja (logFC > 0) i blau si baixa (logFC < 0)
                      col_barras <- ifelse(v_raw > 0, "#e74c3c", "#3498db")
                      
                      # Es dibuixa un punt per gen a la seva posició cromosòmmica
                      circos.genomicPoints(
                        region,
                        value,
                        col = col_barras,
                        pch = 16,
                        cex = 0.5
                      )}
                  )}
                
                # ------- TRACK 3: logFC del contrast 2 ---------
                
                # logFC contrast 2 (només gens significatius del contrast 2)
                if (nrow(bed_c2) > 0) {
                  circlize::circos.genomicTrackPlotRegion(
                    
                    bed_c2[, c("chr", "start", "end", "val2")],
                    
                    ylim = range(c(bed_c2$val2, 0), na.rm = TRUE),
                    
                    panel.fun = function(region, value, ...) {
                      
                      v_raw <- as.numeric(value[[1]])
                      
                      # Groc si el gen puja (logFC > 0) i verd si baixa (logFC < 0)
                      col_barras <- ifelse(v_raw > 0, "#f1c40f", "#2ecc71")
                      
                      circos.genomicPoints(
                        region,
                        value,
                        col = col_barras,
                        pch = 16,
                        cex = 0.5
                      )}
                  )}
                
                # Titol centrat al plot (ex: "Circos: FFXvsWT vs GEMvsWT")
                title(paste("Circos:", c1, "vs", c2), cex.main = 1.5)
                
                legend("topleft", 
                       legend = c("Density", 
                                  paste(c1, "up"), paste(c1, "down"), 
                                  paste(c2, "up"), paste(c2, "down")),
                       col = c("#A84DA2","#e74c3c","#3498db","#f1c40f","#2ecc71"),
                       pch = c(15, 16, 16, 16, 16),
                       pt.cex = 1.5,
                       bty = "n",
                       cex = 1)
              })
            }
  )
  
  datos_reset <- reactiveVal(NULL)
  
  #1. CARREGA I PROCESSAMENT_______________________________________________________________________________________
  
  # Primer reactive() --> Lectura pura del fitxer
  dde_raw <- reactive({
    req(input$dde_file)
    readRDS(input$dde_file$datapath)
  })
  #Segon reactive() --> filtrar/normalitzar
  dde2 <- reactive({
    req(dde_raw())
    obj <- dde_raw()
    
    # Validació de format
    if (!inherits(obj, "DeeDeeExperiment")) {
      showNotification("Error: This file is not a DeeDeeExperiment type object", closeButton = TRUE, type = "error")
      return(NULL)
    }
    
    #validació: existeix la columna Cond?
    if (!"Cond" %in% colnames(colData(obj))) {
      showNotification("Error: Column 'Cond' not found in the metadata", duration = 30, type = "error")
      return(NULL)
    }
    
    #FILTRATGE
    counts <- assay(obj, "counts")
    keep <- filterByExpr(counts, group = colData(obj)[["Cond"]]) #Mira quins gens tenen prou expressió per ser estadísticament útils segons la teva columna de condicions ("Cond")
    obj <- obj[keep, ]#ens quedem amb les files (gens) que han passat el filtre i eliminem la resta
    countsF <- assay(obj, "counts")
    
    #NORMALITZACIÓ
    countsTMM <- normTMM(countsF, log = TRUE) 
    assay(obj, "countsTMM") <- countsTMM
    
    #CÀLCUL DE LA PCA
    pca <- prcomp(t(countsTMM), scale. = TRUE)
    n <- ncol(countsTMM)
    pca_scores_makePCA <- sweep(pca$x, 2,pca$sdev * sqrt(n), FUN = "/")
    
    #guardem pca a obj perque isee la trobi
    reducedDims(obj)$PCA <- pca_scores_makePCA[, 1:2]
    
    return(obj)
  })
  
  #2. OUTPUTS DE LA UI________________________________________________________________________________________________
  
  # Variable per agrupar (metadata) --> cluster_var
  output$cluster_var <- renderUI({
    req(dde2())
    dde <- dde2()
    if (input$tabs %in% c("basic_tab","df_tab")) {
      cols_fil <- colnames(colData(dde))
      selectizeInput("cluster_var","Clustering by:",
                     choices = cols_fil,
                     selected = cols_fil[1])
    }
  })
  
  # Input númerico genes a graficar --> num_genes
  output$num_genes <- renderUI({
    req(dde2(), input$tabs)
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
  
  
  #Selecció de mostres --> mostres
  output$mostres <- renderUI({
    req(dde2())
    req(input$tabs)
    
    if (input$tabs %in% c("basic_tab","df_tab")) {
      coldata <- as.data.frame(colData(dde2()))
      if (input$tabs == "df_tab") {
        req(input$contrast)
        grups <- unlist(strsplit(input$contrast, "vs"))
        grups_finals <- unique(c(grups, "WT")) #per eliminar duplicat WT 
        choices_mostres <- rownames(coldata[coldata$Cond %in% grups_finals, ])
        label_text <- paste("Select samples for selected contrast/s:")
      } else {
        choices_mostres <- rownames(coldata)
        label_text <- "Select samples to visualize"
      }
      if (length(choices_mostres) == 0) choices_mostres <- rownames(coldata)
      selectizeInput(
        inputId = "mostres",
        label = label_text,
        choices = choices_mostres,
        selected = choices_mostres,
        multiple = TRUE,
        options = list(
          placeholder = 'Select sample...',
          plugins = list('remove_button')
        ))
    } else { #help_tab
      return(NULL)
    }
  })#render
  
  # Desplegable selección genes a graficar
  output$gens <- renderUI({
    req(dde2())
    req(input$num_genes)
    req(top_genes())
    
    if (input$tabs %in% c("basic_tab","df_tab")) {
      selectizeInput(
        inputId = "gens",
        label = "Select gens to visualize:",
        choices = top_genes(),
        multiple = TRUE,
        options = list(
          placeholder = 'Select genes...',
          plugins = list('remove_button')
        )
      )}
  })
  
  output$padj <- renderUI({
    req(input$tabs)
    if (input$tabs == "df_tab") {
      numericInput("padj", "P adjusted Value Threshold:", value = 0.05, min = 0, max = 1, step = 0.01)
    }
  })
  
  output$logFC <- renderUI({
    req(input$tabs)
    if (input$tabs == "df_tab") {
      numericInput("logFC", "Log2 Fold Change Threshold:", value = 1, min = 0, max = 5, step = 0.1)
    }
  })
  
  output$chr_selector <- renderUI({
    req(dde2(), input$tabs == "g_tab")
    
    cromosomes <- paste0("chr", c(1:22, "X", "Y"))
    
    selectizeInput(
      inputId = "chr_sel",
      label = "Zoom to chromosome:",
      choices = cromosomes,
      selected = "chr1",
      multiple = FALSE
    )
  })
  
  #3. REACTIUS ISEE_UI________________________________________________________________________________________________________
  
  # Colors i nivells
  color_palette <- reactive({
    req(dde2(), input$cluster_var)
    levs <- unique(as.character(colData(dde2())[[input$cluster_var]]))
    pal <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))
    setNames(pal[seq_along(levs)], levs)
  })
  
  #Selecció de mostres --> quines mostres (cols) conservem després de la selecció de l'usuari
  selected_samples <- reactive({
    req(dde2())
    
    samples <- colnames(dde2())
    sel <- input$mostres
    
    if (is.null(sel) || length(sel) == 0) {
      return(samples)
    }
    intersect(sel, samples)
  })
  
  #Funció per calcular els gens variables
  get_top_variable_genes <- function(obj, n) {
    counts_mat <- assay(obj, "countsTMM")
    gene_vars <- apply(counts_mat, 1, var, na.rm = TRUE)
    
    keep <- !is.na(gene_vars) & gene_vars > 1e-10
    gene_vars <- gene_vars[keep]
    
    if (length(gene_vars) == 0) return(character(0))
    n_final <- min(n, length(gene_vars))
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_final]
    return(top_genes)
  }
  
  selected_samples <- reactive({
    req(dde2(), input$mostres)
    intersect(input$mostres, colnames(dde2()))
  })
  
  #quan selected samples canvia --> dde_subset es genera
  dde_subset <- reactive({
    req(dde2(), selected_samples())
    dde2()[, selected_samples()]
  })
  
  #com dde_subset canvia top_genes es torna a calcular
  top_genes <- reactive({
    req(dde_subset(), input$num_genes)
    get_top_variable_genes(dde_subset(), input$num_genes)
  })
  
  #Selecció de gens --> quins gens (files) conservem després de la selecció de l'usuari
  selected_gens <- reactive({
    req(top_genes(), dde2())
    
    gens_top <- top_genes()
    sel_gens <- input$gens
    
    if (is.null(sel_gens) || length(sel_gens) == 0) {
      return(gens_top)
    }
    intersect(sel_gens, rownames(dde2()))
  })
  
  #Objecte final filtrat = dde_filtrat
  dde_filtrat <- reactive({
    req(dde_subset(), selected_gens())
    
    dde_subset()[selected_gens(), , drop = FALSE]
  })
  
  # Función per a la selecció de columnes (sce --> dde)
  get_selected_dde <- function(dde, columns = NULL) {
    sel_global <- colnames(dde)
    if (!is.null(columns) && length(columns) > 0) {
      sel_local <- unique(unlist(columns))
      sel_local <- sel_local[sel_local %in% sel_global]
      sel <- sel_local
    } else {
      sel <- sel_global
    }
    dde[, sel, drop = FALSE]
  }
  
  cluster_by <- reactive({
    req(input$cluster_var)
    return(input$cluster_var)
    #cluster_by <- input$cluster_var
  })
  
  # Función QC
  QC_fun <- function(dde, rows = NULL, columns = NULL) {
    selected_obj <- get_selected_dde(dde, columns)
    samples_actuals <- as.character(colnames(selected_obj))
    # sample.totals <- apply(counts(dde), 2, sum) 
    # sample_order <- colnames(dde)
    
    all_counts <- assay(isolate(dde2()), "counts")[, samples_actuals, drop = FALSE]
    sample.totals <- colSums(all_counts)
    
    cluster_raw <- as.character(colData(dde)[[input$cluster_var]])
    cluster_levels <- unique(cluster_raw)
    
    if (all(!is.na(suppressWarnings(as.numeric(cluster_levels))))) {
      cluster_levels <- cluster_levels[order(as.numeric(cluster_levels))]
    } else {
      cluster_levels <- sort(cluster_levels)
    }
    
    sample.totals.df <- data.frame(
      sample = factor(samples_actuals, levels = samples_actuals), #sample_order
      total = as.numeric(sample.totals) / 1e6,
      cluster = factor(cluster_raw, levels = cluster_levels)
    )
    
    cond_levels <- levels(sample.totals.df$cluster)
    
    gg_default_palette <- function(n) {
      hues <- seq(15, 375, length = (n + 1))
      hcl(h = hues, l = 65, c = 100)[seq_len(n)]
    }
    
    palette <- gg_default_palette(length(cond_levels))
    color_palette <- setNames(palette, cond_levels)
    
    label_colors <- color_palette[sample.totals.df$cluster]
    
    remove_grid <- ncol(dde) > 50
    #elements a la llegenda:
    num_grups <- length(unique(sample.totals.df$cluster))
    m_sel <- length(input$mostres)
    
    #mida_text <- if(num_grups > 40) 4 else if(num_grups > 20) 7 else 9
    mida_quadrat <- if(num_grups > 40) 2.5 else 3 #if(num_grups > 10) 3 else 5
    mida_titol <- if(num_grups > 20) 8 else 10
    mida_text <- if(num_grups > 40) 5 else if(num_grups > 20) 7 else 9
    mida_eix_x <- if(m_sel > 50) 5 else if(m_sel > 30) 7 else 9
    grid_x <- if(remove_grid) element_blank() else element_line()
    p <- ggplot(data=sample.totals.df, aes(x=sample, y=total)) + 
      geom_bar(aes(fill = total), stat = "identity") + #fill= total
      geom_point(aes(colour = cluster), y = -Inf, alpha = 0) +
      scale_colour_manual(values = color_palette, name = input$cluster_var) +
      guides(
        colour = guide_legend(
          override.aes = list(
            shape = 15,
            size = mida_quadrat, #5
            alpha = 1
          )
        )
      ) +
      theme_bw() +
      theme(
        panel.grid.major.x = grid_x,
        panel.grid.minor.x = grid_x,
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust=1,
          colour = label_colors,
          size = mida_eix_x),
        legend.text = element_text(size = mida_text),
        legend.title = element_text(size = mida_titol)
      ) +
      ylab("Million reads") +
      xlab(NULL) +
      labs(fill = NULL) 
    
    return(p)   
  }
  
  # oneCluster_iSEE
  oneCluster_iSEE <- function(estimates, conditions = NULL,
                              distance="correlation", method="ward.D2",
                              title=NULL, ...) {
    ordered_levels <- function(x) {
      x <- as.character(x)
      lev <- unique(x)
      if (all(!is.na(suppressWarnings(as.numeric(lev))))) {
        lev[order(as.numeric(lev))]
      } else {
        sort(lev)
      }
    }
    
    labels <- colnames(estimates)
    n_mostres <- length(labels)
    
    if (n_mostres < 2) {
      plot.new()
      title(main = paste(title, "\n(Select at least 2 samples)")) #abans: Must have n >= 2 objects selected
      return(recordPlot())
    }
    
    mida_labels <- if(n_mostres > 50) 0.4 else if(n_mostres > 20) 0.6 else 0.8
    
    parameters <- setParameters(labels)
    use.cor <- "pairwise.complete.obs"
    
    mat_dist <- if (distance == "correlation") {
      as.dist(1 - cor(estimates, use = use.cor))
    } else {
      dist(t(estimates))
    }
    
    # Si la matriu de distància té NAs (per variància 0), els convertim a 0 o un valor neutre
    if (any(is.na(mat_dist))) {
      mat_dist[is.na(mat_dist)] <- 0
    }
    
    clust <- hclust(mat_dist, method = method)
    
    #calcul distància
    if (distance == "correlation") {
      clust <- hclust(as.dist(1 - cor(estimates, use = use.cor)), method = method)
      xlab <- paste("Distance: Correlation / Linkage:", method, sep = "-") #abans: Correlation
    } else {
      clust <- hclust(dist(t(estimates)), method = method)
      xlab <- paste("Distance: Euclidean / Linkage:", method, sep = "-") #abans: Euclidean
    }
    
    cond_levels <- ordered_levels(conditions)
    conditions <- factor(conditions, levels = cond_levels)
    # conditions <- as.character(conditions)
    # cond_levels <- sort(unique(conditions))
    
    #paleta de colors estàndard (de ggplot2) si no passen una específica
    gg_color_hue <- function(n) { #abans: gg_default_palette
      hues <- seq(15, 375, length = (n + 1))
      hcl(h = hues, l = 65, c = 100)[seq_len(n)]
    }
    
    color_map <- setNames(gg_color_hue(length(cond_levels)), cond_levels)
    #palette <- gg_default_palette(length(cond_levels)) # ^
    #color_palette <- setNames(palette, cond_levels)   #  |: ara es fa en una sola linia
    sample_colors <- color_map[conditions]
    #colors <- color_palette[conditions]
    
    layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))       # Capa gràfic amb llegenda lateral
    #Dendograma
    par(mar = c(7, 4, 4, 1)) #5                                   # Capa dendrograma
    clust_col <- colorCluster(clust, sample_colors,ce = mida_labels) #0.8
    plot(clust_col, main = title, xlab = xlab, sub = "")
    # Capa leyenda
    par(mar = c(7, 0, 4, 1)) #5
    plot.new()
    legend("center", #llegenda automàticament centrada sempre
           legend = cond_levels,
           fill = color_map, #alinea text i quadrat
           title = "Groups",
           cex = mida_labels, #0.8
           bty = "n") #traiem la caixa negra del voltant de la llegenda
    
    layout(1)
    return(recordPlot())
  }
  
  # Función dendrograma
  dendro_fun <- function(dde, cluster_by, genes_use, rows = NULL, columns = NULL) {
    
    dde <- dde[genes_use, , drop = FALSE]
    
    keep <- filterByExpr(
      assay(dde, "counts"),
      group = colData(dde)[[cluster_by]]
    )
    countsF <- assay(dde, "counts")[keep, , drop = FALSE]
    countsTMM <- normTMM(countsF, log = TRUE)
    
    if (!is.null(columns) && length(columns) > 0) {
      sel <- intersect(unique(unlist(columns)), colnames(countsTMM))
    } else {
      sel <- colnames(countsTMM)
    }
    
    m <- countsTMM[, sel, drop = FALSE]
    #linies de seguretat
    gene_vars <- apply(m, 1, var, na.rm = TRUE)
    #mantenim els gens que realment varien (var > 0)
    m <- m[which(gene_vars > 0 & !is.na(gene_vars)), , drop = FALSE]
    
    cond <- colData(dde)[[cluster_by]][match(sel, colnames(dde))]
    
    p <- oneCluster_iSEE(
      estimates = m,
      distance = "euclidean",
      method = "ward.D2",
      conditions = cond
    )
    return(p)
  }
  
  #4. LLANÇAMENT DE L'ISEE + PANELLS DE LA *MAIN TAB*________________________________________________________________________________
  
  output$isee_ui <- renderUI({
    req(dde_subset(), input$mostres, input$cluster_var, input$num_genes, dde_filtrat(), cluster_by())
    
    
    #Calculem variancia sobre el subset de mostres seleccionades
    mat1 <- assay(dde_subset(), "countsTMM")
    vars1 <- apply(mat1, 1, var, na.rm = TRUE)
    gens_valids <- names(vars1)[!is.na(vars1) & vars1 > 1e-08]
    
    if (length(gens_valids) < 10) {
      return(h4("Not enough variable genes found for these samples.", style="text-align:center; color:red;"))
    }
    
    dde_display <- dde_filtrat()
    current_cluster_var <- cluster_by()#input$cluster_var
    #ens assegurem que la variable sigui factor
    colData(dde_display)[[current_cluster_var]] <- as.factor(colData(dde_display)[[current_cluster_var]])
    
    # dde3 <- dde3[gens_valids, mostres_sub]
    # gens_finals <- intersect(gens_unio, gens_valids)
    
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
    
    initial_panels <- list()
    initial_panels[["QCPlot1"]] <- QC_plot
    
    #Configuració del panell del Dendrograma
    dendro_plot <- createCustomPlot(
      function(dde, rows, columns) {
        
        dendro_fun(
          dde = dde, #iSEE passa automàticament l'objecte ACTUAL --> manté actualitzats els gràfics segons els canvis a altres panells
          cluster_by = current_cluster_var,
          genes_use = rownames(dde), #els gens que volem fer servir son els noms de les files de l'objecte actual
          rows = rows,
          columns = columns
        )
      },
      
      restrict = NULL,
      className = "DendroPlot1",
      fullName = "Dendrograma"
    )(
      PanelHeight = 400L,
      PanelWidth = 6L,
      ColumnSelectionDynamicSource = TRUE,
      ColumnSelectionSource = "ReducedDimensionPlot1",
      RowSelectionRestrict = FALSE,
      ColumnSelectionRestrict = TRUE,
      SelectionHistory = list()
    )
    
    # ----------- PANELLS ----------------
    
    initial_panels[["DendroPlot1"]] <- dendro_plot
    # label <- if (length(input$mostres) > 30) FALSE else TRUE  
    if (length(input$mostres) >= 2) {
      initial_panels[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "PCA", XAxis = 1L, YAxis = 2L,
                                                       ColorByColumnData = current_cluster_var, ColorByFeatureNameAssay = "counts",
                                                       ColorBy = "Column data", ColorBySampleNameColor = "#FF0000",
                                                       SizeByColumnData = NA_character_, TooltipColumnData = character(0),
                                                       FacetRowBy = "None", FacetColumnBy = "None",
                                                       ColorByDefaultColor = "#FFFFFF",
                                                       ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
                                                       ColorBySampleSource = "---",
                                                       ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
                                                       VisualBoxOpen = FALSE, VisualChoices = c("Color", "Size",
                                                                                                "Text"), ContourAdd = FALSE, ContourColor = "#0000FF", FixAspectRatio = FALSE,
                                                       ViolinAdd = TRUE, PointSize = 3, PointAlpha = 1, Downsample = FALSE,
                                                       DownsampleResolution = 200, CustomLabels = FALSE, #FALSE
                                                       FontSize = 1, LegendPointSize = 2, LegendPosition = "Right",
                                                       HoverInfo = TRUE, LabelCenters = FALSE, #TRUE
                                                       LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                                         c(2L, 20L, 0L)), class = c("package_version", "numeric_version"
                                                         ))), PanelId = c(ReducedDimensionPlot = 1L), PanelHeight = 400L,
                                                       PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
                                                       ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
                                                       ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
                                                       ColumnSelectionRestrict = FALSE, SelectionHistory = list())
    }
    
    if (length(input$mostres) <2) {
      return(tagList(
        br(),
        h4("At least 2 samples must be selected.", 
           style = "color: #d9534f; text-align: center; font-weight: bold; padding: 20px; border: 1px solid #ebccd1; background-color: #f2dede; border-radius: 4px;")
      ))
    } else {
      initial_panels[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "counts", CustomRows = TRUE,
                                                     CustomRowsText = selected_gens(), CapRowSelection = length(selected_gens()), ClusterRows = TRUE,
                                                     ClusterRowsDistance = "correlation", ClusterRowsMethod = "ward.D2",
                                                     DataBoxOpen = FALSE, ColumnData = input$cluster_var,#input$cluster_var,
                                                     RowData = character(0), CustomBounds = FALSE, LowerBound = -2L,
                                                     UpperBound = 2L, AssayCenterRows = TRUE, AssayScaleRows = TRUE,
                                                     DivergentColormap = "blue < white < red", ShowDimNames = c("Rows", "Columns"),
                                                     LegendPosition = "Right", LegendDirection = "Vertical",
                                                     VisualBoxOpen = FALSE,  NamesRowFontSize = 6, NamesColumnFontSize = 5,
                                                     ShowColumnSelection = FALSE, OrderColumnSelection = TRUE,
                                                     VersionInfo = list(iSEE = structure(list(c(2L, 20L, 0L)), class = c("package_version",
                                                                                                                         "numeric_version"))), PanelId = c(ComplexHeatmapPlot = 1L),
                                                     PanelHeight = 500L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
                                                     RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = TRUE,
                                                     ColumnSelectionSource = "ReducedDimensionPlot1",
                                                     RowSelectionRestrict = FALSE, ColumnSelectionRestrict = TRUE,
                                                     SelectionHistory = list())
    }
    
    initial_panels[["SampleAssayPlot1"]] <- new("SampleAssayPlot", Assay = "countsTMM", XAxis = "None", XAxisRowData = "",
                                                XAxisSampleSource = "---", XAxisSampleDynamicSource = FALSE,
                                                YAxisSampleSource = "---", YAxisSampleDynamicSource = FALSE,
                                                FacetRowByRowData = NA_character_, FacetColumnByRowData = NA_character_,
                                                ColorByRowData = "", ColorBySampleNameAssay = "counts", ColorByFeatureNameColor = "#FF0000",
                                                ShapeByRowData = NA_character_, SizeByRowData = NA_character_,
                                                TooltipRowData = character(0), FacetRowBy = "None", FacetColumnBy = "None",
                                                ColorBy = "None", ColorByDefaultColor = "#000000",
                                                ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
                                                ColorBySampleSource = "---",
                                                ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
                                                VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
                                                FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1,
                                                PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
                                                CustomLabels = FALSE, FontSize = 1,
                                                LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
                                                LabelCenters = FALSE, LabelCentersBy = NA_character_, LabelCentersColor = "black",
                                                VersionInfo = list(iSEE = structure(list(c(2L, 20L, 0L)), class = c("package_version",
                                                                                                                    "numeric_version"))), PanelId = c(SampleAssayPlot = 1L),
                                                PanelHeight = 400L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
                                                RowSelectionSource = "---", ColumnSelectionSource = "---",
                                                DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                                RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                                SelectionHistory = list())
    
    #Output final de l'iSEE
    output$isee_ui <- renderUI({
      req(dde_filtrat())
      dde_display <- isolate(dde_filtrat()) #aïllem la versió actual de l'objecte filtrat
      
      iSEE(
        dde_display,
        initial = initial_panels,
        appTitle = "Basic exploration of the data"
      )
    })
  }) #output$isee_ui
  
  #5. REACTIUS ISEE_UI2 _____________________________________________________________________________________________________________
  
  #Input per escollir el contrast
  output$contrast <- renderUI({
    req(dde2())
    req(input$tabs)
    if (input$tabs == "df_tab") { #si es canvia de tab
      
      contrastos <- names(dde2()@dea) #FFXvsWT i GEMvsWT
      
      #Retornem el selector només en aquest cas
      selectizeInput("contrast", "Select Contrast for DEA:", 
                     choices = contrastos, 
                     selected = contrastos[1],
                     multiple = TRUE,
                     options = list(
                       placeholder = 'Choose at least one contrast',
                       plugins = list("remove_button")))
    } else {
      #Si estem a 'basic_tab' o 'help_tab', no el mostra
      return(NULL)
    }
  })
  
  #quins gens surten al desplegable
  output$gens <- renderUI({
    req(dde2(), input$tabs)
    
    if (input$tabs == "df_tab") {
      #req(deg_from_dea())
      choices_gen <- deg_from_dea()
      label_text <- "Select DEGs to visualize:"
    } 
    
    if (input$tabs %in% c("basic_tab","df_tab")) {
      #req(top_genes())
      choices_gen <- top_genes()
      label_text <- "Select top variable genes to visualize:"
    } else {
      return(NULL)
    }
    
    if (is.null(choices_gen) || (length(choices_gen) == 0)) {
      choices_gen <- "No genes available"
    }
    
    selectizeInput(
      inputId = "gens", 
      label = label_text,
      choices = choices_gen, 
      multiple = TRUE, 
      options = list(
        placeholder = 'Select genes...',
        plugins = list('remove_button')
      )
    )
  })
  
  # Reactiu per filtrar els gens significatius a partir de l'slot dea
  deg_from_dea <- reactive({
    req(dde2())
    req(input$contrast)
    
    obj <- dde2()
    dea_list <- obj@dea 
    pval_cut <- if(!is.null(input$padj)) input$padj else 0.05
    lfc_cut <- if(!is.null(input$logFC)) input$logFC else 1.0
    
    contrast_actiu <- input$contrast[1]
    if (!contrast_actiu %in% names(dea_list)) return(character(0)) #si el contrast no es a la llista retorna null
    
    res_dea <- as.data.frame(dea_list[[contrast_actiu]]) #dde@dea[["FFXvsWT"]]
    
    col_p <- grep("padj$", colnames(res_dea), value = TRUE)[1]
    col_logfc <- grep("log2FoldChange$", colnames(res_dea), value = TRUE)[1]
    
    if (is.na(col_p) || is.null(col_p)) {
      message("ERROR: No s'ha trobat cap columna que acabi en 'padj'")
      return(character(0))
    }
    
    sig_rows <- res_dea[res_dea[[col_p]] < pval_cut & !is.na(res_dea[[col_p]]) & abs(res_dea[[col_logfc]]) >= lfc_cut, , drop = FALSE]
    
    if (nrow(sig_rows) == 0) return(character(0)) #si no hi ha gens retorna null
    
    sig_rows <- sig_rows[order(abs(sig_rows[[col_logfc]]), decreasing = TRUE), ]
    return(rownames(sig_rows))
  })
  
  #6. LLANÇAMENT DE L'ISEE + PANELLS DE LA *DEA/FEA TAB*________________________________________________________________________________
  
  output$isee_ui2 <- renderUI({
    req(dde2(), input$contrast, input$tabs)
    
    dde3 <- as(dde2(), "SummarizedExperiment")
    #definim rownames (amb el as(, se), es perden)
    noms_gens <- rownames(dde2())
    rownames(dde3) <- noms_gens
    
    #per la rowDataTable
    rd <- as.data.frame(rowData(dde2()))
    
    cols_interes <- c("Geneid", "Symbol", "Chr", "Start", "End", "Strand", "length", "Description")
    existing_info <- unlist(lapply(cols_interes, function(x) {
      grep(paste0("^", x, "$"), colnames(rd), ignore.case = TRUE, value = TRUE)
    }))
    
    # cols_interes <- c("Geneid", "Chr", "Start", "End", "Strand", "Description")
    # existing_info <- intersect(cols_interes, colnames(rd))
    
    # ---------- VOLCANO -----------
    contrast_principal <- input$contrast[1]
    
    col_x <- paste0(contrast_principal, "_log2FoldChange")
    col_y <- paste0(contrast_principal, "_padj")
    
    # Validar que les columnes existeixin
    if (!col_x %in% colnames(rd) || !col_y %in% colnames(rd)) {
      return(tagList(h4("Columns not found for this contrast.", style="color:red;")))
    }
    
    #fem servir dades originals (rd)
    val_logfc <- as.numeric(rd[[col_x]])
    val_padj <- as.numeric(rd[[col_y]])
    
    #thresholds
    pval_cut <- if(!is.null(input$padj)) input$padj else 0.05
    lfc_cut <- if(!is.null(input$logFC)) input$logFC else 1.0
    
    status <- rep("NS", nrow(rd))
    status[is.na(status)] <- "NS"
    status[val_logfc > lfc_cut & val_padj < pval_cut & !is.na(val_padj)] <- "Upregulated"
    status[val_logfc < -lfc_cut & val_padj < pval_cut & !is.na(val_padj)] <- "Downregulated"
    
    original_rd <- as.data.frame(rowData(dde2()))
    
    #per al venn
    obj_original <- dde2()
    tots_noms <- names(obj_original@dea)
    cols_venn <- grep(paste0("^(", paste(tots_noms, collapse="|"), ")_(padj|log2FoldChange)"), 
                      colnames(rd), value = TRUE)
    final_rd <- rd[, unique(c(existing_info, cols_venn)), drop = FALSE] 
    #new_rd$negLog10Padj <- -log10(val_padj)
    
    #creo nou df amb existing_info i li afegeixo les noves cols
    final_rd$logFC <- as.numeric(val_logfc)
    final_rd$padj <- as.numeric(val_padj)
    final_rd$PValue <- final_rd$padj #Perquè el Volcano agafi la col que toca
    final_rd$Significance <- factor(status, levels = c("Upregulated", "Downregulated", "NS"))
    
    final_rd[[col_x]] <- val_logfc
    final_rd[[col_y]] <- val_padj
    
    ordre_final <- c(existing_info, col_x, col_y, "Significance") 
    resta <- setdiff(colnames(final_rd), ordre_final)
    final_rd <- final_rd[, c(ordre_final, resta), drop = FALSE]
    
    # for(col in existing_info) {
    #   final_rd[[col]] <- rd[[col]]
    # }
    
    #rownames(final_rd)  <- noms_gens
    rowData(dde3) <- DataFrame(final_rd,check.names = FALSE)
    
    # ---------- HEATMAP ---------
    selected_cluster_var <- input$cluster_var
    if (is.null(selected_cluster_var) || selected_cluster_var == "") {
      selected_cluster_var <- colnames(colData(dde3))[1]
    }
    
    gens_actuals <- deg_from_dea()
    req(length(gens_actuals) > 0)
    
    gens_unio <- if (!is.null(input$gens) && length(input$gens) > 0) {
      input$gens
    } else {
      gens_actuals
    }
    
    
    gens_unio <- unique(na.omit(gens_unio))
    
    if (is.null(gens_actuals) || length(gens_actuals) == 0) {
      return(tagList(
        br(),
        h4("No significant genes (FDR < 0.05) found for contrast.", 
           style = "color: #888; text-align: left;")
      ))
    }
    
    mostres_sub <- if (!is.null(input$mostres) && length(input$mostres) > 0) {
      input$mostres
    } else {
      req(input$contrast)
      grups <- unlist(strsplit(input$contrast, "vs"))
      rownames(colData(dde3))[colData(dde3)$Cond %in% grups]
    }
    print(mostres_sub)
    
    #validació
    if (length(mostres_sub) == 0) mostres_sub <- rownames(colData(dde3))
    dde3 <- dde3[, mostres_sub]
    
    if (length(mostres_sub) < 2) {
      return(tagList(
        br(),
        h4("At least 2 samples must be selected.", 
           style = "color: #d9534f; text-align: center; font-weight: bold; padding: 20px; border: 1px solid #ebccd1; background-color: #f2dede; border-radius: 4px;")
      ))
    }
    
    gens_unio <- intersect(gens_unio, rownames(dde3))
    #eliminem gens amb variància 0, sino heatmap peta
    mat <- assay(dde3, "countsTMM")[gens_unio, mostres_sub, drop = FALSE]
    vars <- apply(mat, 1, var, na.rm = TRUE)
    gens_valids <- names(vars)[!is.na(vars) & vars > 0]
    
    
    # dde3 <- dde3[gens_valids, mostres_sub]
    # gens_finals <- intersect(gens_unio, gens_valids)
    
    # -------- Funció Venn Diagram-------
    venn_fun <- createCustomPlot(
      function(dde, rows, columns) {
        pval_cut <- if(!is.null(input$padj)) input$padj else 0.05
        lfc_cut <- if(!is.null(input$logFC)) input$logFC else 1.0
        rd <- as.data.frame(rowData(dde))
        
        #en cas q el format de cols sigui "_padj"
        cols_padj <- grep(".+padj$", colnames(rd), value = TRUE)
        #".+" obliga a que hi hagi algun text abans de "padj" ("FFX_")
        
        if (length(cols_padj) < 2) {
          totes_p <- grep("padj", colnames(rd), value = TRUE)
          cols_padj <- totes_p[totes_p != "padj"]
          return(ggplot() +
                   annotate("text", x=0, y=0, label=paste("Select at least 2 contrasts to see intersections.\nAvailable columns:",
                                                          paste(colnames(rd)[grep("padj", colnames(rd))], collapse=", "))) +
                   theme_void())
        }
        
        lists <- list()
        
        for (col in cols_padj) {
          contrast_name <- gsub("(.original_object|.padj|_padj)", "", col)
          print(contrast_name)
          #logFC corresponent a aquest contrast
          col_fc <- grep(paste0(contrast_name, ".*log2FoldChange"), colnames(rd), value = TRUE)[1]
          #filtrem si existeix lfc
          if (!is.na(col_fc)) {
            gens_sig <- rownames(rd)[which(rd[[col]] < pval_cut & abs(rd[[col_fc]]) > lfc_cut)]
          } else {
            #i no hi ha lfc, filtrem només per p-adj
            gens_sig <- rownames(rd)[which(rd[[col]] < pval_cut)]
          }
          
          if (length(gens_sig) > 0) {
            lists[[contrast_name]] <- gens_sig
          }
        }
        if (length(lists) < 2) {
          return(ggplot() + 
                   annotate("text", x=0, y=0, label="There are not enough significant genes found.") + 
                   theme_void())
        }
        
        ggvenn::ggvenn(lists[1:min(3, length(lists))], 
                       fill_color = c("#00AFBB", "#E7B800", "#FC4E07", "#7E4E90")[1:length(lists)], 
                       stroke_size = 1,
                       set_name_size = 6,
                       text_size = 6
        ) + labs(title = "Significant genes")
      },
      restrict = NULL,
      className = "VennDiagram1",
      fullName = "Venn Diagram"
    )
    
    # ------------ Funció FEA DotPlot ----------------
    fea_dotplot_fun <- createCustomPlot(
      function(dde, rows, columns) {
        req(dde2())
        obj <- dde2()
        
        active_contrast <- input$contrast[1]
        fea_list <- obj@fea[[active_contrast]]$original_object
        
        if (is.null(fea_list)) {
          return(ggplot() + annotate("text", x=0, y=0, label="No GSEA data") + theme_void())
        }
        
        df_fea <- as.data.frame(fea_list)
        
        if (nrow(df_fea) == 0) {
          return(ggplot() + annotate("text", x=0, y=0, label="No enrichment found") + theme_void())
        }
        
        #top 20 termes per p val ajustat
        df_sig <- df_fea[df_fea$padj < 0.05 & !is.na(df_fea$padj), ]
        
        if (nrow(df_sig) == 0) {
          return(ggplot() + annotate("text", x=0, y=0, label="No significant pathways found") + theme_void())
        }
        #ordenem pel valor absolut del NES
        top_fea <- head(df_sig[order(abs(df_sig$NES), decreasing = TRUE), ], 30)
        
        #es treu el prefix HALLMARK_
        top_fea$pathway <- gsub("HALLMARK_", "", top_fea$pathway)
        top_fea$pathway <- gsub("_", " ", top_fea$pathway)
        
        library(ggplot2) #adaptat a gsea
        ggplot(top_fea, aes(x = NES, y = reorder(pathway, NES))) +
          geom_point(aes(size = size, color = padj)) +
          scale_color_gradient(low = "red", high = "blue") +
          theme_bw() +
          labs(title = paste("GSEA: Enriched Hallmark Pathways,", active_contrast), x = "Normalized Enrichment Score (NES)", y = "ID") +
          theme(axis.text.y = element_text(size = 8))
      },
      className = "FEADotPlot",
      fullName = "FEA Dot Plot"
    )
    
    # ---------- Funció FEA BarPlot ------------
    fea_barplot_fun <- createCustomPlot(
      function(dde, rows, columns) {
        req(dde2())
        obj <- dde2()
        
        active_contrast2 <- input$contrast[1]
        fea_list2 <- obj@fea[[active_contrast2]]$original_object
        
        if (is.null(fea_list2)) {
          return(ggplot() + annotate("text", x=0, y=0, label="No GSEA data") + theme_void())
        }
        df_fea2 <- as.data.frame(fea_list2)
        
        #top 20
        df_sig2 <- df_fea2[df_fea2$padj < 0.05 & !is.na(df_fea2$padj), ]
        
        if (nrow(df_sig2) == 0) {
          return(ggplot() + annotate("text", x=0, y=0, label="No significant pathways found") + theme_void())
        }
        #ordenem pel valor absolut del NES
        top_fea2 <- head(df_sig2[order(abs(df_sig2$NES), decreasing = TRUE), ], 30)
        
        #traiem hallmark
        top_fea2$pathway <- gsub("HALLMARK_", "", top_fea2$pathway)
        top_fea2$pathway <- gsub("_", " ", top_fea2$pathway)
        
        #creem columna de color segons si el NES és positiu o negatiu
        #top_fea2$Direction <- ifelse(top_fea2$NES > 0, "Enriched in GEM", "Enriched in WT")
        
        library(ggplot2)
        ggplot(top_fea2, aes(x = NES, y = reorder(pathway, NES), fill = padj)) +
          geom_bar(stat = "identity") +
          scale_fill_gradient(low = "red", high = "blue") +
          theme_bw() +
          labs(title = paste("GSEA: Hallmark Normalized Enrichment Score,", active_contrast2),
               x = "Normalized Enrichment Score (NES)",
               y = "ID") +
          theme(axis.text.y = element_text(size = 8),
                legend.position = "bottom") +
          geom_vline(xintercept = 0, linetype = "solid", color = "black")
      },
      className = "FEABarPlot",
      fullName = "FEA Bar Plot"
    )
    gens_heatmap <- intersect(gens_unio, gens_valids)
    
    # ------------- PANELLS --------------
    isolate({
      initial_panels2 <- list()
      
      
      heat_cols <- unique(c("Cond", selected_cluster_var)) 
      dims <- if (length(gens_heatmap) < 100) c("Rows", "Columns") else "Columns"
      initial_panels2[["ComplexHeatmapPlot2"]] <- new("ComplexHeatmapPlot", Assay = "countsTMM", 
                                                      CustomRows = FALSE,
                                                      CustomRowsText = paste(gens_unio, collapse = "\n"), #gens_heatmap
                                                      CapRowSelection = length(gens_unio),
                                                      
                                                      ClusterRows = TRUE,ClusterRowsDistance = "correlation", ClusterRowsMethod = "ward.D2",
                                                      OrderColumnSelection = TRUE,
                                                      DataBoxOpen = FALSE, 
                                                      ColumnData = heat_cols,
                                                      RowData = character(0), CustomBounds = FALSE,LowerBound = -2L,
                                                      UpperBound = 2L, AssayCenterRows = TRUE, AssayScaleRows = TRUE,
                                                      DivergentColormap = "blue < white < red", ShowDimNames = dims,
                                                      LegendPosition = "Right", LegendDirection = "Vertical",
                                                      VisualBoxOpen = TRUE,  NamesRowFontSize = 6, NamesColumnFontSize = 8,
                                                      ShowColumnSelection = FALSE,
                                                      VersionInfo = list(iSEE = structure(list(c(2L, 20L, 0L)), class = c("package_version",
                                                                                                                          "numeric_version"))), PanelId = c(ComplexHeatmapPlot = 1L),
                                                      PanelHeight = 400L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
                                                      RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = TRUE,
                                                      ColumnSelectionSource = "---",
                                                      ColumnSelectionRestrict = TRUE,
                                                      SelectionHistory = list())
      
      
      initial_panels2[["VolcanoPlot1"]] <- new("VolcanoPlot",
                                               XAxis = "Row data", 
                                               XAxisRowData = "logFC",
                                               YAxis = "PValue",
                                               PValueThreshold = pval_cut, 
                                               LogFCThreshold = lfc_cut,
                                               VisualBoxOpen = TRUE, VisualChoices = "Color", 
                                               ColorBy = "Row data",
                                               ColorByRowData = "Significance",
                                               #RowSelectionSource = "ComplexHeatmapPlot1",
                                               PanelWidth = 6L,
                                               PanelHeight = 400L
      )
      
      if (!is.null(input$contrast) && length(input$contrast) >= 2 && length(input$contrast) <= 4) {
        obj_venn <- venn_fun()
        
        slot(obj_venn, "DataBoxOpen") <- FALSE
        slot(obj_venn, "SelectionBoxOpen") <- FALSE
        
        slot(obj_venn, "PanelId") <- 1L 
        slot(obj_venn, "PanelWidth") <- 6L
        slot(obj_venn, "PanelHeight") <- 500L
        
        initial_panels2[["VennDiagram1"]] <- obj_venn
      }
      
      
      gen_seleccionat <- if (length(gens_unio) > 0) as.character(gens_unio[1]) else "" 
      #ordre+cols que es mostren
      tots_els_contrastos <- names(dde2()@dea)
      altres_contrastos <- setdiff(tots_els_contrastos, contrast_principal)
      patro_amagar <- paste0("^(", paste(altres_contrastos, collapse = "|"), ")_")
      cols_altres <- grep(patro_amagar, colnames(final_rd), value = TRUE)
      cols_no_visibles <- unique(c(cols_altres, "logFC", "padj", "PValue"))
      
      initial_panels2[["RowDataTable1"]] <- new ("RowDataTable", 
                                                 Selected = gen_seleccionat,
                                                 Search = "", 
                                                 HiddenColumns = intersect(cols_no_visibles, colnames(final_rd)), 
                                                 PanelWidth = 12L,
                                                 PanelHeight = 400L)
      
      
      # gens_string <- paste(gens_a_mostrar, collapse = "\n")
      # initial_panels2[["AggregatedDotPlot1"]] <- new("AggregatedDotPlot",
      #                                                Assay = "countsTMM",
      #                                                ColumnDataLabel= selected_cluster_var,
      #                                                CustomRows = TRUE,
      #                                                CustomRowsText = gens_string,
      #                                                PanelHeight = 400L,
      #                                                PanelWidth = 6L
      # )
      
      if (!is.null(input$contrast)) {
        obj_bar <- fea_barplot_fun()
        
        slot(obj_bar, "DataBoxOpen") <- FALSE
        slot(obj_bar, "SelectionBoxOpen") <- FALSE
        slot(obj_bar, "PanelId") <- 1L
        slot(obj_bar, "PanelWidth") <- 6L  
        slot(obj_bar, "PanelHeight") <- 500L
        
        initial_panels2[["FEABarPlot"]] <- obj_bar
      }
      
      if (!is.null(input$contrast)) {
        obj_fea <- fea_dotplot_fun()
        
        slot(obj_fea, "DataBoxOpen") <- FALSE
        slot(obj_fea, "SelectionBoxOpen") <- FALSE
        
        slot(obj_fea, "PanelId") <- 1L
        slot(obj_fea, "PanelWidth") <- 6L
        slot(obj_fea, "PanelHeight") <- 500L
        
        initial_panels2[["FEADotPlot"]] <- obj_fea 
      }
      
      
      ecm <- ExperimentColorMap(
        all_discrete = list(
          rowData = function(n) {
            if (n == 3) {
              return(c("Upregulated" = "red", "Downregulated" = "blue", "NS" = "grey"))
            } else {
              return(viridis::viridis(n))
            }
          },
          colData = function(n) { rainbow(n) },
          assays = function(n) { viridis::viridis(n) }
        )
      )
      
      
      isee_obj <- iSEE(dde3, initial = initial_panels2, colormap = ecm,
                       appTitle = "Visualization of DEA and FEA results")
      tags$div(style = "height: 800px;", isee_obj)
    }) # tanca isolate
    
  }) # tanca renderUI
  
  # 7. TAB GENOMIC EXPLORATION + PANELLS______________________________________________________________________
  
  output$isee_ui3 <- renderUI({
    req(dde2())
    
    dde_display <- dde2()
    
    # llegim el cromosoma seleccionat 
    chr_sel <- chr_sel_reactive() 
    
    initial_panels <- list()
    
    # ------------- PANELLS ---------------
    
    # Karyoplot
    initial_panels[["KaryoPlot"]] <- new("KaryoPlot",
                                         PanelHeight = 400L,
                                         PanelWidth  = 6L,
                                         RowSelectionSource = "RowDataTable1",
                                         ZoomChr = FALSE)
    
    # Zoom - reacciona al chr_selector del sidebar
    initial_panels[["KaryoPlot2"]] <- new("KaryoPlot",
                                          PanelHeight = 400L,
                                          PanelWidth  = 6L,
                                          RowSelectionSource = "RowDataTable1",
                                          ZoomChr = TRUE)
    
    # GRanges table
    initial_panels[["RowDataTable1"]] <- new("RowDataTable",
                                             Search = chr_sel,
                                             HiddenColumns = c("FFXvsWT_log2FoldChange", "FFXvsWT_pvalue", 
                                                               "FFXvsWT_padj", "GEMvsWT_log2FoldChange", "GEMvsWT_pvalue", 
                                                               "GEMvsWT_padj", "Chr", "Start", "End", "Strand", "Geneid"
                                             ), VersionInfo = list(iSEE = structure(list(c(2L, 20L, 0L
                                             )), class = c("package_version", "numeric_version"))), PanelId = c(RowDataTable = 1L), 
                                             PanelHeight = 400L, PanelWidth = 12L, SelectionBoxOpen = FALSE, 
                                             RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                             DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                             RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                             SelectionHistory = list())
    
    
    # Circus Plot
    initial_panels[["CircosPlot1"]] <- new("CircosPlot",
                                           PanelHeight = 400L,
                                           PanelWidth  = 6L,
                                           RowSelectionSource = "RowDataTable1",
                                           Contrast1 = "FFXvsWT", 
                                           Contrast2 = "GEMvsWT")
    
    
    iSEE(
      dde_display,
      initial = initial_panels,
      appTitle = "Genomic Exploration")
  })
  
} #server
shinyApp(ui, server)
