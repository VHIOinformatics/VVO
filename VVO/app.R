library(shiny)
library(shinyjs)
library(shinydashboard)
library(rintrojs)

library(iSEE)
library(iSEEu)
library(iSEEde) # per a les subclasses (VolcanoPlot)

library(S4Vectors)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(DeeDeeExperiment)

library(BasicPlots) # Funcions propies de VHIOinformatics
library(BasicFunctions) # Funcions propies de VHIOinformatics

library(RColorBrewer)
library(ggplot2)
library(edgeR)
library(ggvenn)
library(ggrepel)
library(karyoploteR)
library(circlize)

options(shiny.maxRequestSize = 50 * 1024^2)

# classe panell KaryoPlot

source("panels/KaryoPlot.R")

# classe panell CircusPlot

source("panels/CircusPlot.R")

# ---------- UI ---------------

ui <- dashboardPage(
  
  dashboardHeader(title = "VHIO's VISUAL OMICS (VVO)", titleWidth = 240),
  
  dashboardSidebar(
    useShinyjs(),
    width = 240,
    sidebarMenu(
      id = "tabs",
      menuItem("Basic Exploration", tabName = "basic_tab", icon = icon("eye")),
      menuItem("DEA and FEA results", tabName = "df_tab", icon = icon("chart-line")),
      menuItem("Help", tabName = "help_tab", icon = icon("question-circle"))
    ),
    tags$hr(style = "border-top: 2px solid white; margin-top:4px; margin-bottom:4px;"),
    
    # Icona d'ajuda per a les dades requerides
    div(
      style = "display: flex; align-items: center; gap: 8px; margin-bottom: 10px;",
      tags$label("Input data:", style = "font-weight: bold; flex: 1; margin: 0;"),
      actionLink("help_required", icon("circle-question"), title = "Go to Help", style = "color: #0066cc;")
    ),
    
    fileInput("dde_file", "Upload .rds data:(DeeDeeExperiment)", accept = ".rds"),
    
    # Icona d'ajuda per als controls del sidebar
    div(
      style = "display: flex; align-items: center; gap: 8px; margin-top: 15px; margin-bottom: 10px;",
      tags$label("Parameters:", style = "font-weight: bold; flex: 1; margin: 0;"),
      actionLink("help_controls", icon("circle-question"), title = "Go to Help", style = "color: #0066cc;")
    ),
    
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
      tabItem(tabName = "help_tab",
              HTML('
    <h3 style="color:#2c3e50;">User Guide - VHIO Visual Omics (VVO)</h3>
    <p style="color:#555;">VVO lets you explore bulk RNA-seq data through two linked modules: <b>Basic Exploration</b> (QC and clustering) and <b>DEA and FEA results</b> (differential expression and functional enrichment), built on <a href="https://bioconductor.org/packages/iSEE" target="_blank">iSEE</a>.</p>
    <hr>

    <h4 id="section-required" style="color:#2c3e50;">1. Required data format</h4>
    <div style="background-color:#f9f9f9; padding:10px;">
      <ul>
        <li>Upload a single <code>.RDS</code> file containing a <b><code>DeeDeeExperiment</code></b> object (extends <code>SummarizedExperiment</code>/<code>SingleCellExperiment</code>, storing DEA results in <code>@dea</code> and enrichment results in <code>@fea</code>).</li>
        <li>The object must have a <code>Cond</code> column in <code>colData</code>, identifying each sample\'s experimental group.</li>
        <li>The <code>assays</code> slot must include raw <code>counts</code>. VVO filters low-expression genes, normalizes with TMM, and computes a PCA automatically.</li>
        <li>Contrasts shown in the sidebar come from the names of <code>@dea</code>. For each contrast, VVO reads <code>&lt;contrast&gt;_padj</code> and <code>&lt;contrast&gt;_log2FoldChange</code> columns in <code>rowData</code>.</li>
        <li>If <code>@fea</code> has GSEA results for a contrast, the FEA plots will show them.</li>
      </ul>
    </div>

    <h4 id="section-genomic" style="color:#2c3e50;">2. With or without genomic position (rowRanges)</h4>
    <div style="background-color:#f9f9f9; padding:10px;">
      <ul>
        <li>VVO checks automatically whether the object has valid genomic coordinates (<code>rowRanges</code>) for each gene.</li>
        <li><b>With genomic position:</b> three extra panels appear on the DEA/FEA tab - <b>KaryoPlot</b> (whole genome), <b>KaryoPlot (zoom)</b> (one chromosome, chosen in the sidebar), and <b>CircosPlot</b> (circular view, up to two contrasts compared at once).</li>
        <li><b>Without genomic position:</b> those three panels are hidden, everything else works the same. No action needed from you.</li>
        <li>If you expect the genomic panels and don\'t see them, check that <code>rowRanges()</code> is not empty and has valid <code>start</code>/<code>end</code> values.</li>
      </ul>
    </div>

    <h4 id="section-controls" style="color:#2c3e50;">3. Sidebar controls</h4>
    <div style="background-color:#f9f9f9; padding:10px;">
      <ul>
        <li><b>Upload .rds data:</b> your <code>DeeDeeExperiment</code> file.</li>
        <li><b>Select Contrast for DEA</b> (DEA/FEA tab): pick one or more contrasts. 2-4 contrasts enable the Venn diagram; the first one drives the volcano plot, FEA plots and gene tables.</li>
        <li><b>Clustering by:</b> <code>colData</code> column used to color samples and group the dendrogram.</li>
        <li><b>Number of genes to display</b> (Basic Exploration only): how many top variable genes feed the heatmap and gene selector.</li>
        <li><b>Select samples:</b> which samples to include. On DEA/FEA, pre-filtered to the selected contrast(s).</li>
        <li><b>Select genes:</b> restrict plots to specific genes.</li>
        <li><b>P adjusted Value / Log2 Fold Change thresholds</b> (DEA/FEA tab): significance cutoffs used everywhere (volcano, Venn, karyotype/circos).</li>
        <li><b>Zoom to chromosome</b> (only if genomic position available): chromosome shown in the zoomed KaryoPlot.</li>
      </ul>
    </div>

    <h4 id="section-interaction" style="color:#2c3e50;">4. How panels react to each other</h4>
    <div style="background-color:#f9f9f9; padding:10px;">
      <p>Panels are linked: selecting something in one panel (a point, a row, a sample) updates other panels automatically. In Basic Exploration:</p>
      <ul>
        <li>Selecting points in the <b>PCA</b> plot updates <b>Library Size</b>, <b>Dendrograma</b> and the <b>Heatmap</b> to show only those samples.</li>
      </ul>
      <p>In DEA and FEA results:</p>
      <ul>
        <li>Clicking a gene (row) in the <b>RowDataTable</b> highlights that gene in red on the <b>KaryoPlot</b> / <b>KaryoPlot (zoom)</b> (only with genomic position).</li>
        <li>Filtering for Upregulated or Downregulated genes in the <b>RowDataTable</b> highlights that genes on the <b>VolcanoPlot</b>.</li>
        <li>Selecting samples in the <b>ColumnDataTable</b> updates the <b>Heatmap</b> to show only those samples.</li>
        <li>The <b>CircosPlot</b> does not react to row/column selections — it always shows the DEGs of the two contrasts assigned to it (<code>Contrast1</code>/<code>Contrast2</code>), based on the padj/logFC thresholds.</li>
      </ul>
      <p style="color:#888;">Tip: to see a specific gene highlighted on the genome, search or click it in the RowDataTable first.</p>
    </div>

    <h4 style="color:#2c3e50;">5. Basic Exploration tab</h4>
    <div style="background-color:#f9f9f9; padding:10px;">
      <ul>
        <li><b>Library Size:</b> total reads (millions) per sample, colored by group. Quick QC check.</li>
        <li><b>Dendrograma:</b> hierarchical clustering of samples (Euclidean distance, Ward.D2), colored by group.</li>
        <li><b>PCA:</b> principal component plot of samples, colored by group.</li>
        <li><b>Heatmap:</b> expression of the selected genes, scaled by row, samples annotated by group.</li>
        <li><b>SampleAssayPlot:</b> expression values per sample for individual genes.</li>
      </ul>
      <p style="color:#888;">At least 2 samples are required.</p>
    </div>

    <h4 style="color:#2c3e50;">6. DEA and FEA results tab</h4>
    <div style="background-color:#f9f9f9; padding:10px;">
      <ul>
        <li><b>Volcano Plot:</b> logFC vs. p-value for the first selected contrast, colored by significance.</li>
        <li><b>RowDataTable:</b> gene-level table with annotation and DEA stats.</li>
        <li><b>ColumnDataTable:</b> sample metadata table.</li>
        <li><b>Heatmap:</b> expression of the DEGs (or selected genes), annotated by group.</li>
        <li><b>Venn Diagram:</b> overlap of significant genes across 2–4 selected contrasts.</li>
        <li><b>FEA Dot Plot / Bar Plot:</b> top significant Hallmark pathways (GSEA) for the first contrast, by NES and padj.</li>
        <li><b>KaryoPlot / KaryoPlot (zoom) / CircosPlot</b> — only if genomic position is available: show where the DEGs are located on the genome.</li>
      </ul>
      <p style="color:#888;">At least 2 samples and a contrast with significant genes are required.</p>
    </div>
    
      </ul>
    </div>
  ')
      )
    )
  )
)

# detectar genoma
get_genome <- function(dde) {
  n_chrs <- length(unique(as.character(seqnames(rowRanges(dde)))))
  if(n_chrs < 23) "mm10" else "hg38"
}

# --------------- SERVER ---------------------

server <- function(input, output, session) {
  
  # Icones d'ajuda navegació a la secció corresponent del Help tab
  observeEvent(input$help_required, {
    updateTabItems(session, "tabs", selected = "help_tab")
    shinyjs::delay(300, shinyjs::runjs("document.getElementById('section-required').scrollIntoView({behavior: 'smooth'});"))
  })
  
  observeEvent(input$help_controls, {
    updateTabItems(session, "tabs", selected = "help_tab")
    shinyjs::delay(300, shinyjs::runjs("document.getElementById('section-controls').scrollIntoView({behavior: 'smooth'});"))
  })
  
  observeEvent(input$help_interaction, {
    updateTabItems(session, "tabs", selected = "help_tab")
    shinyjs::delay(300, shinyjs::runjs("document.getElementById('section-interaction').scrollIntoView({behavior: 'smooth'});"))
  })
  
  # Valor reactiu per guardar el cromosoma seleccionat.
  chr_sel_reactive <- reactiveVal("chr1") 
  
  # Observador que escolta canvis a input$chr_sel.
  observe({
    req(input$chr_sel)
    chr_sel_reactive(input$chr_sel)
  })
  
  
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
    keep <- filterByExpr(counts, group = colData(obj)[["Cond"]]) 
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
    n_pcs <- min(4, ncol(pca_scores_makePCA))
    reducedDims(obj)$PCA <- pca_scores_makePCA[, 1:n_pcs] 
    
    return(obj)
  })
  
  # Comprova si l'objecte té rowRanges vàlids (per mostrar/amagar la tab genòmica)
  te_rowRanges <- reactive({
    req(dde2())
    obj <- dde2()
    
    tryCatch({
      rr <- rowRanges(obj) # agafem el slot rowRanges de l'objecte
      length(rr) > 0 && !all(is.na(start(rr))) # comprova que el GRanges no estigui buit
    }, error = function(e) {
      FALSE  # si peta per qualsevol motiu, assumim que no té rowRanges 
    })
  })
  
  
  #2. OUTPUTS DE LA UI________________________________________________________________________________________________
  
  # Variable per agrupar (metadata) --> cluster_var
  output$cluster_var <- renderUI({
    req(dde2())
    dde <- dde2()
    if (input$tabs %in% c("basic_tab","df_tab")) {
      cols_fil <- colnames(colData(dde))
      default_sel <- if ("Cond" %in% cols_fil) "Cond" else cols_fil[1]
      selectizeInput("cluster_var","Clustering by:",
                     choices = cols_fil,
                     selected = default_sel)
    }
  })
  
  # Input númerico genes a graficar --> num_genes
  output$num_genes <- renderUI({
    req(dde2(), input$tabs)
    if (input$tabs == "basic_tab") {
      numericInput(
        "num_genes",
        "Number of genes to display (min 2 - max 50):",
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
        # suporta tant "vs" com ".vs." com a separador
        if (any(grepl(".vs.", input$contrast, fixed = TRUE))) split_by <- ".vs." else split_by <- "vs"
        grups <- unlist(strsplit(input$contrast, split_by, fixed = TRUE))
        grups_finals <- unique(grups) #per eliminar duplicats
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
    req(dde2(), input$tabs == "df_tab", te_rowRanges())
    
    cromosomes <- if(genome_sp()=="mm10") paste0("chr",c(1:19,"X","Y")) else paste0("chr",c(1:22,"X","Y"))
    
    selectizeInput(
      inputId = "chr_sel",
      label = "Zoom to chromosome:",
      choices = cromosomes,
      selected = "chr1",
      multiple = FALSE
    )
  })
  
  #3. REACTIUS COMPARTITS________________________________________________________________________________________________________
  
  # Detecció de tipus de genoma (humà o ratolí)
  genome_sp <- reactive({
    req(dde2())
    get_genome(dde2())
  })
  
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
  
  # ----------- Funció QC ----------------
  
  QC_fun <- function(dde, rows = NULL, columns = NULL) {
    selected_obj <- get_selected_dde(dde, columns)
    samples_actuals <- as.character(colnames(selected_obj))
    # sample.totals <- apply(counts(dde), 2, sum) 
    # sample_order <- colnames(dde)
    
    all_counts <- assay(isolate(dde2()), "counts")[, samples_actuals, drop = FALSE]
    sample.totals <- colSums(all_counts)
    
    cluster_raw <- as.character(colData(selected_obj)[[input$cluster_var]])
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
    # elements a la llegenda:
    num_grups <- length(unique(sample.totals.df$cluster))
    m_sel <- length(input$mostres)
    
    # mida_text <- if(num_grups > 40) 4 else if(num_grups > 20) 7 else 9
    mida_quadrat <- if(num_grups > 40) 2.5 else 3 #if(num_grups > 10) 3 else 5
    mida_titol <- if(num_grups > 20) 8 else 10
    mida_text <- if(num_grups > 40) 5 else if(num_grups > 20) 7 else 9
    mida_eix_x <- if(m_sel > 50) 5 else if(m_sel > 30) 7 else 9
    grid_x <- if(remove_grid) element_blank() else element_line()
    num_grups <- length(unique(sample.totals.df$cluster))
    mostra_llegenda <- num_grups <= 20   # NUEVO: umbral a partir del cual se oculta
    
    p <- ggplot(data=sample.totals.df, aes(x=sample, y=total)) + 
      geom_bar(aes(fill = total), stat = "identity") +
      geom_point(aes(colour = cluster), y = -Inf, alpha = 0) +
      scale_colour_manual(values = color_palette, name = input$cluster_var) +
      guides(
        colour = if (mostra_llegenda) {
          guide_legend(override.aes = list(shape = 15, size = mida_quadrat, alpha = 1))
        } else {
          "none"     # oculta llegenda si hi ha massa grups
        }
      ) +
      theme_bw() +
      theme(
        panel.grid.major.x = grid_x,
        panel.grid.minor.x = grid_x,
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = label_colors, size = mida_eix_x),
        legend.text = element_text(size = mida_text),
        legend.title = element_text(size = mida_titol)
      ) +
      ylab("Million reads") +
      xlab(NULL) +
      labs(fill = NULL) 
    
    return(p)   
  }
  
  # ------------- Dendograma i helper functions --------------
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
    mida_llegenda <- 0.8
    
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
    
    if (is.null(conditions) || length(conditions) == 0) {
      conditions <- rep("Group1", ncol(estimates))
    }
    
    cond_levels <- ordered_levels(conditions)
    conditions <- factor(conditions, levels = cond_levels)
    # conditions <- as.character(conditions)
    # cond_levels <- sort(unique(conditions))
    
    mostra_llegenda <- length(cond_levels) <= 20 
    
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
    
    if (mostra_llegenda) {
      layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
    } else {
      layout(matrix(1, nrow = 1))   # si no mostra llegenda surt el plot a tot el panell
    }
    
    par(mar = c(7, 4, 4, 1))
    clust_col <- colorCluster(clust, sample_colors, ce = mida_labels)
    plot(clust_col, main = title, xlab = xlab, sub = "")
    
    if (mostra_llegenda) {
      par(mar = c(7, 0, 4, 1))
      plot.new()
      legend("center",
             legend = cond_levels, fill = color_map, title = "Groups",
             cex = mida_llegenda, bty = "n")
    }
    
    layout(1)
    return(recordPlot())
  }
  
  # Función dendrograma
  dendro_fun <- function(dde, cluster_by, genes_use, rows = NULL, columns = NULL) {
    
    if (is.null(cluster_by) || !(cluster_by %in% colnames(colData(dde)))) {
      cluster_by <- colnames(colData(dde))[1]  # fallback a la primera columna
    }
    
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
  
  #4. BASIC EXPLORATION TAB + Panells ________________________________________________________________________________
  
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
    
    # d'aquesta manera treiem el boc Data Parameters que estava buit
    setMethod(".hideInterface", "QCPlot1", function(x, field) {
      if (field %in% c("DataBoxOpen")) TRUE  # amaga aquest slot
      else callNextMethod() # la resta, comportament per defecte
    })
    
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
      fullName = "Dendrogram"
    )(
      PanelHeight = 400L,
      PanelWidth = 6L,
      ColumnSelectionDynamicSource = TRUE,
      ColumnSelectionSource = "ReducedDimensionPlot1",
      RowSelectionRestrict = FALSE,
      ColumnSelectionRestrict = TRUE,
      SelectionHistory = list()
    )
    
    # d'aquesta manera treiem el boc Data Parameters que estava buit
    setMethod(".hideInterface", "DendroPlot1", function(x, field) {
      if (field %in% c("DataBoxOpen")) TRUE  # amaga aquest slot
      else callNextMethod() # la resta, comportament per defecte
    })
    
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
    iSEE(
      dde_display,
      initial = initial_panels,
      appTitle = "Basic exploration of the data"
    )
    
  }) #output$isee_ui
  
  #5. REACTIUS DEA/FEA TAB _____________________________________________________________________________________________________________
  
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
      message("ERROR: There is no column that finish with 'padj'")
      return(character(0))
    }
    
    sig_rows <- res_dea[res_dea[[col_p]] < pval_cut & !is.na(res_dea[[col_p]]) & abs(res_dea[[col_logfc]]) >= lfc_cut, , drop = FALSE]
    
    if (nrow(sig_rows) == 0) return(character(0)) #si no hi ha gens retorna null
    
    sig_rows <- sig_rows[order(abs(sig_rows[[col_logfc]]), decreasing = TRUE), ]
    return(rownames(sig_rows))
  })
  
  #6. DEA/FEA TAB + Panells________________________________________________________________________________
  
  # ---------- Karyoplot -----------
  
  # Renderitza el panel 
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
                
                #  reactivitat als inputs del sidebar
                input$contrast
                input$padj
                input$logFC
                
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
                
                # Definim aquí les variables que necessita el plot (accés a input des del servidor)
                contrasts_sel   <- input$contrast
                n_con           <- length(contrasts_sel)
                pval_cut        <- if (!is.null(input$padj))  input$padj  else 0.05
                lfc_cut         <- if (!is.null(input$logFC)) input$logFC else 1.0
                rd              <- as.data.frame(SummarizedExperiment::rowData(se))
                contrast_colors <- c("#00AFBB", "#7E4E90", "#E7B800", "#FC4E07")
                
                # Dibuixem el cariotip base (global o zoom segons el panel)
                if (is_zoom) {
                  chr_sel <- chr_sel_reactive()
                  pp <- karyoploteR::getDefaultPlotParams(plot.type = 1)
                  pp$topmargin      <- 0
                  pp$bottommargin   <- 150
                  pp$ideogramheight <- 25
                  pp$leftmargin     <- 0.15
                  kp <- karyoploteR::plotKaryotype(
                    genome = genome_sp(), chromosomes = chr_sel,
                    plot.type = 1, plot.params = pp
                  )
                } else {
                  kp <- karyoploteR::plotKaryotype(genome = genome_sp(), plot.type = 1)
                }
                
                # Dividim l'eix Y en bandes iguals: 1 per grisos + 1 per cada contrast
                # Gris molt prim, contrastos ocupen la resta de l'eix Y
                grey_h <- 0.05  # altura fixa petita per al gris
                band_h <- if (n_con > 0) (1 - grey_h) / n_con else (1 - grey_h)
                
                # Banda 0: gris prim al fons
                karyoploteR::kpRect(kp, data = gr,
                                    y0     = 0,
                                    y1     = grey_h,
                                    col    = "#aaaaaa22",
                                    border = "#aaaaaa88",
                                    lwd    = 0.3)   # molt fi
                
                # Bandes de contrastos: comencen on acaba el gris
                if (n_con > 0) {
                  legend_labels <- c(); legend_cols <- c()
                  for (i in seq_along(contrasts_sel)) {
                    con    <- contrasts_sel[[i]]
                    col_p  <- paste0(con, "_padj")
                    col_fc <- paste0(con, "_log2FoldChange")
                    if (!col_p %in% colnames(rd) || !col_fc %in% colnames(rd)) next
                    sig    <- rownames(rd)[!is.na(rd[[col_p]]) & rd[[col_p]] < pval_cut & abs(rd[[col_fc]]) >= lfc_cut]
                    gr_sig <- gr[intersect(sig, names(gr))]
                    col_i  <- contrast_colors[((i - 1) %% length(contrast_colors)) + 1]
                    
                    y0_i <- grey_h + (i - 1) * band_h
                    y1_i <- grey_h + i * band_h
                    
                    if (length(gr_sig) > 0) {
                      karyoploteR::kpRect(kp, data = gr_sig,
                                          y0     = y0_i,
                                          y1     = y1_i,
                                          col    = paste0(col_i, "44"),
                                          border = col_i,
                                          lwd    = 1)
                    }
                    legend_labels <- c(legend_labels, paste0(con, "\n(n=", length(sig), ")"))
                    legend_cols   <- c(legend_cols, col_i)
                  }
                  
                  legend_title <- paste0("DEGs padj<", pval_cut, "\n|logFC|>", lfc_cut)
                  
                  if (is_zoom) {
                    legend(x      = "bottomright",
                           inset  = c(0.02, 0),
                           legend = legend_labels, fill = legend_cols,
                           title  = legend_title,
                           bty    = "n", cex = 0.65, xpd = NA,
                           y.intersp = 1.8)
                  } else {
                    legend(x      = "right",
                           inset  = c(0.02, 0),
                           legend = legend_labels, fill = legend_cols,
                           title  = legend_title,
                           bty    = "n", cex = 0.65, xpd = NA,
                           y.intersp = 1.8)
                  }
                }
                
                # Gen seleccionat a la taula → vermell per sobre de tot (igual que abans)
                if (length(sel_rows) > 0 && sel_rows %in% names(gr)) {
                  karyoploteR::kpRect(kp, data = gr[sel_rows],
                                      y0 = 0, y1 = 0.95,
                                      col = "red", border = "red", lwd = 3)
                }
              })
            })
  
  # --------- CircosPlot -------------
  
  setMethod(".renderOutput", "CircosPlot",
            function(x, se, ..., output, pObjects, rObjects) {
              
              panel_name <- .getEncodedName(x)
              
              output[[panel_name]] <- renderPlot({
                
                # Forcem reactivitat quan s'actualitza el panel
                force(rObjects[[paste0(panel_name, "_INTERNAL_output_update")]])
                
                # Llegim els contrastos dels slots (inicialitzats des de input$contrast a app.R)
                c1 <- pObjects$memory[[panel_name]]@Contrast1
                c2 <- pObjects$memory[[panel_name]]@Contrast2
                contrasts_sel <- unique(c(c1, c2))  # deduplicar si només hi ha un contrast seleccionat
                
                # Llegim els thresholds del sidebar 
                input$padj
                input$logFC
                pval_cut <- if (!is.null(input$padj))  input$padj  else 0.05
                lfc_cut  <- if (!is.null(input$logFC)) input$logFC else 1.0
                
                # Extraiem posicions genòmiques i dades de rowData
                gr       <- rowRanges(se)
                rd       <- as.data.frame(rowData(se))
                gr_chr   <- as.character(seqnames(gr))
                gr_start <- start(gr)
                gr_end   <- end(gr)
                gene_nms <- rownames(se)
                
                # Filtrem a cromosomes estàndards (descarta scaffolds i patches)
                std_chrs   <- paste0("chr", c(1:22, "X", "Y"))
                valid_mask <- gr_chr %in% std_chrs
                validate(need(sum(valid_mask) > 0, "No valid chromosomal data"))
                
                # Paleta de colors: dos colors per contrast (up i down)
                up_colors   <- c("#e74c3c", "#f39c12", "#9b59b6", "#1abc9c")
                down_colors <- c("#3498db", "#2ecc71", "#e91e63", "#ff9800")
                
                # Filtrem DEGs per cada contrast separant up i down
                sig_list <- list()
                for (i in seq_along(contrasts_sel)) {
                  con    <- contrasts_sel[[i]]
                  col_p  <- paste0(con, "_padj")
                  col_fc <- paste0(con, "_log2FoldChange")
                  # Saltem si les columnes no existeixen a rowData
                  if (!col_p %in% colnames(rd) || !col_fc %in% colnames(rd)) next
                  is_sig <- !is.na(rd[[col_p]]) & rd[[col_p]] < pval_cut & abs(rd[[col_fc]]) >= lfc_cut
                  sig_list[[con]] <- list(
                    up     = rownames(rd)[is_sig & rd[[col_fc]] >  0],  # logFC positiu = upregulated
                    down   = rownames(rd)[is_sig & rd[[col_fc]] <= 0],  # logFC negatiu = downregulated
                    fc_col = col_fc
                  )
                }
                validate(need(
                  any(sapply(sig_list, function(s) length(s$up) + length(s$down)) > 0),
                  "No significant genes found"
                ))
                
                # Inicialitzem el plot circos
                circlize::circos.clear()
                circlize::circos.par(start.degree = 90, gap.degree = 2)  # chr1 comença a dalt
                circlize::circos.initializeWithIdeogram(species = genome_sp(), plotType = c("ideogram", "labels"))
                
                # Vectors per construir la llegenda al final
                legend_labels <- c()
                legend_cols   <- c()
                
                # Un track circos per cada contrast
                for (i in seq_along(contrasts_sel)) {
                  con <- contrasts_sel[[i]]
                  if (!con %in% names(sig_list)) next
                  
                  sig_up   <- sig_list[[con]]$up
                  sig_down <- sig_list[[con]]$down
                  col_fc   <- sig_list[[con]]$fc_col
                  all_sig  <- c(sig_up, sig_down)
                  if (length(all_sig) == 0) next
                  
                  # Colors d'aquest contrast (ciclant si hi ha més contrastos que colors)
                  up_col   <- up_colors[((i - 1) %% length(up_colors)) + 1]
                  down_col <- down_colors[((i - 1) %% length(down_colors)) + 1]
                  
                  # Construïm el BED dels DEGs d'aquest contrast amb el logFC com a valor Y
                  sig_mask <- gene_nms %in% all_sig & valid_mask
                  bed_con <- data.frame(
                    chr   = gr_chr[sig_mask],
                    start = gr_start[sig_mask],
                    end   = gr_end[sig_mask],
                    logfc = as.numeric(rd[[col_fc]][sig_mask])
                  )
                  # Eliminem files amb NA o cromosomes no estàndards
                  bed_con <- bed_con[bed_con$chr %in% std_chrs & !is.na(bed_con$logfc), ]
                  if (nrow(bed_con) == 0) next
                  
                  # Rang de l'eix Y: cobreix el logFC real + marge flexible
                  lfc_range <- range(bed_con$logfc, na.rm = TRUE)
                  margin <- 0.5
                  lfc_ylim <- c(lfc_range[1] - margin, lfc_range[2] + margin)
                  
                  # Capturem colors en variables locals per al closure de panel.fun
                  # (evita que el loop sobreescrigui el valor capturat)
                  .up_col   <- up_col
                  .down_col <- down_col
                  
                  circlize::circos.genomicTrackPlotRegion(
                    bed_con,
                    ylim         = lfc_ylim,
                    track.height = 0.25,
                    bg.border    = "grey80",
                    bg.col       = "grey97",
                    panel.fun = function(region, value, ...) {
                      lfc_vals <- as.numeric(value[[1]])
                      # Color per punt: up si logFC > 0, down si logFC <= 0
                      col_vec  <- ifelse(lfc_vals > 0, .up_col, .down_col)
                      circlize::circos.genomicPoints(region, value, col = col_vec, pch = 16, cex = 0.5)
                      # Línia de referència a logFC = 0
                      circlize::circos.lines(CELL_META$cell.xlim, c(0, 0), col = "grey50", lty = 2, lwd = 0.5)
                    }
                  )
                  
                  # Afegim entrada a la llegenda per up i down d'aquest contrast
                  legend_labels <- c(legend_labels,
                                     paste0(con, " up (n=",   length(sig_up),   ")"),
                                     paste0(con, " down (n=", length(sig_down), ")"))
                  legend_cols <- c(legend_cols, up_col, down_col)
                }
                
                # Títol i llegenda final
                title(paste("Circos:", paste(contrasts_sel, collapse = " vs ")), cex.main = 1.2, line = -2)
                legend("bottomleft",
                       legend = legend_labels,
                       col    = legend_cols,
                       pch    = 16,
                       pt.cex = 1.5,
                       bty    = "n",
                       title  = paste0("padj<", pval_cut, "  |logFC|>", lfc_cut),
                       cex    = 0.9)
              })
            })
  
  # amagar data parameters del circosPlot
  setMethod(".hideInterface", "CircosPlot", function(x, field) {
    if (field %in% c("SelectionBoxOpen", "DataBoxOpen")) TRUE else callNextMethod()
  })
  
  output$isee_ui2 <- renderUI({
    req(dde2(), input$contrast, input$tabs)
    
    dde3 <- as(dde2(), "SummarizedExperiment")
    #definim rownames (amb el as(, se), es perden)
    noms_gens <- rownames(dde2())
    rownames(dde3) <- noms_gens
    
    # conversió a rowRanges
    has_row_ranges <- te_rowRanges()
    if (has_row_ranges) {
      rr <- rowRanges(dde2())[rownames(dde3), ]
      dde3_rse <- SummarizedExperiment(
        assays    = assays(dde3),
        rowRanges = rr,
        colData   = colData(dde3)
      )
      mcols(rowRanges(dde3_rse)) <- rowData(dde3)
      rownames(dde3_rse) <- rownames(dde3)
      dde3 <- dde3_rse
    }
    
    #per la rowDataTable
    rd <- as.data.frame(rowData(dde2()))
    
    cols_interes <- c("Geneid", "Symbol", "GeneName", "Chr", "Start", "End", "Strand", "length", "Description")
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
    
    # forçam a factor per a que no surti la llegenda continua
    colData(dde3)[[selected_cluster_var]] <- as.factor(colData(dde3)[[selected_cluster_var]])
    if ("Cond" %in% colnames(colData(dde3))) {
      colData(dde3)$Cond <- as.factor(colData(dde3)$Cond)
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
        
        req(input$contrast)
        sel_contrasts <- input$contrast
        rd <- as.data.frame(rowData(dde))
        lists <- list()
        for (con in sel_contrasts) {
          col_p <- paste0(con, "_padj")
          col_fc <- paste0(con, "_log2FoldChange")
          #comprovem que les columnes existeixen al rowData
          if (col_p %in% colnames(rd) && col_fc %in% colnames(rd)) {
            gens_sig <- rownames(rd)[which(rd[[col_p]] < pval_cut & 
                                             !is.na(rd[[col_p]]) & 
                                             abs(rd[[col_fc]]) >= lfc_cut)]
            #Només afegim a la llista si hi ha algun gen
            lists[[con]] <- gens_sig
          }
        }
        
        if (length(lists) < 2) {
          return(ggplot() + 
                   annotate("text", x=0, y=0, 
                            label="Select at least 2 contrasts with significant genes\nto see intersections.") + 
                   theme_void())
        }
        #limit de 4
        ggvenn::ggvenn(lists[1:min(4, length(lists))], 
                       fill_color = c("#00AFBB", "#E7B800", "#FC4E07", "#7E4E90")[1:length(lists)], 
                       stroke_size = 0.5,
                       set_name_size = 4,
                       text_size = 4
        ) + labs(title = paste("DEGs intersection (P adjusted Value <", pval_cut, "& | logFC | >", lfc_cut, ")")) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      },
      restrict = NULL,
      className = "VennDiagram1",
      fullName = "Venn Diagram"
    )
    
    # d'aquesta manera treiem el boc Data Parameters que estava buit
    setMethod(".hideInterface", "VennDiagram1", function(x, field) {
      if (field %in% c("DataBoxOpen")) TRUE  # amaga aquest slot
      else callNextMethod() # la resta, comportament per defecte
    })
    
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
    
    # d'aquesta manera treiem el boc Data Parameters que estava buit
    setMethod(".hideInterface", "FEADotPlot", function(x, field) {
      if (field %in% c("DataBoxOpen")) TRUE  # amaga aquest slot
      else callNextMethod() # la resta, comportament per defecte
    })
    
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
    
    # d'aquesta manera treiem el boc Data Parameters que estava buit
    setMethod(".hideInterface", "FEABarPlot", function(x, field) {
      if (field %in% c("DataBoxOpen")) TRUE  # amaga aquest slot
      else callNextMethod() # la resta, comportament per defecte
    })
    
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
                                                      ColumnSelectionSource = "ColumnDataTable1", # linked to ColumnDataTable1
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
                                               RowSelectionSource = "RowDataTable1", 
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
      
      initial_panels2[["ColumnDataTable1"]] <- new("ColumnDataTable",
                                                   Search = "",
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
      
      if (has_row_ranges) {
        chr_sel <- chr_sel_reactive()
        
        
        initial_panels2[["KaryoPlot"]] <- new("KaryoPlot",
                                              PanelHeight = 400L, PanelWidth = 6L,
                                              RowSelectionSource = "RowDataTable1", ZoomChr = FALSE)
        
        initial_panels2[["KaryoPlot2"]] <- new("KaryoPlot",
                                               PanelHeight = 400L, PanelWidth = 6L,
                                               RowSelectionSource = "RowDataTable1", ZoomChr = TRUE)
        
        initial_panels2[["CircosPlot1"]] <- new("CircosPlot",
                                                PanelHeight = 400L, PanelWidth = 6L,
                                                RowSelectionSource = "RowDataTable1",
                                                Contrast1 = input$contrast[1],
                                                Contrast2 = if (length(input$contrast) >= 2) input$contrast[2] else input$contrast[1])
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
  
} #server
shinyApp(ui, server)