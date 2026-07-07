library(shiny)
library(shinydashboard)
library(iSEE)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(fromparse2onco)
library(DT)

options(shiny.host = "0.0.0.0")
options(shiny.port = 8181)
options(shiny.maxRequestSize = 50 * 1024^2)

source("panels/OncoPlot.R")

# Converteix la matriu de mutacions en un SummarizedExperiment per passar a iSEE
build_se <- function(mat, vc_legend = NULL, col_data = NULL, tcga_annot = NULL) {
  
  # Si hi ha llegenda de variants, substitueix els codis numèrics de la matriu pels noms de les mutacions
  if (!is.null(vc_legend)) {
    mat <- apply(mat, c(1,2), function(x) {
      label <- vc_legend[as.character(x)]
      if (is.na(label)) "" else label
    })
  }
  
  # substitueix NA i "NA" per buit 
  mat[is.na(mat) | mat == "NA"] <- ""
  
  # Elimina columnes (mostres) amb nom buit o NA
  valid_cols <- !is.na(colnames(mat)) & colnames(mat) != "" & colnames(mat) != "NA"
  mat <- mat[, valid_cols, drop = FALSE]
  
  # Elimina files (gens) amb nom buit o NA 
  valid_rows <- !is.na(rownames(mat)) & rownames(mat) != "" & rownames(mat) != "NA"
  mat <- mat[valid_rows, , drop = FALSE]
  
  # Crea el SummarizedExperiment amb la matriu de mutacions com a assay
  se <- SummarizedExperiment(assays = list(mutations = mat))
  
  # Si hi ha metadades de mostres (colData), les afegeix al SE ordenades per mostra
  if (!is.null(col_data)) {
    col_data_ordered <- col_data[match(colnames(mat), col_data[["Tumor_Sample_Barcode"]]), , drop = FALSE]
    col_data_df <- as.data.frame(col_data_ordered[, !names(col_data_ordered) %in% "Tumor_Sample_Barcode", drop = FALSE])
    rownames(col_data_df) <- colnames(mat)
    colData(se) <- DataFrame(col_data_df)
  }
  
  # Si hi ha fitxer de referència, l'afegim com a columna a rowData perquè es vegi a la RowDataTable
  if (!is.null(tcga_annot) && all(c("Gene", "Freq") %in% names(tcga_annot))) {
    idx <- match(rownames(se), tcga_annot$Gene)
    rowData(se)$"Reference cohort (%)" <- ifelse(is.na(idx), "0%", tcga_annot$Freq[idx])
  }
  
  # si hi ha file tcga, ficar-lo a metadata del SummarizedExperiment
  if (!is.null(tcga_annot)) {
    metadata(se)$tcga_annotation <- tcga_annot
  }
  
  se
}
# ---------- UI ---------------

ui <- dashboardPage(
  
  dashboardHeader(title = "VHIO's VISUAL OMICS (VVO) - Mutations", titleWidth = 240),
  
  dashboardSidebar(
    width = 240,
    sidebarMenu(
      id = "tabs",
      menuItem("Mutation Analysis", tabName = "mut_tab", icon = icon("th")),
      menuItem("Help", tabName = "help_tab", icon = icon("question-circle"))
    ),
    tags$hr(style = "border-top: 2px solid white; margin-top:4px; margin-bottom:4px;"),
    
    # Mode de variant calling: paired (tumor + control) o tumor only
    radioButtons("tumor_only", "Variant calling mode:",
                 choices = c("Paired" = "FALSE", "Tumor only" = "TRUE"),
                 selected = "FALSE"),
    
    # Càrrega del fitxer Excel amb les variants (output de parseVCF)
    fileInput("excel_file",
              "Upload Excel (.xlsx):",
              accept = c(".xlsx")),
    
    # Selector dinàmic de valors de la columna FILTER (s'omple quan es carrega el fitxer)
    uiOutput("filter_column_ui"),
    
    # VAF mínim a la mostra tumoral
    sliderInput("VAF_tumor",
                "Min VAF tumor:",
                min = 0, max = 1, value = 0, step = 0.05),
    
    # VAF màxim al control
    uiOutput("VAF_control_ui"),  #  només apareix en mode Paired
    
    # Cobertura mínima total a la posició de la variant
    numericInput("total_tumor_reads",
                 "Min total reads tumor:",
                 value = 0, min = 0, step = 1),
    
    # Nombre mínim de reads que suporten la variant
    numericInput("alt_tumor_reads",
                 "Min alternative reads tumor:",
                 value = 0, min = 0, step = 1),
      
    # Nivell d'impacte funcional de la variant segons SnpEff
    checkboxGroupInput("annott",
                       "Annotation impact:",
                       choices  = c("HIGH", "MODERATE", "MODIFIER"),
                       selected = c("HIGH", "MODERATE", "MODIFIER")),
    
    tags$hr(style = "border-top: 1px solid white; margin-top:4px; margin-bottom:4px;"),
    
    # Càrrega opcional de metadades de mostres (colData)
    fileInput("coldata_file",
              "Upload sample metadata (.xlsx):",
              accept = c(".xlsx")),
    
    # Selector dinàmic de quines columnes del colData es volen mostrar (màxim 6)
    uiOutput("coldata_columns_ui"),
    
    # Càrrega de annotació TCGA 
    fileInput("tcga_annotation_file", "Reference cohort annotation file (cBioPortal format)", 
              accept = c(".txt", ".tsv"))
  ),
  
  dashboardBody(
    includeCSS(system.file(package = "iSEE", "www", "iSEE.css")),
    
    tags$head(tags$style(HTML("
      iframe.shiny-frame { height: 1200px !important; }
    "))),
    
    tabItems(
      tabItem(tabName = "mut_tab",  uiOutput("isee_ui")),
      tabItem(tabName = "help_tab",
              HTML('<h3 style="color:#2c3e50;">User Guide</h3><hr>
              <h4 style="color:#2c3e50;">Coming soon</h4>')
      )
    )
  )
)

# --------------- SERVER ---------------------

server <- function(input, output, session) {
  
  # llegeix el fitxer Excel i el converteix a format MAF
  maf_loaded <- reactive({
    req(input$excel_file)
    tumor_only <- as.logical(input$tumor_only)
    fromParse2MAF(path_to_parse = input$excel_file$datapath, tumor_only = tumor_only)
  })
  
  # Llegeix el colData carregat
  coldata_raw <- reactive({
    req(input$coldata_file)
    as.data.frame(readxl::read_xlsx(input$coldata_file$datapath))
  })
  
  # Selector dinàmic de filter_column, per defecte selecciona PASS si existeix
  output$filter_column_ui <- renderUI({
    req(maf_loaded())
    vals <- unique(na.omit(maf_loaded()$FILTER))
    selectInput("filter_column",
                "Filter column values:",
                choices  = vals,
                selected = if ("PASS" %in% vals) "PASS" else vals[1],
                multiple = TRUE)
  })
  
  # Slider de VAF_control, només en mode paired
  output$VAF_control_ui <- renderUI({
    req(input$tumor_only == "FALSE")
    sliderInput("VAF_control",
                "Max VAF control:",
                min = 0, max = 1, value = 0, step = 0.05)
  })
  
  # desplegable per les variables de colData
  output$coldata_columns_ui <- renderUI({
    req(coldata_raw())
    
    # Totes les columnes excepte l'identificador de mostra
    vals <- setdiff(names(coldata_raw()), "Tumor_Sample_Barcode")
    
    # Per defecte, les 4 primeres
    default_sel <- vals[seq_len(min(4, length(vals)))]
    
    selectizeInput("coldata_columns",
                   "Sample metadata columns to show (max 6):",
                   choices  = vals,
                   selected = default_sel,
                   multiple = TRUE,
                   options  = list(maxItems = 6))
  })
  
  # aplica els filtres i construeix el SummarizedExperiment
  se_loaded <- reactive({
    req(maf_loaded())
    tumor_only <- as.logical(input$tumor_only)
    input$tcga_annotation_file
    input$coldata_file
    
    # Recull els valors dels filtres amb valors per defecte si encara no s'han inicialitzat
    filter_col <- if (!is.null(input$filter_column)) input$filter_column else "PASS"
    
    # Aplica els filtres de qualitat al MAF
    filtered_df <- filterMAF(maf_loaded(),
                             tumor_only = tumor_only,
                             filter_column = filter_col,
                             VAF_tumor = if (!is.null(input$VAF_tumor)) input$VAF_tumor else 0,
                             VAF_control = if (!is.null(input$VAF_control)) input$VAF_control else 0,
                             total_tumor_reads = if (!is.null(input$total_tumor_reads)) input$total_tumor_reads else 0,
                             alt_tumor_reads = if (!is.null(input$alt_tumor_reads)) input$alt_tumor_reads else 0,
                             annott = if (!is.null(input$annott)) input$annott else c("HIGH", "MODERATE", "MODIFIER"))
    # Prepara la matriu d'oncoplots a partir del MAF filtrat
    result <- prepareForOncoplot(filtered_df, save_matrix = FALSE, save_tmb = FALSE)
    
    # Carrega el colData si s'ha proporcionat, filtrant només les columnes seleccionades
    col_data <- NULL
    if (!is.null(input$coldata_file)) {
      cols_sel <- if (!is.null(input$coldata_columns)) input$coldata_columns else character(0)
      col_data <- coldata_raw()[, c("Tumor_Sample_Barcode", cols_sel), drop = FALSE]
    }
    
    # Carrega la tcga annotation si s'ha proporcionat
    tcga_annot <- NULL
    if (!is.null(input$tcga_annotation_file)) {
      tcga_annot <- read.delim(input$tcga_annotation_file$datapath,
                               header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE)
    }
    # Construeix i retorna el SummarizedExperiment final  
    build_se(result$oncomatrix, vc_legend = result$vc_legend, col_data = col_data, tcga_annot = tcga_annot)
  })

  # --------- Panells iSEE -----------
  output$isee_ui <- renderUI({
    req(se_loaded())
    se    <- se_loaded()
    
    panels <- list(
      new("OncoPlot",
          MutationAssay      = "mutations",
          TopNGenes          = 5,
          RowSelectionSource = "RowDataTable1",
          ColumnSelectionSource = "ColumnDataTable1",
          PanelWidth         = 12L,
          PanelHeight        = 600L),
      
      new("RowDataTable",
          PanelWidth  = 6L,
          PanelHeight = 600L)
    )
    
    # Afegeix ColumnDataTable només si hi ha colData carregat
    if (!is.null(input$coldata_file) && ncol(colData(se)) > 0) {
      panels <- c(panels, list(
        new("ColumnDataTable",
            PanelWidth  = 6L,
            PanelHeight = 600L)
      ))
    }
    iSEE(se, initial = panels)
  })
  
} # server

shinyApp(ui, server)
