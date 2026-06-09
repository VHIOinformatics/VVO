library(shiny)
library(shinydashboard)
library(iSEE)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(fromparse2onco)

# Permet pujar fitxers de fins a 50MB
options(shiny.maxRequestSize = 50 * 1024^2)

# Carrega la classe OncoPlot i els seus mètodes
source("panels/OncoPlot.R")


# Funció per construir el SE a partir del txt
build_se <- function(path_to_txt) {
  mat     <- read.table(path_to_txt, sep = "\t", header = FALSE,
                        check.names = FALSE, fill = TRUE)
  samples <- as.character(mat[1, -1])   # fila 1 sense la primera cel·la
  genes   <- as.character(mat[-1, 1])   # columna 1 des de fila 2
  mat     <- as.matrix(mat[-1, -1])
  rownames(mat) <- genes
  colnames(mat) <- samples
  mat[is.na(mat) | mat == "NA"] <- ""
  
  # Treu columnes amb nom NA o buit
  valid_cols <- !is.na(colnames(mat)) & colnames(mat) != "" & colnames(mat) != "NA"
  mat <- mat[, valid_cols, drop = FALSE]
  
  # Treu files amb nom NA o buit
  valid_rows <- !is.na(rownames(mat)) & rownames(mat) != "" & rownames(mat) != "NA"
  mat <- mat[valid_rows, , drop = FALSE]
  
  SummarizedExperiment(assays = list(mutations = mat))
}

# ---------- UI ---------------

ui <- dashboardPage(
  
  dashboardHeader(title = "VHIO's VISUAL OMICS (VVO) - Mutations", titleWidth = 240),
  
  dashboardSidebar(
    width = 240,
    sidebarMenu(
      id = "tabs",
      menuItem("Mutation Analysis", tabName = "mut_tab", icon = icon("th")),
      menuItem("Help",              tabName = "help_tab", icon = icon("question-circle"))
    ),
    tags$hr(style = "border-top: 2px solid white; margin-top:4px; margin-bottom:4px;"),
    
    fileInput("excel_file",
              "Upload .txt:",
              accept = c(".txt")),
    tags$hr(style = "border-top: 1px solid white; margin-top:4px; margin-bottom:4px;"),
    
    uiOutput("ui_top_n"),
    uiOutput("status_box")
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
  
  se_loaded <- reactive({
    req(input$excel_file)
    withProgress(message = "Building SummarizedExperiment...", value = 0.5, {
      build_se(input$excel_file$datapath) 
    })
  })
  
  output$status_box <- renderUI({
    if (is.null(input$excel_file)) {
      tags$div(class = "status-wait", icon("clock"), " Waiting for file...")
    } else {
      se <- se_loaded()
      if (is.null(se)) {
        tags$div(class = "status-err", icon("times-circle"), " Error loading data")
      } else {
        tags$div(class = "status-ok",
                 icon("check-circle"),
                 sprintf(" %d genes x %d samples", nrow(se), ncol(se))
        )
      }
    }
  })
  
  # Controls del sidebar
  # Selector numèric: top N gens més mutats a mostrar
  output$ui_top_n <- renderUI({
    req(se_loaded())
    se <- se_loaded()
    numericInput("top_n", "Top N mutated genes:",
                 value = min(25, nrow(se)), min = 2, max = nrow(se))
  })
  

  # -------- Tab Mutation Analysis ----------
  
  output$isee_ui <- renderUI({
    req(se_loaded())
    se    <- se_loaded()
    top_n <- if (!is.null(input$top_n) && !is.na(input$top_n) && input$top_n >= 2) {
      as.integer(input$top_n)
    } else {
      min(25L, nrow(se))
    }
    
    
    # ---------- PANELLS ------------
    iSEE(se, initial = list(
      new("OncoPlot",
          MutationAssay      = "mutations",
          TopNGenes          = as.integer(top_n),
          RowSelectionSource = "RowDataTable1",
          PanelWidth         = 12L,
          PanelHeight        = 600L),
      
      new("RowDataTable",
          PanelWidth  = 12L,
          PanelHeight = 600L)
    ))
  })
  
} # server

shinyApp(ui, server)
