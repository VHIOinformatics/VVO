library(shiny)
library(shinydashboard)
library(iSEE)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(fromparse2onco)
library(DT)
library(waiter)
library(shinyjs)

options(shiny.maxRequestSize = 50 * 1024^2)

source("panels/OncoPlot.R")

# Converteix la matriu de mutacions en un SummarizedExperiment per passar a iSEE
build_se <- function(mat, vc_legend = NULL, col_data = NULL, tcga_annot = NULL, tcga_annot_name = NULL) {
  
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
    metadata(se)$tcga_annotation_name <- tcga_annot_name
  }
  
  se
}
# ---------- UI ---------------

ui <- dashboardPage(
  
  dashboardHeader(title = "VHIO's VISUAL OMICS (VVO) - Mutations", titleWidth = 240),
  
  dashboardSidebar(
    shinyjs::useShinyjs(),
    width = 240,
    sidebarMenu(
      id = "tabs",
      menuItem("Mutation Analysis", tabName = "mut_tab", icon = icon("th")),
      menuItem("Help", tabName = "help_tab", icon = icon("question-circle"))
    ),
    tags$hr(style = "border-top: 2px solid white; margin-top:4px; margin-bottom:4px;"),
    
    # Mode de variant calling: paired (tumor + control) o tumor only
    div(
      style = "display: flex; align-items: center; gap: 8px; margin-bottom: 10px;",
      tags$label("Required inputs:", style = "font-weight: bold; flex: 1; margin: 0;"),
      actionLink("help_required", icon("circle-question"), title = "Go to Help", style = "color: #0066cc;")
    ),
    radioButtons("tumor_only", "Variant calling mode:",
                 choices = c("Paired" = "FALSE", "Tumor only" = "TRUE"),
                 selected = "FALSE"),
    
    # Càrrega del fitxer Excel amb les variants (output de parseVCF)
    fileInput("excel_file",
              "Upload Excel (.xlsx):",
              accept = c(".xlsx")),
    
    tags$hr(style = "border-top: 2px solid white; margin-top:4px; margin-bottom:4px;"),
    
    # Filtres
    div(
      style = "display: flex; align-items: center; gap: 8px; margin-bottom: 10px;",
      tags$label("Filters:", style = "font-weight: bold; flex: 1; margin: 0;"),
      actionLink("help_filters", icon("circle-question"), title = "Go to Help", style = "color: #0066cc;")
    ),
    
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
    # CGI
    uiOutput("cgi_ui"),
    
    # oncoKB
    uiOutput("oncokb_ui"),
    
    # Selector dinàmic dels tipus de mutació a mostrar
    uiOutput("mutation_types_ui"),
    
    tags$hr(style = "border-top: 1px solid white; margin-top:4px; margin-bottom:4px;"),
    
    # Opcional inputs
    div(
      style = "display: flex; align-items: center; gap: 8px; margin-bottom: 10px;",
      tags$label("Optional inputs:", style = "font-weight: bold; flex: 1; margin: 0;"),
      actionLink("help_optional", icon("circle-question"), title = "Go to Help", style = "color: #0066cc;")
    ),
    
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
    # símbol de càrrega 
    use_waiter(),
    
    includeCSS(system.file(package = "iSEE", "www", "iSEE.css")),
    
    tags$head(
      tags$style(HTML("
        iframe.shiny-frame { height: 1200px !important; }
      ")),
      tags$script(HTML("
        function amagaDownloadISEE() {
          var ifr = document.querySelector('iframe.shiny-frame');
          if (!ifr) return;
          try {
            var doc = ifr.contentDocument || ifr.contentWindow.document;
            if (!doc || !doc.head) return;
            if (doc.getElementById('hide-isee-download')) return;  
            var st = doc.createElement('style');
            st.id = 'hide-isee-download';
            st.textContent = 'li.dropdown:has(i.fa-download[aria-label=\"download icon\"]) { display: none !important; }';
            doc.head.appendChild(st);
          } catch (e) {}
        }
        setInterval(amagaDownloadISEE, 800);
      "))
    ),
    
    tabItems(
      tabItem(tabName = "mut_tab",  uiOutput("isee_ui")),
      tabItem(tabName = "help_tab",
              HTML('
    <h3 style="color:#2c3e50;">User Guide - VHIO Visual Omics (VVO) - Mutations</h3>
    <p style="color:#555;">This app visualizes somatic mutations as an OncoPrint, from the output file of the function <code>parseVCF</code>.</p>
    <hr>

    <h4 id="section-required" style="color:#2c3e50;">1. Required data</h4>
    <div style="background-color:#f9f9f9; padding:10px;">
      <ul>
        <li><b>Upload Excel (.xlsx):</b> the variant table exported by <code>parseVCF</code>, converted internally with <code>fromParse2MAF</code>.</li>
        <li><b>Variant calling mode:</b> choose <b>Paired</b> (tumor + control) or <b>Tumor only</b> before/after uploading. Ypu will get an error if the mode you picked doesn\'t match the file.</li>
      </ul>
    </div>

    <h4 id="section-filters" style="color:#2c3e50;">2. Filters (sidebar)</h4>
    <div style="background-color:#f9f9f9; padding:10px;">
      <p>All filters are applied together by <code>filterMAF</code> before building the OncoPrint.</p>
      <ul>
        <li><b>Filter column values:</b> which <code>FILTER</code> values to keep (defaults to <code>PASS</code> if present).</li>
        <li><b>Min VAF tumor:</b> minimum variant allele frequency in the tumor sample.</li>
        <li><b>Max VAF control</b> (Paired mode only): maximum VAF allowed in the control sample.</li>
        <li><b>Min total reads tumor / Min alternative reads tumor:</b> coverage and support-read thresholds at the variant position.</li>
        <li><b>Annotation impact:</b> SnpEff impact levels to include (HIGH / MODERATE / MODIFIER).</li>
        <li><b>CGI Oncogenic Summary</b> (only shown if the column exists in your file): filter by CGI oncogenicity categories.</li>
        <li><b>OncoKB</b> (only shown if the column exists): filter by Likely Oncogenic / Oncogenic.</li>
        <li><b>Mutation types to show:</b> restricts to the variant classes actually present in your file AND drawable by the OncoPrint (Missense, Nonsense, Frame_Shift_Del/Ins, In_Frame_Del/Ins, Splice_Site, Translation_Start_Site, Nonstop_Mutation).</li>
      </ul>
    </div>

    <h4 id="section-optional" style="color:#2c3e50;">3. Optional inputs</h4>
    <div style="background-color:#f9f9f9; padding:10px;">
      <ul>
        <li><b>Upload sample metadata (.xlsx):</b> a table with a <code>Tumor_Sample_Barcode</code> column plus any clinical/technical variables. Adds a bottom annotation to the OncoPrint and enables the <b>ColumnDataTable</b> panel. This column (<code>Tumor_Sample_Barcode</code>) needs to match the names of the column <code>Tumor</code> in the output file from <code>parseVCF</code>.</li>
        <li><b>Sample metadata columns to show (max 6):</b> pick which metadata columns are drawn as annotation tracks.</li>
        <li><b>Reference cohort annotation file (cBioPortal format, .txt/.tsv):</b> a gene-level mutation frequency reference (columns <code>Gene</code>, <code>Freq</code>). Adds a <code>Reference cohort (%)</code> column to the gene table, a left annotation with those percentages, and a <b>Gene order</b> option in the OncoPlot controls. The name of this annotation in the plot depends on the name of the file you are uploading.</li>
      </ul>
    </div>

    <h4 style="color:#2c3e50;">4. OncoPlot panel - controls</h4>
    <div style="background-color:#f9f9f9; padding:10px;">
      <ul>
        <li><b>Show genes mutated in at least (%):</b> minimum mutation frequency (as % of shown samples) for a gene to appear. Maximum 75 genes are plotted even if more pass the threshold.</li>
        <li><b>Add specific genes:</b> manually force genes into the plot even if they don\'t reach the % threshold (and even if they fall outside the top 75).</li>
        <li><b>Gene order</b> (only visible if a reference cohort file is loaded): sort genes by <b>Mutation frequency</b> (in your samples) or by <b>Reference cohort frequency</b> (external reference).</li>
        <li><b>ColData colors:</b> one block per metadata column, hidden behind a checkbox, ypu can change the color of every variable. </li>
        <li><b>Gene name font size / Sample name font size:</b> text size for row and column labels.</li>
        <li><b>Download PNG:</b> exports the exact OncoPrint currently drawn.</li>
      </ul>
    </div>

    <h4 style="color:#2c3e50;">5. How panels react to each other</h4>
    <div style="background-color:#f9f9f9; padding:10px;">
      <ul>
        <li><b>ColumnDataTable → OncoPlot:</b> typing in a search filter in the <b>ColumnDataTable</b> restricts which samples are drawn in the OncoPrint.</li>
        <li>Any sidebar filter change (VAF, reads, impact, mutation types, etc.) re-runs <code>filterMAF</code> and rebuilds the OncoPrint from scratch.</li>
      </ul>
    </div>

    <h4 style="color:#2c3e50;">6. Troubleshooting</h4>
    <div style="background-color:#d5d2d4; padding:10px;">
      <ul>
        <li><b>"Wrong variant calling mode":</b> the selected mode (Paired/Tumor only) doesn\'t match the uploaded file — switch the mode.</li>
        <li><b>"No genes exceed this minimum mutation percentage":</b> lower the % threshold or add genes manually.</li>
        <li><b>ColData annotation missing:</b> check that your metadata file has a <code>Tumor_Sample_Barcode</code> column matching the sample names in the variant file.</li>
        <li><b>Gene order / Reference cohort column missing:</b> no reference cohort file was uploaded, or it lacks <code>Gene</code>/<code>Freq</code> columns.</li>
      </ul>
    </div>
  ')
      )
    )
  )
)

# --------------- SERVER ---------------------

server <- function(input, output, session) {
  
  # llegeix el fitxer Excel i el converteix a format MAF
  maf_loaded <- reactive({
    req(input$excel_file)
    
    waiter_show(
      html = tagList(
        spin_loaders(id = 1, color = "#333"),
        h4("Loading file...", style = "color: #333; margin-top: 15px;")
      ),
      color = "transparent"
    )
    
    on.exit(waiter_hide(), add = TRUE)
    
    tumor_only <- as.logical(input$tumor_only)
    
    # Llegeix NOMÉS els noms de columna (n_max = 0, no carrega dades)
    cols <- names(readxl::read_excel(input$excel_file$datapath, n_max = 0))
    
    # Un fitxer Paired té columnes Control_*; un Tumor only no en té
    has_control <- any(grepl("^Control_", cols, ignore.case = TRUE))
    
    # Comprova que el mode seleccionat coincideix amb el fitxer pujat
    if (tumor_only && has_control) {
      showNotification(
        "Wrong variant calling mode: you selected 'Tumor only', but this file was generated in 'Paired' mode. Please switch the variant calling mode.",
        type = "error", duration = 10
      )
      req(FALSE)
    }
    if (!tumor_only && !has_control) {
      showNotification(
        "Wrong variant calling mode: you selected 'Paired', but this file was generated in 'Tumor Only' mode. Please switch the variant calling mode.",
        type = "error", duration = 10
      )
      req(FALSE)
    }
    
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
  
  # Tipus de mutació que l'oncoplot realment és capaç de dibuixar
  mutation_types_available <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
                                "Translation_Start_Site", "Nonsense_Mutation",
                                "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins",
                                "Missense_Mutation")
  
  # desplegable per a triar el tipus de mutació que vol que surti al plot
  output$mutation_types_ui <- renderUI({
    req(maf_loaded())
    
    # Només mostra com a opcions els tipus que a la vegada existeixen al fitxer I es poden dibuixar
    vals <- intersect(mutation_types_available, unique(na.omit(maf_loaded()$Variant_Classification)))
    
    selectizeInput("mutation_types",
                   "Mutation types to show:",
                   choices  = vals,
                   selected = vals,
                   multiple = TRUE)
  })
  
  output$cgi_ui <- renderUI({
    req(maf_loaded())
    # Si la columna no existeix, no es mostra res
    if (!"CGI-Oncogenic Prediction" %in% names(maf_loaded())) return(NULL)
    # Si existeix, mostra directament les opcions (com annott)
    checkboxGroupInput("cgi_list", "CGI Oncogenic Summary:",
                       choices  = c("Oncogenic (predicted by BoostDM)",
                                    "Oncogenic (predicted by RulesDM)",
                                    "Oncogenic (annotated and predicted by BoostDM)",
                                    "Potentially Oncogenic (predicted by RulesDM)"),
                       selected = c("Oncogenic (predicted by BoostDM)",
                                    "Oncogenic (predicted by RulesDM)",
                                    "Oncogenic (annotated and predicted by BoostDM)",
                                    "Potentially Oncogenic (predicted by RulesDM)")
    )
  })
  
  output$oncokb_ui <- renderUI({
    req(maf_loaded())
    # Si la columna no existeix, no es mostra res
    if (!"OncoKB" %in% names(maf_loaded())) return(NULL)
    # Si existeix, mostra directament les opcions (com annott)
    checkboxGroupInput("oncokb_list", "OncoKB:",
                       choices  = c("Likely Oncogenic", "Oncogenic"),
                       selected = c("Likely Oncogenic", "Oncogenic"))
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
    
    input$cgi_list
    input$oncokb_list
    input$mutation_types
    input$filter_column
    input$VAF_tumor
    input$VAF_control
    input$total_tumor_reads
    input$alt_tumor_reads
    input$annott
    
    waiter_show(
      html = tagList(
        spin_loaders(id = 1, color = "#333"), 
        h4("Processing...", style = "color: #333; margin-top: 15px;")
      ),
      color = "transparent"
    )
    
    on.exit(waiter_hide(), add = TRUE)
    
    tumor_only <- as.logical(input$tumor_only)
    input$tcga_annotation_file
    input$coldata_file
    
    # Recull els valors dels filtres amb valors per defecte si encara no s'han inicialitzat
    filter_col <- if (!is.null(input$filter_column)) input$filter_column else "PASS"
    
    # Aplica els filtres de qualitat al MAF
    filtered_df <- tryCatch({
      filterMAF(maf_loaded(),
                filter_column = filter_col,
                VAF_tumor = if (!is.null(input$VAF_tumor)) input$VAF_tumor else 0,
                VAF_control = if (!is.null(input$VAF_control)) input$VAF_control else 0,
                total_tumor_reads = if (!is.null(input$total_tumor_reads)) input$total_tumor_reads else 0,
                alt_tumor_reads = if (!is.null(input$alt_tumor_reads)) input$alt_tumor_reads else 0,
                annott = if (!is.null(input$annott)) input$annott else c("HIGH", "MODERATE", "MODIFIER"),
                cgi    = !is.null(input$cgi_list) && length(input$cgi_list) > 0 && 'CGI-Oncogenic Prediction' %in% colnames(maf_loaded()),
                cgi_list = if (!is.null(input$cgi_list)) input$cgi_list else character(0),
                oncokb = !is.null(input$oncokb_list) && length(input$oncokb_list) > 0,
                oncokb_list = if (!is.null(input$oncokb_list)) input$oncokb_list else character(0))
    }, error = function(e) {
      maf_loaded()  # si falla, retorna el MAF sin filtrar
    })
    
    # Filtra pels tipus de mutació seleccionats (ABANS de prepareForOncoplot)
    if (!is.null(input$mutation_types) && length(input$mutation_types) > 0) {
      filtered_df <- filtered_df[filtered_df$Variant_Classification %in% input$mutation_types, ]
    }
    
    # Prepara la matriu d'oncoplots a partir del MAF filtrat
    result <- prepareForOncoplot(filtered_df, save_matrix = FALSE, save_tmb = FALSE)
    
    # Carrega el colData si s'ha proporcionat, filtrant només les columnes seleccionades
    col_data <- NULL
    if (!is.null(input$coldata_file)) {
      cols_sel <- if (!is.null(input$coldata_columns)) input$coldata_columns else character(0)
      col_data <- coldata_raw()[, c("Tumor_Sample_Barcode", cols_sel), drop = FALSE]
      
      # Comprova que Tumor_Sample_Barcode del colData coincideix amb les mostres del fitxer de variants
      samples_matrix  <- colnames(result$oncomatrix)
      samples_coldata <- col_data$Tumor_Sample_Barcode
      n_matched <- sum(samples_matrix %in% samples_coldata)
      n_total   <- length(samples_matrix)
      
      if (n_matched == 0) {
        showNotification(
          "Sample metadata: no sample names in 'Tumor_Sample_Barcode' match the uploaded variant file. Check that the metadata file corresponds to this dataset.",
          type = "error", duration = 15
        )
      } else if (n_matched < n_total) {
        showNotification(
          sprintf("Sample metadata: %d out of %d samples did not match 'Tumor_Sample_Barcode' and will appear with no metadata (NA).",
                  n_total - n_matched, n_total),
          type = "warning", duration = 15
        )
      }
    }
    
    # Carrega la tcga annotation si s'ha proporcionat
    tcga_annot <- NULL
    tcga_annot_name <- NULL
    if (!is.null(input$tcga_annotation_file)) {
      tcga_annot <- read.delim(input$tcga_annotation_file$datapath,
                               header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE)
      tcga_annot_name <- tools::file_path_sans_ext(input$tcga_annotation_file$name)  # <-- NOU
    }
    # Construeix i retorna el SummarizedExperiment final  
    build_se(result$oncomatrix, vc_legend = result$vc_legend, col_data = col_data,
             tcga_annot = tcga_annot, tcga_annot_name = tcga_annot_name)
  })
  
  # Observe events per a conectar el help tab amb els botons
  observeEvent(input$help_required, {
    updateTabItems(session, "tabs", selected = "help_tab")
    shinyjs::delay(300, shinyjs::runjs("document.getElementById('section-required').scrollIntoView({behavior: 'smooth'});"))
  })
  
  observeEvent(input$help_filters, {
    updateTabItems(session, "tabs", selected = "help_tab")
    shinyjs::delay(300, shinyjs::runjs("document.getElementById('section-filters').scrollIntoView({behavior: 'smooth'});"))
  })
  
  observeEvent(input$help_optional, {
    updateTabItems(session, "tabs", selected = "help_tab")
    shinyjs::delay(300, shinyjs::runjs("document.getElementById('section-optional').scrollIntoView({behavior: 'smooth'});"))
  })
  
  # --------- Panells iSEE -----------
  output$isee_ui <- renderUI({
    req(se_loaded())
    se    <- se_loaded()
    
    panels <- list(
      new("OncoPlot",
          MutationAssay      = "mutations",
          TopNGenes          = 5,
          ColumnSelectionSource = "ColumnDataTable1",
          PanelWidth         = 12L,
          PanelHeight        = 600L)
    )
    
    # Afegeix ColumnDataTable només si hi ha colData carregat
    if (!is.null(input$coldata_file) && ncol(colData(se)) > 0) {
      panels <- c(panels, list(
        new("ColumnDataTable",
            PanelWidth  = 12L,
            PanelHeight = 600L)
      ))
    }
    iSEE(se, initial = panels, appTitle = "Somatic Mutation Analysis")
  })
  
} # server

shinyApp(ui, server)