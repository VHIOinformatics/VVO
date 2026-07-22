# ----------- DEFINICIÓ DE LA CLASSE CircusPlot ----------------

# Definició de la Classe
setClass("CircosPlot",
         contains = "Panel",
         slots = c(RowSelectionSource = "character",
                   Contrast1 = "character",
                   Contrast2 = "character",
                   PadjCut   = "numeric",   
                   LfcCut    = "numeric"),
         prototype = list(PadjCut = 0.05, LfcCut = 1.0))

# Constructor amb valors per defecte
CircosPlot <- function(RowSelectionSource = "RowDataTable1", 
                       Contrast1 = "FFXvsWT", 
                       Contrast2 = "GEMvsWT",
                       PadjCut   = 0.05,
                       LfcCut    = 1.0, ...) {
  new("CircosPlot", 
      RowSelectionSource = RowSelectionSource, 
      Contrast1 = Contrast1, 
      Contrast2 = Contrast2,
      PadjCut   = PadjCut,
      LfcCut    = LfcCut, 
      ...)
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
            
            # Existents: contrastos (des de Data Parameters)
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
              circlize::circos.initializeWithIdeogram(species = "hg38", plotType = c("ideogram", "labels"))
              
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
                
                # Rang de l'eix Y: cobreix el logFC real + marge per veure el threshold
                lfc_range <- range(bed_con$logfc, na.rm = TRUE)
                lfc_ylim  <- c(min(lfc_range[1], -(lfc_cut + 0.5)), max(lfc_range[2], lfc_cut + 0.5))
                
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
