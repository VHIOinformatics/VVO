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

setMethod(".renderOutput", "CircosPlot",
          function(x, se, ..., output, pObjects, rObjects) {
            
            # Agafem el nom únic d'aquest panell
            # iSEE el genera automàticament a partir del nom de la classe + ID
            panel_name <- .getEncodedName(x)
            
            # Creem el renderPlot associat a aquest panell
            output[[panel_name]] <- renderPlot({
              
              # reactivitat
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
              
              
              # ------- TRACK 1: logFC del contrast 1 ---------
              
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
              
              # ------- TRACK 2: logFC del contrast 2 ---------
              
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
              title(paste("Circos:", c1, "vs", c2), cex.main = 1.5, line = -2)
              
              legend("bottomleft", 
                     legend = c(paste(c1, "up"), paste(c1, "down"), 
                                paste(c2, "up"), paste(c2, "down")),
                     col = c("#e74c3c","#3498db","#f1c40f","#2ecc71"),
                     pch = c(15, 16, 16, 16, 16),
                     pt.cex = 1.5,
                     bty = "n",
                     cex = 1)
            })
          }
)