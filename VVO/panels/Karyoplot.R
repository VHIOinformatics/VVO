# ------------ DEFINICIÓ DE LA CLASSE KaryoPlot -------------------

# Classe nova que hereta directament de Panel (classe base d'iSEE)
# slot RowSelectionSource: nom del panel del qual llegirà la selecció
setClass("KaryoPlot",
         contains = "Panel",
         slots = c(RowSelectionSource = "character",
                   ZoomChr = "logical",
                   ZoomChrValue = "character"))

# Constructor amb valor per defecte apuntant a la taula de files
KaryoPlot <- function(RowSelectionSource = "RowDataTable1",ZoomChr = FALSE, ZoomChrValue = "chr1", ...) {
  new("KaryoPlot", RowSelectionSource = RowSelectionSource, ZoomChr = ZoomChr, ZoomChrValue = ZoomChrValue, ...)
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
              
              # llegeix sempre del slot propi
              chr_sel <- pObjects$memory[[panel_name]]@ZoomChrValue
              
              if (x@ZoomChr) {
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