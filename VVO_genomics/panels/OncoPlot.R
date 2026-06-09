# panels/OncoPlot.R

# ----- CLASSE ONCOPLOT ------

# nova classe heretant de Panel
setClass("OncoPlot",
         contains = "Panel",
         slots = c(
           MutationAssay = "character",
           TopNGenes     = "integer",
           RowFontSize    = "numeric",   # mida lletra gens
           ColFontSize    = "numeric"    # mida lletra mostres
         ),
         prototype = prototype(
           MutationAssay      = "mutations",   # slot propi d'OncoPlot
           TopNGenes          = 20L,           # slot propi d'OncoPlot
           RowSelectionSource = "RowDataTable1", # slot heretat de Panel
           RowFontSize = 8,
           ColFontSize = 6
         )
)

# Nom del panell
setMethod(".fullName",   "OncoPlot", function(x) "OncoPrint")

# Color del panell
setMethod(".panelColor", "OncoPlot", function(x) "#E64B35")

# Defineix el tipus d'output: un plotOutput de l'alçada definida al slot PanelHeight
setMethod(".defineOutput", "OncoPlot", function(x) {
  plotOutput(.getEncodedName(x), height = paste0(x@PanelHeight, "px"))
})

# Controls que apareixen dins el box "Data Parameters" del panell
setMethod(".defineDataInterface", "OncoPlot", function(x, se, select_info) {
  panel_name <- .getEncodedName(x)
  list(
    numericInput(
      paste0(panel_name, "_RowFontSize"),
      label = "Gene name font size:",
      value = x@RowFontSize, min = 4, max = 20, step = 1
    ),
    numericInput(
      paste0(panel_name, "_ColFontSize"),
      label = "Sample name font size:",
      value = x@ColFontSize, min = 4, max = 20, step = 1
    )
  )
})

# Observadors que escolten els canvis dels controls i actualitzen els slots
setMethod(".createObservers", "OncoPlot",
          function(x, se, input, session, pObjects, rObjects) {
            callNextMethod()  # manté els observadors per defecte de Panel
            
            panel_name <- .getEncodedName(x)
            
            # Quan canvia RowFontSize -> actualitza el slot i força re-dibuix
            observeEvent(input[[paste0(panel_name, "_RowFontSize")]], {
              pObjects$memory[[panel_name]]@RowFontSize <- input[[paste0(panel_name, "_RowFontSize")]]
              .requestUpdate(panel_name, rObjects)
            }, ignoreInit = TRUE)
            
            # Quan canvia ColFontSize -> actualitza el slot i força re-dibuix
            observeEvent(input[[paste0(panel_name, "_ColFontSize")]], {
              pObjects$memory[[panel_name]]@ColFontSize <- input[[paste0(panel_name, "_ColFontSize")]]
              .requestUpdate(panel_name, rObjects)
            }, ignoreInit = TRUE)
          }
)

# Funció principal: defineix com es renderitza el plot 
setMethod(".renderOutput", "OncoPlot",
          function(x, se, ..., output, pObjects, rObjects) {
            panel_name <- .getEncodedName(x) # nom codificat del panell (OncoPlot1)
            panel_src  <- x@RowSelectionSource # nom del panell de selecció
            
            output[[panel_name]] <- renderPlot({
              
              # Força que el plot es torni a dibuixar quan canvia la selecció a la taula
              force(rObjects[[paste0(panel_name, "_INTERNAL_output_update")]])
              
              # Llegim els slots de font des de la memòria del panell (ja actualitzats pels observadors)
              current <- pObjects$memory[[panel_name]]
              row_fs  <- if (!is.null(current)) current@RowFontSize else x@RowFontSize
              col_fs  <- if (!is.null(current)) current@ColFontSize else x@ColFontSize
              
              
              # Colors per cada tipus de mutació 
              col <- c(Missense_Mutation = "#2a9134",
                       Nonsense_Mutation = "#ffca3a",
                       Nonstop_Mutation = "#000000",
                       Frame_Shift_Del = "blue",
                       Frame_Shift_Ins = "purple",
                       In_Frame_Ins = "lightblue",
                       In_Frame_Del = "plum1",
                       Translation_Start_Site = "#ff0a54",
                       Splice_Site = "darkorange",
                       Multi_Hit = "#dab49d")
              
              # alter_fun: dibuix per cada tipus de mutació
              alter_fun <- list(
                background = alter_graphic("rect", fill = "#CCCCCC"),
                Missense_Mutation = alter_graphic("rect",
                                                  width = 1,
                                                  height = 1,
                                                  fill = col["Missense_Mutation"]),
                Nonsense_Mutation = alter_graphic("rect",
                                                  width = 1,
                                                  height = 1,
                                                  fill = col["Nonsense_Mutation"]),
                Nonstop_Mutation = alter_graphic("rect",
                                                 width = 1,
                                                 height = 1,
                                                 fill = col["Nonstop_Mutation"]),
                Multi_Hit = alter_graphic("rect",
                                          width = 1,
                                          height = 1,
                                          fill = col["Multi_Hit"]),
                Frame_Shift_Del = alter_graphic("rect",
                                                width = 1,
                                                height = 1,
                                                fill = col["Frame_Shift_Del"]),
                Frame_Shift_Ins = alter_graphic("rect",
                                                width = 1,
                                                height = 1,
                                                fill = col["Frame_Shift_Ins"]),
                Translation_Start_Site = alter_graphic("rect",
                                                       width = 1,
                                                       height = 1,
                                                       fill = col["Translation_Start_Site"]),
                Splice_Site = alter_graphic("rect",
                                            width = 1,
                                            height = 1,
                                            fill = col["Splice_Site"]),
                In_Frame_Ins = alter_graphic("rect",
                                             width = 1,
                                             height = 1,
                                             fill = col["In_Frame_Ins"]),
                In_Frame_Del = alter_graphic("rect",
                                             width = 1,
                                             height = 1,
                                             fill = col["In_Frame_Del"]))
              
              # Llegim la matriu de mutacions de l'assay del SE
              mut_mat <- as.matrix(assay(se, x@MutationAssay))
              
              
              # Calculem la freqüència de mutació per gen i ens quedem els top N
              mut_freq <- rowSums(mut_mat != "" & !is.na(mut_mat))
              top      <- names(sort(mut_freq, decreasing = TRUE))[seq_len(min(x@TopNGenes, sum(mut_freq > 0)))]
              mut_mat  <- mut_mat[top, , drop = FALSE]
              
              
              # OncoPrint 
              ComplexHeatmap::oncoPrint(
                mat                     = mut_mat,
                col                     = col,
                column_order            = colnames(mut_mat),
                alter_fun               = alter_fun,
                alter_fun_is_vectorized = FALSE,
                show_row_names          = TRUE,
                row_names_gp            = gpar(fontsize = row_fs, fontface = "italic"),
                row_names_side          = "left",
                show_pct                = TRUE,
                show_column_names       = TRUE,
                column_names_gp         = gpar(fontsize = col_fs, fontface = "italic"),
                pct_side                = "right",
                column_names_side       = "bottom",
                show_heatmap_legend     = TRUE,
                top_annotation = NULL
              )
            })
          }
)


# Necessari per iSEE: retorna contingut buit
setMethod(".generateOutput", "OncoPlot",
          function(x, se, ..., all_memory, all_contents) list(contents = NULL, commands = list())
)