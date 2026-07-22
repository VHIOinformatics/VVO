# panels/OncoPlot.R

# ----- CLASSE ONCOPLOT ------

# nova classe heretant de Panel
setClass("OncoPlot",
         contains = "Panel",
         slots = c(
           MutationAssay = "character",
           TopNGenes = "numeric",
           RowFontSize = "numeric",   # mida lletra gens
           ColFontSize   = "numeric",   # mida lletra mostres
           ColDataColors = "list",      # canvi de colors del colData
           GeneOrder = "character",
           SelectedGenes = "character"
         ),
         prototype = prototype(
           MutationAssay = "mutations",
           TopNGenes = 5,
           RowSelectionSource = "RowDataTable1",
           ColumnSelectionSource = "ColumnDataTable1",
           RowFontSize = 8,
           ColFontSize = 6,
           ColDataColors = list(),  # per defecte llista buida, s'omple quan l'usuari canvia colors
           GeneOrder = "mutation_freq",
           SelectedGenes = character(0)
         )
)

# Nom del panell
setMethod(".fullName", "OncoPlot", function(x) "OncoPrint")

# Color del panell
setMethod(".panelColor", "OncoPlot", function(x) "#E64B35")

# Defineix el tipus d'output: un plotOutput de l'alçada definida al slot PanelHeight
setMethod(".defineOutput", "OncoPlot", function(x) {
  panel_name <- .getEncodedName(x)
  tagList(
    plotOutput(panel_name, height = paste0(x@PanelHeight, "px")),
    downloadButton(paste0(panel_name, "_downloadPNG"), "Download PNG",
                   class = "btn-sm", style = "margin-top:8px;")
  )
})

# Controls que apareixen dins el box "Data Parameters" del panell
setMethod(".defineDataInterface", "OncoPlot", function(x, se, select_info) {
  panel_name <- .getEncodedName(x)
  all_genes  <- sort(rownames(se))
  
  # llegeix i converteix colData a majúscules
  cd <- as.data.frame(lapply(colData(se), function(col) {
    if (is.character(col) || is.factor(col)) toupper(as.character(col)) else col
  }), stringsAsFactors = FALSE, row.names = rownames(colData(se)))
  
  # Detecta quines variables són contínues (numèriques) per no oferir-los selector de color discret
  is_continuous <- vapply(cd, is.numeric, logical(1))
  
  # Per cada variable CATEGÒRICA del colData, crea un desplegable de colors
  color_inputs <- lapply(names(cd)[!is_continuous], function(var) {
    vals <- sort(unique(na.omit(as.character(cd[[var]]))))
    n <- length(vals)
    
    # Paleta de colors per defecte: 3 colors fixos si n<=3, sino Set2 (color Brewer)
    fixed_colors <- c("#87D2E6", "#CB87E6", "#E69C87")
    default_pal <- if (n <= 3) fixed_colors[seq_len(n)] else colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n)
    
    # ID del checkbox que controla si es mostra o no el bloc de colors d'aquesta variable
    toggle_id <- paste0(panel_name, "_show_", var)
    
    tagList(
      # Checkbox amb el nom de la variable: quan es marca, apareixen els textInputs de colors
      checkboxInput(toggle_id, label = var, value = FALSE),
      
      # Bloc de colors que només es mostra quan el checkbox està marcat
      conditionalPanel(
        condition = paste0("input['", toggle_id, "'] === true"),
        
        # controls de color per cada valor
        lapply(seq_along(vals), function(i) {
          input_id <- paste0(panel_name, "_color_", var, "_", i)
          hex_input_id <- paste0(input_id, "_hex")
          
          current_color <- tryCatch(x@ColDataColors[[var]][i], error = function(e) NULL)
          default_color <- if (!is.null(current_color) && !is.na(current_color)) current_color else default_pal[i]
          
          div(style = "display:flex; align-items:center; gap:8px; margin-bottom:4px; padding-left:10px;",
              tags$span(vals[i], style = "min-width:80px; font-size:12px;"),
              
              # Color picker visual
              colourpicker::colourInput(
                input_id, label = NULL, value = default_color,
                palette = "square", showColour = "background", width = "40px"
              ),
              
              # Text input per hex manual
              div(textInput(
                hex_input_id, label = NULL, value = default_color,
                placeholder = "#XXXXXX"),
                style = "width:90px; font-family:monospace; font-size:11px;"
              )
          )
        })
      ),
      tags$hr()
    )
  })
  
  # Per cada variable contínua, crea un selector de color (gradient blanc -> color)
  continuous_color_inputs <- lapply(names(cd)[is_continuous], function(var) {
    toggle_id <- paste0(panel_name, "_show_", var)
    input_id  <- paste0(panel_name, "_color_", var)
    hex_input_id <- paste0(input_id, "_hex")
    
    saved <- tryCatch(x@ColDataColors[[var]], error = function(e) NULL)
    default_color <- if (!is.null(saved) && !is.na(saved)) saved else "#B2182B"
    
    tagList(
      checkboxInput(toggle_id, label = var, value = FALSE),
      conditionalPanel(
        condition = paste0("input['", toggle_id, "'] === true"),
        div(style = "display:flex; align-items:center; gap:8px; margin-bottom:4px; padding-left:10px;",
            tags$span("Color", style = "min-width:80px; font-size:12px;"),
            
            # Color picker visual
            colourpicker::colourInput(
              input_id, label = NULL, value = default_color,
              palette = "square", showColour = "background", width = "40px"
            ),
            
            # Text input per hex manual (igual que categòriques)
            div(textInput(
              hex_input_id, label = NULL, value = default_color,
              placeholder = "#XXXXXX"),
              style = "width:90px; font-family:monospace; font-size:11px;"
            )
        )
      ),
      tags$hr()
    )
  })
  
  # Retorna tots els controls: mida de lletra + secció de colors del colData
  list(
    numericInput(
      paste0(panel_name, "_TopNGenes"),
      label = "Show genes mutated in at least (%):",
      value = x@TopNGenes, min = 2, max = nrow(se), step = 1),
    
    selectizeInput(
      paste0(panel_name, "_SelectedGenes"),
      label = "Add specific genes:",
      choices = all_genes,
      selected = x@SelectedGenes,
      multiple = TRUE,
      options = list(placeholder = "Type a gene name...")
    ),
    
    if (!is.null(metadata(se)$tcga_annotation)) {
      radioButtons(paste0(panel_name, "_GeneOrder"),
                   label = "Gene order:",
                   choices = c("Mutation frequency" = "mutation_freq",
                               "Reference cohort frequency"     = "tcga"),
                   selected = x@GeneOrder)
    },
    
    # Mostra la secció de colors si hi ha qualsevol columna de colData (categòrica o contínua)
    if (ncol(cd) > 0) tagList(tags$strong("ColData colors:"), tags$hr(), color_inputs, continuous_color_inputs),
    
    numericInput(
      paste0(panel_name, "_RowFontSize"),
      label = "Gene name font size:",
      value = x@RowFontSize, min = 4, max = 20, step = 1),
    numericInput(
      paste0(panel_name, "_ColFontSize"),
      label = "Sample name font size:",
      value = x@ColFontSize, min = 4, max = 20, step = 1)
  )
})

# ----- OBSERVERS -----
# Escolten canvis als controls i actualitzen els slots a pObjects$memory, 
# després demanen re-renderitzat amb .requestUpdate
setMethod(".createObservers", "OncoPlot",
          function(x, se, input, session, pObjects, rObjects) {
            callNextMethod()  # manté els observadors per defecte de Panel
            
            panel_name <- .getEncodedName(x)
            panel_src  <- x@ColumnSelectionSource # nom del panell font de la selecció/filtre
            
            # Funció que filtra colData segons els quadres de cerca específics de cada columna
            apply_filters <- function() {
              cd <- as.data.frame(colData(se))
              if (ncol(cd) == 0) return(rownames(cd))  
              
              match_cols <- rep(TRUE, nrow(cd))
              # search_columns: llista amb el text de cerca per cada columna de la taula
              search_cols <- input[[paste0(panel_src, "_search_columns")]]
              
              if (!is.null(search_cols) && length(search_cols) > 0) {
                for (i in seq_along(search_cols)) {
                  txt <- search_cols[[i]]
                  if (!is.null(txt) && txt != "" && i <= ncol(cd)) {
                    match_cols <- match_cols & grepl(txt, as.character(cd[[i]]), ignore.case = TRUE)
                  }
                }
              }
              rownames(cd)[match_cols]
            }
            
            # Observer: filtre de mostres per columna a ColumnDataTable
            # Guarda el resultat a rObjects perquè .renderOutput el pugui llegir
            observeEvent(input[[paste0(panel_src, "_search_columns")]], {
              tryCatch({
                rObjects[[paste0(panel_name, "_filtered_samples")]] <- apply_filters()
                .requestUpdate(panel_name, rObjects)
              }, error = function(e) {
                message("ColData error: ", conditionMessage(e))
              })
            }, ignoreInit = TRUE, ignoreNULL = TRUE)
            
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
            
            # Quan canvi el Top genes -> actualitza el slot i força re-dibuix
            observeEvent(input[[paste0(panel_name, "_TopNGenes")]], {
              pObjects$memory[[panel_name]]@TopNGenes <- input[[paste0(panel_name, "_TopNGenes")]]
              .requestUpdate(panel_name, rObjects)
            }, ignoreInit = TRUE)
            
            # Quan canvien els gens seleccionats manualment -> actualitza el slot i força re-dibuix
            observeEvent(input[[paste0(panel_name, "_SelectedGenes")]], {
              sel <- input[[paste0(panel_name, "_SelectedGenes")]]
              pObjects$memory[[panel_name]]@SelectedGenes <- if (is.null(sel)) character(0) else sel
              .requestUpdate(panel_name, rObjects)
            }, ignoreInit = TRUE, ignoreNULL = FALSE)
            
            # Observer per canviar ordre dels gens depenent de TCGA 
            observeEvent(input[[paste0(panel_name, "_GeneOrder")]], {
              pObjects$memory[[panel_name]]@GeneOrder <- input[[paste0(panel_name, "_GeneOrder")]]
              .requestUpdate(panel_name, rObjects)
            }, ignoreInit = TRUE)
            
            # --- Observers de colors ---
            cd <- as.data.frame(lapply(colData(se), function(col) {
              if (is.character(col) || is.factor(col)) toupper(as.character(col)) else col
            }), stringsAsFactors = FALSE, row.names = rownames(colData(se)))
            is_continuous <- vapply(cd, is.numeric, logical(1))
            
            # VARIABLES CONTÍNUES 
            for (var in names(cd)[is_continuous]) {
              local({
                v <- var
                input_id <- paste0(panel_name, "_color_", v)
                hex_input_id <- paste0(input_id, "_hex")
                
                # Observador 1
                observeEvent(input[[input_id]], {
                  new_colors <- pObjects$memory[[panel_name]]@ColDataColors
                  new_colors[[v]] <- input[[input_id]]
                  pObjects$memory[[panel_name]]@ColDataColors <- new_colors
                  
                  shinyjs::runjs(paste0("$('#", hex_input_id, "').val('", input[[input_id]], "');"))
                  .requestUpdate(panel_name, rObjects)
                }, ignoreInit = TRUE)
                
                # Observador 2
                observeEvent(input[[hex_input_id]], {
                  hex_val <- input[[hex_input_id]]
                  
                  new_colors <- pObjects$memory[[panel_name]]@ColDataColors
                  new_colors[[v]] <- hex_val
                  pObjects$memory[[panel_name]]@ColDataColors <- new_colors
                  
                  colourpicker::updateColourInput(session, input_id, value = hex_val)
                  .requestUpdate(panel_name, rObjects)
                }, ignoreInit = TRUE, ignoreNULL = FALSE)
              })
            }
            
            # VARIABLES CATEGÒRIQUES
            for (var in names(cd)[!is_continuous]) {
              local({
                v <- var
                
                vals <- sort(unique(na.omit(as.character(cd[[v]]))))
                n <- length(vals)
                
                fixed_colors <- c("#87D2E6", "#CB87E6", "#E69C87")
                default_pal <- if (n <= 3) fixed_colors[seq_len(n)] else colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n)
                
                for (i in seq_along(vals)) {
                  local({
                    idx <- i
                    input_id <- paste0(panel_name, "_color_", v, "_", idx)
                    hex_input_id <- paste0(input_id, "_hex")
                    
                    observeEvent(input[[input_id]], {
                      new_colors <- pObjects$memory[[panel_name]]@ColDataColors
                      if (is.null(new_colors[[v]])) {
                        new_colors[[v]] <- setNames(default_pal, vals)
                      }
                      new_colors[[v]][vals[idx]] <- input[[input_id]]
                      pObjects$memory[[panel_name]]@ColDataColors <- new_colors
                      
                      shinyjs::runjs(paste0("$('#", hex_input_id, "').val('", input[[input_id]], "');"))
                      .requestUpdate(panel_name, rObjects)
                    }, ignoreInit = TRUE)
                    
                    observeEvent(input[[hex_input_id]], {
                      hex_val <- input[[hex_input_id]]
                      
                      new_colors <- pObjects$memory[[panel_name]]@ColDataColors
                      if (is.null(new_colors[[v]])) {
                        new_colors[[v]] <- setNames(default_pal, vals)
                      }
                      new_colors[[v]][vals[idx]] <- hex_val
                      pObjects$memory[[panel_name]]@ColDataColors <- new_colors
                      
                      colourpicker::updateColourInput(session, input_id, value = hex_val)
                      .requestUpdate(panel_name, rObjects)
                    }, ignoreInit = TRUE, ignoreNULL = FALSE)
                  })
                }
              })
            }
          }
)


# Funció principal: defineix com es renderitza el plot
setMethod(".renderOutput", "OncoPlot",
          function(x, se, ..., output, pObjects, rObjects) {
            panel_name <- .getEncodedName(x) # nom codificat del panell (OncoPlot1)
            panel_src  <- x@ColumnSelectionSource # nom del panell font de la selecció/filtre
            panel_src_row <- x@RowSelectionSource  # "RowDataTable1"
            
            output[[panel_name]] <- renderPlot({
              
              # Força que el plot es torni a dibuixar quan canvia la selecció a la taula
              force(rObjects[[paste0(panel_name, "_INTERNAL_output_update")]])
              
              force(rObjects[[paste0(panel_src_row, "_INTERNAL_single_select")]])
              
              # llegim els noms de mostra ja filtrats
              filtered_samples <- rObjects[[paste0(panel_name, "_filtered_samples")]]
              
              # Llegim els slots des de la memòria del panell (ja actualitzats pels observadors)
              current <- pObjects$memory[[panel_name]]
              
              # Llegeix el gen seleccionat a la RowDataTable (selecció simple nativa d'iSEE)
              mem_row  <- pObjects$memory[[panel_src_row]]
              sel_gene <- if (!is.null(mem_row) && !is.null(mem_row@Selected) && nchar(mem_row@Selected) > 0) {
                mem_row@Selected
              } else {
                NULL
              }
              
              row_fs  <- if (!is.null(current)) current@RowFontSize else x@RowFontSize
              col_fs  <- if (!is.null(current)) current@ColFontSize else x@ColFontSize
              top_n   <- if (!is.null(current)) current@TopNGenes  else x@TopNGenes 
              
              # Colors per cada tipus de mutació
              col <- c(Missense_Mutation       = "#2a9134",
                       Nonsense_Mutation       = "#ffca3a",
                       Nonstop_Mutation        = "#000000",
                       Frame_Shift_Del         = "blue",
                       Frame_Shift_Ins         = "purple",
                       In_Frame_Ins            = "lightblue",
                       In_Frame_Del            = "plum1",
                       Translation_Start_Site  = "#ff0a54",
                       Splice_Site             = "darkorange",
                       Multi_Hit               = "#dab49d")
              
              # alter_fun: dibuix per cada tipus de mutació
              alter_fun <- c(list(background = alter_graphic("rect", fill = "#CCCCCC")),
                             lapply(col, function(color) {
                               alter_graphic("rect", width = 1, height = 1, fill = color)
                             }))
              
              # Llegim la matriu de mutacions de l'assay del SE
              mut_mat <- as.matrix(assay(se, x@MutationAssay))
              
              # filtra mostres (columnes) segons el filtre actiu a ColumnDataTable
              if (!is.null(filtered_samples) && length(filtered_samples) > 0) {
                keep <- intersect(colnames(mut_mat), filtered_samples)
                if (length(keep) > 0) {
                  mut_mat <- mut_mat[, keep, drop = FALSE]
                }
              }
              
              # Calculem la freqüència de mutació per gen 
              mut_freq   <- rowSums(mut_mat != "" & !is.na(mut_mat))
              pct_mut    <- mut_freq / ncol(mut_mat) * 100
              gene_order <- if (!is.null(current)) current@GeneOrder else "mutation_freq"
              gene_annot <- metadata(se)$tcga_annotation
              
              # Gens que superen el % mínim (top_n)
              genes_pass <- names(pct_mut[pct_mut >= top_n])
              
              # Afegeix els gens triats manualment encara que no arribin al % mínim
              selected_genes <- if (!is.null(current)) current@SelectedGenes else character(0)
              selected_genes <- intersect(selected_genes, rownames(mut_mat))
              genes_pass <- union(genes_pass, selected_genes)
              
              if (!is.null(gene_annot) && "Gene" %in% names(gene_annot) && gene_order == "tcga") {
                gene_annot_sorted <- gene_annot[order(
                  as.numeric(gsub("%", "", gene_annot$Freq)), decreasing = TRUE), ]
                
                genes_ordered     <- gene_annot_sorted$Gene[gene_annot_sorted$Gene %in% genes_pass]
                genes_not_in_tcga <- setdiff(genes_pass, genes_ordered)
                top <- c(genes_ordered, genes_not_in_tcga)
              } else {
                top <- names(sort(pct_mut[genes_pass], decreasing = TRUE))
              }
              
              validate(need(length(top) > 0, "No genes exceed this minimum mutation percentage"))
              
              # Limit de 75 gens a plotar (més no es veuen al oncoplot)
              max_genes <- 75
              if (length(top) > max_genes) {
                top <- top[seq_len(max_genes)]
              }
              # si es selecciona un gen fora dels 75 que tmb surti
              top <- union(top, selected_genes)
              
              # Si el gen coincideix amb el primer gen (selecció automàtica per defecte d'iSEE), no el marquem
              if (!is.null(sel_gene) && length(top) > 0 && identical(sel_gene, top[1])) {
                sel_gene <- NULL
              }
              
              
              mut_mat <- mut_mat[top, , drop = FALSE]
              
              # Left annotation TCGA (si hi ha fitxer carregat)
              if (!is.null(gene_annot) && all(c("Gene", "Freq") %in% names(gene_annot))) {
                idx       <- match(top, gene_annot$Gene)
                tcga_vals <- gene_annot$Freq[idx]
                tcga_vals <- ifelse(is.na(tcga_vals), "0%", tcga_vals)
                
                annot_title <- metadata(se)$tcga_annotation_name
                if (is.null(annot_title) || annot_title == "") annot_title <- "Reference\ncohort"
                
                left_annot <- HeatmapAnnotation(
                  which = "row",
                  Reference = anno_text(tcga_vals, show_name = TRUE),
                  annotation_label = annot_title,
                  annotation_name_side = "bottom",
                  annotation_name_rot  = 0,
                  annotation_name_gp   = gpar(fontsize = 8, fontface = "bold")
                )
              } else {
                left_annot <- NULL
              }
              
              
              # Llegeix colData del SE si existeix i converteix a majúscules tot 
              cd <- as.data.frame(lapply(colData(se), function(col) {
                if (is.character(col) || is.factor(col)) toupper(as.character(col)) else col
              }), stringsAsFactors = FALSE, row.names = rownames(colData(se)))
              
              if (ncol(cd) > 0) {
                cd <- cd[colnames(mut_mat), , drop = FALSE]
                
                # Recupera els colors guardats al slot; si el slot és buit usa llista buida
                fixed_colors <- c("#87D2E6", "#CB87E6", "#E69C87")
                saved_colors <- if (!is.null(current)) current@ColDataColors else list()
                
                # Detecta quines variables del colData són contínues (numèriques)
                is_continuous <- vapply(cd, is.numeric, logical(1))
                
                # Per cada variable categòrica, decideix quin color li correspon a cada valor
                cat_colors <- lapply(names(cd)[!is_continuous], function(var) {
                  
                  # Valors únics que pot prendre la variable (ex: "MALE", "FEMALE")
                  vals <- sort(unique(na.omit(as.character(cd[[var]]))))
                  n <- length(vals)
                  
                  # Si hi ha MÉS de 12 categories, assigna gris a totes
                  if (n > 12) {
                    setNames(rep("#CCCCCC", n), vals)  # Gris per a totes
                  } else {
                    # Paleta per defecte
                    fixed_colors <- c("#87D2E6", "#CB87E6", "#E69C87")
                    default_pal <- if (n <= 3) fixed_colors[seq_len(n)] else colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n)
                    
                    # Si l'usuari ja ha guardat colors per aquesta variable fem servir aquests
                    if (!is.null(saved_colors[[var]]) && length(saved_colors[[var]]) == n) {
                      saved_colors[[var]]
                    } else {
                      setNames(default_pal, vals)
                    }
                  }
                })
                names(cat_colors) <- names(cd)[!is_continuous]
                
                # Per cada variable numèrica, crea un gradient que va de blanc (valor mínim) fins a un color (valor màxim)
                cont_colors <- lapply(names(cd)[is_continuous], function(var) {
                  
                  # Rang de valors de la variable (mínim i màxim), ignorant NA
                  rng <- range(cd[[var]], na.rm = TRUE)
                  
                  # Color guardat per l'usuari per aquesta variable (si n'hi ha)
                  saved <- saved_colors[[var]]
                  
                  # Si hi ha color guardat el fem servir, si no, color vermell per defecte
                  color <- if (!is.null(saved) && !is.na(saved)) saved else "#B2182B"
                  
                  # colorRamp2 crea la funció de gradient
                  circlize::colorRamp2(c(rng[1], rng[2]), c("white", color))
                })
                names(cont_colors) <- names(cd)[is_continuous]
                
                # Uneix els colors categòrics i continus en una sola llista per passar-la a HeatmapAnnotation
                col_data_colors <- c(cat_colors, cont_colors)
                
                # Crea l'anotació inferior del heatmap amb els colors decidits
                tryCatch({
                  bottom_annot <- HeatmapAnnotation(
                    which                = "col",
                    df                   = cd,
                    col                  = col_data_colors,
                    na_col               = "#CCCCCC",
                    annotation_name_side = "left",
                    annotation_name_rot  = 0,
                    annotation_name_gp   = gpar(fontsize = 8, fontface = "bold"),
                    simple_anno_size     = unit(0.2, "cm"))
                }, error = function(e) {
                  message("Warning: Could not map colors to colData levels. Disabling colData annotation.")
                  bottom_annot <<- NULL  # Si falla, no mostra l'anotació
                })
              } else {
                bottom_annot <- NULL
              }
              
              # Calcula l'ordre fora del oncoPrint
              col_order <- order(apply(mut_mat, 2, function(x) length(which(x != ""))), decreasing = TRUE)
              
              # Calcula l'ordre de files abans del oncoPrint
              row_order <- seq_len(nrow(mut_mat))
              
              ht <- ComplexHeatmap::oncoPrint(
                mat                     = mut_mat,
                col                     = col,
                column_order            = col_order,
                remove_empty_columns    = FALSE,
                alter_fun               = alter_fun,
                alter_fun_is_vectorized = FALSE,
                show_row_names          = TRUE,
                row_names_gp = gpar(
                  fontsize = row_fs,
                  fontface = if (!is.null(sel_gene)) ifelse(top == sel_gene, "bold.italic", "italic") else "italic",
                  col      = if (!is.null(sel_gene)) ifelse(top == sel_gene, "red", "black") else "black"
                ),
                row_names_side          = "left",
                row_order               = row_order,
                show_pct                = TRUE,
                show_column_names       = TRUE,
                column_names_gp         = gpar(fontsize = col_fs, fontface = "italic"),
                pct_side                = "right",
                column_names_side       = "bottom",
                show_heatmap_legend     = TRUE,
                bottom_annotation       = bottom_annot,
                left_annotation         = left_annot
              )
              # guarda l'últim plot dibuixat perquè el botó de descàrrega el reutilitzi
              pObjects[[paste0(panel_name, "_last_oncoplot")]] <- ht
              ht
            })
            
            # Handler de descàrrega: redibuixa l'últim oncoPrint guardat
            output[[paste0(panel_name, "_downloadPNG")]] <- downloadHandler(
              filename = function() paste0("oncoplot_", Sys.Date(), ".png"),
              content  = function(file) {
                ht <- pObjects[[paste0(panel_name, "_last_oncoplot")]]
                png(file, width = 14, height = 10, units = "in", res = 200)
                ComplexHeatmap::draw(ht)
                dev.off()
              }
            )
          }
)

# Necessari per iSEE: retorna contingut buit
setMethod(".generateOutput", "OncoPlot",
          function(x, se, ..., all_memory, all_contents) list(contents = NULL, commands = list())
)
