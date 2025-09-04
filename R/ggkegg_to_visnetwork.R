
#' Transform a ggkegg graph to a visNetwork object.
#'
#' @param path_id KEGG pathway ID (e.g., "hsa04110" or "04110").
#' @param organism KEGG organism code (e.g., "hsa" for human, "mmu" for mouse).
#' @param de_results Named list of differential expression results. Each entry should be a list with elements: de_table (data.frame), value_column (character), feature_column (character), threshold (numeric).
#' @return A visNetwork object representing the pathway with colored nodes based on differential expression results.
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @export
ggkegg_to_visnetwork <- function(path_id, organism = "mmu", de_results = NULL) {
  if (!is_valid_pathway(path_id)) {
    stop("Invalid KEGG pathway ID format. Must be like 'hsa04110' or '04110'.")
  }

  # Validate each entry in de_results
  if (!is.null(de_results)) {
    # Check de_results is a named list with required structure
    if (!is.list(de_results) || length(de_results) == 0) {
      stop("de_results must be a non-empty list")
    }

    # Check names and structure of each entry
    invisible(
      lapply(names(de_results), function(name) {
        entry <- de_results[[name]]
        if (!is.list(entry) || !all(c("de_table", "value_column", "feature_column", "threshold") %in% names(entry))) {
          stop(paste0("Each entry in de_results must be a list with elements: data, value_column, feature_column, threshold (problem in '", name, "')"))
        }
        if (!inherits(entry$de_table, "data.frame")) {
          stop(paste0("de_results[['", name, "']]$de_table must be a data frame"))
        }
        if (!is.character(entry$value_column) || !(entry$value_column %in% colnames(entry$de_table))) {
          stop(paste0("de_results[['", name, "']]$value_column must be a column name in de_results[['", name, "']]$de_table"))
        }
        if (!is.numeric(entry$threshold) || length(entry$threshold) != 1) {
          stop(paste0("de_results[['", name, "']]$threshold must be a single numeric value"))
        }
        if (!is.character(entry$feature_column) || !(entry$feature_column %in% c(colnames(entry$de_table), "rownames"))) {
          stop(paste0("de_results[['", name, "']]$feature_column must be a column name in de_results[['", name, "']]$de_table or 'rownames'"))
        }
      })
    )
  }


  if (!is.character(organism) || length(organism) != 1) {
    stop("organism must be a single character string")
  }

  # Pathway info
  pathway_name <- paste0("(", path_id, ") ", get_pathway_name(path_id, organism))
  graph <- ggkegg::pathway(path_id)

  # ggkegg documentation: https://noriakis.github.io/software/ggkegg/

  # 1. Get nodes and edges data frames
  nodes_df <- get_nodes_df(graph, organism, delete_na = TRUE)
  edges_df <- get_edges_df(graph)

  # 2. Color and style nodes and edges
  # --- Color based on de results ---
  node_colors <- color_nodes(nodes_df, de_results)

  # --- Nodes ---
  nodes_df <- style_nodes(nodes_df, node_colors)

  # --- Edges ---
  edges_df <- style_edges(edges_df)

  # --- 3. Scale node coordinates ---
  nodes_df <- scale_dimensions(nodes_df, factor = 2.5)

  # --- 4. Define legend edges ---
  # legend_elements <- create_legend_dataframe(edges_df)

  # --- 5. Build network with legend ---
  return(
    visNetwork::visNetwork(nodes_df, edges_df, width = "100%", main = pathway_name) %>%
      visNetwork::visPhysics(enabled = FALSE) %>%
      visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") 
      # visNetwork::visLegend(addEdges = legend_elements$edges, useGroups = FALSE) %>% # , addNodes = legend_elements$nodes
      # visNetwork::visGroups(groupname = unique(nodes_df$group)) %>%
      # visNetwork::visLayout(randomSeed = 42)
  )
}


# create_legend_dataframe <- function() {
#   legend_edges <- data.frame(
#     from = c(1, 1, 1, 1),
#     to   = c(2, 3, 4, 5),
#     label = c("Activation / Expression / Phosphorylation",
#               "Inhibition / Repression / Dephosphorylation",
#               "Binding / Association / Compound / State change",
#               "Other"),
#     color = c("blue", "red", "green", "gray"),
#     dashes = c(FALSE, TRUE, FALSE, FALSE),
#     arrows = c("to", "to;middle;from", "none", "to"),
#     width = c(3, 3, 1.5, 1.5),
#     font.size = 20
#   )

#   legend_nodes <- data.frame(
#     id = 1:5,
#     label = c("Legend", "Activation", "Inhibition", "Binding", "Other")
#   )

#   return(list(nodes = legend_nodes, edges = legend_edges))
# }


scale_dimensions <- function(nodes_df, factor = 2) {
  # Scale x and y coordinates to make the graph look nicer
  nodes_df$x <- nodes_df$x * factor
  nodes_df$y <- -nodes_df$y * factor # Invert y-axis

  return(nodes_df)
}


style_nodes <- function(nodes_df, node_colors) {
  # Apply given  node dimensons
  nodes_df$widthConstraint <- abs(nodes_df$xmax - nodes_df$xmin)
  nodes_df$heightConstraint <- abs(nodes_df$ymax - nodes_df$ymin)
  nodes_df$shape <- "box"

  # Assign backfround
  color_vec <- rep("white", nrow(nodes_df))
  idx <- match(nodes_df$kegg_ids, names(node_colors))
  color_vec[!is.na(idx)] <- node_colors[idx[!is.na(idx)]]
  nodes_df$color.background <- color_vec
  nodes_df$color.highlight.border <- color_vec

  # Styling
  nodes_df$shadow <- TRUE
  nodes_df$color.border <- "black"
  # nodes_df$color.highlight.background <- "orange"

  #  Ids and labels
  nodes_df$title <- paste0(
    "<p><b>", nodes_df$type, "</b><br/>",
    "ID: ", nodes_df$kegg_ids, "<br/>",
    "<a href='https://www.kegg.jp/entry/", nodes_df$name, "' target='_blank'>Wiew on KEGG</a></p>",
    "Name: ", nodes_df$name, "<br/>"
  )
  # nodes_df$label <- nodes_df$graphics_name
  nodes_df$group <- nodes_df$type

  return(nodes_df)
}

style_edges <- function(edges_df) {
  edges_df$width <- 2
  edges_df$color <- "gray"
  edges_df$arrows <- "middle"

  edges_df <- edges_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      # Set color by main type
      color = dplyr::case_when(
        type == "PPrel" ~ "blue",
        type == "GErel" ~ "green",
        type == "ECrel" ~ "orange",
        type == "PCrel" ~ "purple",
        TRUE ~ "gray"
      ),

      # Set dashes or solid based on subtype
      dashes = dplyr::case_when(
        subtype_name %in% c("inhibition", "repression", "dephosphorylation") ~ TRUE,
        TRUE ~ FALSE
      ),

      # Set arrow type
      arrows = dplyr::case_when(
        subtype_name %in% c("activation", "expression", "phosphorylation") ~ "to",
        subtype_name %in% c("inhibition", "repression", "dephosphorylation") ~ "to;middle;from",
        subtype_name %in% c("binding/association", "compound", "state change") ~ "none",
        TRUE ~ "to"
      ),

      # Optionally, width by importance
      width = dplyr::case_when(
        subtype_name %in% c("activation", "inhibition") ~ 3,
        TRUE ~ 1.5
      )
    ) %>%
    dplyr::ungroup()

  return(edges_df)
}

color_nodes <- function(nodes_df, de_results_list) {
  # For each entry in de_results_list, color nodes accordingly
  colors_mapping <- lapply(names(de_results_list), function(name) {
    # --- Extract DE results ---
    entry <- de_results_list[[name]]
    de_results <- entry$de_table
    value_column <- entry$value_column
    feature_column <- entry$feature_column
    threshold <- entry$threshold

    # --- Row check ---
    if (feature_column == "rownames") {
      de_results$rownames <- rownames(de_results)
    } else if (!(feature_column %in% colnames(de_results))) {
      stop(paste("Column", feature_column, "not found in de_results"))
    }

    # --- Default color for nodes below threshold ---
    node_colors <- "white" # Default color
    elements_to_color <- de_results %>%
      # dplyr::filter(.data[[value_column]] < threshold) %>%
      dplyr::select(dplyr::all_of(c(feature_column, value_column)))

    # --- Override with DE table colors ---
    if (!is.null(de_results) && value_column %in% colnames(de_results)) {
      # Match using either entity column or rownames
      idx <- sapply(nodes_df$kegg_ids, function(id) {
        match(get_ids(id), ids_from_kegg_mappings(de_results[[feature_column]]))
      }, USE.NAMES = TRUE)

      vals <- sapply(idx, function(i) {
        if (all(is.na(i))) {
          return(NA)
        } else {
          return(de_results[[value_column]][i[!is.na(i)][1]])
        }
      })

      if (all(is.na(vals))) {
        message("All values for value_column are NA. No features to map.")
      } else {
        # Use raw values for coloring (no normalization)
        # Use a color palette for values
        palette <- colorRampPalette(c("blue", "white", "red"))
        # Normalize values to [-1, 1] for logFC, or [0, 1] for p-value
        vals_norm <- vals

        if (all(vals_norm >= 0, na.rm = TRUE) && all(vals_norm <= 1, na.rm = TRUE)) {
          # p-value: invert so significant is red, not significant is blue
          vals_norm <- 1 - pmin(vals_norm, 1)
        } else {
          # logFC: clamp to [-2, 2] for color scaling
          vals_norm <- pmax(pmin(vals_norm, 2), -2)
          vals_norm <- (vals_norm + 2) / 4
        }
        node_colors <- ifelse(is.na(vals_norm), node_colors, palette(100)[as.integer(vals_norm * 99) + 1])
      }
    } else {
      message("de_results not being mapped on graph, if provided check that column names are correct")
    }

    return(node_colors)
  })
  colors_merged <- merge_color_lists(colors_mapping)
  return(setNames(as.character(colors_merged), names(colors_merged)))
}

#' Create and save a colorbar PNG for the node color mapping.
#'
#' @param filename Output PNG file path.
#' @param palette Colors used for the colorbar (default: c("blue", "white", "red")).
#' @param n Number of color steps (default: 100).
#' @param value_range Numeric vector of length 2, value range (default: c(-2, 2)).
#' @param label Label for the colorbar (default: "log2FoldChange").
#' @return Invisibly returns the filename.
#'
#' @export
create_colorbar_png <- function(filename = "colorbar.png",
                                palette = c("blue", "white", "red"),
                                n = 100,
                                value_range = c(-2, 2),
                                label = "log2FoldChange") {
  colors <- colorRampPalette(palette)(n)

  png(filename, width = 400, height = 150)
  par(mar = c(3, 3, 2, 1))

  x <- seq(value_range[1], value_range[2], length.out = n)
  y <- c(0, 1)
  z <- matrix(seq(value_range[1], value_range[2], length.out = n - 1), nrow = 1)

  image(
    x = x, y = y, z = z, col = colors,
    axes = FALSE, xlab = label, ylab = "",
    main = paste("Colorbar:", label)
  )

  axis(1,
    at = seq(value_range[1], value_range[2], length.out = 5),
    labels = round(seq(value_range[1], value_range[2], length.out = 5), 2)
  )
  box()
  dev.off()

  invisible(filename)
}



get_nodes_df <- function(graph, organism, delete_na = TRUE) {
  # Convert KEGG IDs to gene symbols for better readability

  graph <- graph %>%
    dplyr::mutate(
      converted_id = ggkegg::convert_id(organism),
      converted_map = ggkegg::convert_id("pathway")
    )

  # Make sure node names are unique
  igraph::V(graph)$name_old <- igraph::V(graph)$name
  igraph::V(graph)$name <- make.unique(igraph::V(graph)$name)

  # Extract nodes + edges from igraph/tidygraph
  nodes_df <- igraph::as_data_frame(graph, what = "vertices")
  nodes_df$kegg_ids <- vapply(nodes_df$name_old, ids_from_nodes_id, FUN.VALUE = character(1))

  nodes_df <- nodes_df %>%
    dplyr::filter(!grepl("undefined", rownames(.)))

  # If nodes_df does not have 'id', set it to 'name'
  if (!("id" %in% colnames(nodes_df))) {
    # nodes_df$id <- make.unique(nodes_df$converted_id, sep = ".")
    nodes_df$label <- make.unique(nodes_df$converted_id, sep = ".")
    nodes_df$id <- nodes_df$name
  }

  # Delete nodes without label
  if (delete_na) {
    nodes_df <- nodes_df %>%
      dplyr::filter(!is.na(.data$graphics_name))
  }

  return(nodes_df)
}

get_edges_df <- function(graph) {
  edges_df <- igraph::as_data_frame(graph, what = "edges") # %>%
  # dplyr::filter(!grepl("undefined", .data$from) & !grepl("undefined", .data$to))

  # Ensure 'from' and 'to' columns exist
  if (!("from" %in% colnames(edges_df) && "to" %in% colnames(edges_df))) {
    colnames(edges_df)[1:2] <- c("from", "to")
  }

  return(edges_df)
}

#
# edges_df <- edges_df %>%
#   dplyr::mutate(
#     arrows = dplyr::case_when(
#       edges_df$subtype_value == "-->" ~ "to",
#       edges_df$subtype_value == "--|" ~ "to",
#       edges_df$subtype_value == "+p" ~ "to",
#       edges_df$subtype_value == "-p" ~ "to",
#       edges_df$subtype_value == "+u" ~ "to",
#       edges_df$subtype_value == "..>" ~ "to",
#       edges_df$subtype_value == "---" ~ NA_character_,
#       edges_df$subtype_value == "-+-" ~ "to",
#       TRUE ~ NA_character_
#     ),
#     dashes = edges_df$subtype_value %in% c("..>", "---"),
#     width = dplyr::case_when(
#       edges_df$subtype_value == "+p" ~ 3,
#       edges_df$subtype_value == "-p" ~ 2,
#       edges_df$subtype_value == "-+-" ~ 4,
#       TRUE ~ 1
#     ),
#     arrowStrikethrough = FALSE,
#     endPointOffset = dplyr::case_when(
#       edges_df$subtype_value == "--|" ~ 20, # fake T-bar inhibition by stopping short
#       TRUE ~ 0
#     ),
#     smooth = dplyr::case_when(
#       edges_df$subtype_value == "-+-" ~ TRUE,
#       TRUE ~ FALSE
#     ),
#     color = dplyr::case_when(
#       edges_df$subtype_value == "-->" ~ "darkgreen", # Activation
#       edges_df$subtype_value == "--|" ~ "red", # Inhibition
#       edges_df$subtype_value == "+p" ~ "blue", # Phosphorylation
#       edges_df$subtype_value == "-p" ~ "orange", # Dephosphorylation
#       edges_df$subtype_value == "+u" ~ "purple", # Ubiquitination
#       edges_df$subtype_value == "..>" ~ "brown", # Indirect effect
#       edges_df$subtype_value == "---" ~ "gray", # Association
#       edges_df$subtype_value == "-+-" ~ "black", # Complex/other
#       TRUE ~ "lightgray" # fallback
#     )
#   )
