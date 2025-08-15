#' Download KEGG pathway KXML files for a specific organism
#'
#' This function retrieves the list of KEGG pathways for a given organism and downloads their KXML files into a specified cache directory.
#'
#' @param organism Character. KEGG organism code (e.g., "hsa" for human, "mmu" for mouse).
#' @param cache_dir Character. Directory to store cached KXML files. Defaults to tempdir().
#' @param overwrite Logical. If TRUE, existing files will be re-downloaded. Defaults to FALSE.
#' @return Character vector of file paths to downloaded KXML files.
#' @export
download_kegg_kxml_files <- function(organism, cache_dir = tempdir(), overwrite = FALSE) {
  if (missing(organism) || !nzchar(organism)) {
    stop("You must provide a valid KEGG organism code.")
  }

  organism <- to_organism_kegg(organism)

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  # Get the list of pathways for the specified organism
  pathway_list_url <- paste0("http://rest.kegg.jp/list/pathway/", organism)
  pathway_list <- tryCatch(
    {
      readLines(pathway_list_url)
    },
    error = function(e) {
      stop("Failed to retrieve pathway list: ", e$message)
    }
  )

  # Extract pathway IDs
  pathway_ids <- sub("^path:", "", sapply(strsplit(pathway_list, "\t"), `[`, 1))

  downloaded_files <- character(length(pathway_ids))

  for (i in seq_along(pathway_ids)) {
    pathway_id <- pathway_ids[i]
    file_name <- paste0(pathway_id, ".xml")
    file_path <- file.path(cache_dir, file_name)

    if (!file.exists(file_path) || overwrite) {
      tryCatch(
        {
          message(sprintf("Downloading %s...", file_name))
          utils::download.file(
            url = paste0("http://rest.kegg.jp/get/", pathway_id, "/kgml"),
            destfile = file_path,
            mode = "wb"
          )
        },
        error = function(e) {
          warning(sprintf("Failed to download %s: %s", pathway_id, e$message))
          downloaded_files[i] <- NA
          next
        }
      )
    } else {
      message(sprintf("%s already exists. Skipping...", file_name))
    }

    downloaded_files[i] <- file_path
  }

  return(downloaded_files)
}

#' Load a KEGG pathway KXML file
#'
#' This function loads a KEGG pathway KXML file from disk and returns its contents as an XML document.
#'
#' @param file_path Character. Path to the KXML file to load.
#' @return An XML document object.
#' @export
load_kegg_kxml_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  xml2::read_xml(file_path)
}

#' List downloaded KEGG pathway KXML files
#'
#' This function lists all KEGG pathway KXML files present in the specified cache directory.
#'
#' @param cache_dir Character. Directory where KXML files are stored. Defaults to tempdir().
#' @return Character vector of file paths to KXML files.
#' @export
list_kegg_kxml_files <- function(cache_dir = tempdir()) {
  if (!dir.exists(cache_dir)) {
    stop("Cache directory does not exist: ", cache_dir)
  }
  list.files(cache_dir, pattern = "\\.xml$", full.names = TRUE)
}

#' Convert KEGG KGML XML to igraph object
#'
#' This function parses a KEGG KGML XML document and converts it to an igraph object.
#'
#' @param xml_doc An XML document object as returned by load_kegg_kxml_file().
#' @return An igraph object representing the pathway graph.
#' @export
kgml_to_igraph <- function(xml_doc) {
  if (!inherits(xml_doc, "xml_document")) {
    stop("Input must be an xml_document object.")
  }
  browser()
  # Note the difference between .// and //
  # // finds anywhere in the document (ignoring the current node)
  # .// finds anywhere beneath the current node

  entries <- xml2::xml_find_all(xml_doc, ".//entry")
  entry_ids <- xml2::xml_attr(entries, "id")
  entry_names <- xml2::xml_attr(entries, "name")
  entry_types <- xml2::xml_attr(entries, "type")
  # Extract graphics information for each entry
  graphics_nodes <- xml2::xml_find_first(entries, ".//graphics")
  graphics_name <- xml2::xml_attr(graphics_nodes, "name")
  graphics_fgcolor <- xml2::xml_attr(graphics_nodes, "fgcolor")
  graphics_bgcolor <- xml2::xml_attr(graphics_nodes, "bgcolor")
  graphics_type <- xml2::xml_attr(graphics_nodes, "type")
  graphics_x <- as.numeric(xml2::xml_attr(graphics_nodes, "x"))
  graphics_y <- as.numeric(xml2::xml_attr(graphics_nodes, "y"))
  graphics_width <- as.numeric(xml2::xml_attr(graphics_nodes, "width"))
  graphics_height <- as.numeric(xml2::xml_attr(graphics_nodes, "height"))

  # Extract component IDs for each entry
  components_list <- lapply(entries, function(entry) {
    comps <- xml2::xml_find_all(entry, ".//component")
    if (length(comps) == 0) {
      return(NA_character_)
    }
    xml2::xml_attr(comps, "id")
  })

  # Store as a list-column in the nodes data.frame
  nodes <- data.frame(
    id = entry_ids,
    name = entry_names,
    type = entry_types,
    component_ids = I(components_list),
    stringsAsFactors = FALSE
  )
  nodes <- data.frame(
    id = entry_ids,
    name = entry_names,
    type = entry_types,
    stringsAsFactors = FALSE
  )

  relations <- xml2::xml_find_all(xml_doc, ".//relation")
  rel_entry1 <- xml2::xml_attr(relations, "entry1")
  rel_entry2 <- xml2::xml_attr(relations, "entry2")
  rel_type <- xml2::xml_attr(relations, "type")

  # Extract component IDs for each entry
  components_list <- lapply(relations, function(realtion) {
    comps <- xml2::xml_find_all(relation, ".//component")
    if (length(comps) == 0) {
      return(NA_character_)
    }
    xml2::xml_attr(comps, "id")
  })
  # Extract subtype info for each relation
  subtypes_list <- lapply(relations, function(rel) {
    subtypes <- xml2::xml_find_all(rel, ".//subtype")
    if (length(subtypes) == 0) {
      return(list(name = NA, value = NA))
    }
    names <- xml2::xml_attr(subtypes, "name")
    values <- xml2::xml_attr(subtypes, "value")
    # Collapse multiple subtypes into a single string (comma-separated)
    list(
      name = paste(names, collapse = ","),
      value = paste(values, collapse = ",")
    )
  })
  #
  rel_subtype_name <- vapply(subtypes_list, function(x) x$name, character(1))
  rel_subtype_value <- vapply(subtypes_list, function(x) x$value, character(1))
  # Extract reactions as directed edges
  reactions <- xml2::xml_find_all(xml_doc, ".//reaction")
  reaction_ids <- xml2::xml_attr(reactions, "id")
  reaction_names <- xml2::xml_attr(reactions, "name")
  reaction_types <- xml2::xml_attr(reactions, "type")

  # For each reaction, extract substrates and products
  reaction_edges <- do.call(rbind, lapply(seq_along(reactions), function(i) {
    rxn <- reactions[[i]]
    rxn_id <- reaction_ids[i]
    substrates <- xml2::xml_find_all(rxn, ".//substrate")
    products <- xml2::xml_find_all(rxn, ".//product")
    substrate_ids <- xml2::xml_attr(substrates, "id")
    product_ids <- xml2::xml_attr(products, "id")
    # Substrate -> reaction node
    sub_edges <- if (length(substrate_ids) > 0) {
      data.frame(
        from = substrate_ids,
        to = rxn_id,
        type = "reaction_substrate",
        subtype_name = NA,
        subtype_value = NA,
        stringsAsFactors = FALSE
      )
    } else NULL
    # Reaction node -> product
    prod_edges <- if (length(product_ids) > 0) {
      data.frame(
        from = rxn_id,
        to = product_ids,
        type = "reaction_product",
        subtype_name = NA,
        subtype_value = NA,
        stringsAsFactors = FALSE
      )
    } else NULL
    rbind(sub_edges, prod_edges)
  }))

  # Add reaction nodes to nodes data.frame
  if (length(reaction_ids) > 0) {
    reaction_nodes <- data.frame(
      id = reaction_ids,
      name = reaction_names,
      type = paste0("reaction:", reaction_types),
      stringsAsFactors = FALSE
    )
    nodes <- rbind(nodes, reaction_nodes)
  }

  # Combine relation edges and reaction edges
  if (!is.null(reaction_edges) && nrow(reaction_edges) > 0) {
    edges <- rbind(
      data.frame(
        from = rel_entry1,
        to = rel_entry2,
        type = rel_type,
        subtype_name = rel_subtype_name,
        subtype_value = rel_subtype_value,
        stringsAsFactors = FALSE
      ),
      reaction_edges
    )
  } else {
    edges <- data.frame(
      from = rel_entry1,
      to = rel_entry2,
      type = rel_type,
      subtype_name = rel_subtype_name,
      subtype_value = rel_subtype_value,
      stringsAsFactors = FALSE
    )
  }
  edges <- data.frame(
    from = rel_entry1,
    to = rel_entry2,
    type = rel_type,
    subtype_name = rel_subtype_name,
    subtype_value = rel_subtype_value,
    stringsAsFactors = FALSE
  )
  # Use extract_entry_positions to get node positions if available
  positions <- extract_entry_positions(xml_doc)
  # Merge positions into nodes data.frame
  nodes <- merge(nodes, positions[, c("id", "x", "y")], by = "id", all.x = TRUE)
  library(igraph)
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  return(g)
}

#' Plot igraph object interactively using plotly
#'
#' This function creates an interactive visualization of an igraph object using plotly.
#'
#' @param g An igraph object.
#' @return A plotly object.
#' @export
plot_igraph_interactive <- function(g, de_results = NULL, color_column = "pvalue") {
  if (!inherits(g, "igraph")) {
    stop("Input must be an igraph object.")
  }

  nodes <- igraph::V(g)
  edges <- igraph::as_data_frame(g, what = "edges")
  layout <- cbind(nodes$x, nodes$y)

  edge_shapes <- lapply(1:nrow(edges), function(i) {
    from <- which(nodes$name == edges$from[i])
    to <- which(nodes$name == edges$to[i])
    list(
      type = "line",
      x0 = layout[from, 1],
      y0 = layout[from, 2],
      x1 = layout[to, 1],
      y1 = layout[to, 2],
      line = list(color = "#888", width = 1)
    )
  })

  # Determine node shapes based on prefix before ":"
  node_prefix <- sapply(strsplit(nodes$name, ":"), `[`, 1)
  entity <- sapply(strsplit(nodes$name, ":"), `[`, 2)
  node_colors <- rep("grey", length(entity))

  if (!is.null(de_results) && color_column %in% colnames(de_results)) {
    # Match entity to rownames in de_results
    idx <- match(entity, rownames(de_results))
    vals <- de_results[idx, color_column]
    # Normalize values for color mapping (e.g., pvalue: lower is more significant)
    norm_vals <- (vals - min(vals, na.rm = TRUE)) / (max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE))
    # Use RColorBrewer for color scale
    library(RColorBrewer)
    pal <- colorRampPalette(brewer.pal(9, "RdYlBu"))(100)
    color_idx <- round(norm_vals * 99) + 1
    node_colors <- ifelse(is.na(norm_vals), "grey", pal[color_idx])
  }

  # Map prefixes to plotly shapes
  shape_map <- c(
    "cpd" = "circle",
    "mmu" = "square",
    "ko" = "triangle-up",
    "path" = "star"
  )

  node_shapes <- shape_map[node_prefix]
  node_shapes[is.na(node_shapes)] <- "circle" # default shape

  # Add legend group and hover text
  legend_map <- c(
    "cpd" = "Compound",
    "mmu" = "Gene",
    "ko" = "KO",
    "path" = "Pathway"
  )

  node_legend <- legend_map[node_prefix]
  hover_text <- paste0("Entity: ", entity, "<br>Type: ", node_legend)

  # Create plotly scatter plot
  plt <- plot_ly(
    x = layout[, 1],
    y = layout[, 2],
    type = "scatter",
    mode = "markers",
    hoverinfo = "text",
    text = hover_text,
    marker = list(
      size = 15,
      color = vals,
      symbol = node_shapes,
      colorbar = list(
        title = color_column,
        len = 0.5,
        thickness = 20
      ),
      colorscale = "RdYlBu",
      showscale = TRUE
    ),
    showlegend = TRUE,
    legendgroup = node_legend
  )

  # If coloring by de_results, add color legend
  if (!is.null(de_results) && color_column %in% colnames(de_results)) {
    # Add a colorbar legend for the color_column
    plt <- colorbar(plt)
  }

  plt %>%
    layout(
      shapes = edge_shapes,
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      legend = list(title = list(text = "Node Type"))
    )
}

#'
#' This function extracts the (x, y) positions of elements (entries) from a KEGG KGML XML document.
#'
#' @param xml_doc An XML document object as returned by load_kegg_kxml_file().
#' @return A data.frame with columns: id, name, x, y.
#' @export
extract_entry_positions <- function(xml_doc) {
  if (!inherits(xml_doc, "xml_document")) {
    stop("Input must be an xml_document object.")
  }
  entries <- xml2::xml_find_all(xml_doc, ".//entry")
  entry_ids <- xml2::xml_attr(entries, "id")
  entry_names <- xml2::xml_attr(entries, "name")
  graphics_nodes <- xml2::xml_find_first(entries, ".//graphics")
  x <- as.numeric(xml2::xml_attr(graphics_nodes, "x"))
  y <- as.numeric(xml2::xml_attr(graphics_nodes, "y"))
  data.frame(
    id = entry_ids,
    name = entry_names,
    x = x,
    y = y,
    stringsAsFactors = FALSE
  )
}


#' Style igraph nodes and edges for plotting
#'
#' This function generates node and edge styles for igraph plotting, similar to plot_igraph_interactive.
#'
#' @param g An igraph object.
#' @param de_results Optional data.frame for coloring nodes.
#' @param color_column Column name in de_results for coloring.
#' @return A list with node and edge style attributes.
#' @export
style_igraph <- function(g, de_results = NULL, color_column = "pvalue") {
  if (!inherits(g, "igraph")) {
    stop("Input must be an igraph object.")
  }
  nodes <- igraph::V(g)
  edges <- igraph::as_data_frame(g, what = "edges")

  node_prefix <- sapply(strsplit(nodes$name, ":"), `[`, 1)
  entity <- sapply(strsplit(nodes$name, ":"), `[`, 2)
  node_colors <- rep("grey", length(entity))
  vals <- NULL

  if (!is.null(de_results) && color_column %in% colnames(de_results)) {
    idx <- match(entity, rownames(de_results))
    vals <- de_results[idx, color_column]
    norm_vals <- (vals - min(vals, na.rm = TRUE)) / (max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE))
    pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(100)
    color_idx <- round(norm_vals * 99) + 1
    node_colors <- ifelse(is.na(norm_vals), "grey", pal[color_idx])
  }

  # shape_map <- c(
  #   "cpd" = "circle", # circle
  #   "mmu" = "square", # square
  #   "ko" = "triangle", # triangle
  #   "path" = "diamond" # diamond
  # )
  # node_shapes <- shape_map[node_prefix]
  # node_shapes[is.na(node_shapes)] <- 21

  node_labels <- entity
  node_sizes <- rep(15, length(entity))

  edge_colors <- rep("#888888", nrow(edges))
  edge_widths <- rep(1, nrow(edges))

  igraph::V(g)$color <- node_colors
  # igraph::V(g)$shape <- node_shapes
  igraph::V(g)$label <- node_labels
  igraph::V(g)$size <- node_sizes
  igraph::E(g)$color <- edge_colors
  igraph::E(g)$width <- edge_widths
  return(g)
}


#' Style and plot an igraph network interactively with visNetwork
#'
#' @param g An igraph object.
#' @param de_results Optional data.frame for coloring nodes (rownames or 'entity' column should match entity IDs).
#' @param color_column Column in de_results to use for continuous coloring (overrides prefix color if given).
#' @return A visNetwork htmlwidget.
#' @export
style_visnetwork <- function(g, de_results = NULL, color_column = "padj", significance_threshold = 0.05) {
  if (!inherits(g, "igraph")) {
    stop("Input must be an igraph object.")
  }

  # --- Get vertex data without forcing unique rownames ---
  nodes_df <- data.frame(igraph::vertex_attr(g), stringsAsFactors = FALSE)
  edges_df <- igraph::as_data_frame(g, what = "edges")
  rownames(nodes_df) <- NULL

  # Split into prefix and entity
  node_prefix <- sapply(strsplit(nodes_df$name, ":"), `[`, 1)
  entity <- sapply(strsplit(nodes_df$name, ":"), `[`, 2)
  nodes_df$label <- entity
  nodes_df$title <- nodes_df$name

  # --- Color by node prefix ---
  prefix_colors <- c(
    "cpd"  = "lightgrey",
    "mmu"  = "lightgrey",
    "ko"   = "lightgrey",
    "path" = "lightgrey"
  )
  node_colors <- prefix_colors[node_prefix]
  node_colors[is.na(node_colors)] <- "grey"

  # --- Override with DE table colors ---
  if (!is.null(de_results) && color_column %in% colnames(de_results)) {
    # Match using either entity column or rownames
    if ("entity" %in% colnames(de_results)) {
      idx <- match(entity, de_results$entity)
    } else {
      de_results$entity <- rownames(de_results)
      idx <- match(entity, de_results$entity)
    }
    vals <- de_results[idx, color_column]

    if (all(is.na(vals))) {
      warning("All values for color_column are NA; using prefix colors.")
    } else {
      # Use raw values for coloring (no normalization)
      # Two-color palette: blue for not significant (>= 0.05), red for significant (< 0.05)
      node_colors <- ifelse(is.na(vals), node_colors,
        ifelse(vals < significance_threshold, "red", "lightblue")
      )
    }
  }

  # --- Shape by prefix ---
  prefix_shapes <- c(
    "cpd"  = "dot",
    "mmu"  = "square",
    "ko"   = "triangle",
    "path" = "diamond"
  )
  node_shapes <- prefix_shapes[node_prefix]
  node_shapes[is.na(node_shapes)] <- "dot"

  # Ensure unique IDs for visNetwork
  nodes_df$id <- make.unique(nodes_df$name, sep = "_dup")

  # Apply visual attributes
  nodes_df$shape <- node_shapes
  nodes_df$shadow <- TRUE
  nodes_df$color.background <- node_colors
  nodes_df$color.border <- "black"
  nodes_df$color.highlight.background <- "orange"
  nodes_df$color.highlight.border <- node_colors
  nodes_df$size <- 15
  nodes_df$group <- node_prefix # assign group for legend

  # --- Edges ---
  edges_df <- igraph::as_data_frame(g, what = "edges")
  edges_df$width <- 1
  edges_df$color <- "gray"
  edges_df$arrows <- "middle"
  edges_df$smooth <- FALSE
  edges_df$shadow <- FALSE

  # --- Preserve original layout ---
  # layout_mat <- igraph::layout_as_tree(g)  # fallback if no layout stored
  layout_mat <- cbind(nodes_df$x, nodes_df$y)

  if (!is.null(layout_mat) && nrow(layout_mat) == nrow(nodes_df)) {
    # Make layout more disperse by scaling coordinates
    scale_factor <- 2
    nodes_df$x <- layout_mat[, 1] * scale_factor
    nodes_df$y <- -layout_mat[, 2] * scale_factor # flip y for visNetwork coordinate system
  }

  # --- Build interactive plot ---
  vis <- visNetwork::visNetwork(nodes_df, edges_df, width = "100%", height = "500px") %>%
    visNetwork::visGroups(groupname = "cpd", color = "lightgrey", shape = "dot") %>%
    visNetwork::visGroups(groupname = "mmu", color = "lightgrey", shape = "square") %>%
    visNetwork::visGroups(groupname = "ko", color = "lightgrey", shape = "triangle") %>%
    visNetwork::visGroups(groupname = "path", color = "lightgrey", shape = "diamond") %>%
    visNetwork::visOptions(highlightNearest = TRUE, selectedBy = "type") %>%
    visNetwork::visLegend(main = "Node type", position = "right", ncol = 1) %>%
    visNetwork::visPhysics(enabled = FALSE)

  vis

  return(vis)
}
