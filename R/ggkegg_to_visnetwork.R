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
    stop("Invalid KEGG pathway ID format. Must be like 'hsa:04110' or '04110'.")
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
        if (!is.list(entry) || !all(c("de_table", "value_column", "feature_column") %in% names(entry))) {
          stop(paste0("Each entry in de_results must be a list with elements: data, value_column, feature_column (problem in '", name, "')"))
        }
        if (!inherits(entry$de_table, "data.frame")) {
          stop(paste0("de_results[['", name, "']]$de_table must be a data frame"))
        }
        if (!is.character(entry$value_column) || !(entry$value_column %in% colnames(entry$de_table))) {
          stop(paste0("de_results[['", name, "']]$value_column must be a column name in de_results[['", name, "']]$de_table"))
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


  organism <- to_organism_kegg(organism)

  # first call downloads, later calls use cached version
  kgml_file <- download_kgml(path_id)

  # 1. Get nodes and edges data frames
  nodes_df <- parse_kgml_entries(kgml_file)
  edges_df <- parse_kgml_relations(kgml_file)
  nodes_df$kegg_ids <- vapply(nodes_df$name, ids_from_nodes_id, FUN.VALUE = character(1))

  # Pathway info
  pathway_name <- paste0("(", path_id, ") ", get_pathway_name(path_id))
  graph <- ggkegg::pathway(path_id)

  # ggkegg documentation: https://noriakis.github.io/software/ggkegg/

  # 1. Get nodes and edges data frames
  # nodes_df <- get_nodes_df(graph, organism, delete_na = TRUE)
  # edges_df <- get_edges_df(graph)

  # 2. Color and style nodes and edges
  # --- Color based on de results ---
  node_colors <- color_nodes(nodes_df, de_results)

  # --- Nodes ---
  nodes_df <- style_nodes(nodes_df, node_colors, organism)

  # --- Edges ---
  edges_df <- style_edges(edges_df)

  # --- 3. Scale node coordinates ---
  nodes_df <- scale_dimensions(nodes_df, factor = 2.5)

  # --- 4. Define legend edges ---
  # legend_elements <- create_legend_dataframe(edges_df)

  # --- 5. Build network with legend ---
  return(
    visNetwork::visNetwork(nodes_df, edges_df, main = pathway_name) %>%
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
  nodes_df$x <- as.numeric(nodes_df$x) * factor
  nodes_df$y <- - as.numeric(nodes_df$y) * factor # Invert y-axis

  return(nodes_df)
}


style_nodes <- function(nodes_df, node_colors, organism, delete_na = FALSE) {

  # Apply given  node dimensons
  nodes_df <- nodes_df %>%
    dplyr::mutate(
      width = as.numeric(.data$width),
      heigth = as.numeric(.data$height),
      shape = dplyr::case_when(
        .data$graphics_type == "rectangle" ~ "box",
        .data$graphics_type == "circle" ~ "dot",
        .data$graphics_type == "roundrectangle" ~ "box",
        TRUE ~ "box"
      ),
      size = ifelse(.data$graphics_type == "circle", 10, NA),
      widthConstraint = ifelse(.data$graphics_type != "circle", .data$width, NA),
      heightConstraint = ifelse(.data$graphics_type != "circle", .data$height, NA)
        ) %>%
    dplyr::filter(.data$type != "group") %>%
    dplyr::mutate(
      converted_id = convert_id(organism),
      converted_map = convert_id("pathway")
    )


  # Assign backfround
  color_vec <- rep("white", nrow(nodes_df))
  idx <- match(nodes_df$name, names(node_colors))
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
    "<a href='", nodes_df$link, "' target='_blank'>View on KEGG</a></p>",
    "Name: ", nodes_df$name, "<br/>"
  )

  # nodes_df$label <- nodes_df$graphics_name
  nodes_df$group <- nodes_df$type

  # Delete nodes without label
  if (delete_na) {
    nodes_df$id[is.na(nodes_df$id)] <- nodes_df$converted_map
  }


  return(nodes_df)
}

style_edges <- function(edges_df) {
  edges_df$width <- 1.5
  edges_df$color <- "gray"
  edges_df$arrows <- "middle"

  edges_df <- edges_df %>%
    dplyr::mutate(
      # Set color by main type
      color = dplyr::case_when(
        type == "PPrel" ~ "blue",
        type == "GErel" ~ "green",
        type == "ECrel" ~ "orange",
        type == "PCrel" ~ "purple",
        TRUE ~ "gray"
      ),
      label = ifelse(value == "+p", "+p", NA_character_),
    )

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
    # threshold <- entry$threshold

    # --- Row check ---
    if (feature_column == "rownames") {
      de_results$rownames <- rownames(de_results)
    } else if (!(feature_column %in% colnames(de_results))) {
      stop(paste("Column", feature_column, "not found in de_results"))
    }

    # --- Default color ---
    elements_to_color <- de_results %>%
      # dplyr::filter(.data[[value_column]] < threshold) %>%
      dplyr::filter(!is.na(.data[[feature_column]])) %>%
      dplyr::select(dplyr::all_of(c(feature_column, value_column)))
    node_colors <- setNames(rep("white", nrow(nodes_df)), nodes_df$kegg_ids) # Default color)

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

# cache environment to avoid re-downloading
.kegg_cache <- new.env(parent = emptyenv())

#' Convert KEGG IDs to human-readable names
#'
#' @param ids Character vector of KEGG identifiers (e.g. "eco:b4025", "cpd:C00022")
#' @param org KEGG database/organism code (e.g. "eco", "hsa", "compound", "pathway").
#' Defaults to "eco".
#'
#' @return A character vector with mapped descriptions (NA if unknown).
#' @examples
#' \dontrun{
#'   convert_id(c("eco:b4025", "eco:b1779"), org = "eco")
#' }
#' @export
convert_id <- function(ids, org = "eco") {
  if (!exists(".kegg_cache", envir = .GlobalEnv)) {
    assign(".kegg_cache", new.env(parent = emptyenv()), envir = .GlobalEnv)
  }
  cache <- get(".kegg_cache", envir = .GlobalEnv)

  if (!exists(org, envir = cache)) {
    url <- paste0("https://rest.kegg.jp/list/", org)
    resp <- httr::GET(url)
    if (httr::http_error(resp)) stop("Failed to fetch KEGG list for org: ", org)

    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    tbl <- data.table::fread(text = txt, header = FALSE, sep = "\t")

    # KEGG returns 2 or 4 columns
    if (ncol(tbl) == 2) {
      colnames(tbl) <- c("id", "desc")
      tbl <- tbl[, c("id", "desc")]
    } else if (ncol(tbl) >= 4) {
      colnames(tbl)[c(1,4)] <- c("id","desc")
      tbl <- tbl[, c("id","desc")]
    } else {
      stop("Unexpected KEGG table format")
    }

    assign(org, setNames(tbl$desc, tbl$id), envir = cache)
  }

  lookup <- get(org, envir = cache)
  unname(lookup[ids])
}


parse_kgml_entries <- function(file) {
  # read the KGML file
  doc <- read_xml(file)

  # find all entry nodes
  entries <- xml_find_all(doc, ".//entry")
  graphics <- xml_find_all(entries, ".//graphics")

  # extract attributes into a data.frame
  df <- tibble(
    id = xml_attr(entries, "id"),
    name = xml_attr(entries, "name"),
    type = xml_attr(entries, "type"),
    link = xml_attr(entries, "link"),
    reaction = xml_attr(entries, "reaction"),
    graphics_name = xml_attr(graphics, "name"),
    label = xml_attr(graphics, "name"), # for visNetwork
    fgcolor = xml_attr(graphics, "fgcolor"),
    bgcolor = xml_attr(graphics, "bgcolor"),
    graphics_type = xml_attr(graphics, "type"),
    x = xml_attr(graphics, "x"),
    y = xml_attr(graphics, "y"),
    width = xml_attr(graphics, "width"),
    height = xml_attr(graphics, "height")
  )
  #     <entry id="184" name="ko:K15359 ko:K18276" type="ortholog" reaction="rn:R09472"
  # link="https://www.kegg.jp/dbget-bin/www_bget?K15359+K18276">
  # <graphics name="K15359..." fgcolor="#000000" bgcolor="#FFFFFF"
  #      type="rectangle" x="303" y="561" width="46" height="17"/>

  return(df)
}

parse_kgml_relations <- function(file) {
  doc <- read_xml(file)

  rels <- xml_find_all(doc, ".//relation")

  purrr::map_dfr(rels, function(rel) {
    entry1 <- xml_attr(rel, "entry1")
    entry2 <- xml_attr(rel, "entry2")
    type <- xml_attr(rel, "type")
    subnodes <- xml_find_all(rel, ".//subtype")

    if (length(subnodes) == 0) {
      tibble(
        from = entry1,
        to = entry2,
        type = type,
        name = NA_character_,
        value = NA_character_
      )
    } else {
      tibble(
        from = entry1,
        to = entry2,
        type = type,
        subtype = xml_attr(subnodes, "name"),
        value = xml_attr(subnodes, "value")
      )
    }
  })
}

# </entry>
# <entry id="60" name="hsa:991" type="gene"
#     link="https://www.kegg.jp/dbget-bin/www_bget?hsa:991">
#     <graphics name="CDC20, CDC20A, OOMD14, OZEMA14, bA276H19.3, p55CDC" fgcolor="#000000" bgcolor="#BFFFBF"
#          type="rectangle" x="981" y="409" width="46" height="17"/>
# </entry>

#     <relation entry1="64" entry2="66" type="PPrel">
#     <subtype name="activation" value="--&gt;"/>
#     <subtype name="dephosphorylation" value="-p"/>
# </relation>
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


#' Download and cache KEGG KGML files.
#'
#'
#' @param pathway_id KEGG pathway ID (e.g., "hsa04110").
#' @param cache_dir Directory to store cached KGML files (default: "kgml_cache").
#' @return Path to the cached KGML file.
#' @importFrom httr GET content http_error
#' @importFrom xml2 read_xml
download_kgml <- function(pathway_id, cache_dir = "kgml_cache") {
  # Ensure cache dir exists
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  # Local file path
  file_path <- file.path(cache_dir, paste0(pathway_id, ".xml"))

  # If already cached, return path
  if (file.exists(file_path)) {
    message("Using cached file: ", file_path)
    return(file_path)
  }

  # KEGG REST API URL
  url <- paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")

  # Download
  resp <- httr::GET(url)

  if (httr::http_error(resp)) {
    stop("Failed to download KGML for pathway: ", pathway_id)
  }

  writeBin(httr::content(resp, as = "text", encoding = "UTF-8"), file_path)
  message("Downloaded and cached: ", file_path)

  return(file_path)
}

capital = 15000
year = 10000
for (i in 1:10) {
capital = capital*1 + year
}
print(capital)
