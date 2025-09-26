ids_from_nodes_id <- function(kegg_ids) {

  # Remove prefix (e.g., "cpd:", "mmu:", "ko:", "path:")
  separated_elements <- strsplit(kegg_ids, " ")
  ids <- lapply(separated_elements, function(x) sub("^[a-z]+:", "", x))
  ids <- paste(unlist(ids), collapse = ";")

  return(ids)
}
ids_from_kegg_mappings <-function(kegg_ids) {

  # Remove prefix (e.g., "cpd:", "mmu:", "ko:", "path:")
  ids <- lapply(kegg_ids, function(x) sub("^[a-z]+:", "", x))

  return(ids)
}
get_ids <- function(ids_string) {
  ids <- unlist(strsplit(ids_string, ";"))
  return(ids)
}

merge_color_lists <- function(list_of_lists) {
  # collect all unique keys across lists
  all_keys <- unique(unlist(lapply(list_of_lists, names)))

  # for each key, pick the first non-white value
  merged <- setNames(lapply(all_keys, function(k) {

    vals <- unlist(lapply(list_of_lists, function(lst) {lst[[k]]}))
    vals <- vals[!is.na(vals) & vals != "white"]
    if (length(vals) > 0) vals[1] else "white"
  }), all_keys)

  return(merged)
}

# --- 4. Combine nodes and edges using lapply to fill missing columns ---
combine_dfs <- function(df1, df2 = NULL) {

  if (is.null(df2)) {
    return(df1)
  }

  all_cols <- union(names(df1), names(df2))

  fill_na <- function(df) {
    # Which columns are not in the df?
    missing <- setdiff(all_cols, names(df))
    # Set them to NA
    if(length(missing) > 0) {
      df[missing] <- NA
    }
    df[, all_cols]
  }

  # concatenate vertically the two data frames
  rbind(fill_na(df1), fill_na(df2))
}

#' Convert KEGG pathway ID to readable pathway name
#'
#' @param pathway_id A KEGG pathway ID (e.g., "hsa04110", "mmu04110", "04110")
#' @param organism Optional organism code (e.g., "hsa", "mmu"). If not provided,
#'   it will be extracted from the pathway_id if present.
#' @return A character string with the readable pathway name
#' @examples
#' \dontrun{
#'   get_pathway_name("hsa04110")
#'   get_pathway_name("04110", organism = "hsa")
#' }
get_pathway_name <- function(id) {
  res <- KEGGREST::keggGet(id)
  res[[1]]$NAME
}

is_valid_pathway <- function(pathway_id) {
  # Check if pathway_id matches KEGG pathway formats: "hsa04110" or "04110"
  if (!is.character(pathway_id) || length(pathway_id) != 1) {
    return(FALSE)
  }
  grepl("^[a-z]{2,3}\\d{5}$", pathway_id) || grepl("^\\d{5}$", pathway_id)
}

to_organism_kegg <- function(organism) {
  if (missing(organism) || !nzchar(organism)) {
    stop("You must provide a valid KEGG organism code.")
  }

  # Handle common abbreviations
  organism_code <- switch(tolower(organism),
    "hs" = "hsa",
    "mm" = "mmu",
    "human" = "hsa",
    "mouse" = "mmu",
    "hsa" = "hsa",
    "mmu" = "mmu",
    stop(sprintf("Unknown organism abbreviation: %s", organism))
  )

  return(organism_code)
}

# Function to get KEGG compound name
get_kegg_name <- function(id) {
  entry <- KEGGREST::keggGet(id)[[1]]
  return(entry$NAME[1])  # first name
}
