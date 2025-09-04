library(httr2)
library(jsonlite)
library(xml2)
library(stringr)

# Available databses: https://www.uniprot.org/help/id_mapping
#' Map IDs using UniProt ID mapping service
#' @param from_db Source database (e.g., "UniProtKB_AC-ID")
#' @param to_db Target database (e.g., "KEGG")
#' @param ids Vector of IDs to map
#'
#' @importFrom httr2 request req_perform req_body_form resp_body_json resp_body_raw resp_body_string resp_headers resp_status
#' @importFrom jsonlite fromJSON
#' @importFrom xml2 read_xml xml_find_all xml_ns xml_add_child xml_root
#' @importFrom stringr str_match
#'
#' @return Named vector of mapped IDs, names are original IDs
map_ids <- function(from_db, to_db, ids) {
  job_id <- submit_id_mapping(from_db, to_db, ids)
  if (check_id_mapping_results_ready(job_id)) {
    link <- get_id_mapping_results_link(job_id)
    response <- get_id_mapping_results_search(link)
    batch_results <- get_mapping(response)

    return(batch_results)
  } else {
    warning("No results found for the given IDs.")
    return(NULL)
  }
}

# From: https://www.uniprot.org/help/api_idmapping
# ---- HTTP session with retries ----
req_session <- function(url) {
  request(url) |>
    req_retry(max_tries = 5, backoff = ~ runif(1, 0.25, 0.5))
}

check_response <- function(resp) {
  if (resp_status(resp) >= 400) {
    print(resp_body_json(resp))
    stop("HTTP error: ", resp_status(resp))
  }
}

# ---- Submit ID mapping ----
submit_id_mapping <- function(from_db, to_db, ids) {
  resp <- req_session(paste0("https://rest.uniprot.org", "/idmapping/run")) |>
    req_body_form(from = from_db, to = to_db, ids = paste(ids, collapse = ",")) |>
    req_perform()

  check_response(resp)
  jobId <- resp_body_json(resp)$jobId
  return(jobId)
}

# ---- Poll for completion ----
check_id_mapping_results_ready <- function(job_id) {
  POLLING_INTERVAL <- 3

  while (TRUE) {
    resp <- req_session(paste0("https://rest.uniprot.org", "/idmapping/status/", job_id)) |> req_perform()
    check_response(resp)
    j <- resp_body_json(resp)
    if (!is.null(j$jobStatus)) {
      if (j$jobStatus %in% c("NEW", "RUNNING")) {
        message("Retrying in ", POLLING_INTERVAL, "s")
        Sys.sleep(POLLING_INTERVAL) # Try every POLLING_INTERVAL seconds
      } else {
        stop("Job failed: ", j$jobStatus)
      }
    } else {
      return(length(j$results) > 0 || length(j$failedIds) > 0)
    }
  }
}

# ---- Get results link ----
get_id_mapping_results_link <- function(job_id) {
  resp <- req_session(paste0("https://rest.uniprot.org", "/idmapping/details/", job_id)) |> req_perform()
  check_response(resp)
  return(resp_body_json(resp)$redirectURL)
}

# ---- Decode results ----
decode_results <- function(resp, file_format, compressed) {
  if (compressed) {
    raw_data <- resp_body_raw(resp)
    decompressed <- memDecompress(raw_data, type = "gzip")
    text <- rawToChar(decompressed)
    if (file_format == "json") {
      return(fromJSON(text))
    } else if (file_format == "tsv") {
      return(strsplit(text, "\n")[[1]])
    } else if (file_format == "xml") {
      return(text)
    } else {
      return(text)
    }
  } else {
    if (file_format == "json") {
      return(resp_body_json(resp))
    } else if (file_format == "tsv") {
      return(strsplit(resp_body_string(resp), "\n")[[1]])
    } else if (file_format == "xml") {
      return(resp_body_string(resp))
    }
    return(resp_body_string(resp))
  }
}

# ---- Merge XML results ----
merge_xml_results <- function(xml_results) {
  merged_root <- read_xml(xml_results[[1]])
  for (res in xml_results[-1]) {
    root <- read_xml(res)
    entries <- xml_find_all(root, ".//d1:entry", xml_ns(root))
    for (child in entries) {
      xml_add_child(xml_root(merged_root), child)
    }
  }
  return(as.character(merged_root))
}

# ---- Progress printing ----
print_progress_batches <- function(batch_index, size, total) {
  n_fetched <- min((batch_index + 1) * size, total)
  message("Fetched: ", n_fetched, " / ", total)
}

# ---- Get results (search endpoint) ----
get_id_mapping_results_search <- function(url) {
  parsed <- httr2::url_parse(url)
  query <- parsed$query
  file_format <- if (!is.null(query$format)) query$format else "json"
  size <- if (!is.null(query$size)) as.integer(query$size) else 500
  query$size <- size
  compressed <- if (!is.null(query$compressed)) tolower(query$compressed) == "true" else FALSE

  url_full <- httr2::url_build(parsed)

  resp <- req_session(url_full) |> req_perform()
  check_response(resp)

  results <- decode_results(resp, file_format, compressed)
  total <- as.integer(resp_headers(resp)[["x-total-results"]])
  print_progress_batches(0, size, total)

  # batch retrieval if 'Link' headers present
  link_header <- resp_headers(resp)[["Link"]]
  batch_url <- NULL
  if (!is.null(link_header)) {
    m <- str_match(link_header, "<(.+)>; rel=\"next\"")
    if (!is.na(m[1, 2])) batch_url <- m[1, 2]
  }

  i <- 1
  while (!is.null(batch_url)) {
    batch_resp <- req_session(batch_url) |> req_perform()
    batch <- decode_results(batch_resp, file_format, compressed)

    if (file_format == "json") {
      results$results <- c(results$results, batch$results)
      results$failedIds <- c(results$failedIds, batch$failedIds)
    } else {
      results <- c(results, batch)
    }

    print_progress_batches(i, size, total)
    i <- i + 1

    link_header <- resp_headers(batch_resp)[["Link"]]
    if (!is.null(link_header)) {
      m <- str_match(link_header, "<(.+)>; rel=\"next\"")
      batch_url <- if (!is.na(m[1, 2])) m[1, 2] else NULL
    } else {
      batch_url <- NULL
    }
  }

  if (file_format == "xml") {
    return(merge_xml_results(results))
  }
  return(results)
}

# ---- Simplify results ----
get_mapping <- function(response) {

  tos <- sapply(response$results, function(res) res$to)
  names(tos) <- sapply(response$results, function(res) res$from)
  return(tos)
}

# ---- Example usage ----
# data <- map_ids("UniProtKB_AC-ID", "KEGG", c("P05067","P12345"))
# batch_results <- get_mapping(data)
