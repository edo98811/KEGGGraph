#' Ensemble ID → KEGG ID
#' 
#' Map a vector of Ensembl gene IDs to KEGG gene IDs using UniProt as an intermediary
#' 
#' @param ensembl_ids Vector of Ensembl gene IDs
#' @param organism Organism name 
#' 
#' @return Named vector of KEGG gene IDs, names are original Ensembl IDs
#' 
#' @importFrom AnnotationDbi mapIds
#' 
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' 
#' @export
ensembl_to_kegg <- function(ensembl_ids, organism = "human") {
  uniprot_ids <- ensembl_to_uniprot(ensembl_ids, organism = organism)
  kegg_ids <- uniprot_to_kegg_batch(uniprot_ids, organism = organism)
  return(kegg_ids)
}

#' Gene Symbol → KEGG
#' 
#' Map a vector of gene symbols to KEGG gene IDs using UniProt as an intermediary
#' 
#' @param symbols Vector of gene symbols
#' @param organism Organism name
#' 
#' @return Named vector of KEGG gene IDs, names are original gene symbols
#' 
#' @importFrom AnnotationDbi mapIds
#' 
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' 
#' @export
symbol_to_kegg <- function(symbols, organism = "human") {
  uniprot_ids <- symbol_to_uniprot(symbols, organism = organism)
  kegg_ids <- uniprot_to_kegg_batch(uniprot_ids, organism = organism)
  return(kegg_ids)
}


symbol_to_uniprot <- function(symbols, organism = "human") {
  db <- select_database(organism)
  uniprot_ids <- mapIds(db,
    keys = symbols,
    column = "UNIPROT",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  return(uniprot_ids)
}

ensembl_to_uniprot <- function(ensembl_ids, organism = "human") {
  db <- select_database(organism)
  uniprot_ids <- mapIds(db,
    keys = ensembl_ids,
    column = "UNIPROT",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  return(uniprot_ids)
}

uniprot_to_kegg_batch <- function(uniprot_ids, organism = "human", batch_size = 400) {
  organism <- to_organism_kegg(organism)

  all_results <- c()

  for (i in seq(1, length(uniprot_ids), by = batch_size)) {
    batch_ids <- uniprot_ids[i:min(i + batch_size - 1, length(uniprot_ids))]
    mapped_ids <- map_ids("UniProtKB_AC-ID", "KEGG", batch_ids)
    all_results <- c(all_results, mapped_ids)
  }

  # Ensure output order matches input
  all_results <- all_results[uniprot_ids]

  return(all_results)
}

select_database <- function(organism) {
  if (organism %in% c("human", "hsa")) {
    return(org.Hs.eg.db)
  } else if (organism %in% c("mouse", "mmu")) {
    return(org.Mm.eg.db)
  } else {
    stop("Unsupported organism. Please use 'human' or 'mouse'.")
  }
}

#' UniProt ID → KEGG ID
#' Map a vector of UniProt IDs to KEGG gene IDs using UniProt ID mapping service
#' 
#' @param uniprot_ids Vector of UniProt IDs
#' @param organism Organism name
#' 
#' @return Named vector of KEGG gene IDs, names are original UniProt IDs
#'  
#' @export
uniprot_to_kegg <- function(...) {
  return(uniprot_to_kegg_batch(...))
}

