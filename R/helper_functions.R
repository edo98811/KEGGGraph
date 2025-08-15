to_organism_kegg <- function(organism) {
  if (missing(organism) || !nzchar(organism)) {
    stop("You must provide a valid KEGG organism code.")
  }

  # Handle common abbreviations
  organism_code <- switch(tolower(organism),
    "hs" = "hsa",
    "mm" = "mmu",
    tolower(organism)
  )

  # Check if the organism code is valid
  valid_organisms <- c("hsa", "mmu") # Add more as needed
  if (!(organism_code %in% valid_organisms)) {
    stop(sprintf("Invalid KEGG organism code: %s. Valid codes are: %s",
                 organism_code, paste(valid_organisms, collapse = ", ")))
  }


  return(organism_code)
}
