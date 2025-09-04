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
    vals <- unlist(lapply(list_of_lists, function(lst) lst[[k]]))
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
    missing <- setdiff(all_cols, names(df))
    if(length(missing) > 0) {
      df[missing] <- NA
    }
    df[, all_cols]
  }

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
get_pathway_name <- function(pathway_id, organism = NULL) {
  # Input validation
  if (is.null(pathway_id) || !is.character(pathway_id) || length(pathway_id) != 1) {
    stop("pathway_id must be a single character string")
  }

  # Extract organism code if present in pathway_id
  if (is.null(organism)) {
    if (grepl("^[a-z]{2,3}\\d{5}$", pathway_id)) {
      organism <- sub("(^[a-z]{2,3})(\\d{5}$)", "\\1", pathway_id)
      pathway_number <- sub("(^[a-z]{2,3})(\\d{5}$)", "\\2", pathway_id)
    } else if (grepl("^\\d{5}$", pathway_id)) {
      pathway_number <- pathway_id
      organism <- "map" # Use generic map
    } else {
      stop("Invalid pathway_id format. Expected formats: 'hsa04110' or '04110'")
    }
  } else {
    # Remove organism prefix if present in pathway_id
    pathway_number <- sub("^[a-z]{2,3}", "", pathway_id)
    if (!grepl("^\\d{5}$", pathway_number)) {
      stop("Invalid pathway number format. Expected 5 digits.")
    }
  }

  # Predefined mapping of common KEGG pathways
  pathway_names <- list(
    # Metabolism pathways
    "00010" = "Glycolysis / Gluconeogenesis",
    "00020" = "Citrate cycle (TCA cycle)",
    "00030" = "Pentose phosphate pathway",
    "00040" = "Pentose and glucuronate interconversions",
    "00051" = "Fructose and mannose metabolism",
    "00052" = "Galactose metabolism",
    "00053" = "Ascorbate and aldarate metabolism",
    "00071" = "Fatty acid degradation",
    "00072" = "Synthesis and degradation of ketone bodies",
    "00190" = "Oxidative phosphorylation",
    "00230" = "Purine metabolism",
    "00240" = "Pyrimidine metabolism",
    "00250" = "Alanine, aspartate and glutamate metabolism",
    "00260" = "Glycine, serine and threonine metabolism",
    "00270" = "Cysteine and methionine metabolism",
    "00280" = "Valine, leucine and isoleucine degradation",
    "00290" = "Valine, leucine and isoleucine biosynthesis",

    # Signaling pathways
    "04010" = "MAPK signaling pathway",
    "04012" = "ErbB signaling pathway",
    "04014" = "Ras signaling pathway",
    "04015" = "Rap1 signaling pathway",
    "04020" = "Calcium signaling pathway",
    "04024" = "cAMP signaling pathway",
    "04062" = "Chemokine signaling pathway",
    "04064" = "NF-kappa B signaling pathway",
    "04066" = "HIF-1 signaling pathway",
    "04068" = "FoxO signaling pathway",
    "04070" = "Phosphatidylinositol signaling system",
    "04071" = "Sphingolipid signaling pathway",
    "04072" = "Phospholipase D signaling pathway",
    "04110" = "Cell cycle",
    "04115" = "p53 signaling pathway",
    "04120" = "Ubiquitin mediated proteolysis",
    "04136" = "Autophagy",
    "04137" = "Mitophagy",
    "04140" = "Autophagy - animal",
    "04141" = "Protein processing in endoplasmic reticulum",
    "04142" = "Lysosome",
    "04144" = "Endocytosis",
    "04145" = "Phagosome",
    "04146" = "Peroxisome",
    "04150" = "mTOR signaling pathway",
    "04151" = "PI3K-Akt signaling pathway",
    "04152" = "AMPK signaling pathway",
    "04210" = "Apoptosis",
    "04211" = "Longevity regulating pathway",
    "04212" = "Longevity regulating pathway - worm",
    "04213" = "Longevity regulating pathway - multiple species",
    "04214" = "Apoptosis - fly",
    "04215" = "Apoptosis - multiple species",
    "04216" = "Ferroptosis",
    "04217" = "Necroptosis",
    "04218" = "Cellular senescence",

    # Immune system
    "04610" = "Complement and coagulation cascades",
    "04611" = "Platelet activation",
    "04612" = "Antigen processing and presentation",
    "04620" = "Toll-like receptor signaling pathway",
    "04621" = "NOD-like receptor signaling pathway",
    "04622" = "RIG-I-like receptor signaling pathway",
    "04623" = "Cytosolic DNA-sensing pathway",
    "04624" = "Toll and Imd signaling pathway",
    "04625" = "C-type lectin receptor signaling pathway",
    "04630" = "JAK-STAT signaling pathway",
    "04640" = "Hematopoietic cell lineage",
    "04650" = "Natural killer cell mediated cytotoxicity",
    "04660" = "T cell receptor signaling pathway",
    "04662" = "B cell receptor signaling pathway",
    "04664" = "Fc epsilon RI signaling pathway",
    "04666" = "Fc gamma R-mediated phagocytosis",
    "04668" = "TNF signaling pathway",
    "04670" = "Leukocyte transendothelial migration",
    "04672" = "Intestinal immune network for IgA production",

    # Disease pathways
    "05010" = "Alzheimer disease",
    "05012" = "Parkinson disease",
    "05014" = "Amyotrophic lateral sclerosis",
    "05016" = "Huntington disease",
    "05017" = "Spinocerebellar ataxia",
    "05020" = "Prion disease",
    "05022" = "Pathways of neurodegeneration - multiple diseases",
    "05030" = "Cocaine addiction",
    "05031" = "Amphetamine addiction",
    "05032" = "Morphine addiction",
    "05033" = "Nicotine addiction",
    "05034" = "Alcoholism",
    "05200" = "Pathways in cancer",
    "05202" = "Transcriptional misregulation in cancer",
    "05203" = "Viral carcinogenesis",
    "05204" = "Chemical carcinogenesis",
    "05205" = "Proteoglycans in cancer",
    "05206" = "MicroRNAs in cancer",
    "05207" = "Chemical carcinogenesis - receptor activation",
    "05208" = "Chemical carcinogenesis - reactive oxygen species",
    "05210" = "Colorectal cancer",
    "05211" = "Renal cell carcinoma",
    "05212" = "Pancreatic cancer",
    "05213" = "Endometrial cancer",
    "05214" = "Glioma",
    "05215" = "Prostate cancer",
    "05216" = "Thyroid cancer",
    "05217" = "Basal cell carcinoma",
    "05218" = "Melanoma",
    "05219" = "Bladder cancer",
    "05220" = "Chronic myeloid leukemia",
    "05221" = "Acute myeloid leukemia",
    "05222" = "Small cell lung cancer",
    "05223" = "Non-small cell lung cancer",
    "05224" = "Breast cancer",
    "05225" = "Hepatocellular carcinoma",
    "05226" = "Gastric cancer",

    # Endocrine system
    "04910" = "Insulin signaling pathway",
    "04911" = "Insulin secretion",
    "04912" = "GnRH signaling pathway",
    "04913" = "Ovarian steroidogenesis",
    "04914" = "Progesterone-mediated oocyte maturation",
    "04915" = "Estrogen signaling pathway",
    "04916" = "Melanogenesis",
    "04917" = "Prolactin signaling pathway",
    "04918" = "Thyroid hormone synthesis",
    "04919" = "Thyroid hormone signaling pathway",
    "04920" = "Adipocytokine signaling pathway",
    "04921" = "Oxytocin signaling pathway",
    "04922" = "Glucagon signaling pathway",
    "04923" = "Regulation of lipolysis in adipocytes",
    "04924" = "Renin secretion",
    "04925" = "Aldosterone synthesis and secretion",
    "04926" = "Relaxin signaling pathway",
    "04927" = "Cortisol synthesis and secretion",
    "04928" = "Parathyroid hormone synthesis, secretion and action",
    "04929" = "GnRH secretion",
    "04930" = "Type II diabetes mellitus",
    "04931" = "Insulin resistance",
    "04932" = "Non-alcoholic fatty liver disease",
    "04933" = "AGE-RAGE signaling pathway in diabetic complications",
    "04934" = "Cushing syndrome",
    "04935" = "Growth hormone synthesis, secretion and action",
    "04936" = "Alcoholic liver disease"
  )

  # Get pathway name
  pathway_name <- pathway_names[[pathway_number]]

  if (is.null(pathway_name)) {
    # If pathway not found in our mapping, try to format it nicely
    warning(paste("Pathway", pathway_number, "not found in predefined mappings. Using generic format."))
    pathway_name <- paste("KEGG Pathway", pathway_number)
  }

  # Add organism prefix if specified
  if (!is.null(organism) && organism != "map") {
    organism_names <- list(
      "hsa" = "Homo sapiens",
      "mmu" = "Mus musculus",
      "rno" = "Rattus norvegicus",
      "dme" = "Drosophila melanogaster",
      "cel" = "Caenorhabditis elegans",
      "sce" = "Saccharomyces cerevisiae",
      "eco" = "Escherichia coli",
      "bta" = "Bos taurus",
      "ssc" = "Sus scrofa",
      "gga" = "Gallus gallus",
      "dre" = "Danio rerio",
      "ath" = "Arabidopsis thaliana"
    )

    organism_full <- organism_names[[organism]]
    if (is.null(organism_full)) {
      organism_full <- toupper(organism)
    }

    pathway_name <- paste0(pathway_name, " (", organism_full, ")")
  }

  return(pathway_name)
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
