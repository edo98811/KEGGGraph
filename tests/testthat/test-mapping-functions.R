# library(testthat)

# test_that("extract_kegg_ids returns correct KEGG IDs", {
#   mock_data <- list(
#     results = list(
#       list(
#         primaryAccession = "P12345",
#         genes = list(
#           list(geneId = list(type = "EntrezGene", value = "1234")),
#           list(geneId = list(type = "OtherType", value = "5678"))
#         )
#       ),
#       list(
#         primaryAccession = "Q67890",
#         genes = list(
#           list(geneId = list(type = "OtherType", value = "9999"))
#         )
#       ),
#       list(
#         primaryAccession = "A11111",
#         genes = NULL
#       )
#     )
#   )
#   organism_code <- "hsa"
#   result <- extract_kegg_ids(mock_data, organism_code)
#   expect_equal(result[["P12345"]], "hsa:1234")
#   expect_true(is.na(result[["Q67890"]]))
#   expect_true(is.na(result[["A11111"]]))
# })

# test_that("extract_kegg_ids handles empty results", {
#   mock_data <- list(results = list())
#   organism_code <- "eco"
#   result <- extract_kegg_ids(mock_data, organism_code)
#   expect_equal(length(result), 0)
# })

# test_that("extract_kegg_ids handles missing genes field", {
#   mock_data <- list(
#     results = list(
#       list(primaryAccession = "X00001")
#     )
#   )
#   organism_code <- "mmu"
#   result <- extract_kegg_ids(mock_data, organism_code)
#   expect_true(is.na(result[["X00001"]]))
# })

# test_that("symbol_to_kegg returns expected KEGG IDs for known gene symbols (human)", {
#   symbols <- c("TP53", "BRCA1")
#   kegg_ids <- symbol_to_kegg(symbols, organism = "human")
#   expect_true(all(!is.na(kegg_ids)))
#   expect_true(all(grepl("^hsa:", kegg_ids)))
#   expect_equal(names(kegg_ids), symbols)
# })

# test_that("symbol_to_kegg returns expected KEGG IDs for known gene symbols (mouse)", {
#   symbols <- c("Trp53", "Brca1")
#   kegg_ids <- symbol_to_kegg(symbols, organism = "mouse")
#   expect_true(all(!is.na(kegg_ids)))
#   expect_true(all(grepl("^mmu:", kegg_ids)))
#   expect_equal(names(kegg_ids), symbols)
# })

# test_that("ensembl_to_kegg returns expected KEGG IDs for known Ensembl IDs (human)", {
#   ensembl_ids <- c("ENSG00000141510", "ENSG00000012048") # TP53, BRCA1
#   kegg_ids <- ensembl_to_kegg(ensembl_ids, organism = "human")
#   expect_true(all(!is.na(kegg_ids)))
#   expect_true(all(grepl("^hsa:", kegg_ids)))
#   expect_equal(names(kegg_ids), ensembl_ids)
# })

# test_that("ensembl_to_kegg returns expected KEGG IDs for known Ensembl IDs (mouse)", {
#   ensembl_ids <- c("ENSMUSG00000059552", "ENSMUSG00000017146") # Trp53, Brca1
#   kegg_ids <- ensembl_to_kegg(ensembl_ids, organism = "mouse")
#   expect_true(all(!is.na(kegg_ids)))
#   expect_true(all(grepl("^mmu:", kegg_ids)))
#   expect_equal(names(kegg_ids), ensembl_ids)
# })