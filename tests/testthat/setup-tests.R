

test_de <- read.table(test_path("example_data.csv"))

# # Human KEGG pathway IDs
# human_pathways <- c(
#   "hsa04151",  # PI3K–Akt signaling pathway
#   "hsa04115",  # p53 signaling pathway
#   "hsa04014",  # Ras signaling pathway
#   "hsa05200",  # Pathways in cancer
#   "hsa03440"   # Homologous recombination
# )

# # Mouse KEGG pathway IDs
# mouse_pathways <- c(
#   "mmu04151",  # PI3K–Akt signaling pathway
#   "mmu04115",  # p53 signaling pathway
#   "mmu04014",  # Ras signaling pathway
#   "mmu05200",  # Pathways in cancer
#   "mmu03440"   # Homologous recombination
# )

mouse_paths <- c("mmu00760")

fake_de <- data.frame(
  logFC = rnorm(7, mean = 0, sd = 2),
  ensembl_id = c("ENSG00000139618", "ENSG00000157764", "ENSG00000141510", "ENSG00000155657", "ENSG00000171862", 
                 "ENSG00000198763", "ENSG00000142192"),
  uniprot_id = c("P38398", "P04637", "P31749", "Q9Y6K9", "P15056", 
                 "Q9Y243", "P01116"),
  gene_symbol = c("BRCA2", "TP53", "AKT1", "PIK3CA", "EGFR", 
                  "PTEN", "KRAS"),
  gene_symbol_mouse = c("Brca2", "Trp53", "Akt1", "Pik3ca", "Egfr",
                        "Pten", "Kras"),
  stringsAsFactors = FALSE
)

graph <- ggkegg::pathway("mmu00760")

fake_de_ensembl <- fake_de
rownames(fake_de_ensembl) <- fake_de_ensembl$ensembl_id

fake_de_uniprot <- fake_de
rownames(fake_de_uniprot) <- fake_de_uniprot$uniprot_id

fake_de_symbol <- fake_de
rownames(fake_de_symbol) <- fake_de_symbol$gene_symbol