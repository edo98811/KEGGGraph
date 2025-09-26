test_that("ggkegg_to_visnetwork returns a visNetwork object for human pathways", {
  skip("testing mouse")
  for (p in human_paths) {
    result <- ggkegg_to_visnetwork(p, test_de, feature_column = "rownames", organism = "hsa")
    expect_true(inherits(result, "visNetwork"))
    expect_true("nodes" %in% names(result$x))
    expect_true("edges" %in% names(result$x))
  }
})

test_that("ggkegg_to_visnetwork returns a visNetwork object for mouse pathways", {
  # skip("integration test")
  for (p in mouse_paths) {
    result <- ggkegg_to_visnetwork(p, test_de, organism = "mmu")
    expect_true(inherits(result, "visNetwork"))
    expect_true("nodes" %in% names(result$x))
    expect_true("edges" %in% names(result$x))
  }
})

# test_that("test get_nodes_df", {
#   nodes_df <- get_nodes_df(graph, "mmu")
# })

# test_that("test get_edges_df", {
#   nodes_df <- get_edges_df(graph)
# })

# test_that("test color_nodes", {
#   nodes_df <- color_nodes(nodes_df, nodes_color)
# })

# test_that("test style_edges", {
#   nodes_df <- style_edges(nodes_df, de_results, feature_column, value_column, threshold)
# })

# test_that("test style_nodes", {
#   nodes_df <- style_nodes(nodes_df, nodes_color)
# })
