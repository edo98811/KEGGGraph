test_that("merge color lists", {
  # skip("integration test")
  # Test with empty list
  expect_equal(merge_color_lists(list()), list())

  # Test with 3 lists with elements with different names
  l1 <- list(a = "red")
  l2 <- list(b = "blue")
  l3 <- list(c = "green")
  expect_equal(
    merge_color_lists(list(l1, l2, l3)),
    list(a = "red", b = "blue", c = "green")
  )

  # Test with 3 lists with elements with the same names
  l4 <- list(a = "white", b = "yellow")
  l5 <- list(a = "orange", b = "white")
  l6 <- list(a = "white", b = "white")
  expect_equal(
    merge_color_lists(list(l4, l5, l6)),
    list(a = "orange", b = "yellow")
  )
})