# Test data loading and validation

test_that("package data objects are accessible", {
  # Test that example datasets can be loaded
  data("pancreas_sub", package = "SCP", envir = environment())
  expect_true(exists("pancreas_sub"))

  # Check if it's a Seurat object
  expect_s4_class(pancreas_sub, "Seurat")

  # Test palette list
  data("palette_list", package = "SCP", envir = environment())
  expect_true(exists("palette_list"))
  expect_type(palette_list, "list")
  expect_true(length(palette_list) > 0)
})

test_that("check_DataType identifies data types correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  # Create a simple test matrix
  test_matrix <- Matrix::Matrix(c(1.5, 2.3, 3.7, 4.2), nrow = 2, ncol = 2, sparse = TRUE)
  colnames(test_matrix) <- c("cell1", "cell2")
  rownames(test_matrix) <- c("gene1", "gene2")

  # Create a minimal Seurat object
  suppressWarnings({
    test_srt <- Seurat::CreateSeuratObject(counts = test_matrix, project = "test")
  })

  # Test check_DataType function
  result <- check_DataType(test_srt)
  expect_type(result, "list")
  expect_true("datatype" %in% names(result))
})
