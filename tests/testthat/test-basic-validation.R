# Test basic validation functions

test_that("check_R validates R version requirements", {
  # This should not throw an error if R version is compatible
  expect_no_error(check_R(R_version = "4.0.0"))

  # Test with future version (should warn or pass)
  expect_no_error(check_R(R_version = "5.0.0"))
})

test_that("check_srtList validates Seurat list inputs", {
  skip_if_not_installed("Seurat")

  # Create minimal test Seurat objects
  mat1 <- matrix(rpois(100, 5), nrow = 10)
  colnames(mat1) <- paste0("cell", 1:10)
  rownames(mat1) <- paste0("gene", 1:10)

  mat2 <- matrix(rpois(100, 5), nrow = 10)
  colnames(mat2) <- paste0("cell", 11:20)
  rownames(mat2) <- paste0("gene", 1:10)

  suppressWarnings({
    srt1 <- Seurat::CreateSeuratObject(counts = mat1, project = "test1")
    srt2 <- Seurat::CreateSeuratObject(counts = mat2, project = "test2")
  })

  srt_list <- list(srt1 = srt1, srt2 = srt2)

  # Test validation
  result <- check_srtList(srt_list)
  expect_type(result, "list")
  expect_length(result, 2)
})

test_that("package namespace is properly defined", {
  # Test that key functions are exported
  expect_true("CellDimPlot" %in% getNamespaceExports("SCP"))
  expect_true("FeatureDimPlot" %in% getNamespaceExports("SCP"))
  expect_true("RunCellQC" %in% getNamespaceExports("SCP"))
  expect_true("capitalize" %in% getNamespaceExports("SCP"))

  # Test pipe operator is re-exported
  expect_true("%>%" %in% getNamespaceExports("SCP"))
})

test_that("DESCRIPTION file dependencies are valid", {
  # Read DESCRIPTION file
  desc <- read.dcf(system.file("DESCRIPTION", package = "SCP"))

  # Check required fields
  expect_true("Package" %in% colnames(desc))
  expect_equal(desc[,"Package"], "SCP")
  expect_true("Version" %in% colnames(desc))
  expect_true("Imports" %in% colnames(desc))
  expect_true("Suggests" %in% colnames(desc))

  # Check that limma and monocle3 are in Suggests
  suggests <- desc[,"Suggests"]
  expect_true(grepl("limma", suggests))
  expect_true(grepl("monocle3", suggests))
})
