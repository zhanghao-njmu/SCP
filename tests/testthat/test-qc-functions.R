# Test QC Functions

test_that("RunCellQC basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  # Create test Seurat object
  set.seed(123)
  mat <- Matrix::Matrix(rpois(1000, 5), nrow = 100, sparse = TRUE)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:100)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
  })

  # Add mitochondrial genes
  rownames(srt)[1:10] <- paste0("MT-", rownames(srt)[1:10])

  # Test basic QC without doublet detection
  result <- RunCellQC(
    srt = srt,
    db_method = "scDblFinder",
    db_rate = 0.05,
    return_filtered = FALSE
  )

  expect_s4_class(result, "Seurat")
  expect_true("CellQC" %in% colnames(result@meta.data))
  expect_true("nFeature_RNA" %in% colnames(result@meta.data))
  expect_true("nCount_RNA" %in% colnames(result@meta.data))
})

test_that("RunCellQC with filtering", {
  skip_if_not_installed("Seurat")

  mat <- matrix(rpois(500, 5), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
  })

  # Test with return_filtered = TRUE
  result <- RunCellQC(
    srt = srt,
    return_filtered = TRUE,
    min_features = 5,
    max_features = 100
  )

  expect_s4_class(result, "Seurat")
  expect_true(ncol(result) <= ncol(srt))
})

test_that("isOutlier function works", {
  # Test basic outlier detection
  x <- c(1, 2, 3, 100, 2, 3, 4)
  result <- isOutlier(x, nmads = 3)

  expect_type(result, "logical")
  expect_length(result, length(x))
  expect_true(any(result))  # Should detect 100 as outlier
})

test_that("isOutlier handles edge cases", {
  # All same values
  x <- rep(5, 10)
  result <- isOutlier(x)
  expect_type(result, "logical")
  expect_false(any(result))

  # With NA values
  x <- c(1, 2, NA, 3, 100)
  result <- isOutlier(x)
  expect_type(result, "logical")

  # Single value
  x <- 5
  result <- isOutlier(x)
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("doublet detection methods", {
  skip_if_not_installed("Seurat")
  skip("Skipping doublet detection - requires additional packages")

  # These tests would require scDblFinder, scds, etc.
  # Documented for future implementation
})

test_that("CC_GenePrefetch works", {
  skip_if_not_installed("Seurat")

  # Test cell cycle gene prefetch
  result <- CC_GenePrefetch(species = "Homo_sapiens")

  expect_type(result, "list")
  expect_true("s.genes" %in% names(result) || "S.genes" %in% names(result))
  expect_true("g2m.genes" %in% names(result) || "G2M.genes" %in% names(result))
})

test_that("CellScoring basic test", {
  skip_if_not_installed("Seurat")

  mat <- matrix(rpois(500, 5), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  })

  # Test with feature list
  features <- list(test_sig = rownames(srt)[1:5])

  result <- CellScoring(
    srt = srt,
    features = features,
    method = "Seurat"
  )

  expect_s4_class(result, "Seurat")
  expect_true("test_sig" %in% colnames(result@meta.data))
})

test_that("RecoverCounts works", {
  skip_if_not_installed("Seurat")

  mat <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
  })

  # Test count recovery
  result <- RecoverCounts(srt)

  expect_s4_class(result, "Seurat")
})
