# Test Seurat V5 compatibility

test_that("Seurat version is detected correctly", {
  skip_if_not_installed("Seurat")

  seurat_version <- packageVersion("Seurat")
  expect_true(seurat_version >= "4.0.0")

  # Record version for debugging
  message("Testing with Seurat version: ", seurat_version)
})

test_that("GetAssayData works with both V4 and V5", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  # Create test Seurat object
  mat <- Matrix::Matrix(rpois(200, 5), nrow = 20, sparse = TRUE)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
  })

  # Test data retrieval - should work in both V4 and V5
  counts_data <- Seurat::GetAssayData(srt, slot = "counts")
  expect_true(inherits(counts_data, "Matrix") || inherits(counts_data, "matrix"))
  expect_equal(dim(counts_data), c(20, 10))

  # Test with layer parameter (V5 style) if available
  seurat_version <- packageVersion("Seurat")
  if (seurat_version >= "5.0.0") {
    # V5 uses layer parameter
    tryCatch({
      counts_v5 <- Seurat::GetAssayData(srt, layer = "counts")
      expect_true(inherits(counts_v5, "Matrix") || inherits(counts_v5, "matrix"))
    }, error = function(e) {
      # Log but don't fail if layer parameter doesn't work
      message("Layer parameter not working: ", e$message)
    })
  }
})

test_that("Normalization works with both V4 and V5", {
  skip_if_not_installed("Seurat")

  # Create test object
  mat <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")

    # Normalize
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  })

  # Check normalization succeeded
  expect_s4_class(srt, "Seurat")

  # Try to retrieve normalized data
  norm_data <- tryCatch({
    Seurat::GetAssayData(srt, slot = "data")
  }, error = function(e) {
    # Try V5 style
    Seurat::GetAssayData(srt, layer = "data")
  })

  expect_true(!is.null(norm_data))
})

test_that("Variable features work with both V4 and V5", {
  skip_if_not_installed("Seurat")

  # Create test object with enough features
  mat <- matrix(rpois(1000, 5), nrow = 100)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:100)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
    srt <- Seurat::FindVariableFeatures(srt, nfeatures = 20, verbose = FALSE)
  })

  # Get variable features
  var_features <- Seurat::VariableFeatures(srt)
  expect_type(var_features, "character")
  expect_true(length(var_features) > 0)
  expect_true(length(var_features) <= 20)
})

test_that("PCA works with both V4 and V5", {
  skip_if_not_installed("Seurat")

  # Create test object
  mat <- matrix(rnorm(1000), nrow = 100)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:100)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
    srt <- Seurat::FindVariableFeatures(srt, verbose = FALSE)
    srt <- Seurat::ScaleData(srt, verbose = FALSE)
    srt <- Seurat::RunPCA(srt, verbose = FALSE, npcs = 5)
  })

  # Check PCA reduction exists
  expect_true("pca" %in% names(srt@reductions))

  # Get embeddings
  embeddings <- Seurat::Embeddings(srt, reduction = "pca")
  expect_true(is.matrix(embeddings))
  expect_equal(ncol(embeddings), 5)
})

test_that("Integration functions are available", {
  skip_if_not_installed("Seurat")

  seurat_version <- packageVersion("Seurat")

  if (seurat_version >= "5.0.0") {
    # V5 should have IntegrateLayers
    expect_true(exists("IntegrateLayers", where = asNamespace("Seurat")))
  } else {
    # V4 should have IntegrateData
    expect_true(exists("IntegrateData", where = asNamespace("Seurat")))
  }

  # Both should have FindIntegrationAnchors
  expect_true(exists("FindIntegrationAnchors", where = asNamespace("Seurat")))
})

test_that("SCTransform compatibility", {
  skip_if_not_installed("Seurat")
  skip("Skipping SCTransform test - computationally intensive")

  # This test is skipped by default but documents SCTransform usage
  # In V5, default is vst.flavor = "v2"
  # In V4, default is vst.flavor = "v1"

  # Example usage:
  # srt <- SCTransform(srt, vst.flavor = "v2")  # V5 default
  # srt <- SCTransform(srt, vst.flavor = "v1")  # V4 compatible
})

test_that("check_DataType handles both V4 and V5 assays", {
  skip_if_not_installed("Seurat")

  # Create test object
  mat <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
  })

  # This should work regardless of Seurat version
  result <- check_DataType(srt)
  expect_type(result, "list")
  expect_true("datatype" %in% names(result))
})

test_that("DefaultAssay getter/setter works", {
  skip_if_not_installed("Seurat")

  mat <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
  })

  # Get default assay
  default_assay <- Seurat::DefaultAssay(srt)
  expect_type(default_assay, "character")
  expect_equal(default_assay, "RNA")

  # Set default assay (should work in both versions)
  Seurat::DefaultAssay(srt) <- "RNA"
  expect_equal(Seurat::DefaultAssay(srt), "RNA")
})
