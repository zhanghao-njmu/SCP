# Test dimensionality reduction functions

test_that("RunDimReduction basic test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunDimReduction - computationally intensive")

  # Would test multiple reduction methods
})

test_that("RunUMAP2 compatibility", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  mat <- matrix(rnorm(500), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
    srt <- Seurat::FindVariableFeatures(srt, verbose = FALSE, nfeatures = 20)
    srt <- Seurat::ScaleData(srt, verbose = FALSE)
    srt <- Seurat::RunPCA(srt, npcs = 5, verbose = FALSE)
  })

  # Test RunUMAP2
  result <- RunUMAP2(srt, dims = 1:5, verbose = FALSE)

  expect_s4_class(result, "Seurat")
  expect_true("umap" %in% names(result@reductions))
})

test_that("RunPHATE test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunPHATE - requires Python/reticulate")
})

test_that("RunPaCMAP test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunPaCMAP - requires additional package")
})

test_that("RunTriMap test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunTriMap - requires additional package")
})

test_that("RunDM diffusion map test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunDM - requires destiny package")
})

test_that("RunMDS multidimensional scaling", {
  skip_if_not_installed("Seurat")

  mat <- matrix(rnorm(200), nrow = 20)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  })

  # Test MDS
  result <- RunMDS(srt, dims = 2)

  expect_s4_class(result, "Seurat")
  expect_true("mds" %in% names(result@reductions))
})

test_that("RunNMF non-negative matrix factorization", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunNMF - requires NMF package")
})

test_that("RunGLMPCA test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunGLMPCA - requires glmpca package")
})

test_that("RunFR force-directed layout", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunFR - requires igraph")
})

test_that("RunLargeVis test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunLargeVis - requires additional package")
})

test_that("RunHarmony2 test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunHarmony2 - requires harmony package")
})
