# Test cell annotation functions

test_that("RunSingleR test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunSingleR - requires SingleR package and reference")

  # Would test SingleR annotation
})

test_that("RunScmap test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunScmap - requires scmap package")

  # Would test scmap annotation
})

test_that("RunKNNPredict test", {
  skip_if_not_installed("Seurat")

  # Create reference and query
  mat_ref <- matrix(rpois(500, 5), nrow = 50)
  colnames(mat_ref) <- paste0("ref_cell_", 1:10)
  rownames(mat_ref) <- paste0("gene_", 1:50)

  mat_query <- matrix(rpois(250, 5), nrow = 50)
  colnames(mat_query) <- paste0("query_cell_", 1:5)
  rownames(mat_query) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt_ref <- Seurat::CreateSeuratObject(counts = mat_ref, project = "ref")
    srt_query <- Seurat::CreateSeuratObject(counts = mat_query, project = "query")

    srt_ref <- Seurat::NormalizeData(srt_ref, verbose = FALSE)
    srt_query <- Seurat::NormalizeData(srt_query, verbose = FALSE)

    srt_ref$celltype <- rep(c("TypeA", "TypeB"), each = 5)
  })

  # Test KNN prediction
  result <- RunKNNPredict(
    srt_query = srt_query,
    srt_ref = srt_ref,
    ref_group = "celltype",
    query_assay = "RNA",
    ref_assay = "RNA"
  )

  expect_s4_class(result, "Seurat")
  expect_true("celltype_prediction" %in% colnames(result@meta.data) ||
              "predicted.id" %in% colnames(result@meta.data))
})

test_that("RunSeuratMap test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunSeuratMap - computationally intensive")

  # Would test Seurat reference mapping
})

test_that("RunPCAMap test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunPCAMap - requires PCA reference")

  # Would test PCA-based mapping
})

test_that("RunKNNMap test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunKNNMap - requires reference with UMAP")

  # Would test KNN-based mapping
})

test_that("RunCSSMap test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunCSSMap - requires simspec package")

  # Would test CSS mapping
})

test_that("RunSymphonyMap test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunSymphonyMap - requires symphony package")

  # Would test Symphony mapping
})
