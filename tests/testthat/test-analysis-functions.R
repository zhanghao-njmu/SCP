# Test differential expression and analysis functions

test_that("RunDEtest basic functionality", {
  skip_if_not_installed("Seurat")

  mat <- matrix(rpois(500, 5), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
    srt$group <- rep(c("A", "B"), each = 5)
  })

  # Test basic DE
  result <- RunDEtest(
    srt = srt,
    group_by = "group",
    fc.threshold = 1.2,
    only.pos = FALSE
  )

  expect_type(result, "list")
  expect_true(length(result) > 0)
})

test_that("FindExpressedMarkers test", {
  skip_if_not_installed("Seurat")

  mat <- matrix(rpois(500, 5), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
    srt$group <- rep(c("A", "B"), each = 5)
  })

  # Test expressed markers
  result <- FindExpressedMarkers(
    srt = srt,
    group_by = "group",
    min.pct = 0.1
  )

  expect_type(result, "list")
})

test_that("RunEnrichment basic test", {
  skip_if_not_installed("clusterProfiler")
  skip("Skipping RunEnrichment - requires database connection")

  # This would test enrichment analysis
  # Requires actual gene lists and databases
})

test_that("RunGSEA basic test", {
  skip_if_not_installed("clusterProfiler")
  skip("Skipping RunGSEA - requires database connection")

  # This would test GSEA analysis
})

test_that("RunDynamicFeatures test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunDynamicFeatures - requires trajectory data")

  # Would test dynamic feature detection along trajectories
})

test_that("RunDynamicEnrichment test", {
  skip_if_not_installed("Seurat")
  skip("Skipping RunDynamicEnrichment - requires trajectory data")

  # Would test enrichment along trajectories
})

test_that("GeneConvert basic test", {
  skip("Skipping GeneConvert - requires database connection")

  # This would test gene ID conversion
  # Requires biomaRt connection
})

test_that("AnnotateFeatures test", {
  skip_if_not_installed("Seurat")
  skip("Skipping AnnotateFeatures - requires database")

  # Would test feature annotation
})
