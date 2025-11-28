# Enhanced tests for statistical plotting functions

test_that("FeatureStatPlot handles different plot types", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  # Create test object
  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- factor(rep(c("A", "B"), each = 10))

  # Test violin plot
  expect_no_error(
    p1 <- FeatureStatPlot(
      srt = srt,
      stat.by = "gene1",
      group.by = "group",
      plot_type = "violin"
    )
  )
  expect_s3_class(p1, "gg")

  # Test box plot
  expect_no_error(
    p2 <- FeatureStatPlot(
      srt = srt,
      stat.by = "gene1",
      group.by = "group",
      plot_type = "box"
    )
  )
  expect_s3_class(p2, "gg")

  # Test bar plot
  expect_no_error(
    p3 <- FeatureStatPlot(
      srt = srt,
      stat.by = "gene1",
      group.by = "group",
      plot_type = "bar"
    )
  )
  expect_s3_class(p3, "gg")
})

test_that("FeatureStatPlot handles multiple features", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- factor(rep(c("A", "B"), each = 10))

  # Test with multiple features
  expect_no_error(
    p <- FeatureStatPlot(
      srt = srt,
      stat.by = c("gene1", "gene2", "gene3"),
      group.by = "group"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("CellStatPlot works with metadata", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- factor(rep(c("A", "B"), each = 10))
  srt$nFeature_RNA <- colSums(counts > 0)
  srt$nCount_RNA <- colSums(counts)

  # Test with QC metric
  expect_no_error(
    p <- CellStatPlot(
      srt = srt,
      stat.by = "nFeature_RNA",
      group.by = "group"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("StatPlot handles various statistical measures", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$batch <- factor(rep(c("batch1", "batch2"), each = 10))
  srt$condition <- factor(rep(c("WT", "KO", "WT", "KO"), each = 5))

  # Test with different grouping
  expect_no_error(
    p1 <- StatPlot(
      srt = srt,
      stat.by = "gene1",
      group.by = "batch"
    )
  )
  expect_s3_class(p1, "gg")

  # Test with split.by
  expect_no_error(
    p2 <- StatPlot(
      srt = srt,
      stat.by = "gene1",
      group.by = "batch",
      split.by = "condition"
    )
  )
  expect_s3_class(p2, "gg")
})

test_that("FeatureCorPlot creates correlation plots", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)

  # Test correlation between two features
  expect_no_error(
    p <- FeatureCorPlot(
      srt = srt,
      features = c("gene1", "gene2")
    )
  )
  expect_s3_class(p, "gg")
})

test_that("FeatureCorPlot handles multiple feature pairs", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(400, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:40)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- factor(rep(c("A", "B"), each = 20))

  # Test with group coloring
  expect_no_error(
    p <- FeatureCorPlot(
      srt = srt,
      features = c("gene1", "gene2"),
      group.by = "group"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("Statistical plots handle edge cases", {
  skip_if_not_installed("Seurat")

  # Small dataset
  counts <- matrix(rpois(30, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:3)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- factor(c("A", "A", "B"))

  # Should handle small n
  expect_no_error(
    FeatureStatPlot(
      srt = srt,
      stat.by = "gene1",
      group.by = "group"
    )
  )
})

test_that("Statistical plots validate inputs", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)

  # Test with missing feature
  expect_error(
    FeatureStatPlot(
      srt = srt,
      stat.by = "nonexistent_gene",
      group.by = "group"
    )
  )

  # Test with missing metadata
  expect_error(
    CellStatPlot(
      srt = srt,
      stat.by = "nFeature_RNA",
      group.by = "nonexistent_column"
    )
  )
})

test_that("Statistical plots work with normalized data", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B"), each = 10))

  # Test with normalized data
  expect_no_error(
    p <- FeatureStatPlot(
      srt = srt,
      stat.by = "gene1",
      group.by = "group",
      slot = "data"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("CellDensityPlot creates density visualizations", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- factor(rep(c("A", "B"), each = 10))

  # Add dimensionality reduction
  srt[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "PC_",
    assay = "RNA"
  )

  # Test density plot
  expect_no_error(
    p <- CellDensityPlot(
      srt = srt,
      reduction = "pca"
    )
  )
  expect_s3_class(p, "gg")
})
