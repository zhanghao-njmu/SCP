# Enhanced tests for panel functions and dimensionality reduction plots

test_that("panel_fix adjusts plot panels correctly", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  # Create test plots
  p1 <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
    ggplot2::geom_point()

  p2 <- ggplot2::ggplot(mtcars, ggplot2::aes(x = hp, y = mpg)) +
    ggplot2::geom_point()

  # Combine plots
  combined <- p1 + p2

  # Test panel_fix
  expect_no_error(fixed <- panel_fix(combined, nrow = 1, ncol = 2))
  expect_s3_class(fixed, "patchwork")
})

test_that("panel_fix handles single plots", {
  skip_if_not_installed("ggplot2")

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
    ggplot2::geom_point()

  # Should handle single plot
  expect_no_error(fixed <- panel_fix(p))
  expect_s3_class(fixed, "gg")
})

test_that("panel_fix_overall adjusts panel dimensions", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  p1 <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
    ggplot2::geom_point()

  p2 <- ggplot2::ggplot(mtcars, ggplot2::aes(x = hp, y = qsec)) +
    ggplot2::geom_point()

  combined <- p1 + p2

  # Test overall panel fix
  expect_no_error(
    fixed <- panel_fix_overall(
      combined,
      width = 10,
      height = 5
    )
  )
})

test_that("CellDimPlot handles different reduction types", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- factor(rep(c("TypeA", "TypeB"), each = 10))

  # Add PCA
  srt[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "PC_",
    assay = "RNA"
  )

  # Add UMAP
  srt[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "UMAP_",
    assay = "RNA"
  )

  # Test with PCA
  expect_no_error(
    p1 <- CellDimPlot(
      srt = srt,
      group.by = "celltype",
      reduction = "pca"
    )
  )
  expect_s3_class(p1, "gg")

  # Test with UMAP
  expect_no_error(
    p2 <- CellDimPlot(
      srt = srt,
      group.by = "celltype",
      reduction = "umap"
    )
  )
  expect_s3_class(p2, "gg")
})

test_that("CellDimPlot supports custom colors", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- factor(rep(c("A", "B", "C", "D"), each = 5))

  srt[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "PC_",
    assay = "RNA"
  )

  # Test with custom palette
  expect_no_error(
    p <- CellDimPlot(
      srt = srt,
      group.by = "group",
      reduction = "pca",
      palette = "Set1"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("CellDimPlot handles split.by parameter", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(400, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:40)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- factor(rep(c("TypeA", "TypeB"), times = 20))
  srt$condition <- factor(rep(c("WT", "KO"), each = 20))

  srt[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(80), ncol = 2),
    key = "PC_",
    assay = "RNA"
  )

  # Test with split.by
  expect_no_error(
    p <- CellDimPlot(
      srt = srt,
      group.by = "celltype",
      split.by = "condition",
      reduction = "pca"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("FeatureDimPlot displays feature expression", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)

  srt[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "PC_",
    assay = "RNA"
  )

  # Test feature plot
  expect_no_error(
    p <- FeatureDimPlot(
      srt = srt,
      features = "gene1",
      reduction = "pca"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("FeatureDimPlot handles multiple features", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)

  srt[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "PC_",
    assay = "RNA"
  )

  # Test with multiple features
  expect_no_error(
    p <- FeatureDimPlot(
      srt = srt,
      features = c("gene1", "gene2", "gene3"),
      reduction = "pca"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("FeatureDimPlot supports different color scales", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)

  srt[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "PC_",
    assay = "RNA"
  )

  # Test with custom color scale
  expect_no_error(
    p <- FeatureDimPlot(
      srt = srt,
      features = "gene1",
      reduction = "pca",
      palette = "viridis"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("DimPlots handle labels and highlights", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$cluster <- factor(rep(c("C1", "C2", "C3", "C4"), each = 5))

  srt[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "PC_",
    assay = "RNA"
  )

  # Test with labels
  expect_no_error(
    p1 <- CellDimPlot(
      srt = srt,
      group.by = "cluster",
      reduction = "pca",
      label = TRUE
    )
  )
  expect_s3_class(p1, "gg")

  # Test with highlights
  expect_no_error(
    p2 <- CellDimPlot(
      srt = srt,
      group.by = "cluster",
      reduction = "pca",
      highlight = c("C1", "C2")
    )
  )
  expect_s3_class(p2, "gg")
})

test_that("DimPlots validate inputs", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)

  srt[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "PC_",
    assay = "RNA"
  )

  # Test with missing group.by
  expect_error(
    CellDimPlot(
      srt = srt,
      group.by = "missing_column",
      reduction = "pca"
    )
  )

  # Test with missing reduction
  expect_error(
    CellDimPlot(
      srt = srt,
      group.by = "RNA_snn_res.0.5",
      reduction = "missing_reduction"
    )
  )
})

test_that("DimPlots handle cell subsetting", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(400, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:40)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- factor(rep(c("A", "B", "C", "D"), each = 10))

  srt[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(80), ncol = 2),
    key = "PC_",
    assay = "RNA"
  )

  # Test with cell subset
  cells_subset <- colnames(srt)[1:20]

  expect_no_error(
    p <- CellDimPlot(
      srt = srt,
      group.by = "group",
      reduction = "pca",
      cells = cells_subset
    )
  )
  expect_s3_class(p, "gg")
})
