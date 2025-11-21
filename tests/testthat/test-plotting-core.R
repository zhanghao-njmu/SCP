# Test core plotting functions

test_that("CellDimPlot basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  mat <- matrix(rnorm(500), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
    srt <- Seurat::FindVariableFeatures(srt, verbose = FALSE, nfeatures = 20)
    srt <- Seurat::ScaleData(srt, verbose = FALSE)
    srt <- Seurat::RunPCA(srt, npcs = 5, verbose = FALSE)

    # Add UMAP coordinates manually
    umap_coords <- matrix(rnorm(20), ncol = 2)
    rownames(umap_coords) <- colnames(srt)
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
    srt[["umap"]] <- Seurat::CreateDimReducObject(
      embeddings = umap_coords,
      key = "UMAP_",
      assay = "RNA"
    )

    srt$group <- rep(c("A", "B"), each = 5)
  })

  # Test basic plotting
  p <- CellDimPlot(
    srt = srt,
    group.by = "group",
    reduction = "umap"
  )

  expect_s3_class(p, "ggplot")
})

test_that("FeatureDimPlot basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  mat <- matrix(rnorm(500), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)

    # Add UMAP coordinates
    umap_coords <- matrix(rnorm(20), ncol = 2)
    rownames(umap_coords) <- colnames(srt)
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
    srt[["umap"]] <- Seurat::CreateDimReducObject(
      embeddings = umap_coords,
      key = "UMAP_",
      assay = "RNA"
    )
  })

  # Test feature plotting
  p <- FeatureDimPlot(
    srt = srt,
    features = rownames(srt)[1:2],
    reduction = "umap"
  )

  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("CellStatPlot basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  mat <- matrix(rpois(500, 5), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt$group <- rep(c("A", "B"), each = 5)
  })

  # Test stat plotting
  p <- CellStatPlot(
    srt = srt,
    stat.by = c("nFeature_RNA", "nCount_RNA"),
    group.by = "group"
  )

  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("FeatureStatPlot basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  mat <- matrix(rpois(500, 5), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
    srt$group <- rep(c("A", "B"), each = 5)
  })

  # Test feature stat plotting
  p <- FeatureStatPlot(
    srt = srt,
    features = rownames(srt)[1:3],
    group.by = "group"
  )

  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("StatPlot basic functionality", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    x = rep(c("A", "B", "C"), each = 10),
    y = rnorm(30)
  )

  # Test basic stat plot
  p <- StatPlot(
    data = df,
    stat.by = "y",
    group.by = "x",
    plot.by = "x"
  )

  expect_s3_class(p, "ggplot")
})

test_that("VolcanoPlot basic functionality", {
  skip_if_not_installed("ggplot2")
  skip("Skipping VolcanoPlot - requires DE results")

  # Would test volcano plot with DE results
})

test_that("GroupHeatmap basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")
  skip("Skipping GroupHeatmap - computationally intensive")

  # Would test heatmap generation
})

test_that("FeatureHeatmap basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")
  skip("Skipping FeatureHeatmap - computationally intensive")

  # Would test feature heatmap
})

test_that("CellCorHeatmap basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")
  skip("Skipping CellCorHeatmap - computationally intensive")

  # Would test cell correlation heatmap
})

test_that("DynamicHeatmap basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")
  skip("Skipping DynamicHeatmap - requires trajectory data")

  # Would test dynamic heatmap
})

test_that("EnrichmentPlot basic functionality", {
  skip_if_not_installed("ggplot2")
  skip("Skipping EnrichmentPlot - requires enrichment results")

  # Would test enrichment visualization
})

test_that("GSEAPlot basic functionality", {
  skip_if_not_installed("ggplot2")
  skip("Skipping GSEAPlot - requires GSEA results")

  # Would test GSEA visualization
})

test_that("CellDensityPlot basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  mat <- matrix(rnorm(500), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")

    # Add UMAP coordinates
    umap_coords <- matrix(rnorm(20), ncol = 2)
    rownames(umap_coords) <- colnames(srt)
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
    srt[["umap"]] <- Seurat::CreateDimReducObject(
      embeddings = umap_coords,
      key = "UMAP_",
      assay = "RNA"
    )

    srt$group <- rep(c("A", "B"), each = 5)
  })

  # Test density plot
  p <- CellDensityPlot(
    srt = srt,
    group.by = "group",
    reduction = "umap"
  )

  expect_s3_class(p, "ggplot")
})

test_that("FeatureCorPlot basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  mat <- matrix(rnorm(500), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  })

  # Test feature correlation plot
  p <- FeatureCorPlot(
    srt = srt,
    features = rownames(srt)[1:5]
  )

  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CellDimPlot3D basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("plotly")

  mat <- matrix(rnorm(500), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")

    # Add 3D coordinates
    coords_3d <- matrix(rnorm(30), ncol = 3)
    rownames(coords_3d) <- colnames(srt)
    colnames(coords_3d) <- c("UMAP_1", "UMAP_2", "UMAP_3")
    srt[["umap"]] <- Seurat::CreateDimReducObject(
      embeddings = coords_3d,
      key = "UMAP_",
      assay = "RNA"
    )

    srt$group <- rep(c("A", "B"), each = 5)
  })

  # Test 3D plotting
  p <- CellDimPlot3D(
    srt = srt,
    group.by = "group",
    reduction = "umap"
  )

  expect_true(inherits(p, "plotly"))
})

test_that("FeatureDimPlot3D basic functionality", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("plotly")

  mat <- matrix(rnorm(500), nrow = 50)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:50)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)

    # Add 3D coordinates
    coords_3d <- matrix(rnorm(30), ncol = 3)
    rownames(coords_3d) <- colnames(srt)
    colnames(coords_3d) <- c("UMAP_1", "UMAP_2", "UMAP_3")
    srt[["umap"]] <- Seurat::CreateDimReducObject(
      embeddings = coords_3d,
      key = "UMAP_",
      assay = "RNA"
    )
  })

  # Test 3D feature plotting
  p <- FeatureDimPlot3D(
    srt = srt,
    features = rownames(srt)[1],
    reduction = "umap"
  )

  expect_true(inherits(p, "plotly"))
})
