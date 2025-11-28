# Enhanced tests for heatmap functions

test_that("GroupHeatmap creates basic heatmap", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")

  # Create test object
  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B"), each = 10))

  # Test basic heatmap
  expect_no_error(
    ht <- GroupHeatmap(
      srt = srt,
      features = c("gene1", "gene2", "gene3"),
      group.by = "group"
    )
  )
})

test_that("GroupHeatmap handles multiple groups", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")

  counts <- matrix(rpois(400, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:40)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B", "C", "D"), each = 10))

  expect_no_error(
    ht <- GroupHeatmap(
      srt = srt,
      features = rownames(srt)[1:5],
      group.by = "group"
    )
  )
})

test_that("GroupHeatmap supports scaling options", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B"), each = 10))

  # Test with row scaling
  expect_no_error(
    ht <- GroupHeatmap(
      srt = srt,
      features = c("gene1", "gene2", "gene3"),
      group.by = "group",
      scale = "row"
    )
  )

  # Test with column scaling
  expect_no_error(
    ht <- GroupHeatmap(
      srt = srt,
      features = c("gene1", "gene2", "gene3"),
      group.by = "group",
      scale = "column"
    )
  )
})

test_that("FeatureHeatmap creates feature heatmaps", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$celltype <- factor(rep(c("Type1", "Type2"), each = 10))

  # Test feature heatmap
  expect_no_error(
    ht <- FeatureHeatmap(
      srt = srt,
      features = c("gene1", "gene2", "gene3", "gene4"),
      group.by = "celltype"
    )
  )
})

test_that("FeatureHeatmap handles large feature sets", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")

  counts <- matrix(rpois(500, 5), nrow = 25)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:25)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B"), each = 10))

  # Test with many features
  expect_no_error(
    ht <- FeatureHeatmap(
      srt = srt,
      features = rownames(srt)[1:15],
      group.by = "group"
    )
  )
})

test_that("FeatureHeatmap supports clustering", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B"), each = 10))

  # Test with row clustering
  expect_no_error(
    ht <- FeatureHeatmap(
      srt = srt,
      features = rownames(srt)[1:5],
      group.by = "group",
      cluster_rows = TRUE
    )
  )

  # Test with column clustering
  expect_no_error(
    ht <- FeatureHeatmap(
      srt = srt,
      features = rownames(srt)[1:5],
      group.by = "group",
      cluster_columns = TRUE
    )
  )
})

test_that("CellCorHeatmap creates correlation heatmaps", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)

  # Test correlation heatmap
  expect_no_error(
    ht <- CellCorHeatmap(
      srt = srt,
      features = c("gene1", "gene2", "gene3", "gene4")
    )
  )
})

test_that("CellCorHeatmap handles different correlation methods", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)

  # Test Pearson correlation
  expect_no_error(
    ht1 <- CellCorHeatmap(
      srt = srt,
      features = c("gene1", "gene2", "gene3"),
      cor_method = "pearson"
    )
  )

  # Test Spearman correlation
  expect_no_error(
    ht2 <- CellCorHeatmap(
      srt = srt,
      features = c("gene1", "gene2", "gene3"),
      cor_method = "spearman"
    )
  )
})

test_that("Heatmaps handle subset of cells", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")

  counts <- matrix(rpois(400, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:40)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B", "C", "D"), each = 10))

  # Test with cell subset
  cells_subset <- colnames(srt)[1:20]

  expect_no_error(
    ht <- GroupHeatmap(
      srt = srt,
      features = c("gene1", "gene2", "gene3"),
      group.by = "group",
      cells = cells_subset
    )
  )
})

test_that("Heatmaps validate inputs", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B"), each = 10))

  # Test with missing features
  expect_error(
    GroupHeatmap(
      srt = srt,
      features = c("missing_gene1", "missing_gene2"),
      group.by = "group"
    )
  )

  # Test with missing group.by
  expect_error(
    GroupHeatmap(
      srt = srt,
      features = c("gene1", "gene2"),
      group.by = "missing_column"
    )
  )
})

test_that("Heatmaps handle edge cases", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")

  # Very small dataset
  counts <- matrix(rpois(30, 5), nrow = 5)
  colnames(counts) <- paste0("cell", 1:6)
  rownames(counts) <- paste0("gene", 1:5)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B"), each = 3))

  # Test with minimal data
  expect_no_error(
    ht <- GroupHeatmap(
      srt = srt,
      features = c("gene1", "gene2"),
      group.by = "group"
    )
  )
})

test_that("Heatmaps support custom color schemes", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B"), each = 10))

  # Test with custom colors
  expect_no_error(
    ht <- GroupHeatmap(
      srt = srt,
      features = c("gene1", "gene2", "gene3"),
      group.by = "group",
      palette = "RdBu"
    )
  )
})
