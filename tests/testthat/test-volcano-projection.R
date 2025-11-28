# Tests for VolcanoPlot and ProjectionPlot

test_that("VolcanoPlot creates basic volcano plot", {
  skip_if_not_installed("ggplot2")

  # Create test DE results
  de_results <- data.frame(
    gene = paste0("gene", 1:100),
    avg_log2FC = rnorm(100, mean = 0, sd = 1.5),
    p_val_adj = runif(100, min = 0, max = 1),
    stringsAsFactors = FALSE
  )

  # Calculate -log10 p-value
  de_results$neg_log10_p <- -log10(de_results$p_val_adj + 1e-300)

  # Test basic plot
  expect_no_error(
    p <- VolcanoPlot(
      data = de_results,
      x = "avg_log2FC",
      y = "neg_log10_p"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("VolcanoPlot handles significance thresholds", {
  skip_if_not_installed("ggplot2")

  de_results <- data.frame(
    gene = paste0("gene", 1:100),
    log2FC = rnorm(100, mean = 0, sd = 2),
    pvalue = runif(100, min = 0, max = 0.1),
    padj = p.adjust(runif(100, min = 0, max = 0.1), method = "BH"),
    stringsAsFactors = FALSE
  )

  # Test with custom thresholds
  expect_no_error(
    p <- VolcanoPlot(
      data = de_results,
      x = "log2FC",
      y = "pvalue",
      fc_threshold = 1.0,
      p_threshold = 0.05
    )
  )
  expect_s3_class(p, "gg")
})

test_that("VolcanoPlot highlights specific genes", {
  skip_if_not_installed("ggplot2")

  de_results <- data.frame(
    gene = paste0("gene", 1:50),
    logFC = rnorm(50, mean = 0, sd = 1.5),
    pval = runif(50, min = 0, max = 0.1),
    stringsAsFactors = FALSE
  )

  # Highlight top genes
  highlight_genes <- c("gene1", "gene5", "gene10")

  expect_no_error(
    p <- VolcanoPlot(
      data = de_results,
      x = "logFC",
      y = "pval",
      highlight = highlight_genes
    )
  )
  expect_s3_class(p, "gg")
})

test_that("VolcanoPlot handles different color schemes", {
  skip_if_not_installed("ggplot2")

  de_results <- data.frame(
    gene = paste0("gene", 1:100),
    lfc = rnorm(100, mean = 0, sd = 1),
    adj_p = runif(100, min = 0, max = 1),
    stringsAsFactors = FALSE
  )

  # Test with custom colors
  expect_no_error(
    p <- VolcanoPlot(
      data = de_results,
      x = "lfc",
      y = "adj_p",
      color_up = "red",
      color_down = "blue",
      color_ns = "grey"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("VolcanoPlot validates input data", {
  skip_if_not_installed("ggplot2")

  # Test with missing columns
  bad_data <- data.frame(gene = paste0("gene", 1:10))

  expect_error(
    VolcanoPlot(
      data = bad_data,
      x = "missing_column",
      y = "another_missing"
    )
  )

  # Test with empty data
  expect_error(
    VolcanoPlot(
      data = data.frame(),
      x = "x",
      y = "y"
    )
  )
})

test_that("VolcanoPlot handles edge cases", {
  skip_if_not_installed("ggplot2")

  # Very small dataset
  small_data <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    fc = c(0.5, -0.5, 0),
    p = c(0.01, 0.05, 0.9)
  )

  expect_no_error(
    p <- VolcanoPlot(
      data = small_data,
      x = "fc",
      y = "p"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("ProjectionPlot creates projection visualizations", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  # Create test objects
  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  query <- Seurat::CreateSeuratObject(counts = counts)
  ref <- Seurat::CreateSeuratObject(counts = counts)

  # Add embeddings
  query[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "PC_",
    assay = "RNA"
  )

  ref[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "PC_",
    assay = "RNA"
  )

  # Add cell types
  query$celltype <- factor(rep(c("TypeA", "TypeB"), each = 10))
  ref$celltype <- factor(rep(c("TypeA", "TypeB"), each = 10))

  # Test projection plot
  expect_no_error(
    p <- ProjectionPlot(
      srt_query = query,
      srt_ref = ref,
      query_group = "celltype",
      ref_group = "celltype"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("ProjectionPlot handles different reductions", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  query <- Seurat::CreateSeuratObject(counts = counts)
  ref <- Seurat::CreateSeuratObject(counts = counts)

  # Add UMAP reduction
  query[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "UMAP_",
    assay = "RNA"
  )

  ref[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rnorm(40), ncol = 2),
    key = "UMAP_",
    assay = "RNA"
  )

  query$group <- "query"
  ref$group <- "ref"

  expect_no_error(
    p <- ProjectionPlot(
      srt_query = query,
      srt_ref = ref,
      reduction = "umap"
    )
  )
  expect_s3_class(p, "gg")
})

test_that("ProjectionPlot validates inputs", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  query <- Seurat::CreateSeuratObject(counts = counts)
  ref <- Seurat::CreateSeuratObject(counts = counts)

  # Test with missing reduction
  expect_error(
    ProjectionPlot(
      srt_query = query,
      srt_ref = ref,
      reduction = "missing_reduction"
    )
  )
})

test_that("Volcano and projection plots handle labels", {
  skip_if_not_installed("ggplot2")

  # VolcanoPlot with labels
  de_results <- data.frame(
    gene = paste0("gene", 1:20),
    lfc = rnorm(20),
    pval = runif(20, 0, 0.1)
  )

  expect_no_error(
    p <- VolcanoPlot(
      data = de_results,
      x = "lfc",
      y = "pval",
      label = "gene",
      label_threshold = 0.05
    )
  )
  expect_s3_class(p, "gg")
})
