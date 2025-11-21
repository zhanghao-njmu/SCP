# Test workflow functions

test_that("Standard_SCP basic workflow", {
  skip_if_not_installed("Seurat")
  skip("Skipping Standard_SCP - computationally intensive")

  # This is a complex integration test
  # Would require significant computation time
})

test_that("check_DataType identifies data types", {
  skip_if_not_installed("Seurat")

  # Create test matrices
  # Raw counts (integers)
  counts <- matrix(rpois(200, 5), nrow = 20)
  result_counts <- check_DataType(data = counts)
  expect_type(result_counts, "list")
  expect_true("datatype" %in% names(result_counts))

  # Normalized data (floats)
  norm_data <- log1p(counts)
  result_norm <- check_DataType(data = norm_data)
  expect_type(result_norm, "list")

  # Scaled data
  scaled_data <- scale(norm_data)
  result_scaled <- check_DataType(data = scaled_data)
  expect_type(result_scaled, "list")
})

test_that("SrtAppend combines objects", {
  skip_if_not_installed("Seurat")

  mat1 <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat1) <- paste0("cell_", 1:10)
  rownames(mat1) <- paste0("gene_", 1:20)

  mat2 <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat2) <- paste0("cell_", 11:20)
  rownames(mat2) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt1 <- Seurat::CreateSeuratObject(counts = mat1, project = "test1")
    srt2 <- Seurat::CreateSeuratObject(counts = mat2, project = "test2")
  })

  result <- SrtAppend(srt1, srt2)

  expect_s4_class(result, "Seurat")
  expect_equal(ncol(result), ncol(srt1) + ncol(srt2))
  expect_equal(nrow(result), nrow(srt1))
})

test_that("SrtReorder reorders correctly", {
  skip_if_not_installed("Seurat")

  mat <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
  })

  # Reorder cells
  new_order <- rev(colnames(srt))
  result <- SrtReorder(srt, cells = new_order)

  expect_s4_class(result, "Seurat")
  expect_equal(colnames(result), new_order)
})

test_that("check_srtList validates list", {
  skip_if_not_installed("Seurat")

  mat1 <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat1) <- paste0("cell_", 1:10)
  rownames(mat1) <- paste0("gene_", 1:20)

  mat2 <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat2) <- paste0("cell_", 11:20)
  rownames(mat2) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt1 <- Seurat::CreateSeuratObject(counts = mat1, project = "test1")
    srt2 <- Seurat::CreateSeuratObject(counts = mat2, project = "test2")
  })

  srt_list <- list(srt1 = srt1, srt2 = srt2)
  result <- check_srtList(srt_list)

  expect_type(result, "list")
  expect_length(result, 2)
})

test_that("check_srtMerge validates merge", {
  skip_if_not_installed("Seurat")

  mat1 <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat1) <- paste0("cell_", 1:10)
  rownames(mat1) <- paste0("gene_", 1:20)

  mat2 <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat2) <- paste0("cell_", 11:20)
  rownames(mat2) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt1 <- Seurat::CreateSeuratObject(counts = mat1, project = "test1")
    srt2 <- Seurat::CreateSeuratObject(counts = mat2, project = "test2")
  })

  srt_list <- list(srt1 = srt1, srt2 = srt2)
  result <- check_srtMerge(srt_list)

  expect_type(result, "list")
})

test_that("DefaultReduction gets reduction", {
  skip_if_not_installed("Seurat")

  mat <- matrix(rnorm(200), nrow = 20)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
    srt <- Seurat::FindVariableFeatures(srt, verbose = FALSE)
    srt <- Seurat::ScaleData(srt, verbose = FALSE)
    srt <- Seurat::RunPCA(srt, npcs = 3, verbose = FALSE)
  })

  reduction <- DefaultReduction(srt)

  expect_type(reduction, "character")
  expect_true(reduction %in% names(srt@reductions))
})

test_that("RenameClusters renames correctly", {
  skip_if_not_installed("Seurat")

  mat <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
    srt$cluster <- rep(c("0", "1"), each = 5)
  })

  new_names <- c("0" = "TypeA", "1" = "TypeB")
  result <- RenameClusters(srt, new_names = new_names, cluster_col = "cluster")

  expect_s4_class(result, "Seurat")
  expect_true(all(result$cluster %in% c("TypeA", "TypeB")))
})

test_that("RenameFeatures renames genes", {
  skip_if_not_installed("Seurat")

  mat <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
  })

  new_names <- setNames(paste0("NEW_", rownames(srt)), rownames(srt))
  result <- RenameFeatures(srt, new_names = new_names)

  expect_s4_class(result, "Seurat")
  expect_true(all(grepl("^NEW_", rownames(result))))
})
