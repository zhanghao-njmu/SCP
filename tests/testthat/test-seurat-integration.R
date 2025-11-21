# Test Seurat integration functions

test_that("SrtAppend merges Seurat objects correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  # Create two simple test Seurat objects
  mat1 <- Matrix::Matrix(rpois(200, 5), nrow = 20, sparse = TRUE)
  colnames(mat1) <- paste0("cell_", 1:10)
  rownames(mat1) <- paste0("gene_", 1:20)

  mat2 <- Matrix::Matrix(rpois(200, 5), nrow = 20, sparse = TRUE)
  colnames(mat2) <- paste0("cell_", 11:20)
  rownames(mat2) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt1 <- Seurat::CreateSeuratObject(counts = mat1, project = "test1")
    srt2 <- Seurat::CreateSeuratObject(counts = mat2, project = "test2")
  })

  # Append objects
  result <- SrtAppend(srt1, srt2)

  # Check result
  expect_s4_class(result, "Seurat")
  expect_equal(ncol(result), ncol(srt1) + ncol(srt2))
  expect_equal(nrow(result), nrow(srt1))
})

test_that("check_srtMerge validates mergeable objects", {
  skip_if_not_installed("Seurat")

  # Create compatible test objects
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

  # Should not throw error for compatible objects
  result <- check_srtMerge(srt_list)
  expect_type(result, "list")
})

test_that("SrtReorder reorders cells correctly", {
  skip_if_not_installed("Seurat")

  # Create test object
  mat <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")
  })

  # Create new order
  new_order <- rev(colnames(srt))

  # Reorder
  result <- SrtReorder(srt, cells = new_order)

  # Check result
  expect_s4_class(result, "Seurat")
  expect_equal(colnames(result), new_order)
  expect_equal(ncol(result), ncol(srt))
})

test_that("DefaultReduction identifies reduction correctly", {
  skip_if_not_installed("Seurat")

  # Create test object
  mat <- matrix(rpois(200, 5), nrow = 20)
  colnames(mat) <- paste0("cell_", 1:10)
  rownames(mat) <- paste0("gene_", 1:20)

  suppressWarnings({
    srt <- Seurat::CreateSeuratObject(counts = mat, project = "test")

    # Add PCA reduction
    pca_data <- matrix(rnorm(30), nrow = 10, ncol = 3)
    rownames(pca_data) <- colnames(srt)
    colnames(pca_data) <- paste0("PC_", 1:3)
    srt[["pca"]] <- Seurat::CreateDimReducObject(
      embeddings = pca_data,
      key = "PC_",
      assay = "RNA"
    )
  })

  # Get default reduction
  reduction <- DefaultReduction(srt)
  expect_type(reduction, "character")
  expect_true(reduction %in% names(srt@reductions))
})
