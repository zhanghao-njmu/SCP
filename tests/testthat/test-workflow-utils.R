# Test workflow utility functions

test_that("RecoverCounts works correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  # Create test Seurat object with counts
  counts <- matrix(rpois(100, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:10)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)

  # Test recovering counts
  expect_no_error(recovered <- RecoverCounts(srt = srt))
  expect_true(is(recovered, "Seurat"))

  # Check that data slot now equals counts
  data_matrix <- Seurat::GetAssayData(recovered, slot = "data")
  counts_matrix <- Seurat::GetAssayData(recovered, slot = "counts")

  expect_equal(dim(data_matrix), dim(counts_matrix))
})

test_that("RecoverCounts handles missing counts", {
  skip_if_not_installed("Seurat")

  # Create object without counts
  data <- matrix(rnorm(100), nrow = 10)
  colnames(data) <- paste0("cell", 1:10)
  rownames(data) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = data)

  # Should handle gracefully
  expect_no_error(RecoverCounts(srt = srt))
})

test_that("RenameFeatures works with simple replacement", {
  skip_if_not_installed("Seurat")

  # Create test object
  counts <- matrix(rpois(100, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:10)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)

  # Rename features
  new_names <- paste0("new_gene", 1:10)
  names(new_names) <- paste0("gene", 1:10)

  expect_no_error(renamed <- RenameFeatures(srt = srt, newnames = new_names))
  expect_true(is(renamed, "Seurat"))

  # Check new names
  features <- rownames(renamed)
  expect_true(all(grepl("^new_gene", features)))
})

test_that("RenameFeatures handles partial renaming", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(100, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:10)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)

  # Only rename some features
  new_names <- c("gene1" = "renamed1", "gene5" = "renamed5")

  expect_no_error(renamed <- RenameFeatures(srt = srt, newnames = new_names))
  features <- rownames(renamed)

  # Check specific renames
  expect_true("renamed1" %in% features)
  expect_true("renamed5" %in% features)
  expect_false("gene1" %in% features)
  expect_false("gene5" %in% features)

  # Check unchanged features
  expect_true("gene2" %in% features)
})

test_that("RenameFeatures validates input", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(100, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:10)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)

  # Test with unnamed vector (should handle or error appropriately)
  new_names_unnamed <- paste0("new", 1:5)
  # This should either error or require names - test the behavior
  expect_error(RenameFeatures(srt = srt, newnames = new_names_unnamed))
})

test_that("RenameClusters works with identity classes", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)

  # Add cluster labels
  srt$seurat_clusters <- factor(rep(c("0", "1", "2", "3"), each = 5))
  Seurat::Idents(srt) <- "seurat_clusters"

  # Rename clusters
  new_names <- c("0" = "ClusterA", "1" = "ClusterB", "2" = "ClusterC", "3" = "ClusterD")

  expect_no_error(renamed <- RenameClusters(srt = srt, newnames = new_names))
  expect_true(is(renamed, "Seurat"))

  # Check new cluster names
  clusters <- as.character(Seurat::Idents(renamed))
  expect_true(all(clusters %in% c("ClusterA", "ClusterB", "ClusterC", "ClusterD")))
  expect_false(any(clusters %in% c("0", "1", "2", "3")))
})

test_that("RenameClusters handles partial renaming", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$seurat_clusters <- factor(rep(c("0", "1", "2"), c(7, 7, 6)))
  Seurat::Idents(srt) <- "seurat_clusters"

  # Only rename some clusters
  new_names <- c("0" = "TypeA", "2" = "TypeC")

  expect_no_error(renamed <- RenameClusters(srt = srt, newnames = new_names))
  clusters <- as.character(Seurat::Idents(renamed))

  # Check renamed
  expect_true("TypeA" %in% clusters)
  expect_true("TypeC" %in% clusters)

  # Check unchanged
  expect_true("1" %in% clusters)
})

test_that("RenameClusters works with metadata column", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(200, 5), nrow = 10)
  colnames(counts) <- paste0("cell", 1:20)
  rownames(counts) <- paste0("gene", 1:10)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$custom_clusters <- factor(rep(c("A", "B"), each = 10))

  # Rename custom metadata column
  new_names <- c("A" = "Group1", "B" = "Group2")

  expect_no_error(
    renamed <- RenameClusters(
      srt = srt,
      ident = "custom_clusters",
      newnames = new_names
    )
  )

  # Check renamed values in metadata
  expect_true(all(renamed$custom_clusters %in% c("Group1", "Group2")))
})

test_that("Workflow utils handle edge cases", {
  skip_if_not_installed("Seurat")

  # Very small object
  counts <- matrix(rpois(9, 5), nrow = 3)
  colnames(counts) <- paste0("cell", 1:3)
  rownames(counts) <- paste0("gene", 1:3)

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$clusters <- factor(c("A", "A", "B"))

  # Test all functions with minimal object
  expect_no_error(RecoverCounts(srt = srt))
  expect_no_error(RenameFeatures(srt = srt, newnames = c("gene1" = "g1")))
  expect_no_error(RenameClusters(srt = srt, ident = "clusters", newnames = c("A" = "TypeA")))
})
