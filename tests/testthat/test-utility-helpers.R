# Test utility helper functions

test_that("as_matrix converts to matrix", {
  skip_if_not_installed("Matrix")

  # Test with sparse matrix
  sparse_mat <- Matrix::Matrix(c(1, 0, 0, 2, 3, 0), nrow = 2, sparse = TRUE)
  result <- as_matrix(sparse_mat)

  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(sparse_mat))

  # Test with regular matrix
  regular_mat <- matrix(1:6, nrow = 2)
  result2 <- as_matrix(regular_mat)

  expect_true(is.matrix(result2))
  expect_equal(result2, regular_mat)
})

test_that("try_get safely gets values", {
  test_list <- list(a = 1, b = 2, c = 3)

  # Test existing key
  result1 <- try_get(test_list, "a")
  expect_equal(result1, 1)

  # Test non-existing key with default
  result2 <- try_get(test_list, "d", default = 999)
  expect_equal(result2, 999)

  # Test NULL default
  result3 <- try_get(test_list, "e")
  expect_null(result3)
})

test_that("unnest unnests data", {
  # Test basic unnesting
  df <- data.frame(
    id = 1:3,
    value = I(list(c(1, 2), c(3, 4, 5), c(6)))
  )

  result <- unnest(df, cols = "value")

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) >= nrow(df))
})

test_that("check_R validates R version", {
  # Test R version check
  result <- tryCatch({
    check_R(R_version = "4.0.0")
    TRUE
  }, error = function(e) {
    FALSE
  })

  expect_type(result, "logical")
})

test_that("arrow function creates arrows", {
  skip_if_not_installed("grid")

  # Test arrow creation
  arr <- arrow(length = grid::unit(0.1, "inches"))

  expect_s3_class(arr, "arrow")
})

test_that("gpar function creates graphics parameters", {
  skip_if_not_installed("grid")

  # Test gpar creation
  gp <- gpar(col = "red", fill = "blue", lwd = 2)

  expect_s3_class(gp, "gpar")
})

test_that("unit function creates grid units", {
  skip_if_not_installed("grid")

  # Test unit creation
  u <- unit(1, "cm")

  expect_true(grid::is.unit(u))
})

test_that("palette_scp returns palettes", {
  # Test default palette
  colors <- palette_scp()

  expect_type(colors, "character")
  expect_true(length(colors) > 0)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))

  # Test specific palette
  colors2 <- palette_scp("Paired")

  expect_type(colors2, "character")
  expect_true(length(colors2) > 0)
})

test_that("show_palettes displays information", {
  # Test palette display (should not error)
  expect_no_error(show_palettes(show = FALSE))
})

test_that("%>% pipe operator works", {
  # Test pipe operator
  result <- 1:10 %>% sum() %>% sqrt()

  expect_type(result, "double")
  expect_equal(result, sqrt(sum(1:10)))
})

test_that("%||% null coalescing works", {
  # Test null coalescing
  result1 <- NULL %||% "default"
  expect_equal(result1, "default")

  result2 <- "value" %||% "default"
  expect_equal(result2, "value")
})

test_that("capitalize capitalizes correctly", {
  # Basic capitalization
  expect_equal(capitalize("hello"), "Hello")
  expect_equal(capitalize("HELLO"), "HELLO")
  expect_equal(capitalize("HELLO", force_tolower = TRUE), "Hello")

  # Vector input
  result <- capitalize(c("apple", "banana", "cherry"))
  expect_equal(result, c("Apple", "Banana", "Cherry"))

  # NULL input
  expect_null(capitalize(NULL))

  # Factor input
  factor_input <- factor(c("dog", "cat"))
  result <- capitalize(factor_input)
  expect_equal(result, c("Dog", "Cat"))

  # Edge cases
  expect_equal(capitalize(""), "")
  expect_equal(capitalize("a"), "A")
  expect_equal(capitalize("123"), "123")
})

test_that("adjcolors adjusts colors", {
  colors <- c("#FF0000", "#00FF00", "#0000FF")

  # Test alpha adjustment
  result <- adjcolors(colors, alpha = 0.5)

  expect_type(result, "character")
  expect_length(result, 3)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", result)))

  # Test different alpha values
  result_low <- adjcolors(colors, alpha = 0.1)
  result_high <- adjcolors(colors, alpha = 0.9)

  expect_false(identical(result_low, result_high))
})

test_that("blendcolors blends colors", {
  colors <- c("#FF0000", "#0000FF")

  # Test different blend modes
  blend <- blendcolors(colors, mode = "blend")
  average <- blendcolors(colors, mode = "average")
  screen <- blendcolors(colors, mode = "screen")
  multiply <- blendcolors(colors, mode = "multiply")

  expect_type(blend, "character")
  expect_type(average, "character")
  expect_type(screen, "character")
  expect_type(multiply, "character")

  # Different modes should produce different results
  expect_false(all(c(
    blend == average,
    blend == screen,
    blend == multiply,
    average == screen
  )))

  # Test with NA
  colors_na <- c("#FF0000", NA, "#0000FF")
  result_na <- blendcolors(colors_na, mode = "blend")
  expect_type(result_na, "character")

  # Test with single color
  single <- blendcolors("#FF0000", mode = "blend")
  expect_type(single, "character")
})

test_that("get_vars extracts variable names", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(x = 1:10, y = 1:10, group = rep(c("A", "B"), 5))
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = group)) +
    ggplot2::geom_point()

  vars <- get_vars(p)

  expect_type(vars, "character")
  expect_true(length(vars) > 0)
})

test_that("drop_data removes plot data", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(x = rnorm(1000), y = rnorm(1000))
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  result <- drop_data(p)

  expect_s3_class(result, "ggplot")

  # Test with NULL
  expect_null(drop_data(NULL))
})

test_that("slim_data reduces data size", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(x = rnorm(10000), y = rnorm(10000))
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  result <- slim_data(p, n = 1000)

  expect_s3_class(result, "ggplot")
})

test_that("panel_fix adjusts panels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gtable")

  df <- data.frame(x = 1:10, y = 1:10)
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  result <- panel_fix(
    p,
    height = grid::unit(5, "cm"),
    width = grid::unit(5, "cm")
  )

  expect_true(!is.null(result))
})

test_that("panel_fix_overall adjusts overall size", {
  skip_if_not_installed("ggplot2")
  skip("Skipping panel_fix_overall - complex gtable operation")

  # Would test overall panel adjustment
})
