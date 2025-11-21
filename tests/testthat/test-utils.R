# Test utility functions

test_that("capitalize works correctly", {
  # Test basic capitalization
  expect_equal(capitalize("hello"), "Hello")
  expect_equal(capitalize("HELLO", force_tolower = TRUE), "Hello")

  # Test vector input
  result <- capitalize(c("hello", "world"))
  expect_equal(result, c("Hello", "World"))

  # Test NULL input
  expect_null(capitalize(NULL))

  # Test factor input
  factor_input <- factor(c("apple", "banana"))
  result <- capitalize(factor_input)
  expect_equal(result, c("Apple", "Banana"))

  # Test error for non-character input
  expect_error(capitalize(123), "x must be the type of character")
})

test_that("adjcolors adjusts colors with alpha", {
  # Test basic color adjustment
  colors <- c("#FF0000", "#00FF00", "#0000FF")
  alpha <- 0.5

  result <- adjcolors(colors, alpha)

  # Result should be a character vector of the same length
  expect_type(result, "character")
  expect_length(result, 3)

  # Result should contain valid hex colors
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", result)))
})

test_that("blendcolors blends colors correctly", {
  # Test with valid colors
  colors <- c("#FF0000", "#0000FF")

  # Test blend mode
  result_blend <- blendcolors(colors, mode = "blend")
  expect_type(result_blend, "character")
  expect_length(result_blend, 1)

  # Test average mode
  result_avg <- blendcolors(colors, mode = "average")
  expect_type(result_avg, "character")
  expect_length(result_avg, 1)

  # Test screen mode
  result_screen <- blendcolors(colors, mode = "screen")
  expect_type(result_screen, "character")
  expect_length(result_screen, 1)

  # Test multiply mode
  result_multiply <- blendcolors(colors, mode = "multiply")
  expect_type(result_multiply, "character")
  expect_length(result_multiply, 1)

  # Test with NA values (should be filtered out)
  colors_with_na <- c("#FF0000", NA, "#0000FF")
  result_na <- blendcolors(colors_with_na, mode = "blend")
  expect_type(result_na, "character")
})
