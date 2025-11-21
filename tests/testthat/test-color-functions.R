# Test color manipulation functions

test_that("palette_scp returns valid color palettes", {
  skip_if_not_installed("SCP")

  # Test with default palette
  colors <- palette_scp()
  expect_type(colors, "character")
  expect_true(length(colors) > 0)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))

  # Test with specific palette name
  colors_set2 <- palette_scp("Set2")
  expect_type(colors_set2, "character")
  expect_true(length(colors_set2) > 0)
})

test_that("show_palettes displays palette information", {
  skip_if_not_installed("SCP")

  # Test that function runs without error
  expect_no_error(show_palettes(show = FALSE))
})

test_that("color blending produces consistent results", {
  # Test that same colors always produce same blend
  colors1 <- c("#FF0000", "#0000FF")
  result1 <- blendcolors(colors1, mode = "blend")
  result2 <- blendcolors(colors1, mode = "blend")

  expect_equal(result1, result2)

  # Test different blending modes produce different results
  blend <- blendcolors(colors1, mode = "blend")
  average <- blendcolors(colors1, mode = "average")
  screen <- blendcolors(colors1, mode = "screen")
  multiply <- blendcolors(colors1, mode = "multiply")

  # They should not all be the same
  expect_false(all(c(blend == average, blend == screen, blend == multiply)))
})

test_that("adjcolors handles various alpha values", {
  colors <- c("#FF0000", "#00FF00", "#0000FF")

  # Test with different alpha values
  result_low <- adjcolors(colors, alpha = 0.2)
  result_mid <- adjcolors(colors, alpha = 0.5)
  result_high <- adjcolors(colors, alpha = 0.8)

  # All should return valid hex colors
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", result_low)))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", result_mid)))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", result_high)))

  # Results should be different for different alpha
  expect_false(identical(result_low, result_high))
})

test_that("color functions handle edge cases", {
  # Test with single color
  single <- blendcolors("#FF0000", mode = "blend")
  expect_type(single, "character")
  expect_length(single, 1)

  # Test with empty vector after NA removal
  empty_result <- blendcolors(c(NA, NA), mode = "blend")
  expect_true(is.na(empty_result) || is.null(empty_result))

  # Test adjcolors with alpha = 1
  colors <- c("#FF0000", "#00FF00")
  result_alpha1 <- adjcolors(colors, alpha = 1.0)
  expect_type(result_alpha1, "character")

  # Test adjcolors with alpha = 0
  result_alpha0 <- adjcolors(colors, alpha = 0.0)
  expect_type(result_alpha0, "character")
})
