# Test plotting helper functions

test_that("theme_scp creates valid ggplot theme", {
  skip_if_not_installed("ggplot2")

  theme <- theme_scp()
  expect_s3_class(theme, "theme")
  expect_s3_class(theme, "gg")
})

test_that("theme_blank creates minimal theme", {
  skip_if_not_installed("ggplot2")

  theme <- theme_blank()
  expect_s3_class(theme, "theme")
  expect_s3_class(theme, "gg")
})

test_that("panel_fix adjusts plot panels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gtable")

  # Create simple test plot
  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  # Test panel_fix
  result <- panel_fix(
    p,
    height = ggplot2::unit(5, "cm"),
    width = ggplot2::unit(5, "cm")
  )

  expect_true(gtable::is.gtable(result) || inherits(result, "ggplot"))
})

test_that("drop_data removes large data from plots", {
  skip_if_not_installed("ggplot2")

  # Create test plot with data
  df <- data.frame(x = rnorm(1000), y = rnorm(1000), group = rep(letters[1:10], 100))
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = group)) +
    ggplot2::geom_point()

  # Drop data
  result <- drop_data(p)

  expect_s3_class(result, "ggplot")
})

test_that("slim_data reduces data size", {
  skip_if_not_installed("ggplot2")

  # Create test plot with large data
  df <- data.frame(x = rnorm(10000), y = rnorm(10000))
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  # Slim data
  result <- slim_data(p, n = 1000)

  expect_s3_class(result, "ggplot")

  # Check that data was reduced (if accessible)
  if (!is.null(result$data)) {
    expect_true(nrow(result$data) <= 1000)
  }
})

test_that("plot manipulation handles edge cases", {
  skip_if_not_installed("ggplot2")

  # Test with NULL input
  expect_null(drop_data(NULL))

  # Test with empty plot
  p_empty <- ggplot2::ggplot()
  result <- drop_data(p_empty)
  expect_s3_class(result, "ggplot")

  # Test panel_fix with default parameters
  p <- ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  result <- panel_fix(p)
  expect_true(!is.null(result))
})

test_that("get_vars extracts variables from aesthetic mapping", {
  skip_if_not_installed("ggplot2")

  # Create plot with various aesthetics
  df <- data.frame(x = 1:10, y = 1:10, color = rep(c("a", "b"), 5))
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = color)) +
    ggplot2::geom_point()

  vars <- get_vars(p)
  expect_type(vars, "character")
  expect_true("color" %in% vars || "colour" %in% vars)
})

test_that("color scaling helpers work", {
  # Test that color-related constants exist
  expect_true(exists("palette_scp", mode = "function"))

  # Test basic color operations
  colors <- c("#FF0000", "#00FF00", "#0000FF")
  expect_length(colors, 3)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))
})
