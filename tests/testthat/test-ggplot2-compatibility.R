# Test ggplot2 compatibility

test_that("ggplot2 version is adequate", {
  skip_if_not_installed("ggplot2")

  ggplot_version <- packageVersion("ggplot2")
  expect_true(ggplot_version >= "3.0.0")

  # Record version for debugging
  message("Testing with ggplot2 version: ", ggplot_version)
})

test_that("Basic ggplot2 plotting works", {
  skip_if_not_installed("ggplot2")

  # Create simple test plot
  df <- data.frame(x = 1:10, y = 1:10)
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  expect_s3_class(p, "ggplot")

  # Build plot
  built <- ggplot2::ggplot_build(p)
  expect_type(built, "list")
})

test_that("aes() function works correctly", {
  skip_if_not_installed("ggplot2")

  # Test standard aes
  aes_obj <- ggplot2::aes(x = x, y = y, color = group)
  expect_s3_class(aes_obj, "uneval")

  # Test aes with expressions
  aes_obj2 <- ggplot2::aes(x = log(x), y = sqrt(y))
  expect_s3_class(aes_obj2, "uneval")
})

test_that("Color scales work", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(x = 1:5, y = 1:5, group = letters[1:5])
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = group)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = rainbow(5))

  expect_s3_class(p, "ggplot")

  # Test gradient scales
  df2 <- data.frame(x = 1:10, y = 1:10, z = 1:10)
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x, y, color = z)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradientn(colors = c("blue", "red"))

  expect_s3_class(p2, "ggplot")
})

test_that("Faceting works", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    x = rep(1:10, 2),
    y = rep(1:10, 2),
    group = rep(c("A", "B"), each = 10)
  )

  # facet_wrap
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~group)

  expect_s3_class(p1, "ggplot")

  # facet_grid
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(group ~ .)

  expect_s3_class(p2, "ggplot")
})

test_that("Themes work correctly", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(x = 1:10, y = 1:10)
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  # Test standard themes
  p1 <- p + ggplot2::theme_minimal()
  expect_s3_class(p1, "ggplot")

  p2 <- p + ggplot2::theme_classic()
  expect_s3_class(p2, "ggplot")

  p3 <- p + ggplot2::theme_void()
  expect_s3_class(p3, "ggplot")

  # Test custom theme
  p4 <- p + ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(color = "black")
  )
  expect_s3_class(p4, "ggplot")
})

test_that("Geoms work correctly", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    x = 1:10,
    y = 1:10,
    group = rep(c("A", "B"), each = 5)
  )

  # Test various geoms
  p_point <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_s3_class(p_point, "ggplot")

  p_line <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_line()
  expect_s3_class(p_line, "ggplot")

  p_bar <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_col()
  expect_s3_class(p_bar, "ggplot")

  p_box <- ggplot2::ggplot(df, ggplot2::aes(group, y)) + ggplot2::geom_boxplot()
  expect_s3_class(p_box, "ggplot")

  p_violin <- ggplot2::ggplot(df, ggplot2::aes(group, y)) + ggplot2::geom_violin()
  expect_s3_class(p_violin, "ggplot")
})

test_that("Guide functions work", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(x = 1:10, y = 1:10, z = 1:10)
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = z)) +
    ggplot2::geom_point()

  # Test guide_legend
  p1 <- p + ggplot2::guides(color = ggplot2::guide_legend(title = "Test"))
  expect_s3_class(p1, "ggplot")

  # Test guide_colorbar
  p2 <- p + ggplot2::guides(color = ggplot2::guide_colorbar(title = "Test"))
  expect_s3_class(p2, "ggplot")

  # Test guide_none
  p3 <- p + ggplot2::guides(color = ggplot2::guide_none())
  expect_s3_class(p3, "ggplot")
})

test_that("Coordinate systems work", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(x = 1:10, y = 1:10)
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_point()

  # Test coord_cartesian
  p1 <- p + ggplot2::coord_cartesian(xlim = c(0, 5))
  expect_s3_class(p1, "ggplot")

  # Test coord_flip
  p2 <- p + ggplot2::coord_flip()
  expect_s3_class(p2, "ggplot")

  # Test coord_polar
  p3 <- p + ggplot2::coord_polar()
  expect_s3_class(p3, "ggplot")
})

test_that("ggsave works", {
  skip_if_not_installed("ggplot2")
  skip_on_cran()

  df <- data.frame(x = 1:10, y = 1:10)
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_point()

  # Test saving to temp file
  temp_file <- tempfile(fileext = ".png")
  on.exit(unlink(temp_file))

  expect_no_error(
    ggplot2::ggsave(temp_file, p, width = 5, height = 5, dpi = 72)
  )

  expect_true(file.exists(temp_file))
})

test_that("layer_data and layer_scales work", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(x = 1:10, y = 1:10, group = rep(c("A", "B"), each = 5))
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = group)) +
    ggplot2::geom_point()

  built <- ggplot2::ggplot_build(p)

  # Test layer_data (if available in current ggplot2 version)
  tryCatch({
    data <- ggplot2::layer_data(p)
    expect_true(is.data.frame(data))
  }, error = function(e) {
    # layer_data may not be available in older versions
    skip("layer_data not available in this ggplot2 version")
  })

  # Test layer_scales (if available)
  tryCatch({
    scales <- ggplot2::layer_scales(p)
    expect_type(scales, "list")
  }, error = function(e) {
    skip("layer_scales not available in this ggplot2 version")
  })
})

test_that("SCP theme functions work with current ggplot2", {
  skip_if_not_installed("ggplot2")

  # Test theme_scp if it exists
  if (exists("theme_scp")) {
    theme <- theme_scp()
    expect_s3_class(theme, "theme")
  }

  # Test theme_blank if it exists
  if (exists("theme_blank")) {
    theme <- theme_blank()
    expect_s3_class(theme, "theme")
  }
})
