# Test Sankey and Alluvial plot functions

test_that("make_long creates long format", {
  skip_if_not_installed("ggplot2")

  # Create test data
  df <- data.frame(
    id = 1:10,
    node1 = rep(c("A", "B"), each = 5),
    node2 = rep(c("X", "Y"), 5),
    node3 = sample(c("P", "Q"), 10, replace = TRUE)
  )

  # Test make_long
  result <- make_long(df, node1, node2, node3)

  expect_s3_class(result, "data.frame")
  expect_true("node" %in% colnames(result))
  expect_true("x" %in% colnames(result))
})

test_that("geom_sankey creates layer", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    x = rep(1:2, each = 5),
    node = rep(c("A", "B"), each = 5),
    next_x = rep(2:3, each = 5),
    next_node = rep(c("X", "Y"), 5),
    value = runif(10)
  )

  # Test that geom_sankey can be added to ggplot
  p <- ggplot2::ggplot(df) +
    geom_sankey(
      ggplot2::aes(
        x = x, next_x = next_x,
        node = node, next_node = next_node,
        fill = node,
        value = value
      )
    )

  expect_s3_class(p, "ggplot")
})

test_that("geom_alluvial creates layer", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    x = rep(1:2, each = 5),
    stratum = rep(c("A", "B"), each = 5),
    alluvium = rep(1:5, 2),
    value = runif(10)
  )

  # Test that geom_alluvial can be added to ggplot
  p <- ggplot2::ggplot(df) +
    geom_alluvial(
      ggplot2::aes(
        x = x,
        stratum = stratum,
        alluvium = alluvium,
        fill = stratum
      )
    )

  expect_s3_class(p, "ggplot")
})

test_that("geom_sankey_label creates labels", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    x = 1:3,
    node = c("A", "B", "C"),
    value = c(10, 20, 15)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = value, label = node)) +
    geom_sankey_label()

  expect_s3_class(p, "ggplot")
})

test_that("geom_sankey_text creates text", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    x = 1:3,
    node = c("A", "B", "C"),
    value = c(10, 20, 15)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = value, label = node)) +
    geom_sankey_text()

  expect_s3_class(p, "ggplot")
})

test_that("geom_alluvial_label creates labels", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    x = rep(1:2, each = 3),
    stratum = rep(c("A", "B", "C"), 2),
    alluvium = rep(1:3, 2)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, stratum = stratum, alluvium = alluvium)) +
    geom_alluvial_label()

  expect_s3_class(p, "ggplot")
})

test_that("geom_alluvial_text creates text", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    x = rep(1:2, each = 3),
    stratum = rep(c("A", "B", "C"), 2),
    alluvium = rep(1:3, 2)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, stratum = stratum, alluvium = alluvium)) +
    geom_alluvial_text()

  expect_s3_class(p, "ggplot")
})

test_that("geom_sankey_bump creates bump chart", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    x = rep(1:3, each = 2),
    y = rep(1:2, 3),
    group = rep(c("A", "B"), 3)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, group = group)) +
    geom_sankey_bump()

  expect_s3_class(p, "ggplot")
})

test_that("theme_sankey creates theme", {
  skip_if_not_installed("ggplot2")

  theme <- theme_sankey()

  expect_s3_class(theme, "theme")
  expect_s3_class(theme, "gg")
})

test_that("theme_alluvial creates theme", {
  skip_if_not_installed("ggplot2")

  theme <- theme_alluvial()

  expect_s3_class(theme, "theme")
  expect_s3_class(theme, "gg")
})

test_that("theme_sankey_bump creates theme", {
  skip_if_not_installed("ggplot2")

  theme <- theme_sankey_bump()

  expect_s3_class(theme, "theme")
  expect_s3_class(theme, "gg")
})

test_that("sankey themes work with ggplot", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(x = 1:10, y = rnorm(10))

  p1 <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    theme_sankey()

  expect_s3_class(p1, "ggplot")

  p2 <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    theme_alluvial()

  expect_s3_class(p2, "ggplot")

  p3 <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    theme_sankey_bump()

  expect_s3_class(p3, "ggplot")
})
