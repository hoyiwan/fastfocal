# tests/testthat.R

library(testthat)
library(fastfocal)
library(terra)

# tests/testthat/test-basic.R
test_that("basic functionality works", {
  x <- rast(matrix(1:100, 10, 10))
  result <- fastfocal(x, d = 300, w = "circle", fun = "mean")
  expect_s4_class(result, "SpatRaster")
})
