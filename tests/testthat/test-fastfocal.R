test_that("fastfocal returns identical results to terra::focal for small kernels", {
  r <- terra::rast(matrix(runif(100), 10, 10))
  w <- terra::focalMat(r, d = 1, type = "rectangle")
  
  result_terra <- terra::focal(r, w = w, fun = mean)
  result_fast  <- fastfocal(r, d = 1, w = "rectangle", fun = "mean", engine = "cpp")
  
  expect_s4_class(result_fast, "SpatRaster")
  expect_equal(terra::values(result_fast), terra::values(result_terra), tolerance = 1e-8)
})

test_that("fastfocal with FFT backend matches terra::focal for large kernels", {
  r <- terra::rast(matrix(runif(100), 10, 10))
  w <- terra::focalMat(r, d = 5, type = "rectangle")
  
  result_terra <- terra::focal(r, w = w, fun = mean)
  result_fast  <- fastfocal(r, d = 5, w = "rectangle", fun = "mean", engine = "fft")
  
  expect_s4_class(result_fast, "SpatRaster")
  expect_equal(terra::values(result_fast), terra::values(result_terra), tolerance = 1e-8)
})
