test_that("compute_rt_tol_relative computes something", {
  aver.bin.size <- 3
  min.bins <- 2
  max.bins <- 10
  max.num.segments <- 10
  number_of_samples <- 2
  rt <- c(0, 1, 2, 3, 3.5, 3.7, 3.9, 4, 100, 101, 101.5, 101.2, 105, 108, 108.1, 108.12, 108.3)
  breaks <- c(1, 2, 3, 8, 13, 17)

  actual <- compute_rt_tol_relative(
    breaks,
    max.num.segments,
    aver.bin.size,
    number_of_samples,
    rt,
    min.bins,
    max.bins
  )

  expect_equal(actual, 1.04167, tolerance=0.001)
})
