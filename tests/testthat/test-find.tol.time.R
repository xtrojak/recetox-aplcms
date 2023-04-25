test_that("compute_rt_tol_relative computes something", {
  aver.bin.size <- 200
  min.bins <- 50
  max.bins <- 100
  max.num.segments <- 10000
  number_of_samples <- 1
  rt <- c(100, 101, 101.5, 101.2, 105, 108, 108.3)
  breaks <- c(1, 5, 6)

  actual <- compute_rt_tol_relative(
    breaks,
    max.num.segments,
    aver.bin.size,
    number_of_samples,
    rt,
    min.bins,
    max.bins
  )

  expect_equal(actual, 0.7)
})