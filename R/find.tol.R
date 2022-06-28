find.tol <- function(mz_values,
                     mz_max_diff = 1e-4,
                     aver.bin.size = 4000,
                     min.bins = 50,
                     max.bins = 200,
                     do.plot = TRUE) {
    mz_values <- mz_values[order(mz_values)]
    l <- length(mz_values)
    # pairwise difference divided by average
    da <- (mz_values[2:l] - mz_values[1:(l - 1)]) / ((mz_values[2:l] + mz_values[1:(l - 1)]) / 2)
    # filter values outside of tolerance limit
    da <- da[da < mz_max_diff]
    # number of equally spaced points at which the density is to be estimated
    n <- min(max.bins, max(round(length(da) / aver.bin.size), min.bins))
    # estimate probability density function of da
    des <- density(da, kernel = "gaussian", n = n, bw = mz_max_diff / n * 2, from = 0)
    # the n (-1?) coordinates of the points where the density is estimated
    y <- des$y[des$x > 0]
    # the estimated density values
    x <- des$x[des$x > 0]
    
    # select the upper 75% of the sorted data
    to.use <- da[da > max(da) / 4] - max(da) / 4
    # parameter of the exponential distribution is estimated
    this.rate <- MASS::fitdistr(to.use, "exponential")$estimate
    # values of the exponential distribution are calculated at equally spaced points
    exp.y <- dexp(x, rate = this.rate)
    # add the rest of the data?
    exp.y <- exp.y * sum(y[x > max(da) / 4]) / sum(exp.y[x > max(da) / 4])
    
    # cutoff is selected where the density of the empirical distribution is >1.5 times the density of the exponential distribution
    yy <- cumsum(y > 1.5 * exp.y)
    yi <- seq_along(yy)
    # find last index where y is greater than exp.y
    sel <- min(which(yy < yi)) - 1
    
    if (do.plot) {
        tolerance_plot(x, y, exp.y, sel, main = "find m/z tolerance")
    }
    
    # corresponding coordinate is used as m/z tolerance
    return(x[sel])
}
