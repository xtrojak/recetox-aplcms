find.tol.time <- function(mz,
                          chr,
                          lab,
                          number_of_samples,
                          mz_tol_relative = 2e-5,
                          rt_tol_relative = NA,
                          aver.bin.size = 200,
                          min.bins = 50,
                          max.bins = 100,
                          mz_tol_absolute = 0.01,
                          max.num.segments = 10000,
                          do.plot = TRUE) {
    o <- order(mz)
    mz <- mz[o]
    chr <- chr[o]  # rt
    lab <- lab[o]  # experiment/batch ID
    rm(o)
    
    l <- length(mz)
    
    # indices of m/z where difference from its neighbor is greater than allowed tolerance
    # in result separates m/z to groups with similar values
    breaks <-
        c(0, which((mz[2:l] - mz[1:(l - 1)]) 
                    > min(mz_tol_absolute, 
                          mz_tol_relative * ((mz[2:l] + mz[1:(l - 1)]) / 2)
                         )
                  ), 
          l)
    
    for (i in 2:length(breaks)) {
        # sort rt inside each m/z group
        this.o <- order(chr[(breaks[i - 1] + 1):breaks[i]])
        this.o <- this.o + breaks[i - 1]
        # reorder other m/z, rt and batch ID within group based on rt order
        mz[(breaks[i - 1] + 1):breaks[i]] <- mz[this.o]
        chr[(breaks[i - 1] + 1):breaks[i]] <- chr[this.o]
        lab[(breaks[i - 1] + 1):breaks[i]] <- lab[this.o]
    }
    
    # remove fist and last index
    breaks <- breaks[c(-1, -length(breaks))]
    
    if (is.na(rt_tol_relative)) {
        da <- 0
        # create indices for each groups == s (segments)
        if (length(breaks) > max.num.segments) {
            s <- floor(seq(2, length(breaks), length.out = max.num.segments))
        } else {
            s <- 2:length(breaks)
        }
        
        # compute distance matrix of each group using stats::dist
        for (i in s) {
            this.sel <- (breaks[i - 1] + 1):breaks[i]
            
            if (length(this.sel) <= 3 * number_of_samples) {
                this.d <- as.vector(dist(chr[this.sel]))
                if (length(this.d) > 100)
                    this.d <- sample(this.d, 100)
                da <- c(da, this.d)
            }
        }
        
        ### want to statistically decide the smallest rt which is still significant:
        
        # da is a long vector of distances between rt values (with no particular order?)
        da <- da[!is.na(da)]
        # maximal distance
        uppermost <- max(da)
        # number of equally spaced points at which the density is to be estimated
        n = min(max.bins, max(min.bins, round(length(da) / aver.bin.size)))
        # estimate probability density function of da
        des <- density(da, kernel = "gaussian", n = n,
                       bw = uppermost / n * 2, from = 0)
        # the n (-1?) coordinates of the points where the density is estimated
        y <- des$y[des$x > 0]
        # the estimated density values
        x <- des$x[des$x > 0]
        
        # linear regression using uppermost part of data
        this.l <- lm(y[x > uppermost / 4] ~ x[x > uppermost / 4])
        # compute probability density values (y) using the linear function
        exp.y <- this.l$coef[1] + this.l$coef[2] * x
        
        y2 <- y[1:(length(y) - 1)] # density values without the last one
        y3 <- y[2:(length(y))]     # density values without the first one
        
        # pair-wise copy greater value to the left
        y2[which(y2 < y3)] <- y3[which(y2 < y3)]
        y[1:(length(y) - 1)] <- y2
        
        # cumulative sum - where actual y is greater than calculated exp.y using regression?
        yy <- cumsum(y > 1.5 * exp.y)
        yi <- seq_along(yy)
        # find last index where y is greater than exp.y
        sel <- min(which(yy < yi)) - 1
        # corresponding coordinate is used as rt tolerance
        rt_tol_relative <- x[sel]
        
        if (do.plot) {
          tolerance_plot(x, y, exp.y, sel, main = "find retention time tolerance")
        }
    }
    
    # pair-wise rt differences
    da <- chr[2:l] - chr[1:(l - 1)]
    # find those respecting the computed tolerance
    breaks.2 <- which(da > rt_tol_relative)
    # merge and sort both group delimiters 
    all.breaks <- c(0, unique(c(breaks, breaks.2)), l)
    all.breaks <- all.breaks[order(all.breaks)]
    
    # create list of indices corresponding to a representative from each group 
    # (always the first element)
    grps <- seq_along(mz)
    for (i in 2:length(all.breaks)) {
        grps[(all.breaks[i - 1] + 1):all.breaks[i]] <- i
    }
    
    list(
        mz = mz,
        chr = chr,
        lab = lab,
        grps = grps,
        chr.tol = rt_tol_relative,
        mz.tol = mz_tol_relative
    )
}
