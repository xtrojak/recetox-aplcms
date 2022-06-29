to_attach <- function(pick, number_of_samples, use = "sum") {
    strengths <- rep(0, number_of_samples)
    if (is.null(nrow(pick))) {
        # this is very strange if it can ever happen
        # maybe commas are missing? we want the same as below
        # but also why if there are no rows...
        strengths[pick[6]] <- pick[5]
        return(c(pick[1], pick[2], pick[1], pick[1], strengths))
    } else {
        for (i in seq_along(strengths)) {
            if (use == "sum")
                # sum of all areas from the same sample
                strengths[i] <- sum(pick[pick[, 6] == i, 5])
            if (use == "median")
                # median of all areas from the same sample
                strengths[i] <- median(pick[pick[, 6] == i, 5])
        }
        # average of m/z, average of rt, min of rt, max of rt, sum/median of areas
        return(c(mean(pick[, 1]), mean(pick[, 2]), min(pick[, 1]),
                 max(pick[, 1]), strengths))
    }
}

# returns a list of aligned features and original peak times
feature.align <- function(features,
                          min_occurrence = 2,
                          mz_tol_relative = NA,
                          rt_tol_relative = NA,
                          mz_max_diff = 1e-4,
                          mz_tol_absolute = 0.01,
                          do.plot = TRUE,
                          rt_colname = "pos") {
    if (do.plot) {
        par(mfrow = c(3, 2))
        draw_plot(label = "Feature alignment", cex = 2)
        draw_plot()
    }
    
    number_of_samples <- nrow(summary(features))
    if (number_of_samples > 1) {
        values <- get_feature_values(features, rt_colname)
        mz_values <- values$mz  # vector of m/z values from all samples
        chr <- values$chr       # same for rt
        lab <- values$lab       # sample identifiers
        
        # sort all values by mz, if equal by rt
        o <- order(mz_values, chr)
        mz_values <- mz_values[o]
        chr <- chr[o]
        lab <- lab[o]
        
        # find relative m/z tolerance level
        if (is.na(mz_tol_relative)) {
            mz_tol_relative <- find.tol(mz_values, mz_max_diff = mz_max_diff, do.plot = do.plot)
            if (length(mz_tol_relative) == 0) {
                mz_tol_relative <- 1e-5
                warning(
                    "Automatic tolerance finding failed, 10 ppm was assigned. May need to manually assign alignment mz tolerance level."
                )
            }
        } else if (do.plot) {
            draw_plot(main = "alignment m/z tolerance level given", 
                      label = mz_tol_relative, cex = 1.2)
        }
        
        if (!is.na(rt_tol_relative) && do.plot) {
            draw_plot(main = "retention time \n tolerance level given", 
                      label = rt_tol_relative, cex = 1.2)
        }
        
        # find relative retention time tolerance level
        all.ft <- find.tol.time(mz_values,
                                chr,
                                lab,
                                number_of_samples = number_of_samples,
                                mz_tol_relative = mz_tol_relative,
                                rt_tol_relative = rt_tol_relative,
                                mz_tol_absolute = mz_tol_absolute,
                                do.plot = do.plot)
        rt_tol_relative <- all.ft$chr.tol
        
        message("**** performing feature alignment ****")
        message(paste("m/z tolerance level: ", mz_tol_relative))
        message(paste("time tolerance level:", rt_tol_relative))
        
        # create zero vectors of length number_of_samples + 4 ?
        aligned.ftrs <- pk.times <- rep(0, 4 + number_of_samples)
        mz.sd.rec <- 0
        
        labels <- unique(all.ft$grps)
        
        area <- grps <- mz_values
        
        # grouping the features based on their m/z values (assuming the tolerance level)
        sizes <- c(0, cumsum(sapply(features, nrow)))
        for (i in 1:number_of_samples) {
            this <- features[[i]]
            # select values belonging to this batch/sample
            sel <- which(all.ft$lab == i)
            # all.ft values are already grouped by m/z and rt tolerances *
            that <- cbind(all.ft$mz[sel], all.ft$chr[sel], all.ft$grps[sel])
            # order by m/z then by rt
            this <- this[order(this[, 1], this[, 2]),] # this can go 2 lines higher
            # *but here we sort them again?
            that <- that[order(that[, 1], that[, 2]),]
            # aren't this and that identical on the m/z and rt columns???
            
            # update m/z values with ordered ones
            mz_values[(sizes[i] + 1):sizes[i + 1]] <- this[, 1]
            # assign rt
            chr[(sizes[i] + 1):sizes[i + 1]] <- this[, 2]
            # assign area
            area[(sizes[i] + 1):sizes[i + 1]] <- this[, 5]
            # assign row identifier
            grps[(sizes[i] + 1):sizes[i + 1]] <- that[, 3]
            # assign batch identifier
            lab[(sizes[i] + 1):sizes[i + 1]] <- i
        }
        
        # table with number of values in a group
        ttt <- table(all.ft$grps)
        # count those with minimal occurrence (times 3 ??)
        curr.row <- sum(ttt >= min_occurrence) * 3
        mz.sd.rec <- rep(0, curr.row)
        
        # names creates strings which are again converted to number?
        sel.labels <- as.numeric(names(ttt)[ttt >= min_occurrence])
        
        # retention time alignment
        aligned.ftrs <-
            foreach(i = seq_along(sel.labels), .combine = rbind) %do% {
                if (i %% 100 == 0)
                    gc()
                this.return <- NULL
                # select a group
                sel <- which(grps == sel.labels[i])
                if (length(sel) > 1) {
                    # selected data from the group
                    # WHY 3 times rt??? are columns 3 and 4 even used? can be NA or removed...
                    this <- cbind(mz_values[sel], chr[sel], chr[sel], chr[sel], area[sel], lab[sel])
                    # continue if data is from at least 'min_occurrence' samples
                    if (length(unique(this[, 6])) >= min_occurrence) {
                        # Kernel Density Estimation with target standard deviation 'mz_tol_relative' time median  of m/z
                        this.den <- density(this[, 1], bw = mz_tol_relative * median(this[, 1]))
                        # Finds the peaks and valleys of a smooth curve
                        # for Gaussian distribution (default), isn't this straightforward?
                        # always first, mid, last?
                        turns <- find.turn.point(this.den$y)
                        pks <- this.den$x[turns$pks]
                        vlys <- this.den$x[turns$vlys]
                        for (j in seq_along(pks)) {
                            ### extract m/z selection
                            this.lower <- max(vlys[vlys < pks[j]])
                            this.upper <- min(vlys[vlys > pks[j]])
                            # select data with m/z within lower and upper bound from density estimation
                            this.sel <- which(this[, 1] > this.lower & this[, 1] <= this.upper)
                            that <- this[this.sel, ]
                            if (!is.null(nrow(that))) {
                                # continue if data is still from at least 'min_occurrence' samples
                                if (length(unique(that[, 6])) >= min_occurrence) {
                                    # Kernel Density Estimation, this time for retention time
                                    that.den <- density(that[, 2], bw = rt_tol_relative / 1.414)
                                    # again select statistically significant points
                                    that.turns <- find.turn.point(that.den$y)
                                    that.pks <- that.den$x[that.turns$pks]
                                    that.vlys <- that.den$x[that.turns$vlys]
                                    for (k in seq_along(that.pks)) {
                                        ### extract rt selection
                                        that.lower <- max(that.vlys[that.vlys < that.pks[k]])
                                        that.upper <- min(that.vlys[that.vlys > that.pks[k]])
                                        # select data with rt within lower and upper bound from density estimation
                                        thee <- that[that[, 2] > that.lower & that[, 2] <= that.upper, ]
                                        if (!is.null(nrow(thee))) {
                                            if (length(unique(thee[, 6])) >= min_occurrence) {
                                                # continue if data is still from at least 'min_occurrence' samples
                                                this.return <-
                                                    c(to_attach(thee, number_of_samples, use = "sum"),
                                                      to_attach(thee[, c(1, 2, 3, 4, 2, 6)], number_of_samples, use = "median"),
                                                      sd(thee[, 1], na.rm = TRUE)
                                                    )
                                                # vector of :
                                                # average of m/z, average of rt, min of m/z, max of m/z, sum of areas per sample
                                                # average of m/z, average of rt, min of m/z, max of m/z, median of rt per sample
                                                # and standard deviation of m/z
                                                # such vector is rbind-ed to resulting aligned.ftrs
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (min_occurrence == 1) {
                        thee <- c(mz_values[sel], chr[sel], chr[sel], chr[sel], area[sel], lab[sel])
                        this.return <- c(to_attach(thee, number_of_samples, use = "sum"),
                                         to_attach(thee[c(1, 2, 3, 4, 2, 6)], number_of_samples, use = "median"),
                                         NA
                        )
                    }
                }
                this.return
            }
        
        # select columns: average of m/z, average of rt, min of m/z, max of m/z, median of rt per sample (the second to_attach call)
        pk.times <- aligned.ftrs[, (5 + number_of_samples):(2 * (4 + number_of_samples))]
        # select last column, i.e. standard deviation of m/z
        mz.sd.rec <- aligned.ftrs[, ncol(aligned.ftrs)]
        # select columns: average of m/z, average of rt, min of m/z, max of m/z, sum of areas per sample (the first to_attach call)
        aligned.ftrs <- aligned.ftrs[, 1:(4 + number_of_samples)]
        
        # rename columns on both tables, samples are called "exp_i"
        colnames(aligned.ftrs) <-
            colnames(pk.times) <- c("mz", "time", "mz.min", "mz.max", paste("exp", 1:number_of_samples))
        
        # return both tables and both computed tolerances
        rec <- new("list")
        rec$aligned.ftrs <- aligned.ftrs
        rec$pk.times <- pk.times
        rec$mz.tol <- mz_tol_relative
        rec$chr.tol <- rt_tol_relative
        
        if (do.plot) {
            hist(mz.sd.rec, xlab = "m/z SD", ylab = "Frequency",
                 main = "m/z SD distribution")
            hist(apply(pk.times[, -1:-4], 1, sd, na.rm = TRUE), 
                 xlab = "Retention time SD", ylab = "Frequency",
                 main = "Retention time SD distribution")
        }
        
        return(rec)
    } else {
        message("There is but one experiment.  What are you trying to align?")
        return(0)
    }
}
