find.turn.point <-
function(y)
{
    peaks2<-function (x, ties.method)
    {
        # investigate, this is not clear
        z <- embed(rev(as.vector(c(-Inf, x, -Inf))), dim = 3)
        z <- z[rev(seq(nrow(z))), ]
        v <- max.col(z,ties.method=ties.method) == 2
        v
    }
    msExtrema<-function (x)
    {
        l<-length(x)
        # finds maximal values
        index1 <- peaks2(x, ties.method="first")
        # finds minimal values
        index2 <- peaks2(-x, ties.method="last")
        # remove overlaps? is it even possible?
        index.max <- index1 & !index2
        index.min <- index2 & !index1
        list(index.max = index.max, index.min = index.min)
    }
    
    y <- y[!is.na(y)]
    # if all values are the same
    if (length(unique(y)) == 1) {
        # peak is middle
        pks <- round(length(y)/2)
        # valleys are the first and the last
        vlys <- c(1, length(y))
        x <- new("list")
        x$pks <- pks
        x$vlys <- vlys
        return(x)
    }
    
    b<-msExtrema(y)
    pks<-which(b$index.max)
    vlys<-which(b$index.min)
    # if peak is not on index 1, add index 1 to valleys ??
    if(pks[1] != 1) vlys<-c(1, vlys)
    # the same for the last index
    if(pks[length(pks)] != length(y)) vlys<-c(vlys, length(y))
    # finally if there is only one peak, valleys are on the first and the last index
    # (probably because Gaussian was used)
    if(length(pks) == 1) vlys<-c(1,length(y))
    x <- new("list")
    x$pks <- pks
    x$vlys <- vlys
    return(x)
}
