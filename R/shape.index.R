shape.index <- function(x) {
        total.size <- length(x)
        unique.tss <- unique(x)
        n.unique.tss <- length(unique.tss)
        p.array <- array(NA,c(1,n.unique.tss))
        for (i in 1:n.unique.tss) {
            n.o <- x
            n.i <- length(which(unique.tss[i]==n.o))
            p.sub.i <- n.i/total.size
            p.log.2 <- log2(p.sub.i)
            p.array[1,i] <- (p.sub.i*p.log.2)
        }
        p.total <- sum(p.array)
        SI <- 2+p.total
        return(SI)
    }
