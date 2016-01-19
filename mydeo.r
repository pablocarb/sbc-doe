# Use always multiples of 2

latin <- function(n, nrand = 20) {
    x <- matrix(1:n, n, n)
    x <- t(x)
    for (i in 2:n) {
        x[i, ] <- x[i, c(i:n, 1:(i-1))]
    }
    if (nrand > 0) {
        for (i in 1:nrand) {
            x <- x[sample(n), ]
            x <- x[,sample(n)]
        }
    }
    return(x)
}

# Generate full permutation or latin square
permut <- function(n, ptype='latin') {
    require('combinat')
    if (ptype == 'latin') {
        return(latin(n))
    } else {
        return(permn(n))
    }
}

doe1 <- function(factors, nlevels, timeout=5) {
    require('planor')
    require('R.utils')
    require('combinat')
    construct <- planor.factors(
        factors=factors, 
        nlevels=nlevels
        )
    const.ffd <- list()
    # We check if we have enough resolution for the model 
    HAVERESOLUTION <- TRUE
    # Start with the simplest model: ~S
    res <- 3
    const.total <- prod(construct@fact.info$nlev)
    n <- 1
    const.ffd[[n]] <- NULL
    # Change to verify if we have enough resolution
    units <- 0
    while (units < const.total) {
        # Start with full combinations
        k <- 1
        des <- 0
        units <- const.total
        while (!is.null(des)) {
            des <- NULL
            try(
                evalWithTimeout(
                    des <- regular.design(factors=construct, resolution=res, nunits=const.total/k),
                    timeout=timeout, onTimeout='silent'
                    )
                )
                                        # Simplify S4 structure
            if (!is.null(des)) {
                const.ffd[[n]] <- list(design=slot(des, 'design'), resolution=res)
                units <- dim(const.ffd[[n]]$design)[1]
            }
            k <- k*2
        }
        # Increase to the next model
        n <- n+1
        res <- res+2
    }
    return(const.ffd)
}
