#
#mydeo (c) University of Manchester 2015
#
#mydeo is licensed under the MIT License.
#
#To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.
#
#@author: Pablo Carbonell, SYNBIOCHEM
#@description: Design of experiments routines for combinatorial assembly using R libraries
#
                                        # Use always multiples of 2

repo <- "https://cloud.r-project.org"
if (!require('DoE.base')) {
    if (!require('gmp',quietly=TRUE)) {install.packages('gmp', repos=repo, verbose=FALSE)}
    if (!require('gmp', quietly=TRUE)) {install.packages('numbers', repos=repo, verbose=FALSE)}
    if (!require('crossdes', quietly=TRUE)) {install.packages('crossdes', repos=repo, verbose=FALSE)}
    if (!require('combinat', quietly=TRUE, warn.conflicts = FALSE)) {install.packages('combinat', repos=repo, verbose=FALSE)}
    if (!require('planor', quietly=TRUE)) {install.packages('planor', repos=repo, verbose=FALSE)}
    if (!require('R.utils', quietly=TRUE)) {install.packages('R.utils', repos=repo, verbose=FALSE)}
    if (!require('DoE.base', quietly=TRUE)) {install.packages('DoE.base', repos=repo, verbose=FALSE)}
}

latin.augment <- function(m1, m2) {
    ms <- list()
    for (x in seq(1, max(m1))) {
        ms[[x]] <- m2 + (x-1)*max(m2)
    }
    m4 <- NULL
    for (i in seq(1, ncol(m1))) {
        m3 <- NULL
        for (j in seq(1, ncol(m1))) {
            if (is.null(m3)) {
                m3 <- ms[[ m1[i,j] ]]
            } else {
                m3 <- cbind( m3, ms[[ m1[i,j] ]] )
            }
        }
        if (is.null(m4)) {
            m4 <- m3
        } else {
            m4 <- rbind(m4, m3)
        }
            
    }
    return(m4)
}

latin <- function(ntotal) {
    if (!require('gmp', quietly=TRUE)) {
        require('numbers', quietly=TRUE)
        useGMP <- FALSE
    } else {
        useGMP <- TRUE
    }
    require('crossdes', quietly=TRUE)
    # Find all prime factors
    if (useGMP) {
	if (ntotal > 1) {
        	m <- factorize(ntotal)
        	m <-  as.integer(m)
	} else {
		m <- c(1)
	}
    } else {
        m <- primeFactors(ntotal)
    }
    sq <- list()
    # Join any 2 * 2 as 4 (orthogonal latin squares only for n >= 3)
    fours <- c()
    while ( (length(m) >=2) && (m[2] == 2) ) {
        if (length(m) > 2) {
            m <- m[3:length(m) ]
            fours <- c(fours, 4)
        } else {
            m <- c(4)
        }
    }
    m <- c(fours, m)
    # Compute orthogonal squares of each factor
    for (n in m) {
        if (n > 2) {
            x <- des.MOLS(n,n)
        } else {
            if (n ==2 ) {
                x <- matrix(c(1,2,2,1), nrow=2)
            } else {
                x <- matrix(c(1))
            }
        }
        sq[[n]] <- x
    
    }
    # Augment the matrix with all prime factors
    m1 <- sq[[ m[1] ]]
    if (length(m) == 1) {
        return( m1 )
    } else {
        m2 <- sq[[ m[2] ]]
        mx <- latin.augment(m1, m2)
        if (length(m) > 2) {
            for (i in seq(3, length(m))) {
                m3 <- sq[[ m[i] ]]
                mx <- latin.augment(mx, m3)
            }
        }
    }
             
    return(mx)
}



latin.old <- function(n, nrand = 20) {
    library('crossdes')
    if (n > 2) {
        x <- des.MOLS(n,n)[1:n, 1:n]
    } else {
        if (n ==2 ) {
            x <- matrix(c(1,2,2,1), nrow=2)
        } else {
            x <- matrix(c(1))
        }
    }
    
#    x <- matrix(1:n, n, n)
#    x <- t(x)
#    for (i in 2:n) {
#        x[i, ] <- x[i, c(i:n, 1:(i-1))]
#    }
    if (n > 2) {
        if (nrand > 0) {
            for (i in 1:nrand) {
                x <- x[sample(n), ]
                x <- x[,sample(n)]
            }
        }
        x <- x[, order(x[1,])]
    }
    return(x)
}

# Generate full permutation or latin square
permut <- function(n, ptype='latin', randomize=TRUE) {
    require('combinat', quietly=TRUE, warn.conflicts = FALSE)
    if (ptype == 'latin') {
        if (randomize) {
            return(latin(n))
        } else {
            return(latin(n))
        }
    } else {
        return(permn(n))
    }
}

# Regular factorial design with implicit S model
doe1 <- function(factors, nlevels, timeout=5) {
    print(factors)
    require('planor', quietly=TRUE)
    require('R.utils', quietly=TRUE)
    require('combinat', quietly=TRUE)
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

# Orthogonal arrays
doe2 <- function(factors, nlevels, timeout=5, seed=888) {
    require('DoE.base', quietly=TRUE)
    const.oa <- list()
    des <- oa.design(factor.names=factors, nlevels=nlevels, seed=seed)
    const.oa[[1]] <- list(design=des, resolution=NULL)
    return(const.oa)
}
