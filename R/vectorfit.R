"vectorfit" <-
    function (X, P, permutations = 0, strata, w, ...) 
{
    if (missing(w) || is.null(w)) 
        w <- 1
    if (length(w) == 1) 
        w <- rep(w, nrow(X))
    Xw <- .C("wcentre", x = as.double(X), as.double(w), as.integer(nrow(X)),
             as.integer(ncol(X)), PACKAGE = "vegan")$x
    dim(Xw) <- dim(X)
    P <- as.matrix(P)
    Pw <- .C("wcentre", x = as.double(P), as.double(w), as.integer(nrow(P)),
             as.integer(ncol(P)), PACKAGE = "vegan")$x
    dim(Pw) <- dim(P)
    colnames(Pw) <- colnames(P)
    nc <- ncol(X)
    Q <- qr(Xw)
    H <- qr.fitted(Q, Pw)
    heads <- qr.coef(Q, Pw)
    r <- diag(cor(H, Pw)^2)
    r[is.na(r)] <- 0
    heads <- decostand(heads, "norm", 2)
    heads <- t(heads)
    if (is.null(colnames(X))) 
        colnames(heads) <- paste("Dim", 1:nc, sep = "")
    else colnames(heads) <- colnames(X)
    ## make permutation matrix for all variables handled in the next loop
    nr <- nrow(X)
    if (length(permutations) == 1) {
        if (permutations > 0) {
            arg <- if(missing(strata)) NULL else strata
            permat <- t(replicate(permutations,
                                  permuted.index(nr, strata = arg)))
        }
    } else {
        permat <- as.matrix(permutations)
        if (ncol(permat) != nr)
            stop(gettextf("'permutations' have %d columns, but data have %d rows",
                          ncol(permat), nr))
        permutations <- nrow(permutations)
    }
    if (permutations) {
        ptest <- function(indx, ...) {
            take <- P[indx, , drop = FALSE]
            take <- .C("wcentre", x = as.double(take), as.double(w),
                       as.integer(nrow(take)), as.integer(ncol(take)),
                       PACKAGE = "vegan")$x
            dim(take) <- dim(P)
            Hperm <- qr.fitted(Q, take)
            diag(cor(Hperm, take))^2
        }
        permstore <- sapply(1:permutations, function(indx, ...) ptest(permat[indx,], ...))
        ## Single variable is dropped to a vector, and otherwise
        ## permutations are the matrix columns and variables are rows
        if (!is.matrix(permstore))
            permstore <- matrix(permstore, ncol=permutations)
        permstore <- sweep(permstore, 1, r, ">=")
        validn <- rowSums(is.finite(permstore))
        pvals <- (rowSums(permstore, na.rm = TRUE) + 1)/(validn + 1)
    }
    else pvals <- NULL
    sol <- list(arrows = heads, r = r, permutations = permutations, 
                pvals = pvals)
    if (!missing(strata)) {
        sol$strata <- deparse(substitute(strata))
        sol$stratum.values <- strata
    }
    class(sol) <- "vectorfit"
    sol
}
