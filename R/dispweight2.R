### Veganified version of Eduard Szöcs's dispweight. This version
### reuses other vegan utilities as much as possible and follows vegan
### conventions. Basically, we primarily want to get transformed data
### with information on transformation. Only basically needed
### information are the overdispersion coefficients (D) and p-values;
### the weights can be found from these.
`dispweight2` <-
    function(comm, groups, nperm = 999, nullmodel = "c0_ind", plimit = 0.05,
             ...)
{
    if (missing(groups))
        groups <- rep(1, nrow(comm))
    ## vegan uses 'groups' (not 'group') in other similar
    ## functions. First remove empty levels
    groups <- factor(groups)
    ## Statistic is the sum of squared differences by 'groups'
    means <- apply(comm, 2, function(x) tapply(x, groups, mean))
    ## handle 1-level factors
    if (is.null(dim(means)))
        means <- matrix(means, nrow = 1, ncol = length(means),
                        dimnames=list(levels(groups), names(means)))
    fitted <- means[groups,]
    dhat <- colSums((comm - fitted)^2/fitted, na.rm = TRUE)
    ## Get df for non-zero blocks of species. Completely ignoring
    ## all-zero blocks for species sounds strange to me, but this was
    ## done in the original code, and we follow here. However, this
    ## was not done for significance tests, and only concerns 'D' and
    ## 'weights'.
    nreps <- table(groups)
    div <- colSums(sweep(means > 0, 1, nreps - 1, "*"))
    ## "significance" of overdispersion is assessed from Chi-square
    ## evaluated separately for each species. This means fixing only
    ## marginal totals for species but letting row marginals vary
    ## freely, unlike in standard Chi-square where both margins are
    ## fixed. In vegan this is achieved by nullmodel 'c0_ind'. Instead
    ## of one overall simulation, nullmodel is generated separately
    ## for each of 'groups'
    chisq <- function(x) {
        fitted <- colMeans(x)
        colSums(sweep(x, 2, fitted)^2, na.rm = TRUE) / fitted
    }
    simulated <- matrix(0, nrow = ncol(comm), ncol = nperm)
    for (lev in levels(groups)) {
        nm <- nullmodel(comm[groups == lev,], nullmodel)
        tmp <- apply(simulate(nm, nperm), 3, chisq)
        ok <- !is.na(tmp)
        simulated[ok] <- simulated[ok] + tmp[ok] 
    }
    p <- (rowSums(dhat <= simulated) + 1) / (nperm + 1)
    dhat <- dhat / div
    w <- ifelse(p <= plimit, 1/dhat, 1)
    comm <- sweep(comm, 2, w, "*")
    class(comm) <- c("dispweight", class(comm))
    attr(comm, "D") <- dhat
    attr(comm, "p") <- p
    attr(comm, "weights") <- w
    comm
}

