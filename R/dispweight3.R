`dispweight3` <-
    function(formula, data, plimit = 0.05)
{
    comm <- eval(formula[[2]])
    x <- model.matrix(delete.response(terms(formula)), data = data)
    family <- quasipoisson()
    V <- family$variance
    mods <- lapply(comm, function(y) glm.fit(x, y, family = family))
    y <- sapply(mods, function(x) x$y)
    mu <- sapply(mods, function(x) x$fitted.values)
    wts <- sapply(mods, function(x) x$prior.weights)
    res <- (y-mu) * sqrt(wts) / sqrt(V(mu))
    df <- sapply(mods, function(x) x$df.residual)
    stat <- colSums(res^2)
    p <- pchisq(stat, df, lower.tail = FALSE)
    dhat <- stat/df
    w <- ifelse(p < plimit, 1/dhat, 1)
    w <- ifelse(w > 1, 1, w)
    comm <- sweep(comm, 2, w, "*")
    class(comm) <- c("dispweight", class(comm))
    attr(comm, "D") <- dhat
    attr(comm, "p") <- p
    attr(comm, "weights") <- w
    comm
}
