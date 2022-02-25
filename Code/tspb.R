n <- 100
x <- arima.sim(list(ar = 0.4), n)

ks.ts.pb <- function(x, dist, B = 1000) {
    ## example assuming normal distribution; normal copula
    ## get observed stat
    n <- length(x)
    theta <- fit(x) # theta = (mu, sigma)
    stat  <- ks.test(x, dist, theta)
    ## get lag-1 sample auto-spearman rho
    rho  <-  cor(x[-1], x[-n], method = "spearman")
    r <- iRho(normalCopula(), rho)
    ## parametric bootstrap to get an empirical distribution of stat
    ## set up containers
    u <- double(n)
    for (i in 1:B) {
        ## get one bootstrap sample
        ## 1. the marginal distribution is the fitted dist
        ## 2. there is temporal dependence
        u <- runif(1)
        for (j in 2:n) {
            ## could test other copulas later
            u[j]  <- copula::cCopula(u[j - 1], normalCopula(r))
        }
        x.b <- qnorm(u, mean = theta[1], sd = theta[2])
        ## fit it
        theta.b <- fit(x.b)
        stat.b <- ks.test(x.b, dist, theta.b)
    }
    # return empirical p-value    
}
