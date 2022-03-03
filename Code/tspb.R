library(copula)
library(dgof)

ks.ts.pb <- function(x, dist, B = 1000) {
    ## example assuming normal distribution; normal copula
    ## get observed stat
    n <- length(x)
    stat  <- ks.test(x, dist)$statistic
    ## get lag-1 sample auto-spearman rho
    rho  <-  cor(x[-1], x[-n], method = "spearman")
    r <- iRho(normalCopula(), rho)
    ## parametric bootstrap to get an empirical distribution of stat
    stat.b <- double(n)
    ## set up containers
    u <- double(n)
    for (i in 1:B) {
        ## get one bootstrap sample
        ## 1. the marginal distribution is the fitted dist
        ## 2. there is temporal dependence
        u <- runif(1)
        for (j in 2:n) {
            ## could test other copulas later
            u[j]  <- cCopula(u[j - 1], normalCopula(r))
        }
        x.b <- qnorm(u)
        ## fit it
        mu.b <- mean(x.b)
        sd.b <- sd(x.b)
        stat.b[i] <- ks.test(x.b, dist, mu.b, sd.b)$statistic
    }
    
    p.value <- (sum(stat.b >= stat) + 0.5) / (B + 1)
    
    return(p.value)
    # return empirical p-value    
}

n <- 100
ar = 0.5
x <- arima.sim(list(ar = ar), rand.gen = rnorm, sd = sqrt(1-ar^2), n = n)
ks.ts.pb(x, "pnorm")
