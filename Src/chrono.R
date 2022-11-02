library(nimble)
library(coda)

templeCode <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100) # intercept
    for(j in 1:J){
        beta[j] ~ dnorm(0, sd = 100) # regression coefs
    }
    sigma ~ dunif(0, 100) # prior variance for regression model
    for(j in 1:(J - 1)){
        theta[j] ~ dbeta(a[j], b[j]) # prior for hot-coded vars
    }
    for(n in 1:N){
        x[n, 1] ~ dunif(min = 1, max = 360) # azimuth
        for(j in 2:J){
            x[n, j] ~ dbern(theta[j - 1]) # hot-coded covariates
        }
        temple_age[n] ~ dnorm(beta0 + inprod(beta[1:J], x[n, 1:J]), sd = sigma) # core regression model
    }
})

N <- 200
J <- 6

xmat <- matrix(nrow = N, ncol = J)
xmat[, 1] <- runif(n = N, min = 1, max = 360)
for(j in 2:J){
    xmat[, j] <- rbinom(n = N, size = 1, prob = 0.5)
}
mu <- 800 + xmat %*% rep(4, J)
y <- rnorm(n = N, mu, sd = 10)

templeConsts <- list(a = rep(1, J - 1),
                    b = rep(1, J - 1),
                    N = N,
                    J = J)

templeData <- list(temple_age = y,
                    x = xmat)

templeInits <- list(theta = rep(0.5, J - 1),
                    beta0 = 0,
                    beta = rep(0, J),
                    sigma = 10)

temple <- nimbleModel(code = templeCode,
                name = "temple",
                constants = templeConsts,
                data = templeData,
                inits = templeInits)

temple$getNodeNames()
temple$temple_age

mcmc_out <- nimbleMCMC(model = temple,
                        niter = 20000,
                        summary = T,
                        WAIC = T)

traceplot(mcmc(mcmc_out$samples[-c(1:1000), 1]))
hist(mcmc_out$samples[-c(1:1000), 8])

colnames(mcmc_out$samples)
