# load libraries
library(nimble)
library(coda)
library(tidyr)
library(dplyr)

# pull in raw data
temple_known_dates <- read.csv("./Data/temple_known_dates.csv", as.is = T)
temple_vars <- read.csv("./Data/temple_vars.csv", as.is = T)

# clean up the morphology variable
temple_vars[grep("(east)", temple_vars$Morphology), "Morphology"] <- "horseshoe_east"
temple_vars[grep("(north)", temple_vars$Morphology), "Morphology"] <- "horseshoe_north"
temple_vars[grep("(west)", temple_vars$Morphology), "Morphology"] <- "horseshoe_west"
temple_vars[grep("4causeway", temple_vars$Morphology), "Morphology"] <- "causeway_4"
temple_vars[grep("2causeway", temple_vars$Morphology), "Morphology"] <- "causeway_2"
temple_vars[grep("Square", temple_vars$Morphology), "Morphology"] <- "square"
temple_vars[which(temple_vars$Morphology == ""), "Morphology"] <- NA

# collate: raw is in two separate tables that can be combined for easier use and then only required variables isolated and combined into a single working dataframe
t_dates <- data.frame(id = temple_known_dates$id,
                    year_ce = temple_known_dates$Date.to.use)
t_vars <- data.frame(id = temple_vars$Temple.ID,
                    morph = as.factor(temple_vars$Morphology),
                    azimuth = temple_vars$Azimuth,
                    area = log(temple_vars$Area),
                    trait_1 = temple_vars$Principle.Reservoir,
                    trait_2 = temple_vars$Moat,
                    trait_3 = temple_vars$Sandstone,
                    trait_4 = temple_vars$Pink.Sandstone,
                    trait_5 = temple_vars$Laterite,
                    trait_6 = temple_vars$Brick,
                    trait_7 = temple_vars$Thmaphnom,
                    trait_8 = temple_vars$other)

# merge the two dataframesusing the id column to match while ensuring all of the rows in the vars dataframe are included
temples <- right_join(t_dates, t_vars, by = "id")

# set up a Nimble model
templeCode <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100) # intercept
    for(m in 1:M){
        morph_delta[m] ~ dnorm(0, sd = 100)
    }
    for(j in 1:J){
        beta[j] ~ dnorm(0, sd = 100) # regression coefs
    }
    sigma ~ dunif(0, 150) # prior variance for regression model
    for(j in 1:(J - 2)){
        theta[j] ~ dbeta(a[j], b[j]) # prior for hot-encoded vars
    }
    for(n in 1:N){
        morph[n] ~ dcat(prob = morph_prob[1:M]) # morpho types
        x[n, 1] ~ dunif(min = 1, max = 180) # azimuth
        x[n, 2] ~ dlnorm(meanlog = 8, sdlog = 1.03) # area
        for(j in 3:J){
            x[n, j] ~ dbern(theta[j - 2]) # hot-encoded covariates
        }
        temple_age[n] ~ dnorm(beta0 + morph_delta[morph[n]] + inprod(beta[1:J], x[n, 1:J]), sd = sigma) # core model
    }
})

# subset only temples with known dates
temples_known <- subset(temples, !is.na(year_ce))

# subset only complete cases
temples_complete <- temples_known[complete.cases(temples_known), ]

M <- length(levels(temples$morph))

templeConsts <- list(a = rep(1, 8), # a,b are the parameters of the hot-encoded probability priors and these are naive
                    b = rep(1, 8),
                    morph_prob = rep(1 / M, M), # prob of observing morpho types---this is naive
                    N = nrow(temples_complete),
                    M = M,
                    J = 10)

templeData <- list(temple_age = temples_complete$year_ce,
                    x = temples_complete[, c(4:13)],
                    morph = as.numeric(temples_complete$morph))

templeInits <- list(theta = rep(0.5, 8),
                    morph_delta = rep(0, M),
                    beta0 = 0,
                    beta = rep(0, 10),
                    sigma = 100)

temple <- nimbleModel(code = templeCode,
                name = "temple",
                constants = templeConsts,
                data = templeData,
                inits = templeInits)

mcmc_out <- nimbleMCMC(model = temple,
                        niter = 20000,
                        summary = T,
                        WAIC = T)

cv_config <- configureMCMC(model = temple)

cv_out <- runCrossValidate(MCMCconfiguration = cv_config,
                            k = 10,
                            MCMCcontrol = list(niter = 20000, nburnin = 2000),
                            nCores = 1,
                            nBootReps = NA)

traceplot(mcmc(mcmc_out$samples[, 11]))

hist(mcmc_out$samples[-c(1:1000), 2])

pairs(mcmc_out$samples[-c(1:1000), c(1:9)])

# simple lm for comparison
lm_temples <- lm(year_ce ~
                morph +
                azimuth + 
                area + 
                trait_1 + 
                trait_2 + 
                trait_3 + 
                trait_4 + 
                trait_5 + 
                trait_6 + 
                trait_8, 
                data = temples_complete)