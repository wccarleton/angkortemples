# load libraries
library(nimble)
library(coda)
library(tidyr)
library(dplyr)

# pull in raw data
temple_known_dates <- read.csv("./Data/temple_known_dates.csv", as.is = T)
temple_vars <- read.csv("./Data/temple_vars.csv", as.is = T)

# collate: raw is in two separate tables that can be combined for easier use and then only required variables isolated and combined into a single working dataframe
t_dates <- data.frame(id = temple_known_dates$id,
                    year_ce = temple_known_dates$Date.to.use)
t_vars <- data.frame(id = temple_vars$Temple.ID,
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

# clean up the morphology variable
#temple_vars[grep("(east)", temple_vars$Morphology), "Morphology"] <- "horseshoe_east"
#temple_vars[grep("(north)", temple_vars$Morphology), "Morphology"] <- "horseshoe_north"
#temple_vars[grep("(west)", temple_vars$Morphology), "Morphology"] <- "horseshoe_west"
#temple_vars[grep("4causeway", temple_vars$Morphology), "Morphology"] <- "causeway_4"
#temple_vars[grep("2causeway", temple_vars$Morphology), "Morphology"] <- "causeway_2"
#temple_vars[grep("Square", temple_vars$Morphology), "Morphology"] <- "square"
#temple_vars[which(temple_vars$Morphology == ""), "Morphology"] <- NA

#unique(temple_vars$Morphology)

# merge the two dataframesusing the id column to match while ensuring all of the rows in the vars dataframe are included
temples <- right_join(t_dates, t_vars, by = "id")

# set up a Nimble model
templeCode <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100) # intercept
    for(j in 1:J){
        beta[j] ~ dnorm(0, sd = 100) # regression coefs
    }
    sigma ~ dunif(50, 200) # prior variance for regression model
    for(j in 1:(J - 2)){
        theta[j] ~ dbeta(a[j], b[j]) # prior for hot-encoded vars
    }
    for(n in 1:N){
        x[n, 1] ~ dunif(min = 1, max = 180) # azimuth
        x[n, 2] ~ dlnorm(meanlog = 8, sdlog = 1.03) # area
        for(j in 3:J){
            x[n, j] ~ dbern(theta[j - 2]) # hot-encoded covariates
        }
        temple_age[n] ~ dnorm(beta0 + inprod(beta[1:J], x[n, 1:J]), sd = sigma) # core regression model
    }
})

# have a look again at temples dataframe variables
names(temples)

# note there are 7 hot-encoded vars in the dataframe and 1 other continuous (bounded) covariate (azimuth)

# subset only temples with known dates
temples_known <- subset(temples, !is.na(year_ce))

templeConsts <- list(a = rep(1, 8),
                    b = rep(1, 8),
                    N = nrow(temples_known),
                    J = 10)

templeData <- list(temple_age = temples_known$year_ce,
                    x = temples_known[, c(3:12)])

templeInits <- list(theta = rep(0.5, 8),
                    beta0 = 0,
                    beta = rep(0, 10),
                    sigma = 100)

temple <- nimbleModel(code = templeCode,
                name = "temple",
                constants = templeConsts,
                data = templeData,
                inits = templeInits)

cv_config <- configureMCMC(model = temple)

cv_out <- runCrossValidate(MCMCconfiguration = cv_config,
                            k = nrow(temples_known),
                            MCMCcontrol = list(niter = 20000, nburnin = 2000),
                            nCores = 10,
                            nBootReps = NA)

mcmc_out <- nimbleMCMC(model = temple,
                        niter = 20000,
                        summary = T,
                        WAIC = T)

colnames(mcmc_out$samples)[1:25]

traceplot(mcmc(mcmc_out$samples[, 11]))

hist(mcmc_out$samples[-c(1:1000), 2])

pairs(mcmc_out$samples[-c(1:1000), c(1:9)])
