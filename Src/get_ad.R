# libs
library(nimble)
library(tidyverse)

# args
args <- commandArgs(trailingOnly = TRUE)
left_out <- as.numeric(args)

message(paste("temple index:", left_out))

# data

temples <- read.csv(file = "./Data/temples.csv")
temples$morph <- as.factor(temples$morph)

# isolate desired variables

str_pattern <- paste("morph",
                    "azimuth",
                    "area",
                    "trait_1",
                    "trait_2",
                    "trait_3",
                    "trait_4",
                    "trait_5",
                    "trait_6",
                    "trait_7",
                    "trait_8",
                    sep = "|")

covariate_idx <- grep(str_pattern, colnames(temples))

# Create a column containing only the emprical dates with all other entries
# NA. This will be useful for the code that follows because Nimble will
# automatically impute/predict values for NA entries.

empirically_dated_idx <- which(temples$date_type == "empirical")
temples$date_emp <- NA
temples[empirically_dated_idx, "date_emp"] <- temples[empirically_dated_idx, "date"]

# save the column idx for the new column
date_idx <- grep("^date_emp$", names(temples))

# isolate dated temples (for cross validation comparisons)

temples_dated <- temples[which(!is.na(temples$date_emp)),]

# make it all numeric

x <- mutate(temples_dated[, covariate_idx], across(morph:trait_8, as.numeric))

x <- as.data.frame(x)

# save some useful parameters to pass to the nimble model, like number of
# observations etc.
M <- length(levels(temples$morph)) # number of temple morpho types in total
H <- 8 # number of binary variables
Cont <- 2 # number of continuous variables
J <- H + Cont # total number of variables (covariates/predictors)
N <- nrow(temples_dated) # N obs.

# set up and run a nimble model

# set up a Nimble model
# this model is predicting temple ages with a set a covariates and uses 
# temple morphology as an index variable (rather than one-hot encoding etc)
templeCode <- nimbleCode({
    #beta0 ~ dnorm(1000, sd = 500)
    #sigma0 ~ dunif(0, 500) # prior variance for morpho types (index variable)
    morpho_prob[1:M] ~ ddirch(alpha = d_alpha[1:M])
    for(m in 1:M){
        morpho[m] ~ dnorm(1000, sd = 200)#dnorm(beta0, sd = sigma0)
    }
    for(j in 1:J){
        beta[j] ~ dnorm(0, sd = 500) # regression coefs
    }
    sigma ~ dunif(0, 500) # prior variance for regression model
    for(h in 1:H){
        theta[h] ~ dbeta(a[h], b[h]) # prior for binary vars
    }
    for(n in 1:N){
        x[n, 1] ~ dcat(prob = morpho_prob[1:M])
        x[n, 2] ~ dunif(min = 1, max = 360) # azimuth
        x[n, 3] ~ dlnorm(meanlog = 7.7, sdlog = 1.02) # area
        for(h in 1:H){
            x[n, 3 + h] ~ dbern(theta[h]) # binary covariates
        }
        mu[n] <- morpho[x[n, 1]] + inprod(beta[1:J], x[n, 2:(J + 1)])
        temple_age[n] ~ dnorm(mu[n], sd = sigma) # core model
    }
})

templeConsts <- list(a = rep(1, H), # beta prior
                    b = rep(1, H), # beta prior
                    N = N,
                    M = M,
                    H = H,
                    J = J,
                    d_alpha = rep(1, M)) # parameter vector for Dirichlet prior

# capture true date
left_out_date <- temples_dated$date_emp[left_out]

message(paste("Left out date: ", left_out_date))

# set the left out temple date to NA
temples_dated$date_emp[left_out] <- NA

templeData <- list(temple_age = temples_dated$date_emp,
                    x = x)

templeInits <- list(theta = rep(0.5, H),
                    #beta0 = 1000,
                    #sigma0 = 200,
                    beta = rep(0, J),
                    morpho = rep(0, M),
                    sigma = 100)

templeModel <- nimbleModel(code = templeCode,
                name = "temple",
                constants = templeConsts,
                data = templeData,
                inits = templeInits)

params_to_track <- c(#"beta0",
                    #"sigma0",
                    "morpho", 
                    "morpho_prob", 
                    "beta", 
                    "sigma", 
                    "mu", 
                    "temple_age")

# mcmc config options

# change default block sampling for beta[] and morpho[] to AF_slice
templeModel_c <- compileNimble(templeModel)
temple_mcmc_config <- configureMCMC(templeModel_c)
#temple_mcmc_config$removeSamplers(c("beta", "morpho"))
#temple_mcmc_config$addSampler(target = c("beta", "morpho"), type = "AF_slice")
temple_mcmc_config$setMonitors(params_to_track)

# build mcmc
temple_mcmc <- buildMCMC(temple_mcmc_config)
temple_mcmc_c <- compileNimble(temple_mcmc)

# run mcmc

mcmc_out <- runMCMC(temple_mcmc_c, 
                    niter = 50000,
                    nburnin = 5000)

# extract mean predicted age for the left out temple
left_out_post_col <- grep(paste("mu\\[", left_out, "\\]", sep = ""), 
                        colnames(mcmc_out))

cv_abs_dev <- abs(left_out_date - mean(mcmc_out[, left_out_post_col]))

write.table(cv_abs_dev, 
            file = "Output/cv_abs_devs_2.csv", 
            append = T,
            row.names = F,
            col.names = F)