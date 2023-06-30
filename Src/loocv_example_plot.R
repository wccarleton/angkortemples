# plotting examples from LOO-CV

# libs
library(nimble)
library(tidyverse)

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

# save the column idx for the new column
date_idx <- grep("^date_emp$", names(temples))

# isolate dated temples (for cross validation comparisons)

temples_dated <- temples[which(!is.na(temples$date_emp)), ]

# args
left_out <- sample(1:dim(temples_dated)[1], size = 1)

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
    morpho_prob[1:M] ~ ddirch(alpha = d_alpha[1:M])
    for(m in 1:M){
        morpho[m] ~ dnorm(1000, sd = 200)
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
                    beta = rep(0, J),
                    morpho = rep(0, M),
                    sigma = 100)

templeModel <- nimbleModel(code = templeCode,
                name = "temple",
                constants = templeConsts,
                data = templeData,
                inits = templeInits)

params_to_track <- c("morpho", 
                    "morpho_prob", 
                    "beta", 
                    "sigma", 
                    "mu", 
                    "temple_age")

# mcmc config options

# change default block sampling for beta[] and morpho[] to AF_slice if desired
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

df <- data.frame(date = mcmc_out[, left_out_post_col])
ggplot(data = df) +
    geom_density(mapping = aes(x = date), 
        fill = "steelblue", 
        alpha = 0.75) +
    geom_vline(xintercept = left_out_date) +
    theme_minimal(base_size = 20)

ggsave(filename = "Output/loocv_posterior_example_single.pdf",
        height = 10,
        width = 10,
        units = "cm",
        scale = 2,
        device = "pdf")

df$known_date <- left_out_date

 write.table(df,
    file = "Output/test_loocv_post_pred.csv",
    row.names = F,
    col.names = F,
    sep = ",",
    append = T)

# after re-running the above mannually a few times, bring back in the 
# data from the saved csv file and plot it using the bayesplot package

posterior_samples <- read.csv("Output/test_loocv_post_pred.csv",
                    head = F)

samples_matrix <- matrix(data = posterior_samples[ ,1], 
                    nrow = 45000)

true_dates <- paste(floor(unique(posterior_samples[, 2])))

colnames(samples_matrix) <- true_dates

p <- mcmc_intervals(samples_matrix, 
                    prob = 0.65,
                    prob_outer = 0.95)
p + geom_point(data = true_dates_df,
                mapping = aes(x = true_dates, 
                            y = paste(true_dates))) +
    labs(title = "Posterior Date Estimates", 
        x = "Estimated Years CE",
        y = "True Date") +
    theme_minimal(base_size = 20)

ggsave(filename = "Output/loocv_posterior_example.pdf",
        height = 10,
        width = 15,
        units = "cm",
        device = "pdf")