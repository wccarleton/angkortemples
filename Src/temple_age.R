# load libraries
library(nimble)
library(coda)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# pull in raw data
temple_known_dates <- read.csv("./Data/temple_known_dates.csv", as.is = T)
temple_vars <- read.csv("./Data/temple_vars.csv", as.is = T)
temple_vol <- read.csv("./Data/qry-data.csv", as.is = T)

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
t_vol <- data.frame(id = temple_vol$Temple.ID,
                    volume = temple_vol$TV)

# merge the two dataframesusing the id column to match while ensuring all of the rows in the vars dataframe are included
temples <- right_join(t_dates, t_vars, by = "id") # right_join because not all temples have dates
temples <- left_join(temples, t_vol, by = "id") # add volume data

# row 451 has a bad azimuth entry--364 degrees--so I'm removing it
temples <- temples[-451, ]

# set up a Nimble model
templeCode <- nimbleCode({
    beta0 ~ dnorm(0, sd = 1000)
    morpho_prob[1:M] ~ ddirch(alpha = d_alpha[1:M])
    for(m in 1:M){
        morpho[m] ~ dnorm(beta0, sd = 200)
    }
    for(j in 1:J){
        beta[j] ~ dnorm(0, sd = 100) # regression coefs
    }
    sigma ~ dunif(0, 1000) # prior variance for regression model
    for(h in 1:H){
        theta[h] ~ dbeta(a[h], b[h]) # prior for hot-encoded vars
    }
    for(n in 1:N){
        x[n, 1] ~ dcat(prob = morpho_prob[1:M])
        x[n, 2] ~ dunif(min = 1, max = 360) # azimuth
        x[n, 3] ~ dlnorm(meanlog = 2, sdlog = 0.2) # area
        for(h in 1:H){
            x[n, 3 + h] ~ dbern(theta[h]) # hot-encoded covariates
        }
        mu[n] <- morpho[x[n, 1]] + inprod(beta[1:J], x[n, 2:(J + 1)])
        temple_age[n] ~ dnorm(mu[n], sd = sigma) # core model
    }
})

# subset only temples with known dates
temples_known <- subset(temples, !is.na(year_ce))

# subset only complete cases
temples_complete <- temples_known[complete.cases(temples_known), ]

temples_onehot_morph <- pivot_wider(temples_complete, 
                                    names_from = morph, 
                                    values_from = morph, 
                                    values_fill = 0, 
                                    values_fn = function(x)as.numeric(!is.na(x)))

temples_idx_morph <- temples_complete
temples_idx_morph$morph <- as.numeric(temples_idx_morph$morph)

M <- length(levels(temples$morph))
H <- 8
Cont <- 2
J <- H + Cont
N <- nrow(temples_idx_morph)

templeConsts <- list(a = rep(1, H), # a,b are the parameters of the Beta priors for the hot-encoded data and these are naive
                    b = rep(1, H),
                    N = N, # number of observations (temples)
                    M = M, # number of potential morpho types (more types than are present in the complete.cases data)
                    H = H, # number of binary predictor variables
                    J = J, # total number of predictor variables
                    d_alpha = rep(1, M)) # parameter vector for Dirichlet prior

# don't include temple volumes as a predictor

tv_index <- grep("volume", colnames(temples_idx_morph))

templeData <- list(temple_age = temples_idx_morph$year_ce,
                    x = temples_idx_morph[, -c(1, 2, tv_index)])

templeInits <- list(theta = rep(0.5, H),
                    beta0 = 0,
                    beta = rep(0, J),
                    morpho = rep(0, M),
                    sigma = 100)

temple <- nimbleModel(code = templeCode,
                name = "temple",
                constants = templeConsts,
                data = templeData,
                inits = templeInits)

params_to_track <- c("beta0", 
                    "morpho", 
                    "morpho_prob", 
                    "beta", 
                    "sigma", 
                    "mu", 
                    "temple_age")

mcmc_out <- nimbleMCMC(model = temple,
                        niter = 40000,
                        nburnin = 5000,
                        summary = T,
                        WAIC = T,
                        monitors = params_to_track)

mcmc_out$summary

# look at MAD for the model
idx_mu <- grep("mu",colnames(mcmc_out$samples))
summary(abs(temples_idx_morph$year_ce - colMeans(mcmc_out$samples[, idx_mu])))

# MAE (AAE, MAD) loss function for cv
MADlossFunction <- function(simulatedDataValues, actualDataValues){
  MAD <- mean(abs(simulatedDataValues - actualDataValues))
  return(MAD)
}

cv_config <- configureMCMC(model = temple)

cv_out <- runCrossValidate(MCMCconfiguration = cv_config,
                            k = nrow(temples_idx_morph),
                            lossFunction = MADlossFunction,
                            MCMCcontrol = list(niter = 30000, nburnin = 3000),
                            nCores = 1,
                            nBootReps = NA,
                            silent = T)

cv_mad <- do.call(rbind,cv_out$foldCVinfo)[, 1]

# look at the distribution of cross-validated (leave-one-out) MAD values:
cv_mad_rnd <- round(summary(cv_mad))
quantile_labels <- row.names(cv_mad_mat)
cv_quantile_ann <- mapply(function(x, y){paste(x, y, sep = ": ")},quantile_labels,cv_mad_mat)
cv_plot_ann <- data.frame(x = rep(200, 6), y = seq(20, 15, -1), label = cv_quantile_ann)
cv_mad_df <- data.frame(fold = 1:length(cv_mad), mad = cv_mad)
plt_cv <- ggplot(data = cv_mad_df) +
        geom_histogram(mapping = aes(x = mad), 
            bins = 6,
            color = "grey") +
        labs(title = "Temple Age Prediction Errors",
        y = "Frequency",
        x = "Mean Absolue Deviation") +
        geom_text(data = cv_plot_ann, aes(x = x, y = y, label = label), hjust = 0, size = 10) +
        theme_minimal(base_size = 20) +
        theme(plot.title = element_text(hjust = 0.5))
plt_cv

ggsave(filename = "Output/cv_mad.pdf", 
        device = "pdf")

traceplot(mcmc(mcmc_out$samples[, "morpho_prob[2]"]))

pairs(mcmc_out$samples[, c(11:18)])

# lm for comparison
lm_temples <- lm(year_ce ~
                -1 +
                azimuth + 
                area + 
                trait_1 + 
                trait_2 + 
                trait_3 + 
                trait_4 + 
                trait_5 + 
                trait_6 + 
                trait_8 +
                horseshoe_east +
                causeway_2 +
                square +
                horseshoe_north +
                causeway_4 +
                blob,
                data = temples_onehot_morph)

summary(lm_temples)

# IQR and median temple volumes over time

temple_age_samples <- mcmc_out$samples[, grep("temple_age", colnames(mcmc_out$samples))]
dim(temple_age_samples)

datum <- 750
delta <- 100

idx_range <- floor( (range(temple_age_samples) - datum) / delta) + 1

# how many temples in each period? This only works simply with known dates. Otherwise
# it gives a rough estimate but each new sample of probable dates will potentially shift
# these numbers around a bit as temples move between bins
table(floor( (temple_age_samples[1, ] - datum) / delta) + 1)


nbins <- length(idx_range[1]:idx_range[2])
p <- c(0.25, 0.5, 0.75) #seq(0, 1, 0.25)
volume_summaries <- array(dim = c(nrow(temple_age_samples), length(p), nbins))

period_quantiles <- function(x, 
                            y, 
                            datum, 
                            delta, 
                            included_indeces, 
                            p, 
                            ...) {
    bin_idx <- floor((x - datum) / delta) + 1
    nbins <- length(included_indeces)
    quantiles <- array(dim = c(1, length(p), nbins))
    for(j in 1:nbins) {
        temples_in_bin <- which(bin_idx == included_indeces[j])
        quantiles[, , j] <- quantile(y[temples_in_bin], 
                                    probs = p, 
                                    ...)
    }
    return(quantiles)
}

for(j in 1:nrow(temple_age_samples)){
    volume_summaries[j, , ] <- period_quantiles(x = temple_age_samples[j, ],
                                            y = temples_idx_morph$volume,
                                            datum = datum,
                                            delta = delta,
                                            included_indeces = idx_range[1]:idx_range[2],
                                            p = p,
                                            na.rm = T)
}

head(volume_summaries)
volume_summaries[1, , ]
diff(volume_summaries[1 , c(1, 3), ])

# now run the mcmc again, but include temples with missing data that need to be predicted

temples_predict <- temples
temples_predict$morph <- as.numeric(temples_predict$morph)

M <- length(levels(temples$morph))
H <- 8
Cont <- 2
J <- H + Cont
N <- nrow(temples_predict)

templeConsts <- list(a = rep(1, H),
                    b = rep(1, H),
                    N = N,
                    M = M,
                    H = H,
                    J = J,
                    d_alpha = rep(1, M))

# don't include temple volumes as a predictor

volume_column <- grep("volume", colnames(temples_predict))

templeData <- list(temple_age = temples_predict$year_ce,
                    x = temples_predict[, -c(1, 2, volume_column)])

# Since much of the data now included are NA (missing values) intended to be imputed during mcmc, it would be burdensome to fully initialize
# the model. So, we won't, but it means there will be warnings about out of bounds indeces (the random index variables) and initialization
# when the mcmc is run. These warnings can be safely ignored.
templeInits <- list(theta = rep(0.5, H),
                    beta0 = 0,
                    beta = rep(0, J),
                    morpho = rep(0, M),
                    sigma = 100)

temple <- nimbleModel(code = templeCode,
                name = "temple",
                constants = templeConsts,
                data = templeData,
                inits = templeInits)

params_to_track <- c("temple_age")

mcmc_out_predict <- nimbleMCMC(model = temple,
                                niter = 40000,
                                nburnin = 5000,
                                summary = T,
                                WAIC = T,
                                monitors = params_to_track)

mcmc_out_predict$summary

# IQR and median temple volumes over time

temple_age_samples <- mcmc_out_predict$samples
temples_with_volumes <- which(!is.na(temples_predict$volume))
temple_age_samples <- temple_age_samples[, temples_with_volumes]

temple_volumes <- log(temples_predict[temples_with_volumes, "volume"])

datum <- 800
delta <- 50

idx_range <- floor( (range(temple_age_samples) - datum) / delta) + 1

bins <- idx_range[1]:idx_range[2]
nbins <- length(bins)
p <- seq(0, 1, 0.25)
volume_summaries <- array(dim = c(nrow(temple_age_samples), length(p), nbins))

period_quantiles <- function(x, 
                            y, 
                            datum, 
                            delta,
                            bins, 
                            p, 
                            ...) {
    bin_idx <- floor((x - datum) / delta) + 1
    nbins <- length(bins)
    quantiles <- array(dim = c(1, length(p), nbins))
    for(j in 1:nbins) {
        temples_in_bin <- which(bin_idx == bins[j])
        quantiles[, , j] <- quantile(y[temples_in_bin], 
                                    probs = p, 
                                    ...)
    }
    return(quantiles)
}

for(j in 1:nrow(temple_age_samples)){
    volume_summaries[j, , ] <- period_quantiles(x = temple_age_samples[j, ],
                                            y = temple_volumes,
                                            datum = datum,
                                            delta = delta,
                                            bins = bins,
                                            p = p,
                                            na.rm = T)
}

# focus on IQR and median
bin_start_year_ce <- (bins * delta) + datum

iqr_samples <- t(apply(volume_summaries, 1, function(x)x[4, ] - x[2, ]))
median_samples <- t(apply(volume_summaries, 1, function(x)x[3, ]))

# pivot longer for easier plotting

period_names <- paste("P", 1:nbins, sep = "")
period_labels <- bin_start_year_ce + (0.5 * delta)

# select focal periods
focal_periods <- which(bin_start_year_ce >= 800 & bin_start_year_ce <= 1475)

iqr_samples_df <- as.data.frame(iqr_samples)
names(iqr_samples_df) <- period_names
iqr_samples_long <- pivot_longer(iqr_samples_df, 
                                cols = everything(), 
                                names_to="period", 
                                values_to = "iqr")

median_samples_df <- as.data.frame(median_samples)
names(median_samples_df) <- period_names
med_samples_long <- pivot_longer(median_samples_df, 
                                cols = everything(), 
                                names_to="period", 
                                values_to = "median")

plt_iqr <- ggplot() +
    #geom_jitter(data = subset(iqr_samples_long, period %in% period_names[focal_periods]), 
    #        mapping = aes(x = period, y = iqr),
    #        width = 0.2,
    #        alpha = 0.05) +
    geom_boxplot(data = subset(iqr_samples_long, period %in% period_names[focal_periods]), 
            mapping = aes(x = period, y = iqr),
            fill = "blue",
            colour = "black",#"#a8a8fd",
            outlier.shape = 1,
            alpha = 0.75) +
    scale_x_discrete(labels = period_labels[focal_periods]) +
    labs(title = "Logged Temple Volume IQR and Median per Period",
        y = "IQR(log(volume))",
        x = "Period Midpoints in Years CE") +
    theme_minimal(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))
plt_iqr

ggsave(filename = "Output/AngkorTemples_volumes_iqr.pdf", 
        device = "pdf")

plt_med <- ggplot() +
    #geom_jitter(data = subset(med_samples_long, period %in% period_names[focal_periods]), 
    #        mapping = aes(x = period, y = median),
    #        width = 0.2, 
    #        alpha = 0.05) +
    geom_boxplot(data = subset(med_samples_long, period %in% period_names[focal_periods]),
            mapping = aes(x = period, y = median),
            fill = "red",
            colour = "black",#"#ffa6a6",
            outlier.shape = 1,
            alpha = 0.75) +
    scale_x_discrete(labels = period_labels[focal_periods]) +
    labs(title = "Logged Temple Volume IQR and Median per Period",
        y = "median(log(volume))",
        x = "Period Midpoints in Years CE") +
    theme_minimal(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))
plt_med

ggsave(filename = "Output/AngkorTemples_volumes_median.pdf", 
        device = "pdf")

# have a look at temple counts per period including dating uncertainties
temple_counts <- array(dim = c(nrow(temple_age_samples), nbins))

period_counts <- function(x, 
                        datum, 
                        delta, 
                        bins) {
    bin_idx <- floor((x - datum) / delta) + 1
    nbins <- length(bins)
    counts <- rep(NA, nbins)
    for(j in 1:nbins) {
        counts[j] <- sum(bin_idx == bins[j])
    }
    return(counts)
}

for(j in 1:nrow(temple_age_samples)){
    temple_counts[j, ] <- period_counts(x = temple_age_samples[j, ],
                                        datum = datum,
                                        delta = delta,
                                        bins = bins)
}

count_samples_df = as.data.frame(temple_counts)
colnames(count_samples_df) <- period_names
count_samples_long <- pivot_longer(count_samples_df, 
                                cols = everything(), 
                                names_to="period", 
                                values_to = "count")

plt_count <- ggplot() +
    #geom_jitter(data = subset(count_samples_long, period %in% period_names[focal_periods]), 
    #        mapping = aes(x = period, y = count),
    #        width = 0.2, 
    #        alpha = 0.05) +
    geom_boxplot(data = subset(count_samples_long, period %in% period_names[focal_periods]), 
            mapping = aes(x = period, y = count),
            fill = "#008100",
            colour = "black",#"green",
            outlier.shape = 1,
            alpha = 0.75) +
    scale_x_discrete(labels = period_labels[focal_periods]) +
    labs(title = "Temple Counts per Period",
        y = "count",
        x = "Period Midpoints in Years CE") +
    theme_minimal(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))
plt_count

ggsave(filename = "Output/AngkorTemples_counts.pdf", 
        device = "pdf")

ggarrange(plt_count, plt_med, plt_iqr,
        ncol = 1,
        align = "v")

ggsave(filename = "Output/AngkorTemples.pdf", 
        device = "pdf",
        height = 60,
        width = 40,
        units = "cm")


# Add in climate data for comparison
# Buckley et al. 2010 PDSI reconstruction based on tree-rings
climate_Buckley2010_PDSI <- read.table("Data/bidoup-nui-ba2010pdsi.txt", 
                                    as.is = T, 
                                    skip = 105, 
                                    header = T)
head(climate_Buckley2010_PDSI)
names(climate_Buckley2010_PDSI) <- c("year_ce", "pdsi")

# Zhang et al. 2008 Wanxiang Cave speleothem oxygen isotope record
climate_Zhang2008_d18O <- read.table("Data/wanxiang2008.txt", 
                                    as.is = T, 
                                    skip = 105, 
                                    header = T)
head(climate_Zhang2008_d18O)
names(climate_Zhang2008_d18O) <- c("distance", "year_ce", "d18O")
#flip y axis on d18O
climate_Zhang2008_d18O$d18O_neg <- -(climate_Zhang2008_d18O$d18O)

plt_climate <- ggplot() +
                geom_path(data = climate_Zhang2008_d18O,
                    mapping = aes(x = year_ce, y = scale(d18O_neg)),
                    color = "#5a5ab4",
                    size = 2,
                    alpha = 0.8) +
                geom_path(data = climate_Buckley2010_PDSI,
                    mapping = aes(x = year_ce, y = scale(pdsi)),
                    color = "red") +
                annotate("rect", 
                        xmin = 1340,
                        ymin = -3, 
                        xmax = 1370,
                        ymax = 3,
                        alpha = 0.5, 
                        fill = "yellow") +
                annotate("rect", 
                        xmin = 1400,
                        ymin = -3, 
                        xmax = 1430,
                        ymax = 3,
                        alpha = 0.5, 
                        fill = "yellow") +
                scale_x_continuous(breaks = seq(825, 1475, 50),
                    limits = c(825, 1475)) +
                theme_minimal(base_size = 20) +
                theme(plot.title = element_text(hjust = 0.5))
plt_climate                        

ggarrange(plt_count, plt_climate,
        ncol = 1,
        align = "v")