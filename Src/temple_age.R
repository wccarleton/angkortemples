# load libraries
library(nimble)
library(coda)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(GGally)
library(leaflet)
library(sp)
library(maps)
library(readxl)
library(progress)

# data wrangling

# pull in the data from Excel sheets and CSVs as needed

# get sheet names
sheets <- excel_sheets("./Data/20180416_Urban_Morphology_Tables1-3.xlsx")
sheets

# pull the variables sheet ("Measures")
temple_vars <- read_excel("./Data/20180416_Urban_Morphology_Tables1-3.xlsx", 
                    sheet = sheets[2])

# warning message says there's a bad entry: Expecting numeric in F1047 / 
# R1047C6: got '82.827-17' the cell reference doesn't account for column 
# headers in the spreadsheet, which means the offending entry would be 
# in row 1046 after import into R

temple_vars[1046, ]

# The value can be imputed like any other missing data so we'll leave the entry 
# in with NA for the Azimuth

# new dates from a new SI
sheets <- excel_sheets("./Data/si_tables_updated.xlsx")
temple_known_dates <- read_excel("./Data/si_tables_updated.xlsx", 
                    sheet = sheets[1])

# next, we need to subset the new SI table to extract only the columns with a 
# Klassen ID because those are the ones for which we also have
# predictor variables (the basis for modelled dates)
temple_known_dates <- subset(temple_known_dates, !is.na(Klassen_Temple_ID))

# and again to get rid of the cases where Klassen_Temple_ID == 0
temple_known_dates <- subset(temple_known_dates, Klassen_Temple_ID != 0)

# isolate the relevant columns and rearrange them to my liking
temple_known_dates <- temple_known_dates[, c(4, 1, 13, 14, 21, 22)]

# from the "dating notes" column create a simpler column that indicates the 
# source of the date
date_type_regression_idx <- str_which(temple_known_dates[, 4][[1]], "regression")
date_type_ssl_idx <- str_which(temple_known_dates[, 4][[1]], "graph-based")

temple_known_dates$DateType <- "empirical"
temple_known_dates[date_type_regression_idx, "DateType"] <- "regression"
temple_known_dates[date_type_ssl_idx, "DateType"] <- "graphbased"

# rename cols
col_names <- c("id", "name", "date", "dating_notes", "xlong", "ylat", "date_type")
names(temple_known_dates) <- col_names

# at this point we have 1454 temples....

# remove duplicated rows (duplicated on the basis of the Klassen ID column, 
# but I noticed at least one case where a duplicate also had a different date),
# so this has to be revisted...

dups_idx <- duplicated(temple_known_dates$id)

# 154 duplicates identified. The way R handles this is to identify duplicates
# based on order of appearance from top to bottom and, so, if we drop the
# duplicates, we are just getting rid of any rows with duplicated id's after a 
# given entry---this is a fairly naive approach of course, but without 
# additional information about the duplicates, it's the only way to proceed. 

temple_known_dates <- temple_known_dates[!dups_idx, ]

# clean up the morphology column and simplify the names, removing special characters
temple_vars[grep("(east)", temple_vars$Morphology), "Morphology"] <- "horseshoe_east"
temple_vars[grep("(north)", temple_vars$Morphology), "Morphology"] <- "horseshoe_north"
temple_vars[grep("(west)", temple_vars$Morphology), "Morphology"] <- "horseshoe_west"
temple_vars[grep("4causeway", temple_vars$Morphology), "Morphology"] <- "causeway_4"
temple_vars[grep("2causeway", temple_vars$Morphology), "Morphology"] <- "causeway_2"
temple_vars[grep("Square", temple_vars$Morphology), "Morphology"] <- "square"
temple_vars[which(temple_vars$Morphology == ""), "Morphology"] <- NA

# isolate relevant cols and rename
t_vars <- data.frame(id = temple_vars$`Temple ID`,
                    morph = as.factor(temple_vars$Morphology),
                    azimuth = temple_vars$Azimuth,
                    area = temple_vars$Area,
                    trait_1 = temple_vars$`Principle Reservoir`,
                    trait_2 = temple_vars$Moat,
                    trait_3 = temple_vars$Sandstone,
                    trait_4 = temple_vars$`Pink Sandstone`,
                    trait_5 = temple_vars$Laterite,
                    trait_6 = temple_vars$Brick,
                    trait_7 = temple_vars$Thmaphnom,
                    trait_8 = temple_vars$other)

# join date and variable tables

temples <- left_join(t_vars, temple_known_dates, by = "id")

# now that the temples data bases have been joined, we have many entries with
# NA in the 'date_type' column. This is fine because these are undated temples

# pull in coords

temple_xy <- read.csv("./Data/qry-data.csv", as.is = T)

t_xy <- data.frame(id = temple_xy$Temple.ID,
                    xeast = temple_xy$X,
                    ynorth = temple_xy$Y)

# join again, adding volumes to the main dataframe

temples <- left_join(temples, t_xy, by = "id")

# NOTE there is a bad azimuth entry in row 1239 (> 360), so here we set it to NA

temples[1239, "azimuth"] <- NA

# analysis

# create a vector containg the idx of covariate columns to include in models

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

# the principle dataset from here on is now the "temples" dataframe/tibble
# write it as a csv so that I can switch scripts and clear the R workspace 
# of intermediary variables created above while cleaning the data
write.csv(temples, file = "./Data/temples.csv", row.names = F)

# subset only dated temples in order to compare the two approaches with 
# cross-validation and mean absolute deviation

temples_dated <- temples[which(!is.na(temples$date_emp)), ]

# the models work with numeric data, so all tibble columns have 
# to be converted. at the same time we can just select the relevant covariates
# (i.e., excluding the coordinate columns, which aren't needed until much later)

x <- mutate(temples_dated[, covariate_idx], 
            across(morph:trait_8, as.numeric))

# it also can't be a tibble, so...

x <- as.data.frame(x)

date_emp <- temples_dated$date_emp

temples_dated <- cbind(date_emp, x)

## First run the GSSL

# some additional data wrangling required here to prepare the data for 
# the GSSL analysis/algorithm, so copy the data here to avoid overwriting

temples_dated_gssl <- temples_dated

# next, identify the categorical and continuous data columns because GSSL
# treats these types differently

cont_idx <- grep("azimuth|area", colnames(temples_dated_gssl))
cat_idx <- setdiff(1:ncol(temples_dated_gssl), cont_idx)[-1]

# GSSL process doesn't treat NA values differently from any other potential
# value, so just convert them to a numerical value
na_idx <- which(is.na(temples_dated_gssl[, "morph"]))
temples_dated_gssl[na_idx, "morph"] <- max(temples_dated_gssl[, "morph"], na.rm = T) + 1
categories <- as.matrix(temples_dated_gssl[, cat_idx])
categories[which(is.na(categories))] <- 2
temples_dated_gssl[, cat_idx] <- categories

# rescale azimuth and area to (0,1)
max_azimuth <- max(temples_dated_gssl[, "azimuth"], na.rm = T)
temples_dated_gssl[, "azimuth"] <- temples_dated_gssl[, "azimuth"] / max_azimuth
max_area <- max(temples_dated_gssl[, "area"], na.rm = T)
temples_dated_gssl[, "area"] <- temples_dated_gssl[, "area"] / max_area

# provide values for missing data in the continuous variables as per
# the PLoS paper methods and python code
na_idx <- which(is.na(temples_dated_gssl[, "azimuth"]))
temples_dated_gssl[na_idx, "azimuth"] <- 0.5
na_idx <- which(is.na(temples_dated_gssl[, "area"]))
temples_dated_gssl[na_idx, "area"] <- 0.5

# propagate dates/labels by least-squares

# define the gssl function
propagate_labels <- function(x, cat_idx, cont_idx){
    # I've decided to continue to refer to labels in the code to be consistent
    # with the python code and make it easier for co-authors to follow

    # get the row indeces for dated (labelled) and undated (unlabelled) temples
    labelled_idx <- which(!is.na(x[, 1]))
    unlabelled_idx <- which(is.na(x[, 1]))

    unlabelled <- as.matrix(x[unlabelled_idx, ])
    labelled <- as.matrix(x[labelled_idx, ])

    # set precision manually with a variable, as per the original python code
    eps <- 0.2

    # calculate weights
    # unlabelled x unlabelled similarities

    n <- dim(unlabelled)[1]

    w_uu <- matrix(nrow = n, ncol = n)

    for(j in 1:n){
        v <- matrix(rep(unlabelled[j, ], n), nrow = n, byrow = T)
        d <- unlabelled[, cat_idx] - v[, cat_idx]
        # subtracting here means that any pair of elements with the same value
        # (even NA, which has been set to 2 or 9 above depending on the variable)
        # will yield a 0, so we can then just count zeros...
        # rowSums only works with length(dim) > 1, so check here
        if(is.matrix(d)){
            sim <- rowSums(d == 0)
            a <- (unlabelled[, cont_idx] - v[, cont_idx])**2
            w_uu[, j] <- sim + ( 2 - rowSums(a) )
        }else{
            sim <- sum(d == 0)
            a <- (unlabelled[, cont_idx] - v[, cont_idx])**2
            w_uu[, j] <- sim + ( 2 - sum(a) )
        }
    }

    w_uu <- exp(w_uu / eps)

    # unlabelled x labelled similarities
    m <- dim(labelled)[1]

    w_ul <- matrix(nrow = n, ncol = m)

    for(j in 1:m){
        v <- matrix(rep(labelled[j, ], n), nrow = n, byrow = T)
        d <- unlabelled[, cat_idx] - v[, cat_idx]
        if(is.matrix(d)){
            sim <- rowSums(d == 0)
            a <- (unlabelled[, cont_idx] - v[, cont_idx])**2
            w_ul[, j] <- sim + ( 2 - rowSums(a) )
        }else{
            sim <- sum(d == 0)
            a <- (unlabelled[, cont_idx] - v[, cont_idx])**2
            w_ul[, j] <- sim + ( 2 - sum(a) )
        }
    }

    w_ul <- exp(w_ul / eps)

    # sum the similarity measure per temple to produce a 'total similarity'
    # this will be the diagonal of the graph lapacian and replaces the node
    # order values typical of a simple unweighted graph lapacian matrix

    # propagate labels finally

    label_vector <- labelled[, 1]

    # again the rowSums thing...
    if(length(w_uu) > 1){
        d_u <- rowSums(w_uu) + rowSums(w_ul)
        lsqr_label_predictions <- solve(diag(d_u) - w_uu, w_ul %*% label_vector)
    }else{
        d_u <- w_uu + sum(w_ul)
        lsqr_label_predictions <- solve(d_u - w_uu, w_ul %*% label_vector)
    }

    guesses = round(lsqr_label_predictions, 0)

    # return a list containing the computed arrays

    return(list(labels = guesses,
                w_uu = w_uu, 
                w_ul = w_ul,
                d_u = d_u))
}

# leave-one-out cross validation on labelled set (dated temples). This will
# be comparable to the Bayesian model

# isolate the labelled temple data
n_dated_temples <- dim(temples_dated_gssl)[1]

# create a vector for containing MAD values
predicted_dates <- rep(NA, n_dated_temples)

pb <- progress_bar$new(total = n_dated_temples)

for(k in 1:n_dated_temples){
    pb$tick()
    # save a temple date, then set it to NA
    left_out_date <- temples_dated_gssl[k, 1]
    temples_dated_gssl[k, 1] <- NA
    gssl_dates <- propagate_labels(x = temples_dated_gssl,
                                    cat_idx = cat_idx,
                                    cont_idx = cont_idx)
    predicted_dates[k] <- gssl_dates$labels
    temples_dated_gssl[k, 1] <- left_out_date
}

# calculate absolute deviations for each left-out temple
cv_ad_gssl <- abs(temples_dated_gssl$date_emp - predicted_dates)

# save the deviations for later use
write.table(cv_ad_gssl, 
        file = "Output/cv_abs_devs_gssl.csv",
        row.names = F,
        col.names = F)

## Now run the Nimble/Baysian predictive model

# The data are formatted slightly differently for the Bayesian model. Nimble
# makes use of NA values for imputation/prediction, so unlike the GSSL data
# wrangling steps above, we don't need to change NAs to anything and we 
# won't impute any values at this point. Instead, we will just pass the 
# x dataframe from line 164 and the empirical dates column from the 
# temples_dated dataframe to the nimble functions.

# set up a Nimble model
templeCode <- nimbleCode({
    morpho_prob[1:M] ~ ddirch(alpha = d_alpha[1:M])
    for(m in 1:M){
        morpho[m] ~ dnorm(0, sd = 1000)
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

# save some useful parameters to pass to the nimble model, like number of
# observations etc.
M <- length(levels(temples$morph)) # number of possible temple morpho types
H <- 8 # number of binary variables
Cont <- 2 # number of continuous variables
J <- H + Cont # total number of variables (covariates/predictors)
N <- nrow(temples_dated) # N obs.

templeConsts <- list(a = rep(1, H), # beta prior
                    b = rep(1, H), # beta prior
                    N = N,
                    M = M,
                    H = H,
                    J = J,
                    d_alpha = rep(1, M)) # parameter vector for Dirichlet prior

templeData <- list(temple_age = temples_dated$date_emp,
                    x = temples_dated[, -1])

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
#temple_mcmc_config$removeSamplers(c("beta0","beta","morpho"))
#temple_mcmc_config$addSampler(target = c("beta0", "beta", "morpho"), type = "AF_slice")
temple_mcmc_config$setMonitors(params_to_track)

# build mcmc
temple_mcmc <- buildMCMC(temple_mcmc_config)
temple_mcmc_c <- compileNimble(temple_mcmc)

# run mcmc

mcmc_out <- runMCMC(temple_mcmc_c, 
                    niter = 50000,
                    nburnin = 0)

plot(mcmc_out[,"morpho[4]"], type="l")

# convergence checking
nonparam_cols <- grep("mu|temple_age", colnames(mcmc_out))
g <- geweke.diag(mcmc(mcmc_out[, -nonparam_cols]))

# warn if convergence check fails at 2-sigma, two-tailed
if(any(abs(g$z) > 2)){
    param_names <- colnames(mcmc_out[, -nonparam_cols])
    m1 <- "One or more chains may not have converged."
    m2 <- paste("Check:", paste(param_names[which(abs(g$z) > 2)], collapse = ", "))
    warning(paste(m1, m2, sep = "\n"))
}

# run leave-one-out cv by repeated calls to system2() that 
# then passes a command to bash, launches a new R instance,
# runs mcmc, calculates the absolute deviation for the given left-out sample,
# then writes the absolute deviation to a file. This is done
# rather than using nimble::runCrossValidate() because the latter
# balloons in memory until it gets killed---I didn't have time
# to troubleshoot or figure out if there's a problem to fix
# in the nimble code and it was just easier to use this 
# more awkward, less computationally efficient approach.

# NOTE: the script below (Src/get_ad.R) takes the mean of the posterior 
# predicted temple date and compares that to the empirical date to derive 
# the absolute deviation.

for(j in 1:dim(temples_dated)[1]){
    system2("Rscript", args=c("Src/get_ad.R", paste(j)))
}

cv_ad_bayes <- read.csv("Output/cv_abs_devs_2.csv", head = F)[, 1]
cv_ad_bayes <- data.frame(deviation_years = cv_ad_bayes, model = "bayesian")

cv_ad_gssl <- read.csv("Output/cv_abs_devs_gssl.csv", head = F)[, 1]
cv_ad_gssl <- data.frame(deviation_years = cv_ad_gssl, model = "gssl")

# create a long-format dataframe with the cross-validation-derived absolute
# deviations from each of the two approaches.

ad_both <- rbind(cv_ad_bayes, cv_ad_gssl)
ad_both[, 1] <- round(as.numeric(ad_both[, 1]), 0)
names(ad_both) <- c("deviation_years", "model")

# plot the two for comparison
plt_cv_compared <- ggplot(data = ad_both) +
    geom_histogram(aes(x = deviation_years,
                    y = stat(count / sum(count)), 
                    fill = model)) +
    facet_wrap(~model, nrow = 2) +
    labs(x = "Absolute Deviation\nin Years", 
        y = "Normalized Count",
        title = "Cross Validation Results Compared") +
    theme_minimal(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

plt_cv_compared

ggsave(filename = "Output/cv_ad_compared.pdf", 
        device = "pdf")

# have a look at temple counts per period 
# in one plot, show the series using the GSSL labels, and in the other show the
# series from the Bayesian model including dating uncertainties

x <- mutate(temples[, covariate_idx],
            across(morph:trait_8, as.numeric))
x <- as.data.frame(x)

date_emp <- temples$date_emp

# GSSL
temples_gssl <- cbind(date_emp, x)

# next, identify the categorical and continuous data columns because GSSL
# treats these types differently

cont_idx <- grep("azimuth|area", colnames(temples_gssl))
cat_idx <- setdiff(1:ncol(temples_gssl), cont_idx)[-1]

# GSSL process doesn't treat NA values differently from any other potential
# value, so just convert them to a numerical value
na_idx <- which(is.na(temples_gssl[, "morph"]))
temples_gssl[na_idx, "morph"] <- max(temples_gssl[, "morph"], na.rm = T) + 1
categories <- as.matrix(temples_gssl[, cat_idx])
categories[which(is.na(categories))] <- 2
temples_gssl[, cat_idx] <- categories

# rescale azimuth and area to (0,1)
max_azimuth <- max(temples_gssl[, "azimuth"], na.rm = T)
temples_gssl[, "azimuth"] <- temples_gssl[, "azimuth"] / max_azimuth
max_area <- max(temples_gssl[, "area"], na.rm = T)
temples_gssl[, "area"] <- temples_gssl[, "area"] / max_area

# provide values for missing data in the continuous variables as per
# the PLoS paper methods and python code
na_idx <- which(is.na(temples_gssl[, "azimuth"]))
temples_gssl[na_idx, "azimuth"] <- 0.5
na_idx <- which(is.na(temples_gssl[, "area"]))
temples_gssl[na_idx, "area"] <- 0.5


# now, run the gssl and extract the temple foundation date estimates

gssl_dates <- propagate_labels(x = temples_gssl,
                                cat_idx = cat_idx,
                                cont_idx = cont_idx)

# Bayesian model
# run the model again including all the data, not just the dated temples that 
# were used above to evaluate the model's performance. we only need to change
# the variables related to the new dataset. note that this run will take
# a fair while longer to complete because n is going from 169 to 1300, which
# greatly increases the number of parameters to sample since we're also
# imputing the NAs in each case.

x <- mutate(temples[, covariate_idx], 
            across(morph:trait_8, as.numeric))

x <- as.data.frame(x)

date_emp <- temples$date_emp

temples_predict <- cbind(date_emp, x)

N <- nrow(temples) # N obs.

templeConsts <- list(a = rep(1, H), # beta prior
                    b = rep(1, H), # beta prior
                    N = N,
                    M = M,
                    H = H,
                    J = J,
                    d_alpha = rep(1, M)) # parameter vector for Dirichlet prior

templeData <- list(temple_age = temples_predict$date_emp,
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

# change default block sampling for beta[] and morpho[] to AF_slice
templeModel_c <- compileNimble(templeModel)
temple_mcmc_config <- configureMCMC(templeModel_c)
temple_mcmc_config$setMonitors(params_to_track)

# build mcmc
temple_mcmc <- buildMCMC(temple_mcmc_config)
temple_mcmc_c <- compileNimble(temple_mcmc)

# run mcmc

mcmc_out_predict <- runMCMC(temple_mcmc_c, 
                    niter = 50000)

# Next, take the MCMC samples for the predicted temple foundation dates and 
# bin them (count temple foundations per period) for each MCMC iteration.

temple_age_idx <- grep("temple_age", colnames(mcmc_out))
temple_age_samples <- mcmc_out[, temple_age_idx]

# create a container for the counts---each row will refer to one temporal bin,
# while each column will contain one probable count sequence of temple
# foundation events
temple_counts <- array(dim = c(nrow(temple_age_samples), nbins))

# write a function for counting temple foundations given the age samples, a
# starting time (datum), time bin width (delta), and number of desired bins.
# note that the bins will be defined by their lower date edge.
period_counts <- function(x, 
                        datum, 
                        delta, 
                        nbins) {
    bin_idx <- floor((x - datum) / delta) + 1
    counts <- rep(NA, nbins)
    for(j in 1:nbins) {
        counts[j] <- sum(bin_idx == j)
    }
    return(counts)
}

# plot results of both models

# choose the parameters for the temporal binning
datum <- 700
delta <- 100
nbins <- 7

# loop over the rows (MCMC iterations) of the temple_age_samples matrix and 
# then store the counts
for(j in 1:nrow(temple_age_samples)){
    temple_counts[j, ] <- period_counts(x = temple_age_samples[j, ],
                                        datum = datum,
                                        delta = delta,
                                        nbins = nbins)
}

temple_counts_gssl <- period_counts(x = as.vector(gssl_dates$labels),
                                    datum = datum,
                                    delta = delta,
                                    nbins = nbins)

# for easier plotting, we can make it a dataframe and use ggplot2 tools
periods <- 1:nbins

period_labs <- paste(seq(from = datum + (delta / 2), 
                        by = delta, 
                        length.out = nbins))

count_samples_df = as.data.frame(temple_counts)
colnames(count_samples_df) <- periods
count_samples_long <- pivot_longer(count_samples_df, 
                                cols = everything(), 
                                names_to = "period", 
                                values_to = "bayes_count")

count_samples_long$gssl_count <- rep(temple_counts_gssl, dim(temple_counts)[1])

plt_count <- ggplot() +
    geom_boxplot(data = count_samples_long, 
            mapping = aes(x = period, y = bayes_count),
            fill = "#008100",
            colour = "black",
            outlier.shape = 1,
            alpha = 0.75) +
    geom_boxplot(data = count_samples_long,
            mapping = aes(x = period, y = gssl_count),
            colour = "#a5eca5",
            alpha = 0.75,
            size = 2) +
    scale_x_discrete(labels = period_labs) +
    labs(title = "Temple Counts per Period",
        y = "count",
        x = "Period Midpoints in Years CE") +
    theme_minimal(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))
plt_count

ggsave(filename = "Output/AngkorTemples_counts.pdf", 
        device = "pdf")


# Variable selection

# There are two types of "predictor" variables in the model: 1) morpho types, 
# which are modelled with index variables because they are a set of exclusive 
# categories [each temple must be one and only one morpho type]; and 2) a set of 
# continuous variables and binary variables that are not mutually exclusive 
# like the morpho types. To examine the impact of the morpho types (mutally
# exclusive categories), we can examine the pairwise differences between
# the postorior samples. In the model, the morpho type parameter acts as a 
# intercept would and temples of each category are able to have different
# intercepts. The actual values are in calendar years, so the posterior for 
# a morpho type parameter (e.g., morpho[1], morpho[2],...) represents the
# average foundation year for temples of that type. By comparing them by way
# of contrast distributions we can determine whether any of the morpho
# types were on average older/younger than one or more of the others.

# First, plot all of them as densities on the same timeline

# morpho type columns
morpho_idx <- grep("morpho\\[", colnames(mcmc_out))

posterior_contrasts <- function(x, y, ...){
    usr <- par("usr")
    on.exit(par("usr" = usr))
    diffs <- x - y
    h <- hist(diffs, plot = FALSE)
    #d <- density(diffs)
    xcor <- max(abs(diffs))
    ycor <- max(h$counts)
    #ycor <- range(d$y)
    par(usr = c(-xcor, xcor, 0, ycor))
    rect(-xcor, 0, 0, ycor, col = "#0000ff70")
    rect(0, 0, xcor, ycor, col = "#018a016c")
    #lines(x = d$x, y = d$y)
    hist(diffs, add = TRUE, border = "white", col = "white")
    abline(v = 0)
}

gg_posterior_contrasts <- function(data, mapping, ...){
    x <- rlang::eval_tidy(mapping$x, data)
    y <- rlang::eval_tidy(mapping$y, data)
    d <- data.frame(diffs = x - y)
    xcor <- max(abs(d))
    #xcor <- 10
    h <- hist(d$diffs, plot = F)
    bw <- mean(diff(h$mids))
    ggplot() +
        geom_area(data = NULL, 
            mapping = aes(x = c(-xcor, 0), y = 1),
            fill = "steelblue",
            alpha = 0.75) +
        geom_area(data = NULL, 
            mapping = aes(x = c(0, xcor), y = 1),
            fill = "darkgreen",
            alpha = 0.75) +
        stat_bin(data = d, 
            mapping = aes(x = diffs, y = stat(count / max(count))),
            geom = "step",
            color = "black",
            #fill = "transparent",
            binwidth = bw)
}

#hist(apply(mcmc_out[, morpho_idx[c(1, 6)]], 1,diff))

morpho_type_df <- as.data.frame(mcmc_out[, morpho_idx])

plt_contrasts <- ggpairs(morpho_type_df,
                        lower = list(
                                continuous = gg_posterior_contrasts),
                        upper = "blank",
                        diag = "blank") +
                theme_minimal(base_size = 12) +
                theme(plot.title = element_text(hjust = 0.5))

plt_contrasts

pairs(mcmc_out[, morpho_idx],
    lower.panel = posterior_contrasts,
    upper.panel = NULL)

# plot all densities on single timeline
morpho_type_df_long <- pivot_longer(morpho_type_df,
                                    cols = everything(),
                                    names_to = "morpho_type",
                                    values_to = "date")

morpho_labels_map <- levels(temples$morph)
names(morpho_labels_map) <- paste("morpho[",1:8,"]",sep = "")
morpho_labeller <- labeller(morpho_type = morpho_labels_map)

ggplot(data = subset(morpho_type_df_long, morpho_type != "morpho[7]")) + 
    geom_density(mapping = aes(x = date, fill = morpho_type),
                alpha = 0.75,
                position = "identity") +
    xlim(c(600, 1400)) +
    geom_vline(xintercept = c(800, 1400)) +
    facet_wrap(~ morpho_type, 
        ncol = 1,
        labeller = morpho_labeller) +
    labs(title = "Temple Morpho Type Marginal Date Posteriors",
        x = "Date (Years CE)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = "Output/AngkorTemples_morpho_types.pdf", 
        device = "pdf")

# next, calculate pairwise posterior contrasts and plot those densities---any
# that differ significantly from zero would indicate an average foundation
# date difference between the relevant morpho types.

# Next, we can look at the impact of the other predictor variables on the 
# estimated temple ages.

# There are multiple ways this could ultimately be done, but Nimble provides
# a handy tool called "reversible jump" MCMC whereby essentially the variable
# selection is happening throughout the simulation rather than having to
# manually configure seperate models leaving out different variables
# and then comparing the "fit" manually (e.g., with WAIC). Instead, we can 
# use indicator variables to turn each variable on and off at random 
# while exploring the model likelihood and then at the end examine the posteriors
# for the indicator variables to determine which variables were included
# most frequently overall. To use Nimbles RJMCMC we need to reconfigure the
# MCMC from above.

# reset data to include only dated temples again...

x <- mutate(temples_dated[, covariate_idx],
            across(morph:trait_8, as.numeric))
x <- as.data.frame(x)

N <- nrow(temples_dated) # N obs.

templeConsts <- list(a = rep(1, H), # beta prior
                    b = rep(1, H), # beta prior
                    N = N,
                    M = M,
                    H = H,
                    J = J,
                    d_alpha = rep(1, M)) # parameter vector for Dirichlet prior

templeData <- list(temple_age = temples_dated$date_emp,
                    x = x)

templeInits <- list(theta = rep(0.5, H),
                    beta0 = 1000,
                    sigma0 = 200,
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
                    "sigma")

templeModel_c <- compileNimble(templeModel)
temple_mcmc_config <- configureMCMC(templeModel_c)
temple_mcmc_config$setMonitors(params_to_track)

configureRJ(temple_mcmc_config,
            targetNodes = 'beta',
            priorProb = 0.5,
            control = list(mean = 0, scale = 0.2))

# build mcmc
temple_mcmc <- buildMCMC(temple_mcmc_config)
temple_mcmc_c <- compileNimble(temple_mcmc)

# run mcmc
mcmc_out_var <- runMCMC(temple_mcmc_c, 
                        niter = 50000,
                        nburnin = 5000)

# extract inclusion probabilities
betas <- grep("beta\\[", colnames(mcmc_out_var))
post_inclusion_probs <- colMeans(apply(mcmc_out_var[, betas], 
                                2, 
                                function(x) x != 0))

# make a nice table for export
model_names <- colnames(temples)[covariate_idx][-1]

variable_names <- c("azimuth",
                    "area",
                    "principle reservoir",
                    "moat",
                    "sandstone",
                    "pink sandstone",
                    "laterite",
                    "brick",
                    "thmaphom",
                    "other")

inclusion_prob_df <- data.frame(variable = model_names,
                                description = variable_names,
                                include_prob = round(post_inclusion_probs, 2))

write.csv(inclusion_prob_df, 
            file = "Output/variable_inclusion.csv")

# handy function for excluding zeros
rm_zeros <- function(x){
    zeros <- which(x == 0)
    if(length(zeros) != 0L){
        return(x[-zeros])
    }else{
        return(x)
    }
}

# mapping temple foundation events with chronological uncertainty

predicted_dates_idx <- grep("temple_age", colnames(mcmc_out_predict))
predicted_dates <- mcmc_out_predict[, predicted_dates_idx]

# assumes normally distributed dating uncertainty

opacity <- function(x, at){
    sigma <- sd(x)
    if(sigma == 0){
        return(rep(1, length(at)))
    } else {
        mu <- mean(x)
        d <- dnorm(x = at, 
            mean = mu, 
            sd = sigma)
        return(d / max(d))
    }
}

sample_years <- 800:1400

# create an "opacity matrix"---i.e., get opacties for all temples
opacity_matrix <- apply(predicted_dates, 2, opacity, at = sample_years)

# not all temples have coordinates, so we need to subset the dataframe and the opacity matrix in order to plot
# the following logical vector can be used to select from the dataframe (by rows) and opacity matrix (by columns)

has_coords <- complete.cases(temples[, c("xlong", "ylat")])
xys <- temples[has_coords, c("xlong", "ylat")]

sp_temples <- SpatialPointsDataFrame(coords = xys, data = temples[has_coords, ], 
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# a shiny app

library(shiny)

ui <- fluidPage(
    sliderInput("obstime", 
                "Year CE",
                min = sample_years[1], 
                max = max(sample_years), 
                value = median(sample_years),
                step = 1,
                animate = animationOptions(interval = 100)),
    leafletOutput("AngkorTemples"),
    p(),
)

server <- function(input, output, session) {

    temporal_datum <- reactive({
        sample_years[input$obstime - sample_years[1] + 1]
    })
        
    current_opacity <- reactive({
        as.vector(opacity_matrix[input$obstime - sample_years[1] + 1, has_coords])
    })

    relative_temporal_indicator <- reactive({
        as.vector(temple_age_means >= temporal_datum())
    })

    colour_idx <- reactive({
        relative_temporal_indicator() + 1
    })

    current_colour <- reactive({
        sapply(colour_idx(), function(x)c("blue", "red")[x])
    })

    output$AngkorTemples <- renderLeaflet({
        m <- leaflet()
        m <- addTiles(m)
        m <- addCircles(m, 
                data = sp_temples)
        m
    })

    observe({
        p <- leafletProxy("AngkorTemples")
        p <- clearShapes(p)
        p <- addCircles(p, 
                data = sp_temples, 
                opacity = current_opacity(), 
                color = current_colour())
        p
    })
}

shinyApp(ui, server)
