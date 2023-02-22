library(progress)

# data prep
# pull in CSV
temples <- read.csv("Data/temples.csv")

# add an emprical date column to make further processing easier. any temples
# with modelled dates will have an NA in this column.

temples$emprical_date <- NA
dated_temples <- which(temples$date_type == "empirical")
temples[dated_temples, "emprical_date"] <- temples[dated_temples, "date"]

# copy temples to a new variable b/c the gssl has different formatting and data 
# imputation requirements than the bayesian nimble model

temples_gssl <- temples

# rescale azimuth and area to (0,1)
temples_gssl$azimuth <- temples_gssl$azimuth / max(temples_gssl$azimuth, na.rm = T)
temples_gssl$area <- temples_gssl$area / max(temples$area, na.rm = T)

# select temple data colums and then arrange for conformity with the 
# label propagation function expectations

str_pattern <- paste("emprical_date",
                    "morph",
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

col_idx <- grep(str_pattern, names(temples))

# annoyingly the return from grep is sorted ascending instead of as indicated
# by the str_pattern order, so need to rearrange (or refactor propogate_dates)

ncols <- length(col_idx)
rearranged_col_idx <- c(col_idx[ncols], col_idx[-ncols])

# next, identify the categorical and continuous data columns
cont_idx <- grep("azimuth|area", names(temples[, rearranged_col_idx]))
new_col_idx <- 1:ncol(temples[, rearranged_col_idx])
cat_idx <- setdiff(new_col_idx, 
                    new_col_idx[cont_idx])[-1]

# convert every column to a numeric type and give NA a numeric value that
# depends on the given variable type. This gets around using a large number 
# of logical comparisons, which would like be slower

x <- temples[, rearranged_col_idx]

# NA in temple morphology can be any integer not corresponding to a temple
# type. Start by recasting the morphology column as a factor, which will
# apply an integer level number to every type label represented in the column.
# Then, set NAs to something else (not one of the integer levels)---easy 
# to just add one to the max level and use that number
x$morph <- as.numeric(factor(x$morph))


x$morph[which(is.na(x$morph))] <- max(x$morph, na.rm = T) + 1
y <- as.matrix(x[, cat_idx])
y[which(is.na(y))] <- 2
x[, cat_idx] <- y

# provide values for missing data in the continuous variables as per
# the PLoS paper methods and python code
x$azimuth[which(is.na(x$azimuth))] <- 0.5
x$area[which(is.na(x$area))] <- 0.5

# propagate dates/labels by least-squares

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

gssl_dates <- propagate_labels(x = x,
                            cat_idx = cat_idx,
                            cont_idx = cont_idx)

# leave-one-out cross validation on labelled set (dated temples). This will
# be comparable to the Bayesian model

# isolate the labelled temple data
labelled_temples <- x[which(!is.na(x$emprical_date)), ]

# create a vector for containing MAD values
predicted_dates <- rep(NA, dim(labelled_temples)[1])

n_loo <- dim(labelled_temples)[1]

pb <- progress_bar$new(total = n_loo)

for(k in 1:n_loo){
    pb$tick()
    # save a temple date, then set it to NA
    left_out_date <- labelled_temples[k, 1]
    labelled_temples[k, 1] <- NA
    gssl_dates <- propagate_labels(x = labelled_temples,
                                    cat_idx = cat_idx,
                                    cont_idx = cont_idx)
    predicted_dates[k] <- gssl_dates$labels
    labelled_temples[k, 1] <- left_out_date
}

abs_devs <- abs(labelled_temples$emprical_date - predicted_dates)

write.table(abs_devs, 
        file = "Output/cv_abs_devs_gssl.csv",
        row.names = F,
        col.names = F)