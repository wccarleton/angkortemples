# load useful libraries
library(tidyverse)
library(readxl)

# pull in raw data (note some pre-pre-processing has occurred to produce these csv files from the original xlsx workbook)
temple_known_dates <- read.csv("./Data/temple_known_dates.csv", as.is = T)[-c(105:111), ] # unexplained NAs at the bottom from 105:111 caused be errant values in rows 105:111 and cols AN and AO in the original spreadsheet
temple_vars <- read.csv("./Data/temple_vars.csv", as.is = T)

# change the names of the first (temple id) columns of each dataframe so that they are consistent for further processing
names(temple_known_dates)[1] <- "id"
names(temple_vars)[1] <- "id"

# match up the rows based on the "id" column (left join operation) in order to pair known dates with the other variables for each temple
joined_temples <- left_join(t_dates, t_vars, by = "id") # the left join means we are only keeping temples with known dates

# next, select out onl the complete cases (no NAs in any cell---no missing data), and count them
# note that the "complete.cases()" function returns a logical vector indicating whether a given row in joined_temples is "complete"
# in R, these logical vectors can be summed because R internally recasts the logicals as integers where TRUE == 1 and FALSE == 0
sum(complete.cases(joined_temples))

# just in case, check whether the unscripted step in which I saved excel sheets as individual csv files (loaded on ln 2 and 3 above) explains the missing temples
# to do so, we can use the readxl package from tidyverse

# get sheet names
sheets <- excel_sheets("./Data/20180416_Urban_Morphology_Tables1-3.xlsx")
sheets

# pull the dates sheet ("Known Dates")
temple_known_dates <- read_excel("./Data/20180416_Urban_Morphology_Tables1-3.xlsx", sheet = sheets[1])

# new dates from a new SI
sheets_new <- excel_sheets("./Data/si_tables_updated.xlsx")
temple_known_dates_new <- read_excel("./Data/si_tables_updated.xlsx", sheet = sheets_new[1])

# the new spreadsheet distinguishes between model dated and empirically dated temples with text in a notes column.
# the relevant two modelled date notes to exclude are as follows
exclude_1 <- "Date derived through graph-based semi-supervised machine learning (Klassen et al. 2018)."
exclude_2 <- "Date derived through multiple linear regression (Klassen et al. 2018)."

# next, we need to subset the new SI table to extract only the columns with a Klassen ID because those are the ones for which we also have
# predictor variables (the basis for modelled dates)

temple_known_dates_new_Klassen <- filter(temple_known_dates_new, Klassen_Temple_ID != 0)

# isolate the relevant columns and rearrange them to my liking

temple_known_dates_new_Klassen <- temple_known_dates_new_Klassen[, c(4, 1, 13, 14, 21, 22)]

# convert the "dating notes" column to a simple integer column for easier subsetting

date_type_regression_idx <- str_which(temple_known_dates_new_Klassen[, 4][[1]], "regression")
date_type_ssl_idx <- str_which(temple_known_dates_new_Klassen[, 4][[1]], "graph-based")

temple_known_dates_new_Klassen$DateType <- "empirical"
temple_known_dates_new_Klassen[date_type_regression_idx, "DateType"] <- "regression"
temple_known_dates_new_Klassen[date_type_ssl_idx, "DateType"] <- "graphbased"

# rename cols
col_names <- c("id", "name", "date", "dating_notes", "xlong", "ylat", "date_type")
names(temple_known_dates_new_Klassen) <- col_names

# remove duplicated rows (duplicated on the basis of the Klassen ID column, but I noticed at least one case where a duplicate also had a different date)

dups_idx <- duplicated(temple_known_dates_new_Klassen$id)

temple_known_dates_new_Klassen <- temple_known_dates_new_Klassen[!dups_idx, ]

# pull the variables sheet ("Measures")
temple_vars <- read_excel("./Data/20180416_Urban_Morphology_Tables1-3.xlsx", sheet = sheets[2])

# warning message says there's a bad entry: Expecting numeric in F1047 / R1047C6: got '82.827-17'
# the cell reference doesn't account for column headers in the spreadsheet, which means the offending entry would be 
# in row 1046 after import into R

temple_vars[1046, ]

# the column reference is "F", which is the Azimuth data and from above we can see that it's been converted to NA on import
# searching the original spreadsheet brings up the same row and the bad azimuth entry can be seen in F1047 as the warning message
# indicated. The value can be imputed like any other missing data so we'll leave the entry in with NA for the Azimuth

# clean up the morphology variable
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
                    area = log(temple_vars$Area),
                    trait_1 = temple_vars$`Principle Reservoir`,
                    trait_2 = temple_vars$Moat,
                    trait_3 = temple_vars$Sandstone,
                    trait_4 = temple_vars$`Pink Sandstone`,
                    trait_5 = temple_vars$Laterite,
                    trait_6 = temple_vars$Brick,
                    trait_7 = temple_vars$Thmaphnom,
                    trait_8 = temple_vars$other)

# join date and variable tables

temples <- left_join(temple_known_dates_new_Klassen, t_vars, by = "id")

# pull in volumetric data

temple_vol <- read.csv("./Data/qry-data.csv", as.is = T)

t_vol <- data.frame(id = temple_vol$Temple.ID,
                    vol_temple = temple_vol$TV,
                    vol_reservoir = temple_vol$RV,
                    xeast = temple_vol$X,
                    ynorth = temple_vol$Y)

# join again, adding volumes to the main dataframe

temples <- left_join(temples, t_vol, by = "id")

# the principle dataset from here on is now the "temples" dataframe/tibble
# write it as a csv so that I can switch scripts and clear the R workspace of intermediary variables created above while cleaning the data
write.csv(temples, file = "./Data/temples.csv", row.names = F)
