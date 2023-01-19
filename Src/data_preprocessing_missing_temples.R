# load libraries
library(tidyverse)
library(readxl)

# data wrangling

# pull in the data from Excel sheets and CSVs as needed

# get sheet names
sheets <- excel_sheets("./Data/20180416_Urban_Morphology_Tables1-3.xlsx")
sheets

# pull the variables sheet ("Measures")
temple_vars <- read_excel("./Data/20180416_Urban_Morphology_Tables1-3.xlsx", 
                    sheet = sheets[2])

# new dates from a new SI
sheets <- excel_sheets("./Data/si_tables_updated.xlsx")
temple_known_dates <- read_excel("./Data/si_tables_updated.xlsx", 
                    sheet = sheets[1])

# next, we need to subset the new SI table to extract only the columns with a 
# Klassen ID because those are the ones for which we also have
# predictor variables (the basis for modelled dates)
temple_known_dates <- filter(temple_known_dates, Klassen_Temple_ID != 0)

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

# remove duplicated rows (duplicated on the basis of the Klassen ID column, 
# but I noticed at least one case where a duplicate also had a different date),
# so this has to be revisted...

dups_idx <- duplicated(temple_known_dates$id)

temple_known_dates <- temple_known_dates[!dups_idx, ]

# warning message says there's a bad entry: Expecting numeric in F1047 / 
# R1047C6: got '82.827-17' the cell reference doesn't account for column 
# headers in the spreadsheet, which means the offending entry would be 
# in row 1046 after import into R

temple_vars[1046, ]

# the column reference is "F", which is the Azimuth data and from above we can 
# see that it's been converted to NA on import searching the original 
# spreadsheet brings up the same row and the bad azimuth entry can be seen in 
# F1047 as the warning message indicated. The value can be imputed like any 
# other missing data so we'll leave the entry in with NA for the Azimuth

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

temples <- left_join(temple_known_dates, t_vars, by = "id")

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
# write it as a csv so that I can switch scripts and clear the R workspace 
# of intermediary variables created above while cleaning the data
write.csv(temples, file = "./Data/temples.csv", row.names = F)
