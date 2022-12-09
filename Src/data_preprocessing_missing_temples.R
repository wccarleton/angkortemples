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

# pull the variables sheet ("Measures")
temple_vars <- read_excel("./Data/20180416_Urban_Morphology_Tables1-3.xlsx", sheet = sheets[2])

# and then follow the same steps as before...
names(temple_known_dates)[1] <- "id"
names(temple_vars)[1] <- "id"
joined_temples <- left_join(t_dates, t_vars, by = "id")
sum(complete.cases(joined_temples))
