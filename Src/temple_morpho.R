# libraries
library(ggplot2)
# pull in raw data
temple_known_dates <- read.csv("./data/temple_known_dates.csv", as.is = T)
temple_vars <- read.csv("./data/temple_vars.csv", as.is = T)

temple_vars[grep("(east)", temple_vars$Morphology), "Morphology"] <- "horseshoe_east"
temple_vars[grep("(north)", temple_vars$Morphology), "Morphology"] <- "horseshoe_north"
temple_vars[grep("(west)", temple_vars$Morphology), "Morphology"] <- "horseshoe_west"
temple_vars[grep("4causeway", temple_vars$Morphology), "Morphology"] <- "causeway_4"
temple_vars[grep("2causeway", temple_vars$Morphology), "Morphology"] <- "causeway_2"
temple_vars[grep("Square", temple_vars$Morphology), "Morphology"] <- "square"
temple_vars[which(temple_vars$Morphology == ""), "Morphology"] <- NA

unique(temple_vars$Morphology)

# collate: raw is in two separate tables that can be combined for easier use and then only required variables isolated and combined into a single working dataframe
t_dates <- data.frame(id = temple_known_dates$id,
                    year_ce = temple_known_dates$Date.to.use)
t_morpho <- data.frame(id = temple_vars$Temple.ID,
                    morpho = temple_vars$Morphology)
# merge the two dataframesusing the id column to match while ensuring all of the rows in the vars dataframe are included
temples_morpho <- right_join(t_dates, t_morpho, by = "id")

temples_morpho_known <- subset(temples_morpho, !is.na(year_ce))
p <- ggplot(temples_morpho_known, aes(x = year_ce, fill = morpho)) +
        geom_histogram(alpha = 0.5, position = "identity") +
        theme(axis.text = element_text(size = 20),
                axis.title.x = element_text(size = 40),
                axis.title.y = element_text(size = 40),
                legend.text = element_text(size = 20))
p
unique(t_morpho$morpho)
