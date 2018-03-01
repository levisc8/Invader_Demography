# Code to clean the germination data for each species that we used it for.
# This alleviates the need to hardcode information into demography scripts

rm(list = ls(all = TRUE))
library(dplyr)
library(lubridate)
library(magrittr)

GermRaw <- read.csv('Uncleaned_Data/Germination/Uncleaned_Germ.csv',
                    stringsAsFactors = FALSE)

DataUsed <- c('Kummerowia', 'Ligustrum', 'Perilla', 'Verbascum')
GermRaw$species[GermRaw$species == 'verbascum'] <- 'Verbascum'

GermCleaned <- GermRaw %>%
  select(species:started) %>%
  mutate(started = gsub('/', '-', started),
         Date_Started = mdy(started)) %>%
  select(-started) %>%
  filter(species %in% DataUsed) %>%
  setNames(c('Species', 'Cell_Number', 'Surface_Germ_Prop',
           'Buried_Germ_Prop', 'Date_Started')) %T>%
  write.csv('Germination/Clean_Germ.csv',
            row.names = FALSE)
  
save(GermCleaned, file = 'Germination/Clean_Germ.RData')

## Comments for usage during testing.
# Vertha viability was calculated as 0.511 in our experiment. pass!
# Kummerowia - pass!
# Perilla - pass!

