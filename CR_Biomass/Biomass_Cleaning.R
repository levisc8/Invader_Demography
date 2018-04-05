rm(list = ls())
library(dplyr)
library(ggplot2)
library(lubridate)

# get vector of natives and species we abandoned (Torilis + Cirsium)
natives <- c('Geum', 'Symphoricarpos', 'Torilis', 'Desmodium', 'Cirsium',
             'Teucrium', 'Symcar', 'Torillis')

# Read in raw biomass data and do a little cleaning. This will 
# allow us to aggregate by species rather than species x site

# Run the code prior to first %>% operator to see what raw
# biomass data looks like
Biomass <- read.csv('CR_Biomass/Raw/Biomass4R.csv', 
                    stringsAsFactors = FALSE) %>% 
  arrange(Site, Plot, Tag) %>% 
  mutate(Site2 = lapply(Site, FUN = function(x) stringr::str_split(x, " ")[[1]][1]) %>% 
           unlist() %>%
           tolower()) %>%
  filter(!Site2 %in% tolower(natives))

# Correct a stupid typo in raw data
Biomass$Site2[Biomass$Site2 == 'enonymus'] <- 'euonymus'
# get dates in order using lubridate package
Biomass$Date <- gsub('/','-',Biomass$Date) %>% mdy()

# Create single column for ID'ing plants
Biomass$Tag[Biomass$Tag == ""] <- NA
for(i in 1:dim(Biomass)[1]){
  Biomass$ID[i] <- ifelse(is.na(Biomass$Tag[i]),
                          paste0(Biomass$Plot[i],' - P'),
                          paste0(Biomass$Tag[i],' - T'))
}

save(Biomass, file = 'CR_Biomass/Cleaned/TRC_Comp_Removal_Biomass.rda')
write.csv(Biomass, 'CR_Biomass/Cleaned/TRC_Comp_Removal_Biomass.csv',
          row.names = FALSE)

# Create a table of plot sizes so we can standardize total biomass removed
# by plot area

Areas <- tibble(Species = unique(Biomass$Site2),
                BoundarySize = c(1, # Ailanthus
                                 0.5, # Alliaria
                                 1,  # Carduus
                                 0.5, # draba
                                 1, # Euonymus
                                 0.25, # Kummerowia
                                 1, # Lepidium
                                 0.5, # Lespedeza
                                 1, # Lonicera
                                 0.5, # Thlaspi
                                 0.5, # Perilla
                                 0.5, # Potentilla
                                 1, # Ligustrum
                                 1)) %>% # Verbascum
  mutate(Area = BoundarySize^2)

Biomass <- Biomass %>%
  left_join(., Areas, by = c('Site2' = 'Species')) %>%
  mutate(B_m_2 = Mbt/Area)

# Now, we are going to rearrange some things for a couple weird species.
# First three are our winter annuals. We have data from November
# (when they germinate, but when
# everything else dies) and March/April (when they're about to flower, but everything
# else germinates). Using the March/April data as that is probably more indicative
# of the amount of competition they actually experience over the winter
# Kummerowia was mowed in its first year, so data from that year aren't really
# representative. We collected data in 2015 (2 years after mowing), so using those
# instead.

Winter_annuals <- filter(Biomass, Site2 %in% c('draba', 
                                               'microthlaspi',
                                               'lepidium') &
                           month(Date) == 4 |
                           Site2 == 'kummerowia' & year(Date) == 2015)


Biomass <- filter(Biomass, !Site2 %in% c('draba', 
                                         'microthlaspi',
                                         'lepidium', 
                                         'kummerowia')) %>%
  rbind(Winter_annuals) %>% 
  arrange(Site2, ID)

# Really long pipe chain to calculate means per species at plot level.
Means <- Biomass %>% 
  group_by(Site2, ID, Date) %>% 
  summarise(N = n(),
            Tot = sum(B_m_2, na.rm = TRUE)) %>%
  arrange(Site2, Date) %>%
  mutate(Rep = rank(Date)) %>%
  filter(Rep == 1) %>% 
  group_by(Site2) %>%
 summarise(FirstCut = mean(Tot))

# In subsequent analyses, we want to know about how strong competition is.
# however, smaller species are likely to be more affected by the same 
# amount of biomass than larger species. Thus, we standardize our biomass
# measurements by the size of the focal species in question. Our estimates
# are derived from measuring ~20 individuals/species and calculating the average
# mass of those. For our larger species, we have set a value of 2 kg as a place
# holder. However, I am investigating using allometric equations relating
# DBH and wood density to calculate better values for these data summaries.
# Of course, you can skip this step and just use the raw data too
Standards <- read.csv('CR_Biomass/Raw/BM_Standards.csv',
                  stringsAsFactors = FALSE) %>% 
  select(Species, Standard) %>%
  filter(!Species %in% natives | Species != 'CB') %>% 
  group_by(Species) %>% 
  summarise(Standard = first(Standard))

for(i in unique(Means$Site2)){
  Means$StdBM[Means$Site2 == i] <- Means$FirstCut[Means$Site2 == i] /
                                   Standards$Standard[tolower(Standards$Species) == i]
}

Means$LogCRBM <- log(Means$StdBM)

# Add in full scientific names and change column names to be more
# descriptive
SpecNames <- c(ailanthus = 'Ailanthus_altissima',
               alliaria = 'Alliaria_petiolata',
               carduus = 'Carduus_nutans',
               draba = 'Draba_verna',
               euonymus = 'Euonymus_alatus',
               kummerowia = 'Kummerowia_striata',
               lepidium = 'Lepidium_campestre',
               lespedeza = 'Lespedeza_cuneata',
               privet = 'Ligustrum_obtusifolium',
               lonicera = 'Lonicera_maackii',
               perilla = 'Perilla_frutescens',
               potentilla = 'Potentilla_recta',
               microthlaspi = 'Thlaspi_perfoliatum',
               verbascum = 'Verbascum_thapsus')

Means <- Means %>%
  mutate(Species = SpecNames[order(names(SpecNames))]) %>%
  select(Species, FirstCut:LogCRBM) %>%
  setNames(c('Species', 
             'Unstandardized_Biomass', 
             'Standardized_Biomass', 
             'Log_Standardized_Biomass')) %>%
  arrange(Species)

# Save outputs
save(Means, file = 'CR_Biomass/Cleaned/TRC_Biomass_Means.rda')
write.csv(Means, 'CR_Biomass/Cleaned/TRC_Biomass_means.csv',
          row.names = FALSE)

