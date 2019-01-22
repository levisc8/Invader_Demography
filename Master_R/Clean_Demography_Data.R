# Script to take raw data files with potentially extraneous data 
# and filter it down to necessary information

# rm(list = ls(all = TRUE))
# graphics.off()
library(dplyr)
library(stringr)
library(magrittr)
library(fs)

# Start w/ Ailanthus data. This is sourced from Dryad and will need more 
# cleaning than the others because it contains additional treatments
# that we don't want to give away
AARaw <- read.csv('Uncleaned_Data/Ailanthus/Ailanthus_Raw.csv', 
                  stringsAsFactors = FALSE)
AAGermRaw <- read.csv('Uncleaned_Data/Ailanthus/Ailanthus_Germ.csv',
                      stringsAsFactors = FALSE)

AASDLRaw <- read.csv('Uncleaned_Data/Ailanthus/Ailanthus_Seedlings.csv',
                     stringsAsFactors = FALSE)

# Remove HR and burned plants, select columns and standardize names
AAClean <- AARaw %>%
  select(Site, Plot, Trt, Burn,
         Stage2012, Size, Stage2013, SizeNext, 
         Survival, Quality2013, Seeds, Fecundity) %>%
  setNames(c(
    'Site',
    'Plot',
    'Treatment',
    'Burn',
    'Stage',
    'Size',
    'StageNext',
    'SizeNext',
    'Survival',
    'Quality',
    'Seeds',
    'Reproductive'
  ))

# Turn non-reproductive seed counts into NAs
AAClean$Seeds[AAClean$Reproductive == 0] <- NA_integer_

# Standardize Male/Female
AAClean$Stage[AAClean$Stage == 'Girl'] <- 'Female'
AAClean$StageNext[AAClean$StageNext == 'Girl'] <- 'Female'
AAClean$StageNext[AAClean$StageNext == 'Boy'] <- 'Male'

# Now, germ data
AAGermClean <- AAGermRaw %>%
  filter(Study == 'Germ')  %>%
  mutate(Protected = case_when(
    Protected == 'Yes' ~ 'Y',
    Protected == 'No' ~ 'N'
  )) %>%
  select(-c(ID:Burn, Location)) %>%
  setNames(c(
    'Study',
    'Protected',
    'Total',
    'Germinated',
    "Proportion",
    'Viable',
    'v',
    'g'
  ))
  
# Finally, seedling data
AASdlClean <- AASDLRaw %>%
  mutate(Trt2 = case_when(
    Trt == 'Competitor Removal' ~ 'Comp',
    Trt == 'Herbivore Removal' ~ 'Herb',
    Trt == 'Control' ~ 'Control'
  )) %>%
  select(Burn, Trt2, Trt_Burn, Plot_Value, Year, SDL_no) %>%
  setNames(c(
    'Burn',
    'Treatment',
    'Treatment_X_Burn',
    'Plot',
    'Year',
    'SDL_Count'
  )) %T>%
  write.csv('Ailanthus_IPM/Ailanthus_Seedlings_Clean.csv',
            row.names = FALSE)

save(AASdlClean, file = 'Ailanthus_IPM//Ailanthus_Seedlings_Clean.rda')

write.csv(AAClean, 'Ailanthus_IPM//Ailanthus_Clean.csv',
          row.names = FALSE)

save(AAClean, file = 'Ailanthus_IPM/Ailanthus_Clean.rda')


write.csv(AAGermClean, 'Ailanthus_IPM/Ailanthus_Germ_Clean.csv',
          row.names = FALSE)

save(AAGermClean, file = 'Ailanthus_IPM//Ailanthus_Germ_Clean.rda')


# Alliaria
GMRaw <- read.csv("Uncleaned_Data/Alliaria/GM4R.csv",
                  stringsAsFactors = FALSE) %>%
  filter(Treatment != 'Herb')

GMRaw$Seeds <- NA

GMFec <- read.csv('Uncleaned_Data/Alliaria/GM_Fec_Data.csv',
                  stringsAsFactors = FALSE) %>%
  select(-Plant.Number) %>%
  setNames(c('Plant', 'Height', 'Seeds'))

FecModel <- glm(Seeds ~ Height, 
                data = GMFec,
                family = poisson())

GMRaw$Seeds <- exp(coef(FecModel)[1] + coef(FecModel)[2] * GMRaw$Stem.Ht)


# Calculate number of plots per treatment
Plot.N <- GMRaw %>%
  group_by(Treatment,Density,Burn) %>% 
  summarise(plot.n = length(unique(Plot)))

# calculate seedlings per plot to test for density dependence
sdlPerPlot <- GMRaw %>% 
  group_by(Treatment, Density, Plot) %>%
  summarise(N.SDL = n(),
            N.RA = sum(RA, na.rm = TRUE),
            s2 = N.RA/N.SDL,
            Fec = mean(Seeds,na.rm = TRUE))

# cutting off plots above 200 seedlings/plot
idx <- sdlPerPlot$Plot[sdlPerPlot$N.SDL < 200]

GMClean <- GMRaw %>%
  filter(Plot %in% idx) %>%
  select(-Density) %>%
  setNames(c(
    'Treatment',
    'Burn',
    'Plot',
    'Subquad',
    'Plant',
    'SDLto6NRA',
    'SDLtoFall',
    'Height',
    'Seeds',
    'RA'
  ))%T>% 
  write.csv('Alliaria_MPM/Alliaria_Clean.csv',
            row.names = FALSE) 

GMFec %>%
  write.csv('Alliaria_MPM/Alliaria_Fec_Clean.csv',
            row.names = FALSE)


save(GMClean, file = 'Alliaria_MPM/Alliaria_Clean.rda')
save(GMFec, file = 'Alliaria_MPM/Alliaria_Fec_Clean.rda')  


# Now, Carduus
# rm(list = ls())
CNRaw <- read.csv('Uncleaned_Data/Carduus/tenative.carduus.model.file.csv',
                  stringsAsFactors = FALSE) %>%
  filter(Treatment != 'Herbivore')

CNClean <- CNRaw %>%
  select(-c(Date,
            Tlength,
            Lf10, 
            Therb1014, 
            Swithin, 
            Therb415,
            Seeds,
            Alive615, 
            Notes.1)) %>%
  setNames(c('Treatment', 'Plot', 'Subquad', 
             'Plant', 'Leaf', 
             'Diameter', 'Notes', 'LeafNext', 
             'DiameterNext', 'Survival',
             'Fruit')) 

for(i in seq_len(dim(CNClean)[1])) {
  
  if(CNClean$Survival[i] == 0 |
     is.na(CNClean$Survival[i])) {
    CNClean$LeafNext[i] <- NA_integer_
  }
}

write.csv(CNClean, 'Carduus_MPM/Carduus_Clean.csv',
            row.names = FALSE)
save(CNClean, file = 'Carduus_MPM/Carduus_Clean.rda')


# Draba
# rm(list = ls())
DVRaw <- read.csv('Uncleaned_Data/Draba/Draba4R.csv',
                  stringsAsFactors = FALSE) %>%
  filter(treatment == 'Control' |
           treatment == 'Comp')

DVClean <- DVRaw %>%
  setNames(c('Treatment', 'Plot', 'Plant', 'Survival', 'Fruit')) %T>%
  write.csv('Draba_MPM/Draba_Clean.csv',
            row.names = FALSE)

save(DVClean, file = 'Draba_MPM/Draba_Clean.rda')


# Euonymus - first population data
# rm(list = ls())
EARaw <- read.csv('Uncleaned_Data/Euonymus/Euoala4R.csv', 
                  stringsAsFactors = FALSE) %>%
  filter(Trt != 'Herb')

EAClean <- EARaw %>% 
  select(-c(Quadsite, Cht14, Clen14, 
            circumference.of.stem.at.canopy.base..cm.,
            C.Ht15, Cwid15, Clen15, Fruits15, 
            Therb14)) %>%
  setNames(c('Site', 
             'Plot', 
             'Treatment',
             'Plant',
             'Stage',
             'Height',
             'DBH',
             'Notes',
             'StageNext',
             'HeightNext',
             'DBHNext', 'DBHNext_2',
             'Survival',
             'NotesNext'))

for(i in seq_len(dim(EAClean)[1])) {
  
  if(EAClean$Survival[i] == 0 |
     is.na(EAClean$Survival[i])) {
    EAClean$StageNext[i] <- NA_character_
  }
}
 
write.csv(EAClean, 'Euonymus_IPM/Euonymus_Clean.csv',
            row.names = FALSE)

save(EAClean, file = 'Euonymus_IPM/Euonymus_Clean.rda')

# Reproduction data
EA_RARaw <- read.csv('Uncleaned_Data/Euonymus/RA4R.csv')

EA_RAClean <- EA_RARaw %>%
  select(-(Cht:Cwd)) %>%
  setNames(c(
    'Plant',
    'Height',
    'Fruit'
    )) %T>% 
  write.csv('Euonymus_IPM/Euonymus_RA_Clean.csv', 
            row.names = FALSE)

save(EA_RAClean, file = 'Euonymus_IPM/Euonymus_RA_Clean.rda')


# Kummerowia
# rm(list = ls())
KS14Raw <- read.csv('Uncleaned_Data/Kummerowia/kummerowia AG.csv') %>%
  filter(Treatment != 'Herb')
KS15Raw <- read.csv('Uncleaned_Data/Kummerowia/kstra2015.2 AG.csv') %>%
  filter(Treatment != 'Herb')

KS15Clean <- KS15Raw %>%
  select(Treatment, Plot, Plant, Alive, Notes) %>%
  setNames(c('Treatment', 'Plot', 'Plant', 
             'Survival','Notes')) %T>%
  write.csv('Kummerowia_MPM/Kummerowia_15_Clean.csv',
            row.names = FALSE)

KS14Clean <- KS14Raw %>%
  select(-X) %>%
  setNames(c('Treatment', 'Plot', 'Plant', 'Survival', 'Seeds')) %T>%
  write.csv('Kummerowia_MPM/Kummerowia_14_Clean.csv',
            row.names = FALSE)

save(KS15Clean, file = 'Kummerowia_MPM/Kummerowia_15_Clean.rda')
save(KS14Clean, file = 'Kummerowia_MPM/Kummerowia_14_Clean.rda')


# Lepidium
# rm(list = ls())
LCRaw <- read.csv('Uncleaned_Data/Lepidium/lep4R.csv') %>%
  filter(treatment != 'Herb')

LCClean <- LCRaw %>%
  select(-Seeds) %>%
  setNames(c('Treatment', 'Plot', 'Plant',
             'Survival', 'Fruit')) %T>% 
  write.csv('Lepidium_MPM/Lepidium_Clean.csv',
            row.names = FALSE)

save(LCClean, file = 'Lepidium_MPM/Lepidium_Clean.rda')


# Lespedeza
# rm(list = ls())

LCRaw <- read.csv('Uncleaned_Data/Lespedeza/Lespedeza4R.csv', 
                  stringsAsFactors = FALSE) %>%
  filter(Trt_Burn2013 == 'Competition-Unburned' &
           !is.na(Stage2012) |
           Trt_Burn2013 == 'Control-Unburned' &
           !is.na(Stage2012))

LCClean <- LCRaw %>%
  select(c(Trt2013, Plot, Quadrat, Plant, 
           Stage2012, Size2012, Stage2013,
           Size2013, Survival2013, Resprout2013,
           Fec2013, Seeds2013)) %>%
  setNames(c('Treatment', 'Block', 'Plot', 'Plant',
             'Stage', 'Size', 'StageNext',
             'SizeNext', 'Survival', 'Resprout',
             'Reproductive', 'Seeds'))

for(i in seq_len(dim(LCClean)[1])) {
  
  if(LCClean$Reproductive[i] == 0 |
     is.na(LCClean$Reproductive[i])) {
    LCClean$Seeds[i] <- NA_real_
  }
}


write.csv(LCClean, 'Lespedeza_MPM/Lespedeza_Clean.csv', 
            row.names = FALSE)

save(LCClean, file = 'Lespedeza_MPM/Lespedeza_Clean.rda')


# Ligustrum
# rm(list = ls())

LORaw <- read.csv('Uncleaned_Data/Ligustrum/Ligobt4R-8.8.15.csv') %>%
  filter(Trt != 'Herb') %>%
  filter(Clone == 0 |  is.na(Clone))

LOClean <- LORaw %>%
  select(c(Trt, Quadrat, Plant, Stg14,
           Ht14, Sum.DBH14,
           Notes, Stg15, Ht15,
           Sum.DBH15, Alive2015, 
           Fruits2015)) %>%
  setNames(c('Treatment', 'Plot', 'Plant', 'Stage',
             'Height', 'DBH', 'Notes', 'StageNext',
             'HeightNext', 'DBHNext', 'Survival', 'Seeds'))

for(i in seq_len(dim(LOClean)[1])) {
  
  if(LOClean$Survival[i] == 0 |
     is.na(LOClean$Survival[i])) {
    LOClean$StageNext[i] <- NA_character_
  }
} 

write.csv(LOClean, 'Ligustrum_IPM/Ligustrum_Clean.csv',
          row.names = FALSE)

save(LOClean, file = 'Ligustrum_IPM/Ligustrum_Clean.rda')

LO_RARaw <- read.csv('Uncleaned_Data/Ligustrum/RA4R.csv')

LO_RAClean <- LO_RARaw %>%
  setNames(c('Plant', 'Height', 'Seeds')) %T>%
  write.csv('Ligustrum_IPM/Ligustrum_RA_Clean.csv',
            row.names = FALSE)

save(LO_RAClean, file = 'Ligustrum_IPM/Ligustrum_RA_Clean.rda')

# Lonicera
# rm(list = ls())
LMRaw <- read.csv('Uncleaned_Data/Lonicera/lonmaa1213.csv',
                  stringsAsFactors = FALSE) %>%
  filter(Trt_Burn2013 == 'ContN' |
           Trt_Burn2013 == 'CompN' |
           Trt_Burn2013 == 'AllN')
  
LMClean <- LMRaw %>%
  select(c(Site,
           Plot,
           Trt_Burn2013,
           Size2012, 
           Hght.Jul.2013,
           Stage2012,
           StageJuly2013,
           Fruits2013, 
           Survival2013,
           Resprout2013,
           Fec2013)) %>%
  setNames(c('Site', 'Plot', 'Treatment', 
             'Height', 'HeightNext', 'Stage',
             'StageNext', 'Fruit', 'Survival', 
             'Resprout', 'Reproductive')) 

LMClean$Reproductive[is.na(LMClean$StageNext)] <- NA_integer_
  
write.csv(LMClean, 'Lonicera_IPM/Lonicera_Clean.csv',
            row.names = FALSE)

save(LMClean, file = 'Lonicera_IPM/Lonicera_Clean.rda')


# Perilla
# rm(list = ls())
PFRaw <- read.csv('Uncleaned_Data/Perilla/Perilla4R2013 AG.csv') %>%
  filter(Treatment != 'Herb')

PFClean <- PFRaw %>%
  select(c(Treatment:Height1, Alive, Fruits, Seeds, Comments)) %>%
  setNames(c('Treatment', 'Plot', 'Subquad', 
             'Height', 'Survival', 'Fruit',
             'Seeds', 'Notes')) 

for(i in seq_len(dim(PFClean)[1])) {
  if(PFClean$Survival[i] == 0) {
    PFClean$Fruit[i] <- PFClean$Seeds[i] <- NA_real_
  }
  
}

write.csv(PFClean, 'Perilla_MPM/Perilla_Clean.csv',
          row.names = FALSE)

save(PFClean, file = 'Perilla_MPM/Perilla_Clean.rda')

# Potentilla
# rm(list = ls())
PRRaw <- read.csv('Uncleaned_Data/Potentilla/potrec20154R AG.csv',
                  stringsAsFactors = FALSE) %>%
  filter(Treatment != 'Herb')

PR_SDLRaw <- read.csv("Uncleaned_Data/Potentilla/Potentilla Seedlings 2015 AG.csv") %>%
  filter(Treatment != 'Herb')

PRCleaned <- PRRaw %>%
  select(Treatment, Plot, Plant, Lf14,
         LL14, Stg14, Fruits14, Notes,
         Stg15, Lf15, LL15,
         Fruits15, Alive15, Notes15) %>%
  setNames(c('Treatment', 'Plot', 'Plant',
             'Leaf', 'LongestLeaf', 'Stage',
             'Fruit', 'Notes', 'StageNext', 'LeafNext',
             'LongestLeafNext', 'FruitNext', 'Survival', 'NotesNext')) %T>%
  write.csv('Potentilla_MPM/Potentilla_Clean.csv',
            row.names = FALSE)

save(PRCleaned, file = 'Potentilla_MPM/Potentilla_Clean.rda')

PR_SDLClean <- PR_SDLRaw %>%
  select(-Subquad) %>%
  setNames(c('Plot', 'Treatment', 'Survival')) %T>%
  write.csv('Potentilla_MPM/Potentilla_SDL_Clean.csv',
            row.names = FALSE)

save(PR_SDLClean, file = 'Potentilla_MPM/Potentilla_SDL_Clean.rda')

# Thlaspi
# rm(list = ls())
TPRaw <- read.csv('Uncleaned_Data/Thlaspi/ThlaspiD.csv',
                  stringsAsFactors = FALSE) %>%
  filter(Treatment == 'Control' | 
           Treatment == 'Comp')

TPClean <- TPRaw %>%
  select(Treatment:Fruit) %>%
  setNames(c('Treatment', 'Plot', 'Plant', 'Survival', 'Fruit')) %T>%
  write.csv('Thlaspi_MPM/Thlaspi_Clean.csv',
            row.names = FALSE)

save(TPClean, file = 'Thlaspi_MPM/Thlaspi_Clean.rda')

# Verbascum
# rm(list = ls())
VTRaw <- read.csv('Uncleaned_Data/Verbascum/Herbdemo4R.csv',
                  stringsAsFactors = FALSE)

VTClean <- VTRaw %>%
  select(c(treatment, plot, plant, diameter, notes,
           Alive15, Seeds, Notes15)) %>%
  setNames(c('Treatment', 'Plot', 'Plant', 'Diameter',
             'Notes', 'Survival', 'Seeds', 'NotesNext')) %T>%
  write.csv('Verbascum_MPM/Verbascum_Clean.csv', 
            row.names = FALSE)

save(VTClean, file = 'Verbascum_MPM/Verbascum_Clean.rda')
