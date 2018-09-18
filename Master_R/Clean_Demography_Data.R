# Script to take raw data files with potentially extraneous data 
# and filter it down to necessary information

rm(list = ls(all = TRUE))
graphics.off()
library(dplyr)
library(stringr)
library(magrittr)

# Ailanthus is already on Dryad, so whoever wants it can get it.
# Starting with Alliaria
GMRaw <- read.csv("Uncleaned_Data/Alliaria/GM4R.csv",
                  stringsAsFactors = FALSE) %>%
  filter(Treatment != 'Herb')

GMRaw$Seeds <- NA

GMFec <- read.csv('Uncleaned_Data/Alliaria/GM_Fec_Data.csv',
                  stringsAsFactors = FALSE) %>%
  select(-Plant.Number) %>%
  setNames(c('Plant_Number', 'Stem_Height', 'Seeds'))

FecModel <- glm(Seeds ~ Stem_Height, 
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
  select(-Density) %T>% 
  write.csv('Alliaria/GM_Clean.csv',
            row.names = FALSE) 

GMFec %>%
  write.csv('Alliaria/GM_Fec_Clean.csv',
            row.names = FALSE)


save(GMClean, file = 'Alliaria/GM_Clean.RData')
save(GMFec, file = 'Alliaria/GM_Fec_Clean.RData')  

# Now, Carduus
rm(list = ls())
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
             'Plant_Number', 'Leaf_N_May14', 
             'Plant_Diameter14', 'Notes', 'Leaf_N_Apr15', 
             'Plant_Diameter15', 'Survival',
             'Flower_Number')) 

for(i in seq_len(dim(CNClean)[1])) {
  
  if(CNClean$Survival[i] == 0 |
     is.na(CNClean$Survival[i])) {
    CNClean$Leaf_N_Apr15[i] <- NA_integer_
  }
}

write.csv(CNClean, 'Carduus/CN_Clean.csv',
            row.names = FALSE)
save(CNClean, file = 'Carduus/CN_Clean.RData')

# Draba
rm(list = ls())
DVRaw <- read.csv('Uncleaned_Data/Draba/Draba4R.csv',
                  stringsAsFactors = FALSE) %>%
  filter(treatment == 'Control' |
           treatment == 'Comp')

DVClean <- DVRaw %>%
  setNames(c('Treatment', 'Plot', 'Plant', 'Survival', 'Fruits')) %T>%
  write.csv('Draba/DV_Clean.csv',
            row.names = FALSE)

save(DVClean, file = 'Draba/DV_Clean.RData')

# Euonymus - first population data
rm(list = ls())
EARaw <- read.csv('Uncleaned_Data/Euonymus/Euoala4R.csv', 
                  stringsAsFactors = FALSE) %>%
  filter(Trt != 'Herb')

EAClean <- EARaw %>% 
  select(-c(Quadsite, Cht14, Clen14, 
            circumference.of.stem.at.canopy.base..cm.,
            C.Ht15, Cwid15, Clen15, Fruits15, 
            Therb14)) %>%
  setNames(c('Site', 'Plot', 'Treatment',
             'Plant_Number', 'Stage14', 'Plant_Height14',
             'DBH14', 'Notes_2014', 'Stage15', 'Plant_Height15',
             'DBH15', 'DBH15_2', 'Survival', 'Notes_2015'))

for(i in seq_len(dim(EAClean)[1])) {
  
  if(EAClean$Survival[i] == 0 |
     is.na(EAClean$Survival[i])) {
    EAClean$Stage15[i] <- NA_character_
  }
}
 
write.csv(EAClean, 'Euonymus/IPM/Data/EA_Clean.csv',
            row.names = FALSE)

save(EAClean, file = 'Euonymus/IPM/Data/EA_Clean.RData')

# Reproduction data
EA_RARaw <- read.csv('Uncleaned_Data/Euonymus/RA4R.csv')

EA_RAClean <- EA_RARaw %>%
  select(-(Cht:Cwd)) %T>% 
  write.csv('Euonymus/IPM/Data/EA_RA_Clean.csv', 
            row.names = FALSE)
save(EA_RAClean, file = 'Euonymus/IPM/Data/EA_RA_Clean.RData')

# Kummerowia
rm(list = ls())
KS14Raw <- read.csv('Uncleaned_Data/Kummerowia/kummerowia AG.csv') %>%
  filter(Treatment != 'Herb')
KS15Raw <- read.csv('Uncleaned_Data/Kummerowia/kstra2015.2 AG.csv') %>%
  filter(Treatment != 'Herb')

KS15Clean <- KS15Raw %>%
  select(Treatment, Plot, Plant, Alive, Notes) %>%
  setNames(c('Treatment', 'Plot', 'Plant_Number', 
             'Survival','Notes')) %T>%
  write.csv('Kummerowia/KS_15_Clean.csv',
            row.names = FALSE)

KS14Clean <- KS14Raw %>%
  select(-X) %>%
  setNames(c('Treatment', 'Plot', 'Plant_Number', 'Survival', 'Seeds')) %T>%
  write.csv('Kummerowia/KS_14_Clean.csv',
            row.names = FALSE)

save(KS15Clean, file = 'Kummerowia/KS_15_Clean.RData')
save(KS14Clean, file = 'Kummerowia/KS_14_Clean.RData')

# Lepidium
rm(list = ls())
LCRaw <- read.csv('Uncleaned_Data/Lepidium/lep4R.csv') %>%
  filter(treatment != 'Herb')

LCClean <- LCRaw %>%
  select(-Seeds) %>%
  setNames(c('Treatment', 'Plot', 'Plant_Number',
             'Survival', 'Fruits')) %T>% 
  write.csv('Lepidium/LepCam_Clean.csv',
            row.names = FALSE)

save(LCClean, file = 'Lepidium/LepCam_Clean.RData')


# Lespedeza
rm(list = ls())

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
  setNames(c('Treatment', 'Block', 'Plot', 'Plant_Number',
             'Stage12', 'Size12', 'Stage13',
             'Size13', 'Survival', 'Resprout',
             'Reproductive', 'Seeds'))

for(i in seq_len(dim(LCClean)[1])) {
  
  if(LCClean$Reproductive[i] == 0 |
     is.na(LCClean$Reproductive[i])) {
    LCClean$Seeds[i] <- NA_real_
  }
}


write.csv(LCClean, 'Lespedeza/LesCun_Clean.csv', 
            row.names = FALSE)

save(LCClean, file = 'Lespedeza/LesCun_Clean.RData')
  
# Ligustrum
rm(list = ls())

LORaw <- read.csv('Uncleaned_Data/Ligustrum/Ligobt4R-8.8.15.csv') %>%
  filter(Trt != 'Herb') %>%
  filter(Clone == 0 |  is.na(Clone))

LOClean <- LORaw %>%
  select(c(Trt, Quadrat, Plant, Stg14,
           Ht14, Sum.DBH14,
           Notes, Stg15, Ht15,
           Sum.DBH15, Alive2015, 
           Fruits2015)) %>%
  setNames(c('Treatment', 'Plot', 'Plant_Number', 'Stage14',
             'Plant_Height14', 'DBH14', 'Notes_14', 'Stage15',
             'Plant_Height15', 'DBH15', 'Survival', 'Seeds'))

for(i in seq_len(dim(LOClean)[1])) {
  
  if(LOClean$Survival[i] == 0 |
     is.na(LOClean$Survival[i])) {
    LOClean$Stage15[i] <- NA_character_
  }
} 

write.csv(LOClean, 'Ligustrum/IPM/Data/LO_Clean.csv',
          row.names = FALSE)

save(LOClean, file = 'Ligustrum/IPM/Data/LO_Clean.RData')

LO_RARaw <- read.csv('Uncleaned_Data/Ligustrum/RA4R.csv')

LO_RAClean <- LO_RARaw %>%
  setNames(c('Plant', 'Height', 'Seeds')) %T>%
  write.csv('Ligustrum/IPM/Data/LO_RA_Clean.csv',
            row.names = FALSE)

save(LO_RAClean, file = 'Ligustrum/IPM/Data/LO_RA_Clean.RData')

# Lonicera
rm(list = ls())
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
             'Plant_Height12', 'Plant_Height13', 'Stage12',
             'Stage13', 'Fruit', 'Survival', 
             'Resprout', 'Reproductive')) 

LMClean$Reproductive[is.na(LMClean$Stage13)] <- NA_integer_
  
write.csv(LMClean, 'Lonicera/IPM/Data/LM_Clean.csv',
            row.names = FALSE)

save(LMClean, file = 'Lonicera/IPM/Data/LM_Clean.RData')


# Perilla
rm(list = ls())
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

write.csv(PFClean, 'Perilla/PF_Clean.csv',
          row.names = FALSE)

save(PFClean, file = 'Perilla/PF_Clean.RData')

# Potentilla
rm(list = ls())
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
  setNames(c('Treatment', 'Plot', 'Plant_Number',
             'Leaf_Number14', 'Longest_Leaf14', 'Stage14',
             'Fruit14', 'Notes14', 'Stage15', 'Leaf_Number15',
             'Longest_Leaf15', 'Fruit15', 'Survival', 'Notes15')) %T>%
  write.csv('Potentilla/PR_Clean.csv',
            row.names = FALSE)

save(PRCleaned, file = 'Potentilla/PR_Clean.RData')

PR_SDLClean <- PR_SDLRaw %>%
  select(-Subquad) %>%
  setNames(c('Plot', 'Treatment', 'Survival')) %T>%
  write.csv('Potentilla/PR_SDL_Clean.csv',
            row.names = FALSE)

save(PR_SDLClean, file = 'Potentilla/PR_SDL_Clean.RData')


# Thlaspi
rm(list = ls())
TPRaw <- read.csv('Uncleaned_Data/Thlaspi/ThlaspiD.csv',
                  stringsAsFactors = FALSE) %>%
  filter(Treatment == 'Control' | 
           Treatment == 'Comp')

TPClean <- TPRaw %>%
  select(Treatment:Fruit) %>%
  setNames(c('Treatment', 'Plot', 'Plant', 'Survival', 'Fruit')) %T>%
  write.csv('Thlaspi/TP_Clean.csv',
            row.names = FALSE)

save(TPClean, file = 'Thlaspi/TP_Clean.RData')

# Verbascum
rm(list = ls())
VTRaw <- read.csv('Uncleaned_Data/Verbascum/Herbdemo4R.csv',
                  stringsAsFactors = FALSE)

VTClean <- VTRaw %>%
  select(c(treatment, plot, plant, diameter, notes,
           Alive15, Seeds, Notes15)) %>%
  setNames(c('Treatment', 'Plot', 'Plant_Number', 'Diameter',
             'Notes14', 'Survival', 'Seeds', 'Notes15')) %T>%
  write.csv('Verbascum/VT_Clean.csv', 
            row.names = FALSE)

save(VTClean, file = 'Verbascum/VT_Clean.RData')
