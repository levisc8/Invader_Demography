# Lespedeza Matrix Model for manuscript -------------
# Hopefully this is a bit more rigorous
# than my last pass at this
# Begin: 12/4/17

# Quick notes: I can only use 2012-2013 transition to parameterize
# values for seedlings. We apparently have no data
# seedlings in 2013 :(

rm(list = ls())
graphics.off()

library(dplyr)


PopData <- read.csv('Lespedeza/LesCun_Clean.csv', 
                    stringsAsFactors = FALSE) 

# Delineating stage classes for adults at 250, 251-1000, 1000 > cm.
# This provides adequate sample sizes to capture demographic differences
# between size classes. It's not an IPM, but our data is a little too
# quirky for me to feel comfortable parameterizing one with it.

PopData$StageMatMod <- NA

# Reclassify stages based on logic from above ----------
for(i in 1:dim(PopData)[1]) {
  if(PopData$Stage12[i] == 'SDL') {
    PopData$StageMatMod[i] <- 'SDL'
  } else if(PopData$Stage12[i] == 'SDL2') {
    PopData$StageMatMod[i] <- 'SDL2'
  } else if(PopData$Stage12[i] == 'Adult') {
    if(PopData$Size12[i] > 8 &&
       PopData$Size12[i] <= 250) {
      PopData$StageMatMod[i] <- 'Small'
    } else if(PopData$Size12[i] > 250 &&
              PopData$Size12[i] < 1000) {
      PopData$StageMatMod[i] <- 'Med'
    } else {
      PopData$StageMatMod[i] <- 'Large'
    }
  }
}

# Add stageNext so we know how to calculate transition probabilities --------

PopData$StageMatModNext <- NA
for(i in 1:dim(PopData)[1]) { 
  if(PopData$Survival[i] == 0 |
     is.na(PopData$Survival[i])) { 
    # if it died, stagenext is NA
    PopData$StageMatModNext[i] <- NA
    
  } else if(PopData$Survival[i] == 1 &
            is.na(PopData$Stage13[i]) |
            is.na(PopData$Size13[i])) { 
    if(is.na(PopData$Size13[i])){
    # if it didn't die, but we don't have any other info, stageNext is same Stage
      PopData$StageMatModNext[i] <- PopData$StageMatMod[i]
    } else {
      PopData$StageMatModNext[i] <- PopData$Stage13[i]
    }
    # The rest is pretty self explanatory
  } else if(PopData$Stage13[i] == 'SDL2') {
    PopData$StageMatModNext[i] <- 'SDL2'
  } else if(PopData$Stage13[i] == 'Adult') {
    if(PopData$Size13[i] <= 250) {
      PopData$StageMatModNext[i] <- 'Small'
    } else if(PopData$Size13[i] > 250 &&
              PopData$Size13[i] < 1000) {
      PopData$StageMatModNext[i] <- 'Med'
    } else {
      PopData$StageMatModNext[i] <- 'Large'
    }
  }
}

# Select relevant columns and then filter out the seedling counts
# by quadrat.
MatData <- PopData %>%
  select(Plant_Number, Plot, Quadrat, 
         Treatment, StageMatMod, StageMatModNext, 
         Survival, Reproductive, Seeds) %>%
  filter(!is.na(Plant_Number)) %>%
  setNames(c('Plant_Number', 'Plot', 'Quad', 'Trt',
             'Stage','StageNext',
             'Survival', 'Reproductive', 'Seeds'))
MatData$Seeds[MatData$Reproductive == 0] <- NA

# One plant is incorrectly classified as SDL two years in a row.
# Correcting that below
MatData$StageNext[MatData$Plant_Number == 7563] <- 'SDL2'
# Now, begin extracting vital rates!

# First, fecundity parameters----------
FecParams <- MatData %>%
  filter(!is.na(StageNext)) %>%
  group_by(StageNext, Trt) %>%
  summarise(N = n(),
            Seeds = mean(Seeds, na.rm = TRUE),
            pRepro = mean(Reproductive, na.rm = TRUE)) %>%
  data.frame



# next, survival parameters ---------------
SurvParams <-MatData %>%
  group_by(Stage, Trt) %>%
  summarise(N = n(),
            Surv = mean(Survival, na.rm = TRUE)) %>%
  data.frame()

# Finally, transition parameters. These are bit trickier and require-----------------
# some restructuring of the data before they can be summarized

MatData$Transition <- paste(MatData$Stage, MatData$StageNext, sep = ' - ')
MatData$Transition[grepl('NA', MatData$Transition)] <- NA

for(Trans in unique(MatData$Transition)[-1]) {
  MatData[ ,Trans] <- ifelse(MatData$Transition == Trans, 1, 0)
}

TransParams <- MatData %>%
  group_by(Stage, Trt) %>%
  summarise(N = n(),
            SDL_SDL2 = mean(`SDL - SDL2`, na.rm = TRUE),
            SDL_Small = mean(`SDL - Small`, na.rm = TRUE),
            SDL2_Small = mean(`SDL2 - Small`, na.rm = TRUE),
            Small_Small = mean(`Small - Small`, na.rm = TRUE),
            Small_Med = mean(`Small - Med`, na.rm = TRUE),
            Med_Small = mean(`Med - Small`, na.rm = TRUE),
            Med_Med = mean(`Med - Med`, na.rm = TRUE),
            Med_Large = mean(`Med - Large`, na.rm = TRUE),
            Large_Med = mean(`Large - Med`, na.rm = TRUE),
            Large_Large = mean(`Large - Large`, na.rm = TRUE)) %>%
  data.frame()

# Check to make sure that all rows sum to 1
for(i in 1:dim(TransParams)[1]) {
  print(sum(TransParams[i, 4:13]))
}

# Now, start extracting and storing parameters.-----------
# Starting with CR Fecundity
CR_F_Small <- FecParams[FecParams$StageNext == 'Small' & 
                          FecParams$Trt == 'Competition', 'Seeds'] *
              FecParams[FecParams$StageNext == 'Small' & 
                          FecParams$Trt == 'Competition', 'pRepro']

CR_F_Med <- FecParams[FecParams$StageNext == 'Med' & 
                        FecParams$Trt == 'Competition', 'Seeds'] *
  FecParams[FecParams$StageNext == 'Med' & 
              FecParams$Trt == 'Competition', 'pRepro'] 

CR_F_Large <- FecParams[FecParams$StageNext == 'Large' & 
                        FecParams$Trt == 'Competition', 'Seeds'] *
  FecParams[FecParams$StageNext == 'Large' & 
              FecParams$Trt == 'Competition', 'pRepro'] 

# Control Fecundity
Cont_F_Small <- FecParams[FecParams$StageNext == 'Small' & 
                          FecParams$Trt == 'Control', 'Seeds'] *
  FecParams[FecParams$StageNext == 'Small' & 
              FecParams$Trt == 'Control', 'pRepro']

Cont_F_Med <- FecParams[FecParams$StageNext == 'Med' & 
                        FecParams$Trt == 'Control', 'Seeds'] *
  FecParams[FecParams$StageNext == 'Med' & 
              FecParams$Trt == 'Control', 'pRepro'] 

Cont_F_Large <- FecParams[FecParams$StageNext == 'Large' & 
                          FecParams$Trt == 'Control', 'Seeds'] *
  FecParams[FecParams$StageNext == 'Large' & 
              FecParams$Trt == 'Control', 'pRepro'] 

# Survival parameter extraction --------------

# CR
CR_Surv_SDL <- SurvParams[SurvParams$Stage == 'SDL' & 
                            SurvParams$Trt == 'Competition', 'Surv']

CR_Surv_SDL2 <- SurvParams[SurvParams$Stage == 'SDL2' & 
                            SurvParams$Trt == 'Competition', 'Surv']

CR_Surv_Small <- SurvParams[SurvParams$Stage == 'Small' & 
                            SurvParams$Trt == 'Competition', 'Surv']

CR_Surv_Med <- SurvParams[SurvParams$Stage == 'Med' & 
                            SurvParams$Trt == 'Competition', 'Surv']

CR_Surv_Large <- SurvParams[SurvParams$Stage == 'Large' & 
                            SurvParams$Trt == 'Competition', 'Surv']

# Control
Cont_Surv_SDL <- SurvParams[SurvParams$Stage == 'SDL' & 
                            SurvParams$Trt == 'Control', 'Surv']

Cont_Surv_SDL2 <- SurvParams[SurvParams$Stage == 'SDL2' & 
                             SurvParams$Trt == 'Control', 'Surv']

Cont_Surv_Small <- SurvParams[SurvParams$Stage == 'Small' & 
                              SurvParams$Trt == 'Control', 'Surv']

Cont_Surv_Med <- SurvParams[SurvParams$Stage == 'Med' & 
                            SurvParams$Trt == 'Control', 'Surv']

Cont_Surv_Large <- SurvParams[SurvParams$Stage == 'Large' & 
                              SurvParams$Trt == 'Control', 'Surv']

# Transition probabilities -------------------
# Broken down into stasis, progressive growth, and
# retrogressive growth (shrinkage)
CR_Prog_SDL <- TransParams[TransParams$Stage == 'SDL' &
                             TransParams$Trt == 'Competition', 
                           'SDL_SDL2']
CR_Prog2_SDL <- TransParams[TransParams$Stage == 'SDL' &
                              TransParams$Trt == 'Competition',
                            'SDL_Small']
CR_Prog_SDL2 <- TransParams[TransParams$Stage == 'SDL2' &
                              TransParams$Trt == 'Competition', 
                            'SDL2_Small']
CR_Stasis_Small <- TransParams[TransParams$Stage == 'Small' &
                                 TransParams$Trt == 'Competition', 
                               'Small_Small']
CR_Prog_Small <- TransParams[TransParams$Stage == 'Small' &
                               TransParams$Trt == 'Competition' ,
                             'Small_Med']
CR_Retro_Med <- TransParams[TransParams$Stage == 'Med' &
                              TransParams$Trt == 'Competition', 
                            'Med_Small']
CR_Stasis_Med <- TransParams[TransParams$Stage == 'Med' &
                               TransParams$Trt == 'Competition',
                             'Med_Med']
CR_Prog_Med <- TransParams[TransParams$Stage == 'Med' &
                             TransParams$Trt == 'Competition',
                           'Med_Large']
CR_Stasis_Large <- TransParams[TransParams$Stage == 'Large' &
                                 TransParams$Trt == 'Competition',
                               'Large_Large']
CR_Retro_Large <- TransParams[TransParams$Stage == 'Large' &
                                TransParams$Trt == 'Competition',
                              'Large_Med']


# Control
Cont_Prog_SDL <- TransParams[TransParams$Stage == 'SDL' &
                               TransParams$Trt == 'Control', 
                             'SDL_SDL2']
Cont_Prog2_SDL <- TransParams[TransParams$Stage == 'SDL' &
                                TransParams$Trt == 'Control',
                              'SDL_Small']
Cont_Prog_SDL2 <- TransParams[TransParams$Stage == 'SDL2' &
                                TransParams$Trt == 'Control', 
                              'SDL2_Small']
Cont_Stasis_Small <- TransParams[TransParams$Stage == 'Small' &
                                   TransParams$Trt == 'Control', 
                                 'Small_Small']
Cont_Prog_Small <- TransParams[TransParams$Stage == 'Small' &
                                 TransParams$Trt == 'Control' ,
                               'Small_Med']
Cont_Retro_Med <- TransParams[TransParams$Stage == 'Med' &
                                TransParams$Trt == 'Control', 
                              'Med_Small']
Cont_Stasis_Med <- TransParams[TransParams$Stage == 'Med' &
                                 TransParams$Trt == 'Control',
                               'Med_Med']
Cont_Prog_Med <- TransParams[TransParams$Stage == 'Med' &
                               TransParams$Trt == 'Control',
                             'Med_Large']
Cont_Stasis_Large <- TransParams[TransParams$Stage == 'Large' &
                                   TransParams$Trt == 'Control',
                                 'Large_Large']
Cont_Retro_Large <- TransParams[TransParams$Stage == 'Large' &
                                  TransParams$Trt == 'Control',
                                'Large_Med']

# Fixed parameters from Schutzenhofer et al 2009----------
# using the average of their values for cleistagamous and chasmogamous seeds
# values from Table 1 in paper
germCL <- 0.471
germCH <- 0.883
seedViabilityCL <- 0.357
seedViabilityCH <- 0.297

meanGerm <- mean(c(germCL, germCH))
meanSV <- mean(c(seedViabilityCL, seedViabilityCH))

Emergence <- meanGerm * meanSV

# Create matrices------------

A_CR <- matrix(c(0, 0, 0, CR_F_Small * (1 - meanGerm) * meanSV, CR_F_Med * (1 - meanGerm) * meanSV, CR_F_Large * (1 - meanGerm) * meanSV,
                 Emergence, 0 , 0, CR_F_Small * Emergence, CR_F_Med * Emergence, CR_F_Large * Emergence,
                 0, CR_Prog_SDL * CR_Surv_SDL, 0, 0, 0, 0,
                 0, CR_Prog2_SDL * CR_Surv_SDL, CR_Prog_SDL2 * CR_Surv_SDL2, CR_Surv_Small * CR_Stasis_Small, CR_Surv_Med * CR_Retro_Med, 0,
                 0, 0, 0, CR_Surv_Small * CR_Prog_Small, CR_Surv_Med * CR_Stasis_Med, CR_Surv_Large * CR_Retro_Large,
                 0, 0, 0, 0, CR_Surv_Med * CR_Prog_Med, CR_Surv_Large * CR_Stasis_Large), 
               nrow = 6, byrow = TRUE, 
               dimnames = list(c('SB', 'SDL', 'SDL2', 'Small', 'Med', 'Large'),
                               c('SB', 'SDL', 'SDL2', 'Small', 'Med', 'Large')))

A_Cont <- matrix(c(0, 0, 0, Cont_F_Small * (1 - meanGerm) * meanSV, Cont_F_Med * (1 - meanGerm) * meanSV, Cont_F_Large * (1 - meanGerm) * meanSV,
                 Emergence, 0 , 0, Cont_F_Small * Emergence, Cont_F_Med * Emergence, Cont_F_Large * Emergence,
                 0, Cont_Prog_SDL * Cont_Surv_SDL, 0, 0, 0, 0,
                 0, Cont_Prog2_SDL * Cont_Surv_SDL, Cont_Prog_SDL2 * Cont_Surv_SDL2, Cont_Surv_Small * Cont_Stasis_Small, Cont_Surv_Med * Cont_Retro_Med, 0,
                 0, 0, 0, Cont_Surv_Small * Cont_Prog_Small, Cont_Surv_Med * Cont_Stasis_Med, Cont_Surv_Large * Cont_Retro_Large,
                 0, 0, 0, 0, Cont_Surv_Med * Cont_Prog_Med, Cont_Surv_Large * Cont_Stasis_Large), 
               nrow = 6, byrow = TRUE, 
               dimnames = list(c('SB', 'SDL', 'SDL2', 'Small', 'Med', 'Large'),
                               c('SB', 'SDL', 'SDL2', 'Small', 'Med', 'Large')))

# Extract per-capita growth rates
ev <- eigen(A_Cont)
lmax <- which(Re(ev$values) == max(Re(ev$values)))
lmax <- lmax[1]  
Cont_Lambda <- Re(ev$values[lmax])
Cont_Lambda # Seems more reasonable than previous model

ev <- eigen(A_CR)

lmax <- which(Re(ev$values) == max(Re(ev$values)))
lmax <- lmax[1]  
CR_Lambda  <- Re(ev$values[lmax])
CR_Lambda # ibid

# Store parameter values for plotting later
values <- c(Cont_F_Large, CR_F_Large,
            Cont_F_Med, CR_F_Med,
            Cont_F_Small, CR_F_Small,
            Cont_Prog_SDL, CR_Prog_SDL,
            Cont_Prog2_SDL, CR_Prog2_SDL,
            Cont_Prog_SDL2, CR_Prog_SDL2,
            Cont_Prog_Small, CR_Prog_Small,
            Cont_Stasis_Small, CR_Stasis_Small,
            Cont_Prog_Med, CR_Prog_Med,
            Cont_Stasis_Med, CR_Stasis_Med,
            Cont_Retro_Med, CR_Retro_Med,
            Cont_Stasis_Large, CR_Stasis_Large,
            Cont_Retro_Large, CR_Retro_Large,
            Cont_Surv_SDL, CR_Surv_SDL,
            Cont_Surv_SDL2, CR_Surv_SDL2,
            Cont_Surv_Small, CR_Surv_Small,
            Cont_Surv_Med, CR_Surv_Med,
            Cont_Surv_Large, CR_Surv_Large,
            Cont_Lambda, CR_Lambda)

# set up bootstrapping-------------
nResamp <- 1000
# This nasty looking line of code initializes all of the vectors to store the outputs
# from the bootstrapping loop
eval(parse(text = paste0('boot_', 
                         ls()[grepl('Cont|CR', ls())], 
                         ' <- rep(NA, nResamp)')))

# Remove the erroneously created vectors which are actually matrices. 
# the above step saves a lot of typing, but it's far from perfect!
rm(boot_A_Cont, boot_A_CR)
boot_ES_Lambda <- rep(NA, nResamp)

CR_SDL <- filter(MatData, Stage == 'SDL' & Trt == 'Competition')
CR_SDL2 <- filter(MatData, Stage == 'SDL2' & Trt == 'Competition')
CR_Small <- filter(MatData, Stage == 'Small' & Trt == 'Competition')
CR_Med <- filter(MatData, Stage == 'Med' & Trt == 'Competition')
# This is to deal with fecundity to calculations since we have 
# small numbers of plants that were large at t+1, but are 
# using seeds at t+1 to calculate fecundity
CR_Large <- filter(MatData, Stage == 'Large'  &
                     Trt == 'Competition')
Large_CR_T_1 <- which(CR_Large$StageNext == 'Large') 


Cont_SDL <- filter(MatData, Stage == 'SDL' & Trt == 'Control')
Cont_SDL2 <- filter(MatData, Stage == 'SDL2' & Trt == 'Control')
Cont_Small <- filter(MatData, Stage == 'Small' & Trt == 'Control')
Cont_Med <- filter(MatData, Stage == 'Med' & Trt == 'Control')
Cont_Large <- filter(MatData, Stage == 'Large' &
                       Trt == 'Control')
# add these in in case they aren't resampled due to tiny sample size in this
# stage class x trt combination
Large_Control_T_1 <- which(Cont_Large$StageNext == 'Large')

# Begin bootstrapping loop --------------
for(i in 1:nResamp) {
  # Initialize resamplers for each stage class x trt combo
  x1 <- sample(1:dim(CR_SDL)[1], dim(CR_SDL)[1], replace = TRUE)
  x2 <- sample(1:dim(CR_SDL2)[1], dim(CR_SDL2)[1], replace = TRUE)
  x3 <- sample(1:dim(CR_Small)[1], dim(CR_Small)[1], replace = TRUE)
  x4 <- sample(1:dim(CR_Med)[1], dim(CR_Med)[1], replace = TRUE)
  x5 <- sample(1:dim(CR_Large)[1], dim(CR_Large)[1], replace = TRUE)
  
  x6 <- sample(1:dim(Cont_SDL)[1], dim(Cont_SDL)[1], replace = TRUE)
  x7 <- sample(1:dim(Cont_SDL2)[1], dim(Cont_SDL2)[1], replace = TRUE)
  x8 <- sample(1:dim(Cont_Small)[1], dim(Cont_Small)[1], replace = TRUE)
  x9 <- sample(1:dim(Cont_Med)[1], dim(Cont_Med)[1], replace = TRUE)
  x10 <- sample(1:dim(Cont_Large)[1], dim(Cont_Large)[1], replace = TRUE)
  
  if(all(!Large_Control_T_1 %in% x10)){
    x10[1] <- sample(Large_Control_T_1, 1)
  }
  if(all(!Large_CR_T_1 %in% x5)) {
    x5[1] <- sample(Large_CR_T_1, 1)
  }
  
  # Resample!
  boot_CR_SDL <- CR_SDL[x1, ]
  boot_CR_SDL2 <- CR_SDL2[x2, ]
  boot_CR_Small <- CR_Small[x3, ]
  boot_CR_Med <- CR_Med[x4, ]
  boot_CR_Large <- CR_Large[x5, ]
  
  boot_Cont_SDL <- Cont_SDL[x6, ]
  boot_Cont_SDL2 <- Cont_SDL2[x7, ]
  boot_Cont_Small <- Cont_Small[x8, ]
  boot_Cont_Med <- Cont_Med[x9, ]
  boot_Cont_Large <- Cont_Large[x10, ]
  
  boot_MatData <- rbind(boot_CR_SDL, boot_Cont_SDL,
                        boot_CR_SDL2, boot_Cont_SDL2,
                        boot_CR_Small, boot_Cont_Small,
                        boot_CR_Med, boot_Cont_Med,
                        boot_CR_Large, boot_Cont_Large)
  
  # bootstrapped parameters
  boot_FecParams <- boot_MatData %>%
    filter(!is.na(StageNext)) %>%
    group_by(StageNext, Trt) %>%
    summarise(N = n(),
              Seeds = mean(Seeds, na.rm = TRUE),
              pRepro = mean(Reproductive, na.rm = TRUE)) %>%
    data.frame
  # if a given vital rate can't be computed for an iteration,
  # then impute the observed value
  for(j in 1:dim(boot_FecParams)[1]) {
    for(k in 1:dim(boot_FecParams)[2]){
      if(is.na(boot_FecParams[j, k])){
        boot_FecParams[j, k] <- FecParams[j, k]
      }
    }
  }
  
  boot_SurvParams <- boot_MatData %>%
    group_by(Stage, Trt) %>%
    summarise(N = n(),
              Surv = mean(Survival, na.rm = TRUE)) %>%
    data.frame
  
  # if a given vital rate can't be computed for an iteration,
  # then impute the observed value
  for(j in 1:dim(boot_SurvParams)[1]) {
    for(k in 1:dim(boot_SurvParams)[2]){
      if(is.na(boot_SurvParams[j, k])){
        boot_SurvParams[j, k] <- SurvParams[j, k]
      }
    }
  }
  
  boot_TransParams <- boot_MatData %>%
    group_by(Stage, Trt) %>%
    summarise(N = n(),
              SDL_SDL2 = mean(`SDL - SDL2`, na.rm = TRUE),
              SDL_Small = mean(`SDL - Small`, na.rm = TRUE),
              SDL2_Small = mean(`SDL2 - Small`, na.rm = TRUE),
              Small_Small = mean(`Small - Small`, na.rm = TRUE),
              Small_Med = mean(`Small - Med`, na.rm = TRUE),
              Med_Small = mean(`Med - Small`, na.rm = TRUE),
              Med_Med = mean(`Med - Med`, na.rm = TRUE),
              Med_Large = mean(`Med - Large`, na.rm = TRUE),
              Large_Med = mean(`Large - Med`, na.rm = TRUE),
              Large_Large = mean(`Large - Large`, na.rm = TRUE)) %>%
    data.frame
  
  # if a given vital rate can't be computed for an iteration,
  # then impute the observed value
  for(j in 1:dim(boot_TransParams)[1]) {
    for(k in 1:dim(boot_TransParams)[2]){
      if(is.na(boot_TransParams[j, k])){
        boot_TransParams[j, k] <- TransParams[j, k]
      }
    }
  }
  
  # extracting bootstrapped vital rates -------------------
  boot_CR_F_Small[i] <- CR_F_Small <- boot_FecParams[boot_FecParams$StageNext == 'Small' & 
                                                       boot_FecParams$Trt == 'Competition', 
                                                     'Seeds'] *
    boot_FecParams[boot_FecParams$StageNext == 'Small' & 
                     boot_FecParams$Trt == 'Competition', 
                   'pRepro']
  
  boot_CR_F_Med[i] <- CR_F_Med <- boot_FecParams[boot_FecParams$StageNext == 'Med' & 
                          boot_FecParams$Trt == 'Competition', 'Seeds'] *
    boot_FecParams[boot_FecParams$StageNext == 'Med' & 
                boot_FecParams$Trt == 'Competition', 'pRepro'] 
  
  boot_CR_F_Large[i] <- CR_F_Large <- boot_FecParams[boot_FecParams$StageNext == 'Large' & 
                            boot_FecParams$Trt == 'Competition', 'Seeds'] *
    boot_FecParams[boot_FecParams$StageNext == 'Large' & 
                boot_FecParams$Trt == 'Competition', 'pRepro'] 
  
  # Control Fecundity
  boot_Cont_F_Small[i] <- Cont_F_Small <- boot_FecParams[boot_FecParams$StageNext == 'Small' & 
                              boot_FecParams$Trt == 'Control', 'Seeds'] *
    boot_FecParams[boot_FecParams$StageNext == 'Small' & 
                boot_FecParams$Trt == 'Control', 'pRepro']
  
  boot_Cont_F_Med[i] <- Cont_F_Med <- boot_FecParams[boot_FecParams$StageNext == 'Med' & 
                            boot_FecParams$Trt == 'Control', 'Seeds'] *
    boot_FecParams[boot_FecParams$StageNext == 'Med' & 
                boot_FecParams$Trt == 'Control', 'pRepro'] 
  
  boot_Cont_F_Large[i] <- Cont_F_Large <- boot_FecParams[boot_FecParams$StageNext == 'Large' & 
                              boot_FecParams$Trt == 'Control', 'Seeds'] *
    boot_FecParams[boot_FecParams$StageNext == 'Large' & 
                boot_FecParams$Trt == 'Control', 'pRepro'] 
  
  # Survival parameter extraction --------------
  
  # CR
  boot_CR_Surv_SDL[i] <- CR_Surv_SDL <- boot_SurvParams[boot_SurvParams$Stage == 'SDL' & 
                              boot_SurvParams$Trt == 'Competition', 'Surv']
  
  boot_CR_Surv_SDL2[i] <- CR_Surv_SDL2 <- boot_SurvParams[boot_SurvParams$Stage == 'SDL2' & 
                               boot_SurvParams$Trt == 'Competition', 'Surv']
  
  boot_CR_Surv_Small[i] <- CR_Surv_Small <- boot_SurvParams[boot_SurvParams$Stage == 'Small' & 
                                boot_SurvParams$Trt == 'Competition', 'Surv']
  
  boot_CR_Surv_Med[i] <- CR_Surv_Med <- boot_SurvParams[boot_SurvParams$Stage == 'Med' & 
                              boot_SurvParams$Trt == 'Competition', 'Surv']
  
  boot_CR_Surv_Large[i] <- CR_Surv_Large <- boot_SurvParams[boot_SurvParams$Stage == 'Large' & 
                                boot_SurvParams$Trt == 'Competition', 'Surv']
  
  # Control
  boot_Cont_Surv_SDL[i] <- Cont_Surv_SDL <- boot_SurvParams[boot_SurvParams$Stage == 'SDL' & 
                                boot_SurvParams$Trt == 'Control', 'Surv']
  
  boot_Cont_Surv_SDL2[i] <- Cont_Surv_SDL2 <- boot_SurvParams[boot_SurvParams$Stage == 'SDL2' & 
                                 boot_SurvParams$Trt == 'Control', 'Surv']
  
  boot_Cont_Surv_Small[i] <- Cont_Surv_Small <- boot_SurvParams[boot_SurvParams$Stage == 'Small' & 
                                  boot_SurvParams$Trt == 'Control', 'Surv']
  
  boot_Cont_Surv_Med[i] <- Cont_Surv_Med <- boot_SurvParams[boot_SurvParams$Stage == 'Med' & 
                                boot_SurvParams$Trt == 'Control', 'Surv']
  
  boot_Cont_Surv_Large[i] <- Cont_Surv_Large <- boot_SurvParams[boot_SurvParams$Stage == 'Large' & 
                                  boot_SurvParams$Trt == 'Control', 'Surv']
  
  # Transition probabilities -------------------
  # Broken down into stasis, progressive growth, and
  # retrogressive growth (shrinkage)
  boot_CR_Prog_SDL[i] <- CR_Prog_SDL <- boot_TransParams[boot_TransParams$Stage == 'SDL' &
                               boot_TransParams$Trt == 'Competition', 
                             'SDL_SDL2']
  boot_CR_Prog2_SDL[i] <- CR_Prog2_SDL <- boot_TransParams[boot_TransParams$Stage == 'SDL' &
                                boot_TransParams$Trt == 'Competition',
                              'SDL_Small']
  boot_CR_Prog_SDL2[i] <- CR_Prog_SDL2 <- boot_TransParams[boot_TransParams$Stage == 'SDL2' &
                                boot_TransParams$Trt == 'Competition', 
                              'SDL2_Small']
  boot_CR_Stasis_Small[i] <- CR_Stasis_Small <- boot_TransParams[boot_TransParams$Stage == 'Small' &
                                   boot_TransParams$Trt == 'Competition', 
                                 'Small_Small']
  boot_CR_Prog_Small[i] <- CR_Prog_Small <- boot_TransParams[boot_TransParams$Stage == 'Small' &
                                 boot_TransParams$Trt == 'Competition' ,
                               'Small_Med']
  boot_CR_Retro_Med[i] <- CR_Retro_Med <- boot_TransParams[boot_TransParams$Stage == 'Med' &
                                boot_TransParams$Trt == 'Competition', 
                              'Med_Small']
  boot_CR_Stasis_Med[i] <- CR_Stasis_Med <- boot_TransParams[boot_TransParams$Stage == 'Med' &
                                 boot_TransParams$Trt == 'Competition',
                               'Med_Med']
  boot_CR_Prog_Med[i] <- CR_Prog_Med <- boot_TransParams[boot_TransParams$Stage == 'Med' &
                               boot_TransParams$Trt == 'Competition',
                             'Med_Large']
  boot_CR_Stasis_Large[i] <- CR_Stasis_Large <- boot_TransParams[boot_TransParams$Stage == 'Large' &
                                   boot_TransParams$Trt == 'Competition',
                                 'Large_Large']
  boot_CR_Retro_Large[i] <- CR_Retro_Large <- boot_TransParams[boot_TransParams$Stage == 'Large' &
                                  boot_TransParams$Trt == 'Competition',
                                'Large_Med']
  
  
  # Control
  boot_Cont_Prog_SDL[i] <- Cont_Prog_SDL <- boot_TransParams[boot_TransParams$Stage == 'SDL' &
                                                           boot_TransParams$Trt == 'Control', 
                                                         'SDL_SDL2']
  boot_Cont_Prog2_SDL[i] <- Cont_Prog2_SDL <- boot_TransParams[boot_TransParams$Stage == 'SDL' &
                                                             boot_TransParams$Trt == 'Control',
                                                           'SDL_Small']
  boot_Cont_Prog_SDL2[i] <- Cont_Prog_SDL2 <- boot_TransParams[boot_TransParams$Stage == 'SDL2' &
                                                             boot_TransParams$Trt == 'Control', 
                                                           'SDL2_Small']
  boot_Cont_Stasis_Small[i] <- Cont_Stasis_Small <- boot_TransParams[boot_TransParams$Stage == 'Small' &
                                                                   boot_TransParams$Trt == 'Control', 
                                                                 'Small_Small']
  boot_Cont_Prog_Small[i] <- Cont_Prog_Small <- boot_TransParams[boot_TransParams$Stage == 'Small' &
                                                               boot_TransParams$Trt == 'Control' ,
                                                             'Small_Med']
  boot_Cont_Retro_Med[i] <- Cont_Retro_Med <- boot_TransParams[boot_TransParams$Stage == 'Med' &
                                                             boot_TransParams$Trt == 'Control', 
                                                           'Med_Small']
  boot_Cont_Stasis_Med[i] <- Cont_Stasis_Med <- boot_TransParams[boot_TransParams$Stage == 'Med' &
                                                               boot_TransParams$Trt == 'Control',
                                                             'Med_Med']
  boot_Cont_Prog_Med[i] <- Cont_Prog_Med <- boot_TransParams[boot_TransParams$Stage == 'Med' &
                                                           boot_TransParams$Trt == 'Control',
                                                         'Med_Large']
  boot_Cont_Stasis_Large[i] <- Cont_Stasis_Large <- boot_TransParams[boot_TransParams$Stage == 'Large' &
                                                                   boot_TransParams$Trt == 'Control',
                                                                 'Large_Large']
  boot_Cont_Retro_Large[i] <- Cont_Retro_Large <- boot_TransParams[boot_TransParams$Stage == 'Large' &
                                                                 boot_TransParams$Trt == 'Control',
                                                               'Large_Med']
  
  
  # Create matrices------------
  
  A_CR <- matrix(c(0, 0, 0, CR_F_Small * (1 - meanGerm) * meanSV, CR_F_Med * (1 - meanGerm) * meanSV, CR_F_Large * (1 - meanGerm) * meanSV,
                   Emergence, 0 , 0, CR_F_Small * Emergence, CR_F_Med * Emergence, CR_F_Large * Emergence,
                   0, CR_Prog_SDL * CR_Surv_SDL, 0, 0, 0, 0,
                   0, CR_Prog2_SDL * CR_Surv_SDL, CR_Prog_SDL2 * CR_Surv_SDL2, CR_Surv_Small * CR_Stasis_Small, CR_Surv_Med * CR_Retro_Med, 0,
                   0, 0, 0, CR_Surv_Small * CR_Prog_Small, CR_Surv_Med * CR_Stasis_Med, CR_Surv_Large * CR_Retro_Large,
                   0, 0, 0, 0, CR_Surv_Med * CR_Prog_Med, CR_Surv_Large * CR_Stasis_Large), 
                 nrow = 6, byrow = TRUE, 
                 dimnames = list(c('SB', 'SDL', 'SDL2', 'Small', 'Med', 'Large'),
                                 c('SB', 'SDL', 'SDL2', 'Small', 'Med', 'Large')))
  
  A_Cont <- matrix(c(0, 0, 0, Cont_F_Small * (1 - meanGerm) * meanSV, Cont_F_Med * (1 - meanGerm) * meanSV, Cont_F_Large * (1 - meanGerm) * meanSV,
                     Emergence, 0 , 0, Cont_F_Small * Emergence, Cont_F_Med * Emergence, Cont_F_Large * Emergence,
                     0, Cont_Prog_SDL * Cont_Surv_SDL, 0, 0, 0, 0,
                     0, Cont_Prog2_SDL * Cont_Surv_SDL, Cont_Prog_SDL2 * Cont_Surv_SDL2, Cont_Surv_Small * Cont_Stasis_Small, Cont_Surv_Med * Cont_Retro_Med, 0,
                     0, 0, 0, Cont_Surv_Small * Cont_Prog_Small, Cont_Surv_Med * Cont_Stasis_Med, Cont_Surv_Large * Cont_Retro_Large,
                     0, 0, 0, 0, Cont_Surv_Med * Cont_Prog_Med, Cont_Surv_Large * Cont_Stasis_Large), 
                   nrow = 6, byrow = TRUE, 
                   dimnames = list(c('SB', 'SDL', 'SDL2', 'Small', 'Med', 'Large'),
                                   c('SB', 'SDL', 'SDL2', 'Small', 'Med', 'Large')))
  
  # Extract per-capita growth rates
  ev <- eigen(A_Cont)
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lmax <- lmax[1]  
  boot_Cont_Lambda[i] <- Re(ev$values[lmax])

  ev <- eigen(A_CR)
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lmax <- lmax[1]  
  boot_CR_Lambda[i]  <- Re(ev$values[lmax])
  
  boot_ES_Lambda[i] <- log(boot_CR_Lambda[i] + 0.5) - log(boot_Cont_Lambda[i] + 0.5)
  
}

# Sort all of these vectors ------------
ForSorting <- ls(pattern = 'boot_CR_|boot_Cont_|boot_ES_')

for(i in 1:length(ls())){
  if(length(eval(parse(text = ForSorting[i]))) > 999)
    eval(parse(text = paste(ForSorting[i],
                            '<- sort(',
                            ForSorting[i],
                            ')', 
                            sep = "")))
}

# Extract confidence intervals ------------
UpCI <- c(boot_Cont_F_Large[975], boot_CR_F_Large[975],
          boot_Cont_F_Med[975], boot_CR_F_Med[975],
          boot_Cont_F_Small[975], boot_CR_F_Small[975],
          boot_Cont_Prog_SDL[975], boot_CR_Prog_SDL[975],
          boot_Cont_Prog2_SDL[975], boot_CR_Prog2_SDL[975],
          boot_Cont_Prog_SDL2[975], boot_CR_Prog_SDL2[975],
          boot_Cont_Prog_Small[975], boot_CR_Prog_Small[975],
          boot_Cont_Stasis_Small[975], boot_CR_Stasis_Small[975],
          boot_Cont_Prog_Med[975], boot_CR_Prog_Med[975],
          boot_Cont_Stasis_Med[975], boot_CR_Stasis_Med[975],
          boot_Cont_Retro_Med[975], boot_CR_Retro_Med[975],
          boot_Cont_Stasis_Large[975], boot_CR_Stasis_Large[975],
          boot_Cont_Retro_Large[975], boot_CR_Retro_Large[975],
          boot_Cont_Surv_SDL[975], boot_CR_Surv_SDL[975],
          boot_Cont_Surv_SDL2[975], boot_CR_Surv_SDL2[975],
          boot_Cont_Surv_Small[975], boot_CR_Surv_Small[975],
          boot_Cont_Surv_Med[975], boot_CR_Surv_Med[975],
          boot_Cont_Surv_Large[975], boot_CR_Surv_Large[975],
          boot_Cont_Lambda[975], boot_CR_Lambda[975])

LoCI <- c(boot_Cont_F_Large[25], boot_CR_F_Large[25],
          boot_Cont_F_Med[25], boot_CR_F_Med[25],
          boot_Cont_F_Small[25], boot_CR_F_Small[25],
          boot_Cont_Prog_SDL[25], boot_CR_Prog_SDL[25],
          boot_Cont_Prog2_SDL[25], boot_CR_Prog2_SDL[25],
          boot_Cont_Prog_SDL2[25], boot_CR_Prog_SDL2[25],
          boot_Cont_Prog_Small[25], boot_CR_Prog_Small[25],
          boot_Cont_Stasis_Small[25], boot_CR_Stasis_Small[25],
          boot_Cont_Prog_Med[25], boot_CR_Prog_Med[25],
          boot_Cont_Stasis_Med[25], boot_CR_Stasis_Med[25],
          boot_Cont_Retro_Med[25], boot_CR_Retro_Med[25],
          boot_Cont_Stasis_Large[25], boot_CR_Stasis_Large[25],
          boot_Cont_Retro_Large[25], boot_CR_Retro_Large[25],
          boot_Cont_Surv_SDL[25], boot_CR_Surv_SDL[25],
          boot_Cont_Surv_SDL2[25], boot_CR_Surv_SDL2[25],
          boot_Cont_Surv_Small[25], boot_CR_Surv_Small[25],
          boot_Cont_Surv_Med[25], boot_CR_Surv_Med[25],
          boot_Cont_Surv_Large[25], boot_CR_Surv_Large[25],
          boot_Cont_Lambda[25], boot_CR_Lambda[25])

# Create data frame to store everything -------------------
results <- tibble(Trt = rep(c('Control', 'CR'), 19),
                  VitalRate = factor(c(rep('paste(italic(f[i]))', 6),
                                       rep('paste(italic(p[ij]))', 8),
                                       rep('paste(italic(c[ii]))', 2),
                                       rep('paste(italic(p[ij]))', 2),
                                       rep('paste(italic(c[ii]))', 2),
                                       rep('paste(italic(r[ij]))', 2),
                                       rep('paste(italic(c[ii]))', 2),
                                       rep('paste(italic(r[ij]))', 2),
                                       rep('paste(italic(s[i]))', 10),
                                       rep('lambda',2)),
                                     levels = c('paste(italic(s[i]))',
                                                'paste(italic(c[ii]))',
                                                'paste(italic(p[ij]))',
                                                'paste(italic(r[ij]))',
                                                'paste(italic(f[i]))',
                                                'lambda'),
                                     ordered = TRUE),
                  Stage = c(rep('Large Adult', 2),
                            rep('Medium Adult', 2),
                            rep('Small Adult', 2),
                            rep('Seedling', 4),
                            rep('2nd Year Seedling', 2),
                            rep('Small Adult', 4),
                            rep('Medium Adult', 6),
                            rep('Large Adult', 4),
                            rep('Seedling', 2),
                            rep('2nd Year Seedling', 2),
                            rep('Small Adult', 2),
                            rep('Medium Adult', 2),
                            rep('Large Adult', 2),
                            rep('Full Population', 2)),
                  values = values,
                  UpCI = UpCI,
                  LoCI = LoCI) %>%
  mutate(Trt_Stage = paste(.$Trt, .$Stage,.$VitalRate, sep = " "))

ggplot(data = results, aes(x = Trt)) + 
  geom_point(aes(x = Trt_Stage,
                 y = values,
                 color = Trt,
                 shape = Stage),
             size = 4.5) +
  geom_linerange(aes(x = Trt_Stage,
                     ymin = LoCI,
                     ymax = UpCI,
                     color = Trt),
                 size = 1.25) + 
  facet_wrap(~VitalRate,
             scales = 'free',
             labeller = label_parsed,
             nrow = 3,
             ncol = 2) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = 'black',
                                    fill = NA),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 18,
                                    margin = margin(t = 20,
                                                    l = 0,
                                                    r = 0,
                                                    b = 5)),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  scale_color_manual('Treatment',
                     breaks = c('Control', 'CR'),
                     values = c('black', 'green')) + 
  scale_shape_manual('Stage',
                     breaks = c('Seedling',
                                '2nd Year Seedling',
                                'Small Adult',
                                'Medium Adult',
                                'Large Adult',
                                'Full Population'),
                     values = c(15, 16, 17, 18, 25, 4)) +
  scale_y_continuous('') +
  scale_x_discrete('')

ggsave(filename = 'Lespedeza_VR_Panel.png',
       path = 'Lespedeza/MatModel',
       height = 6,
       width = 9,
       unit = 'in')
