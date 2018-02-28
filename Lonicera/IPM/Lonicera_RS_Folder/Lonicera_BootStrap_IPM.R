# Lonicera IPM for unburned controls and CR treatments
# Includes boot strapping and is intended to run on iDiv's
# RStudio server. COmplete model selection process can
# be found in IPM/R/Lonicera_IPM.R

rm(list = ls())
graphics.off()

library(dplyr)
library(mgcv)
library(brms)

setwd('Lonicera/IPM/Lonicera_RS_Folder') # local de-bugging
# setwd('~/Lonicera_RS_Folder') # RS server
# Read in data and get rid of extraneous info
AllPlants <- read.csv('LM_Clean.csv',
                      stringsAsFactors = FALSE)

AllCR <- filter(AllPlants, Treatment == 'CompN')
AllCont <- filter(AllPlants, Treatment == 'ContN')
AllBig <- filter(AllPlants, Treatment == 'AllN')

# Growth Regressions. Model selection process omitted here
CompGrowLM <- lm(Plant_Height13 ~ Plant_Height12, data = rbind(AllCR, AllBig))
ContGrowLM <- lm(Plant_Height13 ~ Plant_Height12, data = rbind(AllCont, AllBig))

# Survival. Logistic regressions w/ linear and polynomial terms to start off
CompSurvBRM_Quad <- brm(Survival~ Plant_Height12 + I(Plant_Height12^2), 
                        data = rbind(AllCR, AllBig),
                        family = 'bernoulli',
                        control = list(adapt_delta = 0.97,
                                       max_treedepth = 15),
                        cores = getOption('mc.cores', 2L),
                        iter = 4000,
                        chains = 4,
                        warmup = 1000,
                        save_model = 'CR_Surv_BRM_Quad.stan') %>%
  add_waic()

ContSurvBRM_Quad <-  brm(Survival~ Plant_Height12 + I(Plant_Height12^2), 
                         data = rbind(AllCont, AllBig),
                         family = 'bernoulli',
                         control = list(adapt_delta = 0.97,
                                        max_treedepth = 15),
                         cores = getOption('mc.cores', 2L),
                         iter = 4000,
                         chains = 4,
                         warmup = 1000,
                         save_model = 'Cont_Surv_BRM_Quad.stan') %>%
  add_waic()

CompSurvBRM_Lin <- brm(Survival~ Plant_Height12, 
                       data = rbind(AllCR, AllBig),
                       family = 'bernoulli',
                       control = list(adapt_delta = 0.97,
                                      max_treedepth = 15),
                       cores = getOption('mc.cores', 2L),
                       iter = 4000,
                       chains = 4,
                       warmup = 1000,
                       save_model = 'CR_Surv_BRM_Lin.stan') %>%
  add_waic()

ContSurvBRM_Lin <-  brm(Survival~ Plant_Height12, 
                        data = rbind(AllCont, AllBig),
                        family = 'bernoulli',
                        control = list(adapt_delta = 0.97,
                                       max_treedepth = 15),
                        cores = getOption('mc.cores', 2L),
                        iter = 4000,
                        chains = 4,
                        warmup = 1000,
                        save_model = 'Cont_Surv_BRM_Lin.stan') %>%
  add_waic()

# Store chains to resample during bootstrapping
SurvChainsCR_Quad <- CompSurvBRM_Quad$fit@sim$samples[[1]][1:3] %>%
  data.frame()
SurvParamsCR_Quad <- fixef(CompSurvBRM_Quad)[ ,1]
SurvChainsCont_Quad <- ContSurvBRM_Quad$fit@sim$samples[[1]][1:3] %>%
  data.frame()
SurvParamsCont_Quad <- fixef(ContSurvBRM_Quad)[ ,1]


SurvChainsCR_Lin <- CompSurvBRM_Lin$fit@sim$samples[[1]][1:2] %>%
  data.frame()
SurvParamsCR_Lin <- fixef(CompSurvBRM_Lin)[ ,1]
SurvChainsCont_Lin <- ContSurvBRM_Lin$fit@sim$samples[[1]][1:2] %>%
  data.frame()
SurvParamsCont_Lin <- fixef(ContSurvBRM_Lin)[ ,1]



# Reproduction consists of 4 parameters:
# Pr(Repro) = Reproductive ~ size (t+1) 
# Seeds = Seeds ~ size(t+1)
# Sdl size dist = dnorm(mu, sd)
# est prob = constant

# Pr(Repro)
ReproGLM <- glm(Reproductive ~ Plant_Height13, data = AllPlants, family = binomial())

# Regression is underdispersed. Using a quasibinomial instead
ReproGLM <- glm(Reproductive ~ Plant_Height13, data = AllPlants, family = quasibinomial())

# Seeds next. Creating a column for seeds by multiplying fruit# by 
# average seeds/fruit from Amibeth's control treatment.
AllPlants$Seeds <- round(AllPlants$Fruit * 2.79)
SeedGLM <- glm(Seeds ~ Plant_Height13, data = AllPlants, family = quasipoisson())

# Recruit size distribution. 
Seedlings <- filter(AllPlants, Stage13 == 'S')
SdlMean <- mean(Seedlings$Plant_Height13, na.rm = TRUE)
SdlSD <- sd(Seedlings$Plant_Height13, na.rm = TRUE)


# Data hard coded from Rae's 2011-12 population data in unburned controls
# (density x fire x invasive removal experiment)
EstProb <- 0.003562864


FecModels <- list(PRep = ReproGLM,
                  Seeds = SeedGLM,
                  EstProb = EstProb,
                  SizeDist = list(Mean = SdlMean,
                                  SD = SdlSD))  


source('Lonicera_Utility_Functions.R')

# Begin setting up IPM 
MinSize <- 0.9 * min(c(AllPlants$Plant_Height12, AllPlants$Plant_Height13), na.rm = TRUE)
MaxSize <- 1.1 * max(c(AllPlants$Plant_Height12, AllPlants$Plant_Height13), na.rm = TRUE)
nMeshPts <- 500
Boundaries <- MinSize + c(0:nMeshPts) * (MaxSize - MinSize) / nMeshPts
MeshPts <- 0.5 * (Boundaries[1:nMeshPts] + Boundaries[2:(nMeshPts + 1)])
CellSize <- MeshPts[2] - MeshPts[1]

# Construct the P kernel. First, make growth transition matrix and correct for eviction.
# then create survival vector and multiply it in
GMat_Cont <- CellSize* outer(MeshPts, MeshPts, FUN = GrowFun, Model = ContGrowLM)
GMat_Comp <- CellSize * outer(MeshPts, MeshPts, FUN = GrowFun, Model = CompGrowLM)

GMat_Cont <- GMat_Cont/matrix(as.vector(apply(GMat_Cont,
                                              2,
                                              sum)),
                              nrow = nMeshPts,
                              ncol = nMeshPts,
                              byrow = TRUE)

GMat_Comp <- GMat_Comp/matrix(as.vector(apply(GMat_Comp,
                                              2,
                                              sum)),
                              nrow = nMeshPts,
                              ncol = nMeshPts,
                              byrow = TRUE)

SVec_Cont_Quad <- SurvFun(MeshPts, SurvParamsCont_Quad)
SVec_Comp_Quad <- SurvFun(MeshPts, SurvParamsCR_Quad)

SVec_Cont_Lin <- SurvFun(MeshPts, SurvParamsCont_Lin)
SVec_Comp_Lin <- SurvFun(MeshPts, SurvParamsCR_Lin)


PMat_Cont_Quad <- PMat_Comp_Quad <- matrix(0, nMeshPts, nMeshPts)

PMat_Cont_Lin <- PMat_Comp_Lin <- matrix(0, nMeshPts, nMeshPts)

for(i in seq_len(nMeshPts)){
  PMat_Cont_Quad[ ,i] <- GMat_Cont[ ,i] * SVec_Cont_Quad[i]
  PMat_Comp_Quad[ ,i] <- GMat_Comp[ ,i] * SVec_Comp_Quad[i]
  
  PMat_Cont_Lin[ ,i] <- GMat_Cont[ ,i] * SVec_Cont_Lin[i]
  PMat_Comp_Lin[ ,i] <- GMat_Comp[ ,i] * SVec_Comp_Lin[i]
}
# Fecundity kernel is constant for both treatments. Create that, then
# create K kernel and calculate lambdas
FMat <- CellSize * outer(MeshPts, MeshPts, FUN = FecFun, Models = FecModels)

KMat_Cont_Quad <- PMat_Cont_Quad + FMat
KMat_Comp_Quad <- PMat_Comp_Quad + FMat

KMat_Cont_Lin <- PMat_Cont_Lin + FMat
KMat_Comp_Lin <- PMat_Comp_Lin + FMat

ContEV_Quad <- eigen(KMat_Cont_Quad)
CompEV_Quad <- eigen(KMat_Comp_Quad)
ContLambda_Quad_obs <- max(Re(ContEV_Quad$values))
CompLambda_Quad_obs <- max(Re(CompEV_Quad$values))
ContEV_Lin <- eigen(KMat_Cont_Lin)
CompEV_Lin <- eigen(KMat_Comp_Lin)
ContLambda_Lin_obs <- max(Re(ContEV_Lin$values))
CompLambda_Lin_obs <- max(Re(CompEV_Lin$values))

ContLambda_Quad_obs
ContLambda_Lin_obs

CompLambda_Quad_obs
CompLambda_Lin_obs


# Bootstrapping------------
CRData <- filter(AllPlants, Treatment == 'CompN')
ContData <- filter(AllPlants, Treatment == 'ContN')
BigData <- filter(AllPlants, Treatment == 'AllN')

nCR <- dim(CRData)[1]
nCont <- dim(ContData)[1]
nBig <- dim(BigData)[1]
nSurv <- dim(SurvChainsCont_Quad)[1]

nBootSamples <- 1000

OutputValues <- list(GrowSlope_CR = rep(NA, nBootSamples),
                     GrowInt_CR = rep(NA, nBootSamples),
                     GrowSlope_Cont = rep(NA, nBootSamples),
                     GrowInt_Cont = rep(NA, nBootSamples),
                     RepSlope = rep(NA, nBootSamples),
                     RepInt = rep(NA, nBootSamples),
                     SeedSlope = rep(NA, nBootSamples),
                     SeedInt = rep(NA, nBootSamples),
                     RecMean = rep(NA, nBootSamples),
                     RecSD = rep(NA, nBootSamples),
                     Lambda_CR_Quad = rep(NA, nBootSamples),
                     Lambda_CR_Lin = rep(NA, nBootSamples),
                     Lambda_Cont_Quad = rep(NA, nBootSamples),
                     Lambda_Cont_Lin = rep(NA, nBootSamples))

ObservedValues <- list(GrowSlope_CR = coef(CompGrowLM)[2],
                       GrowInt_CR = coef(CompGrowLM)[1],
                       GrowSlope_Cont = coef(ContGrowLM)[2],
                       GrowInt_Cont = coef(ContGrowLM)[1],
                       RepSlope = coef(ReproGLM)[2],
                       RepInt = coef(ReproGLM)[1],
                       SeedSlope = coef(SeedGLM)[2],
                       SeedInt = coef(SeedGLM)[1],
                       RecMean = SdlMean,
                       RecSD = SdlSD,
                       Lambda_CR_Quad = CompLambda_Quad_obs,
                       Lambda_CR_Lin = CompLambda_Lin_obs,
                       Lambda_Cont_Quad = ContLambda_Quad_obs,
                       Lambda_Cont_Lin = ContLambda_Lin_obs)

for(i in seq_len(nBootSamples)) {
  
  # Set up resampling vectors
  CRSampler <- sample(1:nCR, nCR, replace = TRUE)
  ContSampler <- sample(1:nCont, nCont, replace = TRUE)
  BigSampler <- sample(1:nBig, nBig, replace = TRUE)
  SurvSampler <- sample(1:nSurv, 1)
  
  BootCRData <- CRData[CRSampler, ]
  BootContData <- ContData[ContSampler, ]
  BootBigData <- BigData[BigSampler, ]
  
  # Growth regressions
  BootCompGrowLM <- lm(Plant_Height13 ~ Plant_Height12, data = rbind(BootCRData,
                                                   BootBigData))
  BootContGrowLM <- lm(Plant_Height13 ~ Plant_Height12, data = rbind(BootContData,
                                                     BootBigData))
  
  # Survival. These are created by resampling the joint posterior distribution
  # for each of our regressions
  
  BootContSurv_Lin <- SurvChainsCont_Lin[SurvSampler, ] %>% unlist()
  BootCompSurv_Lin <- SurvChainsCR_Lin[SurvSampler, ] %>% unlist()
  BootContSurv_Quad <- SurvChainsCont_Quad[SurvSampler, ] %>% unlist()
  BootCompSurv_Quad <- SurvChainsCR_Quad[SurvSampler, ] %>% unlist()
  
  # Recreate reproduction parameters with pooled bootstrap samples
  AllBootData <- rbind(BootCRData, BootContData, BootBigData)
  
  BootReproGLM <- glm(Reproductive ~ Plant_Height13, 
                      data = AllBootData,
                      family = quasibinomial())
  
  BootSeedGLM <- glm(Seeds ~ Plant_Height13,
                     data = AllBootData,
                     family = quasipoisson())
  
  
  BootSdls <- filter(AllBootData, Stage13 == 'S')
  BootSdlMean <- mean(BootSdls$Plant_Height13, na.rm = TRUE)
  BootSdlSD <- sd(BootSdls$Plant_Height13, na.rm = TRUE)
  
  
  BootFecModels <- list(PRep = BootReproGLM,
                        Seeds = BootSeedGLM,
                        EstProb = EstProb,
                        SizeDist = list(Mean = BootSdlMean,
                                        SD = BootSdlSD)) 
  
  # Now, begin constructing kernels. We don't need to re-assign any of the constants
  # so I'll skip that
  GMat_Cont <- CellSize * outer(MeshPts, MeshPts, FUN = GrowFun, Model = BootContGrowLM)
  GMat_Comp <- CellSize * outer(MeshPts, MeshPts, FUN = GrowFun, Model = BootCompGrowLM)
  
  GMat_Cont <- GMat_Cont/matrix(as.vector(apply(GMat_Cont,
                                                2,
                                                sum)),
                                nrow = nMeshPts,
                                ncol = nMeshPts,
                                byrow = TRUE)
  
  GMat_Comp <- GMat_Comp/matrix(as.vector(apply(GMat_Comp,
                                                2,
                                                sum)),
                                nrow = nMeshPts,
                                ncol = nMeshPts,
                                byrow = TRUE)
  
  SVec_Cont_Quad <- SurvFun(MeshPts, BootContSurv_Quad)
  SVec_Comp_Quad <- SurvFun(MeshPts, BootCompSurv_Quad)
  
  SVec_Cont_Lin <- SurvFun(MeshPts, BootContSurv_Lin)
  SVec_Comp_Lin <- SurvFun(MeshPts, BootCompSurv_Lin)
  
  # P, F, and K matrices
  PMat_Cont_Quad <- PMat_Comp_Quad <- matrix(0, nMeshPts, nMeshPts)
  
  PMat_Cont_Lin <- PMat_Comp_Lin <- matrix(0, nMeshPts, nMeshPts)
  
  
  for(j in seq_len(nMeshPts)){
    PMat_Cont_Quad[ ,j] <- GMat_Cont[ ,j] * SVec_Cont_Quad[j]
    PMat_Comp_Quad[ ,j] <- GMat_Comp[ ,j] * SVec_Comp_Quad[j]
    
    PMat_Cont_Lin[ ,j] <- GMat_Cont[ ,j] * SVec_Cont_Lin[j]
    PMat_Comp_Lin[ ,j] <- GMat_Comp[ ,j] * SVec_Comp_Lin[j]
  }
  
  FMat <- CellSize * outer(MeshPts, MeshPts, FUN = FecFun, Models = FecModels)
  
  KMat_Cont_Quad <- PMat_Cont_Quad + FMat
  KMat_Comp_Quad <- PMat_Comp_Quad + FMat
  
  KMat_Cont_Lin <- PMat_Cont_Lin + FMat
  KMat_Comp_Lin <- PMat_Comp_Lin + FMat
  
  # Extract lambdas
  ContEV_Quad <- eigen(KMat_Cont_Quad)
  CompEV_Quad <- eigen(KMat_Comp_Quad)
  ContLambda_Quad_boot <- max(Re(ContEV_Quad$values))
  CompLambda_Quad_boot <- max(Re(CompEV_Quad$values))
  ContEV_Lin <- eigen(KMat_Cont_Lin)
  CompEV_Lin <- eigen(KMat_Comp_Lin)
  ContLambda_Lin_boot <- max(Re(ContEV_Lin$values))
  CompLambda_Lin_boot <- max(Re(CompEV_Lin$values))
  
  # Store values
  OutputValues$GrowSlope_CR[i] <- coef(BootCompGrowLM)[2]
  OutputValues$GrowInt_CR[i] <- coef(BootCompGrowLM)[1]
  OutputValues$GrowSlope_Cont[i] <- coef(BootContGrowLM)[2]
  OutputValues$GrowInt_Cont[i] <- coef(BootContGrowLM)[1]
  OutputValues$RepSlope[i] <- coef(BootReproGLM)[2]
  OutputValues$RepInt[i] <- coef(BootReproGLM)[1]
  OutputValues$SeedSlope[i] <- coef(BootSeedGLM)[2]
  OutputValues$SeedInt[i] <- coef(BootSeedGLM)[1]
  OutputValues$RecMean[i] <- BootSdlMean
  OutputValues$RecSD[i] <- BootSdlSD
  OutputValues$Lambda_CR_Quad[i] <- CompLambda_Quad_boot
  OutputValues$Lambda_CR_Lin[i] <- CompLambda_Lin_boot
  OutputValues$Lambda_Cont_Quad[i] <- ContLambda_Quad_boot
  OutputValues$Lambda_Cont_Lin[i] <- ContLambda_Lin_boot
  
  if(i %% 100 == 0) {
    message(i / 10, '% of data crunched')
  }
  
}


OutputData <- as.data.frame(OutputValues) %>%
  mutate(Boot_Obs = 'Boot')

AllValues <- as.data.frame(ObservedValues) %>%
  mutate(Boot_Obs = 'Observed') %>%
  rbind(OutputData) 

write.csv(AllValues, 'BootStrap_Output_Lonicera.csv',
          row.names = FALSE)

