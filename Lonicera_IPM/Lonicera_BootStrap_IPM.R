# Lonicera IPM for unburned controls and CR treatments
# Includes boot strapping and is intended to run on iDiv's
# RStudio server. COmplete model selection process can
# be found at 
# https://github.com/levisc8/Invader_Demography/blob/master/Lonicera/IPM/R/Lonicera_IPM.R

# rm(list = ls())
# graphics.off()

library(dplyr)
library(mgcv)
library(brms)
library(ggplot2)
library(tidyr)
library(stringr)

# Read in data and get rid of extraneous info
AllPlants <- read.csv('Lonicera_IPM/Lonicera_Clean.csv',
                      stringsAsFactors = FALSE)

AllCR <- filter(AllPlants, Treatment == 'CompN')
AllCont <- filter(AllPlants, Treatment == 'ContN')
AllBig <- filter(AllPlants, Treatment == 'AllN')

# Growth Regressions. Model selection process omitted here
CompGrowLM <- lm(HeightNext ~ Height, data = rbind(AllCR, AllBig))
ContGrowLM <- lm(HeightNext ~ Height, data = rbind(AllCont, AllBig))

# Survival. Logistic regressions w/ linear and polynomial terms to start off
CompSurvBRM_Quad <- brm(Survival~ Height + I(Height^2), 
                        data = rbind(AllCR, AllBig),
                        family = 'bernoulli',
                        control = list(adapt_delta = 0.97,
                                       max_treedepth = 15),
                        cores = getOption('mc.cores', 2L),
                        iter = 4000,
                        chains = 4,
                        warmup = 1000,
                        save_model = 'Lonicera_IPM/CR_Surv_BRM_Quad.stan') %>%
  add_waic()

ContSurvBRM_Quad <-  brm(Survival~ Height + I(Height^2), 
                         data = rbind(AllCont, AllBig),
                         family = 'bernoulli',
                         control = list(adapt_delta = 0.97,
                                        max_treedepth = 15),
                         cores = getOption('mc.cores', 2L),
                         iter = 4000,
                         chains = 4,
                         warmup = 1000,
                         save_model = 'Lonicera_IPM/Cont_Surv_BRM_Quad.stan') %>%
  add_waic()

CompSurvBRM_Lin <- brm(Survival~ Height, 
                       data = rbind(AllCR, AllBig),
                       family = 'bernoulli',
                       control = list(adapt_delta = 0.97,
                                      max_treedepth = 15),
                       cores = getOption('mc.cores', 2L),
                       iter = 4000,
                       chains = 4,
                       warmup = 1000,
                       save_model = 'Lonicera_IPM/CR_Surv_BRM_Lin.stan') %>%
  add_waic()

ContSurvBRM_Lin <-  brm(Survival~ Height, 
                        data = rbind(AllCont, AllBig),
                        family = 'bernoulli',
                        control = list(adapt_delta = 0.97,
                                       max_treedepth = 15),
                        cores = getOption('mc.cores', 2L),
                        iter = 4000,
                        chains = 4,
                        warmup = 1000,
                        save_model = 'Lonicera_IPM/Cont_Surv_BRM_Lin.stan') %>%
  add_waic()

# Store chains to resample during bootstrapping
thin <- sample(1001:4000, 2000, replace = FALSE)

SurvChainsCR_Quad <- CompSurvBRM_Quad$fit@sim$samples[[1]][1:4] %>%
  data.frame() %>%
  .[thin, ] %>%
  setNames(c(
    'Intercept', 'LinearTerm', 'QuadraticTerm', 'Likelihood'
  )) %>%
  mutate(Treatment = 'CR') %>%
  select(Treatment, Intercept, LinearTerm, QuadraticTerm, Likelihood) 

SurvParamsCR_Quad <- fixef(CompSurvBRM_Quad)[ ,1]

SurvChainsCont_Quad <- ContSurvBRM_Quad$fit@sim$samples[[1]][1:4] %>%
  data.frame() %>%
  .[thin, ] %>%
  setNames(c(
    'Intercept', 'LinearTerm', 'QuadraticTerm', 'Likelihood'
  )) %>%
  mutate(Treatment = 'Control')%>%
  select(Treatment, Intercept, LinearTerm, QuadraticTerm, Likelihood) 

SurvParamsCont_Quad <- fixef(ContSurvBRM_Quad)[ ,1]


SurvChainsCR_Lin <- CompSurvBRM_Lin$fit@sim$samples[[1]][1:3] %>%
  data.frame() %>%
  .[thin, ] %>%
  setNames(c(
    'Intercept', 'LinearTerm', 'Likelihood'
  )) %>%
  mutate(QuadraticTerm = NA_real_,
         Treatment = 'CR')%>%
  select(Treatment, Intercept, LinearTerm, QuadraticTerm, Likelihood) 

SurvParamsCR_Lin <- fixef(CompSurvBRM_Lin)[ ,1]
SurvChainsCont_Lin <- ContSurvBRM_Lin$fit@sim$samples[[1]][1:3] %>%
  data.frame() %>%
  .[thin, ] %>%
  setNames(c(
    'Intercept', 'LinearTerm', 'Likelihood'
  )) %>%
  mutate(QuadraticTerm = NA_real_,
         Treatment = 'Control') %>%
  select(Treatment, Intercept, LinearTerm, QuadraticTerm, Likelihood) 
SurvParamsCont_Lin <- fixef(ContSurvBRM_Lin)[ ,1]

AllSurvChains <- rbind(SurvChainsCont_Lin, SurvChainsCont_Quad,
                       SurvChainsCR_Lin, SurvChainsCR_Quad)

write.csv(AllSurvChains, file = 'Lonicera_IPM/Lonicera_Survival_Model_Chains.csv',
          row.names = FALSE)


# Reproduction consists of 4 parameters:
# Pr(Repro) = Reproductive ~ size (t+1) 
# Seeds = Seeds ~ size(t+1)
# Sdl size dist = dnorm(mu, sd)
# est prob = constant

# Pr(Repro)
ReproGLM <- glm(Reproductive ~ HeightNext, 
                data = AllPlants, 
                family = binomial())

# Regression is underdispersed. Using a quasibinomial instead

ReproGLM <- glm(Reproductive ~ HeightNext,
                data = AllPlants,
                family = quasibinomial())

# Seeds next. Creating a column for seeds by multiplying fruit# by 
# average seeds/fruit from a colleague's control treatment. See manuscript
# for further details

AllPlants$Seeds <- round(AllPlants$Fruit * 2.79)
SeedGLM <- glm(Seeds ~ HeightNext, data = AllPlants, family = quasipoisson())

# Recruit size distribution. 
Seedlings <- filter(AllPlants, StageNext == 'S')
SdlMean <- mean(Seedlings$HeightNext, na.rm = TRUE)
SdlSD <- sd(Seedlings$HeightNext, na.rm = TRUE)


# Data hard coded from Rae's 2011-12 population data in unburned controls
# (density x fire x invasive removal experiment)
EstProb <- 0.003562864


FecModels <- list(PRep = ReproGLM,
                  Seeds = SeedGLM,
                  EstProb = EstProb,
                  SizeDist = list(Mean = SdlMean,
                                  SD = SdlSD))  


source('Lonicera_IPM/Lonicera_Utility_Functions.R')

# Begin setting up IPM 
MinSize <- 0.9 * min(c(AllPlants$Height, AllPlants$HeightNext), na.rm = TRUE)
MaxSize <- 1.1 * max(c(AllPlants$Height, AllPlants$HeightNext), na.rm = TRUE)
nMeshPts <- 500
Boundaries <- MinSize + c(0:nMeshPts) * (MaxSize - MinSize) / nMeshPts
MeshPts <- 0.5 * (Boundaries[1:nMeshPts] + Boundaries[2:(nMeshPts + 1)])
CellSize <- MeshPts[2] - MeshPts[1]

# Construct the P kernel. First, make growth transition matrix and correct for eviction.
# then create survival vector and multiply it in
GMat_Cont <- CellSize * outer(MeshPts, MeshPts, FUN = GrowFun, Model = ContGrowLM)
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

# Extract eigenvalues

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

# The results of the bootstrapping are saved as Lonicera_Bootstrap_Ouput.csv,
# so the lines below will be commented out. However, if you have the computational
# resources to run them, feel free to uncomment them and re-run the procedure.
# The commented code is provided to facilitate understanding how the procedure
# was run. 
# 
# CRData <- filter(AllPlants, Treatment == 'CompN')
# ContData <- filter(AllPlants, Treatment == 'ContN')
# BigData <- filter(AllPlants, Treatment == 'AllN')
# 
# nCR <- dim(CRData)[1]
# nCont <- dim(ContData)[1]
# nBig <- dim(BigData)[1]
# nSurv <- dim(SurvChainsCont_Quad)[1]
# 
# nBootSamples <- 1000
# 
# OutputValues <- list(GrowSlope_CR = rep(NA, nBootSamples),
#                      GrowInt_CR = rep(NA, nBootSamples),
#                      GrowSlope_Cont = rep(NA, nBootSamples),
#                      GrowInt_Cont = rep(NA, nBootSamples),
#                      RepSlope = rep(NA, nBootSamples),
#                      RepInt = rep(NA, nBootSamples),
#                      SeedSlope = rep(NA, nBootSamples),
#                      SeedInt = rep(NA, nBootSamples),
#                      RecMean = rep(NA, nBootSamples),
#                      RecSD = rep(NA, nBootSamples),
#                      Lambda_CR_Quad = rep(NA, nBootSamples),
#                      Lambda_CR_Lin = rep(NA, nBootSamples),
#                      Lambda_Cont_Quad = rep(NA, nBootSamples),
#                      Lambda_Cont_Lin = rep(NA, nBootSamples))
# 
# ObservedValues <- list(GrowSlope_CR = coef(CompGrowLM)[2],
#                        GrowInt_CR = coef(CompGrowLM)[1],
#                        GrowSlope_Cont = coef(ContGrowLM)[2],
#                        GrowInt_Cont = coef(ContGrowLM)[1],
#                        RepSlope = coef(ReproGLM)[2],
#                        RepInt = coef(ReproGLM)[1],
#                        SeedSlope = coef(SeedGLM)[2],
#                        SeedInt = coef(SeedGLM)[1],
#                        RecMean = SdlMean,
#                        RecSD = SdlSD,
#                        Lambda_CR_Quad = CompLambda_Quad_obs,
#                        Lambda_CR_Lin = CompLambda_Lin_obs,
#                        Lambda_Cont_Quad = ContLambda_Quad_obs,
#                        Lambda_Cont_Lin = ContLambda_Lin_obs)
# 
# for(i in seq_len(nBootSamples)) {
  
  # Set up resampling vectors. This generates a vector by sampling 
  # rows of the data frame with replacement, then subsets the treatment
  # data using that vector.
  
  # CRSampler <- sample(1:nCR, nCR, replace = TRUE)
  # ContSampler <- sample(1:nCont, nCont, replace = TRUE)
  # BigSampler <- sample(1:nBig, nBig, replace = TRUE)
  # SurvSampler <- sample(1:nSurv, 1)
  # 
  # BootCRData <- CRData[CRSampler, ]
  # BootContData <- ContData[ContSampler, ]
  # BootBigData <- BigData[BigSampler, ]
  # 
  # # Growth regressions
  # BootCompGrowLM <- lm(HeightNext ~ Height, data = rbind(BootCRData,
  #                                                  BootBigData))
  # BootContGrowLM <- lm(HeightNext ~ Height, data = rbind(BootContData,
  #                                                    BootBigData))
  # 
  # # Survival. These are created by resampling the joint posterior distribution
  # # for each of our regressions rather than re-fitting the model (which would be
  # # extremely time consuming).
  # 
  # BootContSurv_Lin <- SurvChainsCont_Lin[SurvSampler, ] %>% unlist()
  # BootCompSurv_Lin <- SurvChainsCR_Lin[SurvSampler, ] %>% unlist()
  # BootContSurv_Quad <- SurvChainsCont_Quad[SurvSampler, ] %>% unlist()
  # BootCompSurv_Quad <- SurvChainsCR_Quad[SurvSampler, ] %>% unlist()
  # 
  # # Recreate reproduction parameters with pooled bootstrap samples
  # AllBootData <- rbind(BootCRData, BootContData, BootBigData)
  # 
  # BootReproGLM <- glm(Reproductive ~ HeightNext, 
  #                     data = AllBootData,
  #                     family = quasibinomial())
  # 
  # BootSeedGLM <- glm(Seeds ~ HeightNext,
  #                    data = AllBootData,
  #                    family = quasipoisson())
  # 
  # Use the boot strapped data set to generate a recruit size distribution
  
#   BootSdls <- filter(AllBootData, StageNext == 'S')
#   BootSdlMean <- mean(BootSdls$HeightNext, na.rm = TRUE)
#   BootSdlSD <- sd(BootSdls$HeightNext, na.rm = TRUE)
#   
#   
#   BootFecModels <- list(PRep = BootReproGLM,
#                         Seeds = BootSeedGLM,
#                         EstProb = EstProb,
#                         SizeDist = list(Mean = BootSdlMean,
#                                         SD = BootSdlSD)) 
#   
#   # Now, begin constructing kernels. We don't need to re-assign any of the 
#   # constants used to construct the models so I'll skip that
#   GMat_Cont <- CellSize * outer(MeshPts, MeshPts, FUN = GrowFun, Model = BootContGrowLM)
#   GMat_Comp <- CellSize * outer(MeshPts, MeshPts, FUN = GrowFun, Model = BootCompGrowLM)
#   
#   # Correct for eviction in the growth kernels
#   GMat_Cont <- GMat_Cont/matrix(as.vector(apply(GMat_Cont,
#                                                 2,
#                                                 sum)),
#                                 nrow = nMeshPts,
#                                 ncol = nMeshPts,
#                                 byrow = TRUE)
#   
#   GMat_Comp <- GMat_Comp/matrix(as.vector(apply(GMat_Comp,
#                                                 2,
#                                                 sum)),
#                                 nrow = nMeshPts,
#                                 ncol = nMeshPts,
#                                 byrow = TRUE)
#   
#   SVec_Cont_Quad <- SurvFun(MeshPts, BootContSurv_Quad)
#   SVec_Comp_Quad <- SurvFun(MeshPts, BootCompSurv_Quad)
#   
#   SVec_Cont_Lin <- SurvFun(MeshPts, BootContSurv_Lin)
#   SVec_Comp_Lin <- SurvFun(MeshPts, BootCompSurv_Lin)
#   
#   # P, F, and K matrices
#   PMat_Cont_Quad <- PMat_Comp_Quad <- matrix(0, nMeshPts, nMeshPts)
#   
#   PMat_Cont_Lin <- PMat_Comp_Lin <- matrix(0, nMeshPts, nMeshPts)
#   
#   # Generate P Matrix
#   
#   for(j in seq_len(nMeshPts)){
#     PMat_Cont_Quad[ ,j] <- GMat_Cont[ ,j] * SVec_Cont_Quad[j]
#     PMat_Comp_Quad[ ,j] <- GMat_Comp[ ,j] * SVec_Comp_Quad[j]
#     
#     PMat_Cont_Lin[ ,j] <- GMat_Cont[ ,j] * SVec_Cont_Lin[j]
#     PMat_Comp_Lin[ ,j] <- GMat_Comp[ ,j] * SVec_Comp_Lin[j]
#   }
#   
#   # Generate the F Matrix. This is the same for both treatments
#   FMat <- CellSize * outer(MeshPts, MeshPts, FUN = FecFun, Models = FecModels)
#   
#   KMat_Cont_Quad <- PMat_Cont_Quad + FMat
#   KMat_Comp_Quad <- PMat_Comp_Quad + FMat
#   
#   KMat_Cont_Lin <- PMat_Cont_Lin + FMat
#   KMat_Comp_Lin <- PMat_Comp_Lin + FMat
#   
#   # Extract lambdas
#   ContEV_Quad <- eigen(KMat_Cont_Quad)
#   CompEV_Quad <- eigen(KMat_Comp_Quad)
#   ContLambda_Quad_boot <- max(Re(ContEV_Quad$values))
#   CompLambda_Quad_boot <- max(Re(CompEV_Quad$values))
#   ContEV_Lin <- eigen(KMat_Cont_Lin)
#   CompEV_Lin <- eigen(KMat_Comp_Lin)
#   ContLambda_Lin_boot <- max(Re(ContEV_Lin$values))
#   CompLambda_Lin_boot <- max(Re(CompEV_Lin$values))
#   
#   # Store values
#   OutputValues$GrowSlope_CR[i] <- coef(BootCompGrowLM)[2]
#   OutputValues$GrowInt_CR[i] <- coef(BootCompGrowLM)[1]
#   OutputValues$GrowSlope_Cont[i] <- coef(BootContGrowLM)[2]
#   OutputValues$GrowInt_Cont[i] <- coef(BootContGrowLM)[1]
#   OutputValues$RepSlope[i] <- coef(BootReproGLM)[2]
#   OutputValues$RepInt[i] <- coef(BootReproGLM)[1]
#   OutputValues$SeedSlope[i] <- coef(BootSeedGLM)[2]
#   OutputValues$SeedInt[i] <- coef(BootSeedGLM)[1]
#   OutputValues$RecMean[i] <- BootSdlMean
#   OutputValues$RecSD[i] <- BootSdlSD
#   OutputValues$Lambda_CR_Quad[i] <- CompLambda_Quad_boot
#   OutputValues$Lambda_CR_Lin[i] <- CompLambda_Lin_boot
#   OutputValues$Lambda_Cont_Quad[i] <- ContLambda_Quad_boot
#   OutputValues$Lambda_Cont_Lin[i] <- ContLambda_Lin_boot
#   
#   if(i %% 100 == 0) {
#     message(i / 10, '% of data crunched')
#   }
#   
# }

# Save all of the computed values
# OutputData <- as.data.frame(OutputValues) %>%
#   mutate(Boot_Obs = 'Boot')
# 
# AllValues <- as.data.frame(ObservedValues) %>%
#   mutate(Boot_Obs = 'Observed') %>%
#   rbind(OutputData) 
# 
# write.csv(AllValues, 'Lonicera_IPM/BootStrap_Output_Lonicera.csv',
#           row.names = FALSE)
# 

# Now, to make the figures. This next bit does a bunch of re-structuring of the
# data so that it's ggplot2-friendly

AllValues <- read.csv('Lonicera_IPM/Lonicera_Bootstrap_Output.csv',
                      stringsAsFactors = FALSE)

AllData <- AllValues %>%
  gather(key = 'Variable', value = 'Value', -Boot_Obs) %>%
  mutate(Trt = vapply(.$Variable,
                      FUN = function(x) str_split(x, '_')[[1]][2],
                      FUN.VALUE = ''),
         SurvModel = vapply(.$Variable,
                            FUN = function(x) str_split(x, '_')[[1]][3],
                            FUN.VALUE = '')) %>%
  group_by(Variable, Trt, SurvModel) %>%
  arrange(desc(Value)) %>%
  summarise(obs = Value[Boot_Obs == 'Observed'],
            UpCI = Value[25],
            LoCI = Value[975]) %>%
  ungroup() %>%
  mutate(Variable = vapply(.$Variable,
                           FUN = function(x) str_split(x, '_')[[1]][1],
                           FUN.VALUE = ''))

SurvData <- tibble(Variable = rep(c('SurvInt', 'SurvSlope',
                                      'SurvInt', 'SurvSlope', 'SurvSlope2'), 2),
                       Trt = c(rep('Cont', 5),
                               rep('CR', 5)),
                       SurvModel = c(rep('Lin', 2),
                                     rep('Quad', 3),
                                     rep('Lin', 2),
                                     rep('Quad', 3)))


PlotData <- rbind(fixef(ContSurvBRM_Lin),
                  fixef(ContSurvBRM_Quad),
                  fixef(CompSurvBRM_Lin),
                  fixef(CompSurvBRM_Quad)) %>%
  as_tibble() %>%
  .[ ,-2] %>%
  setNames(c('obs', 'LoCI', 'UpCI')) %>%
  select(obs, UpCI, LoCI) %>%
  cbind(SurvData, .) %>%
  rbind(AllData)
write.csv(PlotData, file = 'Lonicera_IPM/Lonicera_Summarized_Output.csv',
          row.names = FALSE)

# Make labels for facet_wrap()
PlotData$Trt[is.na(PlotData$Trt)] <- 'Pooled'
PlotData$Trt[PlotData$Trt == 'Cont'] <- 'Control'
PlotData$Variable[PlotData$Variable == 'GrowInt'] <- "paste(italic(g(y,x)), ' Intercept')"
PlotData$Variable[PlotData$Variable == 'GrowSlope'] <- "paste(italic(g(y,x)), ' Slope')"
PlotData$Variable[PlotData$Variable == 'RecMean'] <- "paste(mu, ' Recruit Size Mean' )"
PlotData$Variable[PlotData$Variable == 'RecSD'] <- "paste(sigma, ' Recruit Size SD')"
PlotData$Variable[PlotData$Variable == 'RepInt'] <- "paste( italic(f[p](x)), ' Intercept')"
PlotData$Variable[PlotData$Variable == 'RepSlope'] <- "paste(italic(f[p](x)), ' Slope')"
PlotData$Variable[PlotData$Variable == 'SeedInt'] <- "paste(italic(f[s](x)), ' Intercept')"
PlotData$Variable[PlotData$Variable == 'SeedSlope'] <- "paste(italic(f[s](x)),' Slope')"
PlotData$Variable[PlotData$Variable == 'SurvInt'] <- "paste(italic(s(x)),' Intercept')"
PlotData$Variable[PlotData$Variable == 'SurvSlope'] <- "paste(italic(s(x)),' Linear Term')"
PlotData$Variable[PlotData$Variable == 'SurvSlope2'] <- "paste(italic(s(x)),' Quadratic Term')"
PlotData$Variable[PlotData$Variable == 'Lambda'] <- 'lambda'

# Create data set for surivival insensitive vital rates
PlotData2 <- filter(PlotData, is.na(SurvModel))

LM_SurvIns_Plot <- ggplot(data = PlotData2,
                          aes(x = Trt)) + 
  geom_point(aes(y = obs, 
                 color = Trt),
             size = 4.5) + 
  geom_linerange(aes(ymin = LoCI,
                     ymax = UpCI,
                     color = Trt),
                 size = 1.25) + 
  facet_wrap(~Variable,
             scales = 'free',
             labeller = label_parsed) +
  theme(panel.background = element_blank(), # Set up plot thematic elements
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.87, 0.1), # legend in top bottom right corner
        legend.title = element_text(size = 22), # make legend text bigger
        legend.text = element_text(size = 17),
        legend.background = element_blank(),
        axis.title.y = element_text(size = 22,
                                    margin = margin(t = 0, # pad the axis title with a bit of whitespace
                                                    r = 20,
                                                    b = 0,
                                                    l = 0)),
        axis.title.x = element_text(size = 22,
                                    margin = margin(t = 20,
                                                    r = 0,
                                                    b = 10,
                                                    l = 0)),
        axis.line = element_line(colour = 'black', # make axis lines a bit bolder
                                 size = 2),
        axis.ticks = element_line(size = 1.2),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = 'white')) +
  scale_color_manual('A\nTreatment',
                     breaks = c('Control', 'CR', 'Pooled'),
                     values = c('black', 'green', 'blue')) + 
  scale_y_continuous('Observed Value + Confidence Intervals') + 
  scale_x_discrete('Treatment (if applicable)') 

LM_SurvIns_Plot

ggsave(filename = 'Lonicera_Survival_Insensitive_Vital_Rates.png',
       device = 'png',
       width = 10,
       height = 8,
       units = 'in',
       path = 'Lonicera_IPM/')

# Create two dummy rows to occupy the right side of the 
# x-axis for SurvSlope^2
DummyRows <- tibble(Variable = rep("paste(italic(s(x)),' Quadratic Term')", 2),
                    Trt = rep(NA, 2),
                    SurvModel = rep(NA, 2),
                    obs = rep(NA, 2),
                    UpCI = rep(NA_real_, 2),
                    LoCI = rep(NA_real_, 2))

PlotData1 <- filter(PlotData, !is.na(SurvModel)) %>%
  rbind(., DummyRows) %>%
  mutate(Dummy = paste(.$Trt, .$SurvModel, sep = '-')) 
PlotData1$Dummy <- gsub('NA-NA', 'Legend', PlotData1$Dummy) 
PlotData1$SurvModel[PlotData1$SurvModel == 'Lin'] <- 'Linear'
PlotData1$SurvModel[PlotData1$SurvModel == 'Quad'] <- 'Quadratic'

# Create horizontal line for lambda = 1
PlotData1$Facet <- relevel(as.factor(PlotData1$Variable), ref = 'lambda')
PlotData1$yints <- NA
PlotData1$yints[PlotData1$Variable == 'lambda'] <- 1

PlotData1 %>%
  as_tibble() %>%
  ggplot(data = .,
         aes(x = Dummy)) + 
  geom_point(aes(y = obs, 
                 shape = SurvModel,
                 color = Trt),
             size = 4.5) + 
  geom_linerange(aes(ymax = UpCI,
                     ymin = LoCI,
                     color = Trt),
                 size = 1.25) +
  geom_hline(data = PlotData1,
             aes(yintercept = yints),
             color = 'grey50',
             alpha = 0.7,
             linetype = 'dashed',
             size = 2, 
             na.rm = TRUE,
             show.legend = FALSE) + 
  facet_wrap(~Facet,
             scales = 'free',
             labeller = label_parsed) + 
  theme(panel.background = element_blank(), # Set up plot thematic elements
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.91, 0.2), # legend in top bottom right corner
        legend.title = element_text(size = 20), # make legend text bigger
        legend.text = element_text(size = 17),
        legend.background = element_blank(),
        axis.title.y = element_text(size = 22,
                                    margin = margin(t = 0, # pad the axis title with a bit of whitespace
                                                    r = 20,
                                                    b = 0,
                                                    l = 0)),
        axis.title.x = element_text(size = 22,
                                    margin = margin(t = 20,
                                                    r = 0,
                                                    b = 10,
                                                    l = 0)),
        axis.line = element_line(colour = 'black', # make axis lines a bit bolder
                                 size = 2),
        axis.ticks = element_line(size = 1.2),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = 'white')) +
  scale_y_continuous('Observed Value + Confidence Intervals') + 
  scale_x_discrete('Treatment and Survival Model') +
  labs(color = 'B\nTreatment', shape = 'Survival Model') +
  scale_color_manual(breaks = c('Control', 'CR'),
                     values = c('black', 'green')) + 
  scale_shape_manual(breaks = c('Quadratic', 'Linear'),
                     values = c(17,19))

ggsave(filename = 'Lonicera_Survival_Sensitive_Vital_Rates.png',
       device = 'png',
       width = 12,
       height = 8,
       units = 'in',
       path = 'Lonicera_IPM/')

