# Lonicera IPM for unburned controls and CR treatments

rm(list = ls())
graphics.off()

library(dplyr)
library(mgcv)
library(brms)

# Read in data and get rid of extraneous info
AllPlants <- read.csv('Lonicera/IPM/Data/LM_Clean.csv',
                      stringsAsFactors = FALSE)


AllCR <- filter(AllPlants, Treatment == 'CompN')
AllCont <- filter(AllPlants, Treatment == 'ContN')
AllBig <- filter(AllPlants, Treatment == 'AllN')

# Growth Regressions. Try Linear fit and splines for comparison
CompGrowLM <- lm(Plant_Height13 ~ Plant_Height12, data = rbind(AllCR, AllBig))
ContGrowLM <- lm(Plant_Height13 ~ Plant_Height12, data = rbind(AllCont, AllBig))
CompGrowGAM <- gam(Plant_Height13 ~ s(Plant_Height12), data = rbind(AllCR, AllBig))
ContGrowGAM <- gam(Plant_Height13 ~ s(Plant_Height12), data = rbind(AllCont, AllBig))

plot(Plant_Height13 ~ Plant_Height12, data = AllPlants, type = 'n')
points(AllCR$Plant_Height12, AllCR$Plant_Height13, col = 'green', pch = 3)
points(AllCont$Plant_Height12, AllCont$Plant_Height13, col = 'black', pch = 1)
points(AllBig$Plant_Height12, AllBig$Plant_Height13, col = 'red', pch = 2)
abline(CompGrowLM, lty = 2, col = 'green')
abline(ContGrowLM, lty = 2, col = 'black')
xx <- seq(0, max(AllPlants$Plant_Height12, na.rm = TRUE) * 1.2, 1)
lines(xx, predict(CompGrowGAM, data.frame(Plant_Height12 = xx), type = 'response'),
      col = 'green', lty = 1)
lines(xx, predict(ContGrowGAM, data.frame(Plant_Height12 = xx), type = 'response'),
      col = 'black', lty = 1)

summary(CompGrowLM)
summary(ContGrowLM)
summary(CompGrowGAM)
summary(ContGrowGAM)

AIC(CompGrowLM, CompGrowGAM)
AIC(ContGrowLM, ContGrowGAM)

# GAMs all have lower AIC, but the CR GAM makes ludicrous predictions.
# The control one actually seems somewhat realistic, but I'd prefer to use
# the linear fit for simplicity's sake.

# Survival. Logistic regressions w/ linear and polynomial terms to start off
CompSurvLin <- glm(Survival ~ Plant_Height12, data = rbind(AllCR, AllBig), family = binomial())
ContSurvLin <- glm(Survival ~ Plant_Height12, data = rbind(AllCont, AllBig), family = binomial())
CompSurvQuad <- glm(Survival ~ Plant_Height12 + I(Plant_Height12^2),
                    data = rbind(AllCR, AllBig), 
                    family = binomial())
ContSurvQuad <- glm(Survival ~ Plant_Height12 + I(Plant_Height12^2),
                    data = rbind(AllCont, AllBig),
                    family = binomial())

plot(Survival ~ Plant_Height12, data = AllPlants)
lines(xx, predict(CompSurvLin, data.frame(Plant_Height12 = xx), type = 'response'),
      col = 'green', lty = 1)
lines(xx, predict(ContSurvLin, data.frame(Plant_Height12 = xx), type = 'response'),
      col = 'black', lty = 1)
lines(xx, predict(CompSurvQuad, data.frame(Plant_Height12 = xx), type = 'response'),
      col = 'green', lty = 2)
lines(xx, predict(ContSurvQuad, data.frame(Plant_Height12 = xx), type = 'response'),
      col = 'black', lty = 2)


# None of these are working particularly well :( Trying brms next
CompSurvBRM_Quad <- brm(Survival ~ Plant_Height12 + I(Plant_Height12^2), 
                        data = rbind(AllCR, AllBig),
                        family = 'bernoulli',
                        control = list(adapt_delta = 0.97,
                                       max_treedepth = 15),
                        cores = getOption('mc.cores', 4L),
                        iter = 4000,
                        chains = 4,
                        warmup = 1000,
                        save_model = 'Lonicera/IPM/R/CR_Surv_BRM_Quad.stan') %>%
  add_waic()

ContSurvBRM_Quad <-  brm(Survival ~ Plant_Height12 + I(Plant_Height12^2), 
                         data = rbind(AllCont, AllBig),
                         family = 'bernoulli',
                         control = list(adapt_delta = 0.97,
                                        max_treedepth = 15),
                         cores = getOption('mc.cores', 4L),
                         iter = 4000,
                         chains = 4,
                         warmup = 1000,
                         save_model = 'Lonicera/IPM/R/Cont_Surv_BRM_Quad.stan') %>%
  add_waic()

CompSurvBRM_Lin <- brm(Survival ~ Plant_Height12, 
                       data = rbind(AllCR, AllBig),
                       family = 'bernoulli',
                       control = list(adapt_delta = 0.97,
                                      max_treedepth = 15),
                       cores = getOption('mc.cores', 4L),
                       iter = 4000,
                       chains = 4,
                       warmup = 1000,
                       save_model = 'Lonicera/IPM/R/CR_Surv_BRM_Lin.stan') %>%
  add_waic()

ContSurvBRM_Lin <-  brm(Survival ~ Plant_Height12, 
                        data = rbind(AllCont, AllBig),
                        family = 'bernoulli',
                        control = list(adapt_delta = 0.97,
                                       max_treedepth = 15),
                        cores = getOption('mc.cores', 4L),
                        iter = 4000,
                        chains = 4,
                        warmup = 1000,
                        save_model = 'Lonicera/IPM/R/Cont_Surv_BRM_Lin.stan') %>%
  add_waic()


plot(Survival ~ Plant_Height12, data = AllPlants, xlim = c(0, 515))
lines(xx, predict(CompSurvBRM_Quad,
                  data.frame(Plant_Height12 = xx), 
                  type = 'response')[ ,1],
      col = 'green')
lines(xx, predict(ContSurvBRM_Quad,
                  data.frame(Plant_Height12 = xx),
                  type = 'response')[ ,1])
lines(xx, predict(CompSurvBRM_Lin,
                  data.frame(Plant_Height12 = xx),
                  type = 'response')[ ,1],
      col = 'green',
      lty = 2)
lines(xx, predict(ContSurvBRM_Lin, 
                  data.frame(Plant_Height12 = xx),
                  type = 'response')[ ,1],
      lty = 2)

# None of these survival curves look particularly good. I think the bayesian 
# regression is somewhat more realistic, but only because there is mortality. 
# It's predictions of mortality rate are high (likely unrealistically so). 
# I'll try a binomial spline next

SurvCRGAM <- gam(Survival ~ s(Plant_Height12),
                 data = rbind(AllCR, AllBig),
                 family = binomial())

SurvContGAM <- gam(Survival ~ s(Plant_Height12),
                 data = rbind(AllCont, AllBig),
                 family = binomial())

plot(Survival ~ Plant_Height12, data = AllPlants,
     xlim = c(0, 515))
lines(xx, predict(SurvCRGAM, 
                  data.frame(Plant_Height12 = xx),
                  type = 'response'),
      col = 'green')
lines(xx, predict(SurvContGAM,
                  data.frame(Plant_Height12 = xx),
                  type = 'response'))

# Still not any better. I think I'll stick with the 
# quadratic bayesian, but this will likely come up in peer
# review

SurvChainsCR_Quad <- CompSurvBRM_Quad$fit@sim$samples[[1]][1:3] %>%
  data.frame()
SurvParamsCR_Quad <- fixef(CompSurvBRM_Quad)[ ,1]
SurvChainsCont_Quad <- ContSurvBRM_Quad$fit@sim$samples[[1]][1:3] %>%
  data.frame()
SurvParamsCont_Quad <- fixef(ContSurvBRM_Quad)[ ,1]


SurvChainsCR_Lin <- CompSurvBRM_Lin$fit@sim$samples[[1]][1:3] %>%
  data.frame()
SurvParamsCR_Lin <- fixef(CompSurvBRM_Lin)[ ,1]
SurvChainsCont_Lin <- ContSurvBRM_Lin$fit@sim$samples[[1]][1:3] %>%
  data.frame()
SurvParamsCont_Lin <- fixef(ContSurvBRM_Lin)[ ,1]

summary(CompSurvBRM_Quad)
summary(ContSurvBRM_Quad)
summary(CompSurvBRM_Lin)
summary(ContSurvBRM_Lin)

plot(CompSurvBRM_Quad)
ggsave(filename = 'Survival_Regression_Chains_CR.png',
       width = 6, height = 6, unit = 'in',
       path = 'Lonicera/IPM/Figures',
       device = 'png')

plot(ContSurvBRM_Quad)
ggsave(filename = 'Survival_Regression_Chains_Cont.png',
       width = 6, height = 6, unit = 'in',
       path = 'Lonicera/IPM/Figures',
       device = 'png')

plot(CompSurvBRM_Lin)
ggsave(filename = 'Survival_Regression_Chains_Comp_Lin.png',
       width = 6, height = 6, unit = 'in',
       path = 'Lonicera/IPM/Figures',
       device = 'png')

plot(ContSurvBRM_Lin)
ggsave(filename = 'Survival_Regression_Chains_Cont_Lin.png',
       width = 6, height = 6, unit = 'in',
       path = 'Lonicera/IPM/Figures',
       device = 'png')

# Reproduction consists of 4 parameters:
# Pr(Repro) = Reproductive~ Plant_Height13 (t+1) 
# Seeds = Seeds ~ Plant_Height13(t+1)
# Sdl size dist = dnorm(mu, sd)
# est prob = constant

# Pr(Repro)
ReproGLM <- glm(Reproductive ~ Plant_Height13, data = AllPlants, family = binomial())
summary(ReproGLM)

# Regression is underdispersed. Using a quasibinomial instead
ReproGLM <- glm(Reproductive ~ Plant_Height13, data = AllPlants, family = quasibinomial())
summary(ReproGLM)

# Seeds next. Creating a column for seeds by multiplying fruit# by 
# average seeds/fruit from Amibeth's control treatment.
AllPlants$Seeds <- round(AllPlants$Fruit * 2.79)
SeedGLM <- glm(Seeds ~ Plant_Height13, data = AllPlants, family = poisson())
summary(SeedGLM)

# Again, overdispersed. using quasipoisson instead

SeedGLM <- glm(Seeds ~ Plant_Height13, data = AllPlants, family = quasipoisson())
summary(SeedGLM)

# Recruit size distribution. I'll see if we have enough per treatment
# to keep these separate. otherwise, combine them all

table(AllPlants$Stage13, AllPlants$Treatment)
table(AllPlants$Stage12, AllPlants$Treatment)

# We only have 17 total at t+1. Not great even for pooling. 
# Definitely going to combine across treatments.
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


source('Lonicera/IPM/R/Lonicera_Utility_Functions.R')

# Begin setting up IPM 
MinSize <- 0.9 * min(c(AllPlants$Plant_Height12, AllPlants$Plant_Height13), na.rm = TRUE)
MaxSize <- 1.1 * max(c(AllPlants$Plant_Height12, AllPlants$Plant_Height13), na.rm = TRUE)
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

ContEV_Quad <- eigen(KMat_Cont_Quad)
CompEV_Quad <- eigen(KMat_Comp_Quad)
ContLambda_Quad_obs <- max(Re(ContEV_Quad$values))
CompLambda_Quad_obs <- max(Re(CompEV_Quad$values))
ContEV_Lin <- eigen(KMat_Cont_Lin)
CompEV_Lin <- eigen(KMat_Comp_Lin)
ContLambda_Lin_obs <- max(Re(ContEV_Lin$values))
CompLambda_Lin_obs <- max(Re(CompEV_Lin$values))

ContLambda_Quad_obs
CompLambda_Quad_obs
ContLambda_Lin_obs
CompLambda_Lin_obs




