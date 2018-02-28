# Euonymus IPM with bootstrapping. This script is intended to run
# on iDiv's RStudio server and/or the EVE cluster at iDiv/UFZ. It will
# likely not run as is on a local machine. I've also eliminated a lot
# the models that were fitted for comparison as those can be found
# in R/Euonymus_Exporation_and_IPM.R

rm(list = ls()) 
graphics.off()

library(dplyr)
library(magrittr)
library(stringr)
library(mgcv)
library(brms)
# setwd('~/EuoAla_RS_Folder')
setwd('Euonymus/IPM/EuoAla_RS_Folder')

source('IPM_Functions_Euonymus.R')  


# load the data
RAs<-read.csv("EA_RA_Clean.csv")
AllPlants1 <- read.csv("EA_Clean.csv") 

AllPlants1$incr<-AllPlants1$Plant_Height15-AllPlants1$Plant_Height14

# Generate data sets for CR and control treatments

AllSmall <- filter(AllPlants1, Treatment != 'All')
ControlGrowData <- filter(AllPlants1, Treatment == 'Control' | 
                            Treatment == 'All')
CrGrowData <- filter(AllPlants1, Treatment == 'Comp' |
                       Treatment == 'All')

# Growth---------------------------------------------------------------------------------------
# Prior investigation indicated that a GAM provided a better description 
# of the growth dynamics

ControlGam <- gam(Plant_Height15 ~ s(Plant_Height14), data = ControlGrowData)
CRGam <- gam(Plant_Height15 ~ s(Plant_Height14), data = CrGrowData)

# Survival across treatments----------------------------------------------------------------------
# As with above, model selection procedures have been omitted in favor of brevity.
# Those procedures can similarly be found in the script intended for smaller machines

BrmSurvQuad_Cont <- brm(Survival ~ Plant_Height14 + I(Plant_Height14^2), 
                          data = ControlGrowData,
                          family = 'bernoulli',
                          iter = 6000,
                          warmup = 2000,
                          cores = getOption('mc.cores', 4L),
                          control = list(adapt_delta = 0.97,
                                         max_treedepth = 15)) 

BrmSurvQuad_CR <- brm(Survival ~ Plant_Height14 + I(Plant_Height14^2), 
                          data = CrGrowData,
                          family = 'bernoulli',
                          iter = 6000,
                          warmup = 2000,
                          cores = getOption('mc.cores', 4L),
                          control = list(adapt_delta = 0.97,
                                         max_treedepth = 15)) 
 
# Store paremeters and chains for later     
SurvParamsCont <- fixef(BrmSurvQuad_Cont)[ ,1] %>%
  setNames(c('Int', 'Lin', 'Quad'))

SurvParamsCR <- fixef(BrmSurvQuad_CR)[ ,1] %>%
  setNames(c('Int', 'Lin', 'Quad'))

thin <- seq(1, 6000, length.out = 1000) %>% round()

CrSurvChains <- BrmSurvQuad_CR$fit@sim$samples[[1]][1:3] %>%
  data.frame() %>%
  .[thin, ]

ContSurvChains <- BrmSurvQuad_Cont$fit@sim$samples[[1]][1:3] %>%
  data.frame() %>%
  .[thin, ]

# Fecundity-------------------------------------------------------------------------------
# first, we will estimate fecundity from our sample of 18 RAs
#subsetting out two outliers with small Hts but lots of fruits
RAs2 <- filter(RAs, Fruit < 2500 & Ht < max(RAs$Ht))

# quasipoisson
fecquasi <- glm(Fruit ~ Ht,
                data = RAs2,
                family = quasipoisson())
summary(fecquasi)

#Pr(Reproductive) - Logistic regression ----------------------------------
AllPlants1$Repro <- ifelse(AllPlants1$Stage15 == "RA", 1, 0)

Regression.Data <- filter(AllPlants1, Survival != "NA")

Repro.Glm <- glm(Repro ~ Plant_Height15,
                 data = Regression.Data,
                 family = binomial())

AllPlants1$Fruits15 <- exp(coefficients(fecquasi)[1] + 
                            coefficients(fecquasi)[2] *AllPlants1$Plant_Height15)
AllPlants1$Fruits15[AllPlants1$Stage15 != "RA"] <- NA
summary(AllPlants1$Fruits15)

# Recruit size distribution----------------------------------------------------------------------
sdls <- filter(AllPlants1, Stage14 == "SDL")
Sdl.mean <- mean(sdls$Plant_Height14, na.rm = TRUE)
Sdl.SD <- sd(sdls$Plant_Height14, na.rm = TRUE)

# Discrete parameters------------------------------------

# Establishment probability from Brand et al 2012. I'm taking the average of all cultivars from both 
# time replications. seed bank germination is the difference between cumulative germination in years 2 and 3
# i.e. cumulative germ yr3 - cumulative germ yr 2. 

G2 <- mean(c(34.8,.3,.3,28.5))/100
G3 <- mean(c((37.8-34.8)/37.8,(1-.3)/1,(.3-.3)/.3,(32.3-28.5)/32.3))/100

Establishment.prob2 <- G2 
Establishment.prob3 <- G3 

# Create fecundity parameters data frame-------------------------------
f.params<-data.frame(prob.repro.int=as.numeric(coefficients(Repro.Glm)[1]),
                     prob.repro.slope=as.numeric(coefficients(Repro.Glm)[2]),
                     recruit.size.mean=Sdl.mean,
                     recruit.size.sd=Sdl.SD,
                     seed.int=as.numeric(coefficients(fecquasi)[1]),
                     seed.slope=coefficients(fecquasi)[2],
                     E2=Establishment.prob2,
                     E3=Establishment.prob3)

# Build P Matrix and F Vectors-----------------------------------
# the size range must extend beyond the limits of the data
min.size <- min(AllPlants1$Plant_Height14, na.rm = TRUE) * .8
max.size <- max(AllPlants1$Plant_Height14, na.rm = TRUE) * 1.2
S <- 500 # Number of cells in matrix  

# matrix variables 
b <-  min.size + c(0:S)*(max.size - min.size)/S  # boundary points of mesh cells
Y <- 0.5* (b[1:S] + b[2:(S+1)])  # mid points of mesh cells 
h <- Y[2]-Y[1]  # cell widths

# NOTES SPECIFIC TO EUONYMUS
# The fecundity matrix for this species is actually more like 2 vectors that we'll wrap around the growth matrix
# adults at time T don't produce seedlings at T+1 because of complex dormancy, so we have a seedbank vector
# for each year (2 and 3) to account for seeds that emerge from each of those. I'll be making those vectors right now
# and then binding them onto the P matrix later.
# Because we conducted post-reproductive censuses at time t, plants must survive in 
# order to reproduce to t+1. Thus, through survival, there is a treatment effect
# in fecundity.

SB2.Go <- SB.Go(Y, params = f.params,
                SB.Year = 2)
SB3.Go <- SB.Go(Y, params = f.params,
                SB.Year = 3)

SB2.Emerge <- SB.Emerge(Y, params = f.params)
SB3.Emerge <- rep(0, S)
Discrete.Corner <- matrix(c(0, 0, 1, 0), 
                          nrow = 2, byrow = TRUE)

# I'm also going to try correcting for eviction using a different
# method that Maria Paniw recommended. This re-scales all of the columns
# based on the column sums.

Gmat_Cont_2 <- h * outer(Y, Y, 
                         FUN = grow.prob,
                         model = ControlGam)
Gmat_CR_2 <- h * outer(Y, Y,
                       FUN = grow.prob,
                       model = CRGam)

Gmat_Cont_2 <- Gmat_Cont_2/matrix(as.vector(apply(Gmat_Cont_2,
                                                  2,
                                                  sum)),
                                  nrow = S,
                                  ncol = S,
                                  byrow = TRUE)

Gmat_CR_2 <- Gmat_CR_2/matrix(as.vector(apply(Gmat_CR_2,
                                              2,
                                              sum)),
                              nrow = S,
                              ncol = S,
                              byrow = TRUE)


# combine growth and survival into P matrix------------------------------------
# 5b. Create a matrix with survival probabilities on the diagonal,
# use matrix multiplication to combine it with the growth matrix

S_Cont <- predict.surv(Y, BrmSurvQuad_Cont)
S_CR <- predict.surv(Y, BrmSurvQuad_CR)

Pmat_Cont <- Pmat_CR <- matrix(0, S, S)

# make a P matrix (growth and survival)
for(i in seq_len(S)){
  Pmat_Cont[ ,i] <- Gmat_Cont_2[ ,i] * S_Cont[i]
  Pmat_CR[ ,i] <- Gmat_CR_2[ ,i] * S_CR[i]
}

# build full K kernel-----------------------------------
Kstart <- rbind(Discrete.Corner, cbind(SB3.Emerge, SB2.Emerge))
K_Cont <- cbind(Kstart, rbind(SB3.Go, SB2.Go, Pmat_Cont))
K_CR <- cbind(Kstart, rbind(SB3.Go, SB2.Go, Pmat_CR))

# Calculate lambdas --------------------------------
eigen_Cont <- eigen(K_Cont)
lambda_Cont_obs <- max(Re(eigen_Cont$values))

eigen_CR <- eigen(K_CR)
lambda_CR_obs <- max(Re(eigen_CR$values))

# Bootstrapping --------------------
# Now, we need to bootstrap everything (hooray!!!).
# I doubt there will be significant differences, but you never
# know. 

# Initialize data storage areas and underlying data sets.
# I am not sure if it makes sense to boot strap the seed production
# data set because we have so few individuals. Re-building
# a regression from the values we predicted above doesn't make
# much sense to me either. Finally, I have no idea how to store
# a GAM in a way that makes sense, so I won't store those
# parameters for now. However, this should be fairly
# quick to re-run, so hopefully we can add that in
# when I've had a chance to explore it some more. 



CRData <- filter(AllPlants1, Treatment == 'Comp')
ContData <- filter(AllPlants1, Treatment == 'Control')
BigData <- filter(AllPlants1, Treatment == 'All')

nCR <- dim(CRData)[1]
nCont <- dim(ContData)[1]
nBig <- dim(BigData)[1]
nBootSamples <- 1000

OutputValues <- list(RepSlope = rep(NA, nBootSamples),
                     RepInt = rep(NA, nBootSamples),
                     RecMean = rep(NA, nBootSamples),
                     RecSD = rep(NA, nBootSamples),
                     Lambda_Cont = rep(NA, nBootSamples),
                     Lambda_CR = rep(NA, nBootSamples))

for(i in seq_len(nBootSamples)) {

  # Set up vectors and resample data sets and survival parameter posterior
  CrSampler <- sample(1:nCR, nCR, replace = TRUE)
  ContSampler <- sample(1:nCont, nCont, replace = TRUE)
  BigSampler <- sample(1:nBig, nBig, replace = TRUE)
  SurvSampler <- sample(1:1000, 1)
  
  BootCRData <- CRData[CrSampler, ]
  BootContData <- ContData[ContSampler, ]
  BootBigData <- BigData[BigSampler, ]
  
  # Growth models with bootstrapped data sets
  BootCrGAM <- gam(Plant_Height15 ~ s(Plant_Height14), data = rbind(BootCRData,
                                                BootBigData))
  BootControlGAM <- gam(Plant_Height15 ~ s(Plant_Height14), data = rbind(BootContData,
                                                  BootBigData))
  
  # Resample the posterior distributions for survival parameters
  BootCRSurvParams <- CrSurvChains[SurvSampler, ] %>% unlist()
  BootContSurvParams <- ContSurvChains[SurvSampler, ] %>% unlist()
  
  # pr(Repro)
  AllBootData <- rbind(BootCRData, BootContData, BootBigData)
  BootReproData <- filter(AllBootData, Survival != 'NA')
  
  BootReproGLM <- glm(Repro ~ Plant_Height15,
                      data = BootReproData,
                      family = binomial())
  
  # Seedling sizes
  BootSdls <- filter(AllBootData, Stage14 == 'SDL')
  BootSdlMean <- mean(BootSdls$Plant_Height14, na.rm = TRUE)
  BootSdlSD <- sd(BootSdls$Plant_Height14, na.rm = TRUE)
  
  
  # Create fecundity vectors and corner
  boot.f.params <- data.frame(prob.repro.int = as.numeric(coefficients(BootReproGLM)[1]),
                              prob.repro.slope = as.numeric(coefficients(BootReproGLM)[2]),
                              recruit.size.mean = BootSdlMean,
                              recruit.size.sd = BootSdlSD,
                              seed.int = as.numeric(coefficients(fecquasi)[1]),
                              seed.slope = coefficients(fecquasi)[2],
                              E2 = Establishment.prob2,
                              E3 = Establishment.prob3)
  
  SB2.Go <- SB.Go(Y, params = boot.f.params,
                  SB.Year = 2)
  SB3.Go <- SB.Go(Y, params = boot.f.params,
                  SB.Year = 3)
  
  SB2.Emerge <- SB.Emerge(Y, params = boot.f.params)
  SB3.Emerge <- rep(0, S)
  Discrete.Corner <- matrix(c(0, 0, 1, 0), 
                            nrow = 2, byrow = TRUE)
  
  # P matrix with bootstrapped data. Build G, correct for eviction,
  # then multiply by survival
  
  Gmat_Cont_2 <- h * outer(Y, Y, 
                           FUN = grow.prob,
                           model = BootControlGAM)
  Gmat_CR_2 <- h * outer(Y, Y,
                         FUN = grow.prob,
                         model = BootCrGAM)
  
  Gmat_Cont_2 <- Gmat_Cont_2/matrix(as.vector(apply(Gmat_Cont_2,
                                                    2,
                                                    sum)),
                                    nrow = S,
                                    ncol = S,
                                    byrow = TRUE)
  
  Gmat_CR_2 <- Gmat_CR_2/matrix(as.vector(apply(Gmat_CR_2,
                                                2,
                                                sum)),
                                nrow = S,
                                ncol = S,
                                byrow = TRUE)
  
  BootSurvCR <- boot_predict_surv(Y, BootCRSurvParams)
  BootSurvCont <- boot_predict_surv(Y, BootContSurvParams)
  
  BootPMatCR <- BootPMatCont <- matrix(0, S, S)
  
  for(x in seq_len(S)) {
    BootPMatCR[ ,x] <- Gmat_CR_2[ ,x] * BootSurvCR[x]
    BootPMatCont[ ,x] <- Gmat_Cont_2[ ,x] * BootSurvCont[x]
  }
  
  
  Kstart <- rbind(Discrete.Corner, cbind(SB3.Emerge, SB2.Emerge))
  K_Cont <- cbind(Kstart, rbind(SB3.Go, SB2.Go, Pmat_Cont))
  K_CR <- cbind(Kstart, rbind(SB3.Go, SB2.Go, Pmat_CR))
  
  # Calculate lambdas --------------------------------
  eigen_Cont <- eigen(K_Cont)
  OutputValues$Lambda_Cont[i] <- max(Re(eigen_Cont$values))
  
  eigen_CR <- eigen(K_CR)
  OutputValues$Lambda_CR[i] <- max(Re(eigen_CR$values))
  
  OutputValues$RepSlope[i] <- boot.f.params$prob.repro.slope
  OutputValues$RepInt[i] <- boot.f.params$prob.repro.int
  OutputValues$RecMean[i] <- BootSdlMean
  OutputValues$RecSD[i] <- BootSdlSD

  if(i %% 100 == 0) {
    message(i/10, '% of data crunched\n')
  }
}

OutputData <- as.data.frame(OutputValues)

write.csv(OutputData, 'BootStrap_Output_Euonymus.csv', row.names = FALSE)

library(gmailr)
sig <- source("gmailr_signature.R")

draft <- mime() %>%
  to(c("levisc8@gmail.com")) %>%
  from("samlevin.rstudio@gmail.com") %>%
  subject("Bootstrapping complete") %>%
  text_body(paste("Script complete, check RS Server for results",
                  sig$value, sep = '\n\n')) %>%
  send_message()


