# I presume that those attempting to reproduce this have downloaded the archive
# and created an Rstudio project with it. Thus, file paths are relative to the 
# top level folder

# Euonymus IPM with bootstrapping. This script is intended to run
# on iDiv's RStudio server and/or the EVE cluster at iDiv/UFZ. It will
# likely not run as is on a local machine. I've also eliminated a lot
# the models that were fitted for comparison as those can be found
# in R/Euonymus_Exporation_and_IPM.R
# 
# rm(list = ls()) 
# graphics.off()

library(dplyr)
library(magrittr)
library(stringr)
library(mgcv)
library(brms)
library(tidyr)
library(ggplot2)
# setwd('~/EuoAla_RS_Folder')

source('Euonymus_IPM/IPM_Functions_Euonymus.R')  


# load the data
RAs<-read.csv("Euonymus_IPM/Euonymus_RA_Clean.csv")
AllPlants1 <- read.csv("Euonymus_IPM/Euonymus_Clean.csv") 

AllPlants1$incr<-AllPlants1$HeightNext-AllPlants1$Height

# Generate data sets for CR and control treatments

AllSmall <- filter(AllPlants1, Treatment != 'All')

ControlGrowData <- filter(AllPlants1, 
                          Treatment == 'Control' | 
                            Treatment == 'All')
CrGrowData <- filter(AllPlants1, 
                     Treatment == 'Comp' |
                       Treatment == 'All')

# Growth------------------------------------------------------------------
# Prior investigation indicated that a GAM provided a better description 
# of the growth dynamics

ControlGam <- gam(HeightNext ~ s(Height), data = ControlGrowData)
CRGam <- gam(HeightNext ~ s(Height), data = CrGrowData)

# Survival across treatments----------------------------------------------
# As with above, model selection procedures have been omitted in favor of brevity

BrmSurvQuad_Cont <- brm(Survival ~ Height + I(Height^2), 
                        data = ControlGrowData,
                        family = 'bernoulli',
                        iter = 6000,
                        warmup = 2000,
                        cores = getOption('mc.cores', 4L),
                        control = list(adapt_delta = 0.97,
                                       max_treedepth = 15),
                        save_model = 'Euonymus_IPM/Cont_Surv_BRM_Quad.stan') 

BrmSurvQuad_CR <- brm(Survival ~ Height + I(Height^2), 
                      data = CrGrowData,
                      family = 'bernoulli',
                      iter = 6000,
                      warmup = 2000,
                      cores = getOption('mc.cores', 4L),
                      control = list(adapt_delta = 0.97,
                                     max_treedepth = 15),
                      save_model = 'Euonymus_IPM/CR_Surv_BRM_Quad.stan') 

# Store paremeters and chains for later     
SurvParamsCont <- fixef(BrmSurvQuad_Cont)[ ,1] %>%
  setNames(c('Int', 'Lin', 'Quad'))

SurvParamsCR <- fixef(BrmSurvQuad_CR)[ ,1] %>%
  setNames(c('Int', 'Lin', 'Quad'))

thin <- seq(2001, 6000, length.out = 2000) %>% round()

CrSurvChains <- BrmSurvQuad_CR$fit@sim$samples[[1]][1:4] %>%
  data.frame() %>%
  setNames(c('Intercept', 'LinearTerm','QuadraticTerm', 'Likelihood')) %>%
  mutate(Treatment = 'CR') %>%
  .[thin, ]

ContSurvChains <- BrmSurvQuad_Cont$fit@sim$samples[[1]][1:4] %>%
  data.frame() %>%
  setNames(c('Intercept', 'LinearTerm','QuadraticTerm', 'Likelihood')) %>%
  mutate(Treatment = 'Control') %>%
  .[thin, ]

survChains <- rbind(CrSurvChains, ContSurvChains) %>% 
  select(Treatment, Intercept, LinearTerm, QuadraticTerm, Likelihood) %T>%
  write.csv(file = 'Euonymus_IPM/Euonymus_Survival_Model_Chains.csv', row.names = FALSE) 

save(survChains, file = 'Euonymus_IPM/Euonymus_Survival_Model_Chains.rda')

# Fecundity-------------------------------------------------------------------------------
# first, we will estimate fecundity from our sample of 18 RAs
#subsetting out two outliers with small Hts but lots of fruits

RAs2 <- filter(RAs, Fruit < 2500 & Height < max(RAs$Height))

# quasipoisson
fecquasi <- glm(Fruit ~ Height,
                data = RAs2,
                family = quasipoisson())
summary(fecquasi)

#Pr(Reproductive) - Logistic regression ----------------------------------
AllPlants1$Repro <- ifelse(AllPlants1$StageNext == "RA", 1, 0)

Regression.Data <- filter(AllPlants1, Survival != "NA")

Repro.Glm <- glm(Repro ~ HeightNext,
                 data = Regression.Data,
                 family = binomial())

AllPlants1$Fruits15 <- exp(coefficients(fecquasi)[1] + 
                             coefficients(fecquasi)[2] *AllPlants1$HeightNext)

AllPlants1$Fruits15[AllPlants1$StageNext != "RA"] <- NA
summary(AllPlants1$Fruits15)

# Recruit size distribution----------------------------------------------------------------------
sdls <- filter(AllPlants1, Stage == "SDL")
Sdl.mean <- mean(sdls$Height, na.rm = TRUE)
Sdl.SD <- sd(sdls$Height, na.rm = TRUE)

# Discrete parameters------------------------------------

# Establishment probability from Brand et al 2012. I'm taking the average 
# of all cultivars from both  time replications. seed bank germination is the 
# difference between cumulative germination in years 2 and 3
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

min.size <- min(AllPlants1$Height, na.rm = TRUE) * .8
max.size <- max(AllPlants1$Height, na.rm = TRUE) * 1.2
S <- 500 # Number of cells in matrix  

# matrix variables 
b <-  min.size + c(0:S)*(max.size - min.size)/S  # boundary points of mesh cells
Y <- 0.5* (b[1:S] + b[2:(S+1)])  # mid points of mesh cells 
h <- Y[2]-Y[1]  # cell widths

# NOTES SPECIFIC TO EUONYMUS
# The fecundity matrix for this species is actually more like 2 vectors
# that we'll wrap around the survival + growth matrix
# adults at time T don't produce seedlings at T+1 because of complex dormancy, 
# so we have a seedbank vector
# for each year (2 and 3) to account for seeds that emerge from each of those.
# I'll be making those vectors right now
# and then binding them onto the P matrix later.

SB2.Go <- SB.Go(Y, params = f.params,
                SB.Year = 2)
SB3.Go <- SB.Go(Y, params = f.params,
                SB.Year = 3)

SB2.Emerge <- SB.Emerge(Y, params = f.params)
SB3.Emerge <- rep(0, S)
Discrete.Corner <- matrix(c(0, 0, 1, 0), 
                          nrow = 2, byrow = TRUE)

# discretize the growth kernel
Gmat_Cont_2 <- h * outer(Y, Y, 
                         FUN = grow.prob,
                         model = ControlGam)
Gmat_CR_2 <- h * outer(Y, Y,
                       FUN = grow.prob,
                       model = CRGam)

# This re-scales all of the columns based on the column sums to correct for 
# evicted indidivduals (Williams et al 2012)

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
# Each entry in the survival vector must be multiplied by the same column in the 
# the growth matrix

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

# Initialize data storage areas and underlying data sets.

CRData <- filter(AllPlants1, Treatment == 'Comp')
ContData <- filter(AllPlants1, Treatment == 'Control')
BigData <- filter(AllPlants1, Treatment == 'All')

# Compute number of individuals in each data set so that we can resample them 
# with replacement in the for loop
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
  
  # Set up vectors and resample data sets and survival parameter posteriors
  # This creates an integer vector by sampling with replacement. When used to 
  # index the raw data sets, it acts as a bootstrapper
  
  CrSampler <- sample(1:nCR, nCR, replace = TRUE)
  ContSampler <- sample(1:nCont, nCont, replace = TRUE)
  BigSampler <- sample(1:nBig, nBig, replace = TRUE)
  
  # This is just 1 integer because we need to sample the posterior distribution
  # by row. All parameters in a given draw are conditional on each other, so it
  # it is incorrect to mix and match values from different rows!
  SurvSampler <- sample(1:1000, 1)  
  
  # The actual bootstrapping
  BootCRData <- CRData[CrSampler, ]
  BootContData <- ContData[ContSampler, ]
  BootBigData <- BigData[BigSampler, ]
  
  # Growth models with bootstrapped data sets
  BootCrGAM <- gam(HeightNext ~ s(Height), data = rbind(BootCRData,
                                                        BootBigData))
  BootControlGAM <- gam(HeightNext ~ s(Height), data = rbind(BootContData,
                                                             BootBigData))
  
  # Resample the posterior distributions for survival parameters
  BootCRSurvParams <- CrSurvChains[SurvSampler, ] %>% unlist()
  BootContSurvParams <- ContSurvChains[SurvSampler, ] %>% unlist()
  
  # pr(Repro)
  AllBootData <- rbind(BootCRData, BootContData, BootBigData)
  BootReproData <- filter(AllBootData, Survival != 'NA')
  
  BootReproGLM <- glm(Repro ~ HeightNext,
                      data = BootReproData,
                      family = binomial())
  
  # Seedling size distribution
  BootSdls <- filter(AllBootData, Stage == 'SDL')
  BootSdlMean <- mean(BootSdls$Height, na.rm = TRUE)
  BootSdlSD <- sd(BootSdls$Height, na.rm = TRUE)
  
  
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
  
  # cheap version of a progress bar 
  if(i %% 100 == 0) {
    message(i/10, '% of data crunched\n')
  }
}

OutputData <- as.data.frame(OutputValues)

write.csv(OutputData, 'BootStrap_Output_Euonymus.csv', row.names = FALSE)

# Build figures. This reshaping is ugly, but leads to prettier figures

EuoAlaData <- OutputData %>%
  gather(key = 'Variable', value = 'Value') %>%
  mutate(Trt = vapply(.$Variable,
                      FUN = function(x) str_split(x, '_')[[1]][2],
                      FUN.VALUE = '')) %>%
  group_by(Variable, Trt) %>%
  arrange(desc(Value)) %>%
  summarise(obs = NA,
            UpCI = Value[25],
            LoCI = Value[975])

EuoAlaData$obs[1] <- lambda_Cont_obs
EuoAlaData$obs[2] <- lambda_CR_obs
EuoAlaData$obs[3] <- Sdl.mean
EuoAlaData$obs[4] <- Sdl.SD
EuoAlaData$obs[5] <- f.params$prob.repro.int
EuoAlaData$obs[6] <- f.params$prob.repro.slope

SurvCIsCR <- fixef(BrmSurvQuad_CR)[ ,-2]
SurvCIsCont <- fixef(BrmSurvQuad_Cont)[ ,-2]

EuoAlaData$Variable <- c('Lambda', 'Lambda',
                         'Recruit Size Mean',
                         'Recruit Size SD',
                         'Repro Intercept',
                         'Repro Slope')


SurvSummary <- tibble(Variable = rep(c('Surv Intercept',
                                       'Surv Linear Term',
                                       'Surv Quadratic Term'), 2),
                      Trt = c(rep('Cont', 3),
                              rep('CR', 3)),
                      obs = c(SurvCIsCont[ ,1],
                              SurvCIsCR[ ,1]),
                      UpCI = c(SurvCIsCont[ ,3],
                               SurvCIsCR[ ,3]),
                      LoCI = c(SurvCIsCont[ ,2],
                               SurvCIsCR[ ,2]))


PlotData <- rbind.data.frame(SurvSummary,
                             EuoAlaData,
                             stringsAsFactors = FALSE)
write.csv(PlotData, 'Euonymus_Summarized_Output.csv', row.names = FALSE)

PlotData$Trt[is.na(PlotData$Trt)] <- 'Pooled'
PlotData$Trt[PlotData$Trt == 'Cont'] <- 'Control'
PlotData$Variable[PlotData$Variable == 'Recruit Size Mean'] <- "paste(mu, ' Recruit Size Mean')"
PlotData$Variable[PlotData$Variable == 'Recruit Size SD'] <- "paste(sigma, ' Recruit Size SD')"
PlotData$Variable[PlotData$Variable == 'Repro Intercept'] <- "paste(italic(f[p](x)), ' Intercept')"
PlotData$Variable[PlotData$Variable == 'Repro Slope'] <- "paste(italic(f[p](x)), ' Slope')"
PlotData$Variable[PlotData$Variable == 'Surv Intercept'] <- "paste(italic(s(x)),' Intercept')"
PlotData$Variable[PlotData$Variable == 'Surv Linear Term'] <- "paste(italic(s(x)),' Linear Term')"
PlotData$Variable[PlotData$Variable == 'Surv Quadratic Term'] <- "paste(italic(s(x)),' Quadratic Term')"
PlotData$Variable[PlotData$Variable == 'Lambda'] <- 'lambda'


VR_Plot <- ggplot(data = PlotData,
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
        legend.position = c(0.85, 0.1), # legend in bottom right corner
        legend.title = element_text(size = 22), # make legend text bigger
        legend.text = element_text(size = 18),
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
  scale_y_continuous('Observed Value + Confidence Intervals') + 
  scale_x_discrete('Treatment (if applicable)') +
  scale_color_manual('Treatment',
                     breaks = c('Control', 'CR', "Pooled"),
                     values = c('black', 'green', 'blue'))

VR_Plot    

ggsave(filename = 'Euonymus_Vital_Rate_Coefficients.png',
       path = 'Euonymus_IPM/',
       height = 8,
       width = 10,
       unit = 'in')

