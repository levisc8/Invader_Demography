#####################################################################################
### Euonymus Exploration and IPM ######################################
#####################################################################################
# notes from 1.12-13.17
# Assumptions checked against Merow et al 2014 Advancing population 
# ecology with integral projection models: a practical guide, Methods in 
# Ecol and Evo

# 1. Explore data behavior.
#   a. I think we should remove the two adults that grew and shrank ~200cm. I checked the raw
#      data and that is not an entry error, so I think it is a field recording error. I am not
#      sure how to correct this otherwise, and it's very weird 
# 2. Fit growth models
#   a. Spline fit, see Rees et al 2014
# 3. Fit survival models
#   a.Again, combine for survival. I think the weirdness with Bunker HR is an aberration
#     and not representative of the actual treatment effect
# 4. fit fecundity models
#   a. Combined for all sites
# 
# rm(list = ls()) 
# cat("\014")

#####################################################################################
### LIBRARIES ###
library(dplyr)
library(magrittr)
library(stringr)
library(mgcv)

source('Euonymus/IPM/R/IPM_Functions_Euonymus.R')

# load the data
RAs <- read.csv("Euonymus/IPM/Data/RA_Clean.csv")
AllPlants1 <- read.csv("Euonymus/IPM/Data/EA_Clean.csv") 

AllPlants1$incr<-AllPlants1$Plant_Height15-AllPlants1$Plant_Height14

AllSmall <- filter(AllPlants1, Treatment != 'All')
ControlGrowData <- filter(AllPlants1, Treatment == 'Control' | 
                            Treatment == 'All')
CrGrowData <- filter(AllPlants1, Treatment == 'Comp' |
                       Treatment == 'All')

# Growth---------------------------------------------------------------------------------------
# trying a linear fit. First, model with all data < 200 cm (because anything larger was not 
# actually treated in our plots) and use Treatment as an interaction 
# term to determine significance. If it is significant, then construct separate models
# for each treatment by combining large plants in with the small ones.

AllLM<- lm(Plant_Height15 ~ Plant_Height14 * Treatment, data = AllSmall)
summary(AllLM)

# Same as above, so separate and do LMs
ControlLM <- lm(Plant_Height15 ~ Plant_Height14, data = ControlGrowData)
summary(ControlLM)

CRLM <- lm(Plant_Height15 ~ Plant_Height14, data = CrGrowData)
summary(CRLM)


# Going to try a spline fit as well since these relationships don't
# really look linear.

xx <- seq(0, max(AllPlants1$Plant_Height14, na.rm = TRUE) + 50, 1)

ControlGam <- gam(Plant_Height15 ~ s(Plant_Height14), data = ControlGrowData)
CRGam <- gam(Plant_Height15 ~ s(Plant_Height14), data = CrGrowData)
summary(ControlGam)
summary(CRGam)

ControlGamPred <- predict(ControlGam, data.frame(Plant_Height14 = xx), type = 'response')
CrGamPred <- predict(CRGam, data.frame(Plant_Height14 = xx), type = 'response')


# plot results
plot(Plant_Height15 ~ Plant_Height14, data = AllPlants1, type = 'n')
points(CrGrowData$Plant_Height14[CrGrowData$Plant_Height14 < 100],
       CrGrowData$Plant_Height15[CrGrowData$Plant_Height14 < 100], 
       col = 'green',
       pch = 4)
points(ControlGrowData$Plant_Height14[ControlGrowData$Plant_Height14 < 100],
       ControlGrowData$Plant_Height15[ControlGrowData$Plant_Height14 < 100],
       col = 'blue',
       pch = 3)
points(ControlGrowData$Plant_Height14[ControlGrowData$Plant_Height14 > 100],
       ControlGrowData$Plant_Height15[ControlGrowData$Plant_Height14 > 100],
       col = 1,
       pch = 1)
lines(xx, ControlGamPred, col = 'blue')
lines(xx, CrGamPred, col = 'green')

# test model
AIC(ControlLM, ControlGam)
AIC(CRLM, CRGam)

# Rae: Did you test the fit for linear vs. B-spline vs. quadratic? I will give you some useful code in email.
# Justification will require that you test the fit as compared to other models. If it does not differ from the
# linear model, I would suggest using the simpler model. In this case, it does apear that a curve would be the 
# best fit.



# Survival across treatments----------------------------------------------------------------------
# Testing each treatment to see if treatments are different or not.


# View logistic regressions for each treatment to evaluate treatment differencesb 
AllSurvGLM <- glm(Survival ~ Plant_Height14 * Treatment, 
                  data = AllSmall, 
                  family = binomial())
AllSurvQuadGLM <- glm(Survival ~ Plant_Height14 * Treatment + I(Plant_Height14^2) * Treatment,
                      data = AllSmall,
                      family = binomial)

summary(AllSurvGLM)
summary(AllSurvQuadGLM)

# There is a marginally significant interaction between Treatment and size
# for survival with no polynomial term, but none when the polynomial
# is included. The AIC score for the non-polynomial model is lower,
# so i will use the results from that model and keep the treatments 
# separate
Sample.Sizes <- AllPlants1 %>% group_by(Site,Treatment) %>% summarise(n = n())
Sample.Sizes

# Create separate models for each treatment
ContSurv <- glm(Survival ~ Plant_Height14, 
                data = ControlGrowData,
                family = binomial())
CompSurv <- glm(Survival ~ Plant_Height14,
                data = CrGrowData,
                family = binomial())


xx <- seq(0, max(AllPlants1$Plant_Height14) + 20, 1)
plot(jitter(Survival, amount = 0.05) ~ Plant_Height14, 
     data = AllPlants1,
     xlab = 'Height 2014',
     ylab = 'Survival')
lines(xx, predict(ContSurv,
                  data.frame(Plant_Height14 = xx),
                  type = 'response'),
      col = 'black',
      lty = 2)
lines(xx, predict(CompSurv,
                  data.frame(Plant_Height14 = xx),
                  type = 'response'),
      col = 'green', lty = 2)


# These models both predict survival probabilities of 1! But we know plants die,
# so I'll try fitting quadratic models next
ContQuadSurv <- glm(Survival ~ Plant_Height14 + I(Plant_Height14^2),
                    data = ControlGrowData,
                    family = binomial())
CompQuadSurv <- glm(Survival ~ Plant_Height14 + I(Plant_Height14^2),
                    data = CrGrowData,
                    family = binomial())

lines(xx, predict(ContQuadSurv,
                  data.frame(Plant_Height14 = xx),
                  type = 'response'),
      col = 'black',
      lty = 1)
lines(xx, predict(CompQuadSurv,
                  data.frame(Plant_Height14 = xx),
                  type = 'response'),
      col = 'green',
      lty = 1)

legend('right',
       c("Control",
         "CR",
         "Control - Polynomial",
         "CR - Polynomial"),
       col = c("black", "green",
               "black", "green"),
       lty = c(2, 2, 1, 1))


# The control model's quadratic term wildly overcompensates for lack of mortality
# by murdering everything above 300 cm. This is bad. The CR model has the opposite problem -
# it overcompensates by setting the intercept too high, dipping early, then converging
# to 1 (yielding the same predictions as the non-quadratic model). Next, I'll try fitting
# linear and quadratic models with brms + Stan. These should weight the quadratic term
# a bit less, hopefully yielding useful predictions
# 3. Fit survival models---------------------------------------------------------------------

library(brms)

BrmSurv_Cont <- brm(Survival ~ Plant_Height14, 
                     data = ControlGrowData,
                     family = 'bernoulli',
                     iter = 4000,
                     warmup = 1000,
                     cores = getOption('mc.cores', 2L),
                     save_model = 'Euonymus/IPM/R/Cont_Linear_BRM.stan',
                     control = list(adapt_delta = 0.97,
                                    max_treedepth = 15)) %>%
  add_waic()

BrmSurvQuad_Cont <- brm(Survival ~ Plant_Height14 + I(Plant_Height14^2), 
                          data = ControlGrowData,
                          family = 'bernoulli',
                          iter = 4000,
                          warmup = 1000,
                          cores = getOption('mc.cores', 2L),
                          save_model = 'Euonymus/IPM/R/Cont_Quadratic_BRM.stan',
                          control = list(adapt_delta = 0.97,
                                         max_treedepth = 15)) %>%
  add_waic()

BrmSurv_CR <- brm(Survival ~ Plant_Height14, 
                     data = CrGrowData,
                     family = 'bernoulli',
                     iter = 4000,
                     warmup = 1000,
                     cores = getOption('mc.cores', 2L),
                     save_model = 'Euonymus/IPM/R/CR_Linear_BRM.stan',
                     control = list(adapt_delta = 0.97,
                                    max_treedepth = 15)) %>%
  add_waic()

BrmSurvQuad_CR <- brm(Survival ~ Plant_Height14 + I(Plant_Height14^2), 
                          data = CrGrowData,
                          family = 'bernoulli',
                          iter = 4000,
                          warmup = 1000,
                          cores = getOption('mc.cores', 2L),
                          save_model = 'Euonymus/IPM/R/CR_Quadratic_BRM.stan',
                          control = list(adapt_delta = 0.97,
                                         max_treedepth = 15)) %>%
  add_waic()

BrmSurv_Cont$waic$waic
BrmSurvQuad_Cont$waic$waic
BrmSurv_CR$waic$waic
BrmSurvQuad_CR$waic$waic

BrmSurv_Cont$waic$waic - BrmSurvQuad_Cont$waic$waic
BrmSurv_CR$waic$waic - BrmSurvQuad_CR$waic$waic


plot(BrmSurv_Cont)
plot(BrmSurvQuad_Cont)
plot(BrmSurv_CR)
plot(BrmSurvQuad_CR)

summary(BrmSurv_Cont)
summary(BrmSurvQuad_Cont)
summary(BrmSurv_CR)
summary(BrmSurvQuad_CR)


# plot the survival curve
par(mfrow = c(1,1))
xx<-seq(0,400, 1)
plot(Survival ~ Plant_Height14, data = AllPlants1)
lines(xx, predict(BrmSurv_Cont, data.frame(Plant_Height14 = xx))[,1],
      col = 'black')
lines(xx, predict(BrmSurvQuad_Cont, data.frame(Plant_Height14 = xx))[,1],
      col = 'black',
      lty = 2)
lines(xx, predict(BrmSurv_CR, data.frame(Plant_Height14 = xx))[,1],
      col = 'green')
lines(xx, predict(BrmSurvQuad_CR, data.frame(Plant_Height14 = xx))[,1],
      col = 'green',
      lty = 2)

legend('bottomright', c('Control', 'Control - Quadratic',
                        'CR', 'CR - Quadratic'),
       col = c('black', 'black',
               'green', 'green'), 
       lty = c(1, 1, 2, 2))
      
SurvParamsCont <- fixef(BrmSurvQuad_Cont)[ ,1] %>%
  setNames(c('Int', 'Lin', 'Quad'))

SurvParamsCR <- fixef(BrmSurvQuad_CR)[ ,1] %>%
  setNames(c('Int', 'Lin', 'Quad'))

# Fecundity-------------------------------------------------------------------------------
# exploratory for now, unsure of how we are going to do this.
# we are missing a recruit size distribution from year two. For now, I am using the initial
# seedling size distribution and seeing how that goes. We can figure out something else
# later. Establishment probability will come from Brand et al 2012 as the average of all 
# germination across multiple habitat types

# first, we will estimate fecundity from our sample of 18 RAs
fec <- glm(Fruit ~ Ht,
           data = RAs,
           family = quasipoisson())
summary(fec)
# normal poisson is very overdispersed, trying lm and splines next
feclm <- lm(Fruit ~ Ht, data = RAs)
summary(feclm)

fec.gam <- gam(Fruit ~ s(Ht), data = RAs)
summary(fec.gam)

#subsetting out two outliers with small Hts but lots of fruits
RAs2 <- filter(RAs, Fruit < 2500 & Ht < max(RAs$Ht))
# poisson, lm, gam
fec <- glm(Fruit ~ Ht,
           data = RAs2,
           family = poisson())
summary(fec)
# still overdispersed
fecquasi <- glm(Fruit ~ Ht,
                data = RAs2,
                family = quasipoisson())
summary(fecquasi)

# not bad, but I'll try the lm and gam spline next
feclm <- lm(Fruit ~ Ht, data = RAs2)
summary(feclm)

fec.gam <- gam(Fruit ~ s(log(Ht)),
               data = RAs2)
summary(fec.gam)

xx <- seq(min(RAs2$Ht), max(RAs2$Ht), 0.1)
plot(Fruit ~ Ht, data = RAs2)
lines(xx, predict(fecquasi,
                  data.frame(Ht = xx),
                  type = 'response'))
lines(xx, predict(fec.gam, data.frame(Ht = xx),
                  type = 'response'),
      col = 2)

anova(fec, fec.gam, test = 'Chi') # keep poisson
anova(fec, feclm, test = 'Chi')
#Pr(Reproductive) - Logistic regression ----------------------------------
AllPlants1$Repro <- ifelse(AllPlants1$Stage15 == "RA", 1, 0)

Regression.Data <- filter(AllPlants1, Survival != "NA")

Repro.Glm <- glm(Repro ~ Plant_Height15,
                 data = Regression.Data,
                 family = binomial())
summary(Repro.Glm)
plot(Repro ~ Plant_Height15,
     data = Regression.Data)
xx <- seq(0, max(AllPlants1$Plant_Height15, na.rm = TRUE), 0.1)
lines(xx, predict(Repro.Glm,
                  data.frame(Plant_Height15 = xx),
                  type = 'response'),
      lty = 2, col = 'red')
AllPlants1$Fruits15 <- exp(coefficients(fecquasi)[1] + 
                            coefficients(fecquasi)[2] *AllPlants1$Plant_Height15)
AllPlants1$Fruits15[AllPlants1$Stage15 != "RA"] <- NA
summary(AllPlants1$Fruits15)




# Recruit size distribution----------------------------------------------------------------------
# Not really sure how to do this with no data from 2015. Using 2014 size distribution for now,
# but there's really no way to separate out the treatments, and I think this is the only place
# we could potentially see treatment effect. However, we may be able to justify this by pointing out
# that there was no difference in growth or survival for small plants between the two censuses,
# so maybe we're ok
sdls <- filter(AllPlants1, Stage14 == "SDL")
hist(sdls$Plant_Height14)
Sdl.mean <- mean(sdls$Plant_Height14, na.rm = TRUE)
Sdl.SD <- sd(sdls$Plant_Height14, na.rm = TRUE)

sdls$Treatment <- droplevels(sdls$Treatment)
sdl.tab <- table(sdls$Treatment, sdls$Survival)
sdl.chisq <- chisq.test(sdl.tab, correct = FALSE)
sdl.chisq
sdl.tab
sdl.chisq$exp

# Create figure panel for each regression
par(mfrow = c(2, 2), 
    mar = c(5,6,4,2) + 0.2)
xx<-seq(0, max(AllPlants1$Plant_Height14, na.rm = TRUE) + 50, 1)
plot(Plant_Height15 ~ Plant_Height14,
     data = AllPlants1, 
     xlab = 'Size (t)', 
     ylab = 'Size (t+1)',
     cex.lab = 1.8)
size2MeanCont <- predict(ControlGam,
                         data.frame(Plant_Height14 = xx), 
                         type = 'response')
size2MeanComp <- predict(CRGam,
                         data.frame(Plant_Height14 = xx),
                         type = 'response')
lines(xx, size2MeanCont, col = 'black', lty = 2)
lines(xx, size2MeanComp, col = 'green', lty = 2)
legend('bottomright',
       c('Control', 'CR'),
       col = c('black', 'green'),
       lty = c(2, 2))

plot(Survival ~ Plant_Height14, data = AllPlants1,
     xlab = 'Size (t)',
     ylab = 'Survival (t+1)',
     cex.lab = 1.8)
lines(xx, predict(BrmSurvQuad_Cont, 
                  data.frame(Plant_Height14 = xx))[,1],
      col = 'black',
      lty = 2)
lines(xx, predict(BrmSurvQuad_CR, 
                  data.frame(Plant_Height14 = xx))[,1],
      col = 'green',
      lty = 2)
legend('bottomright',
       c('Control', 'CR'),
       col = c('black', 'green'),
       lty = c(2, 2))

xx<-seq(150, max(RAs2$Ht), 0.1)
plot(Fruit ~ Ht, 
     data = RAs2,
     xlab = 'Size (t+1)',
     ylab = 'Seed Production (t+1)',
     xlim = c(0, 380),
     cex.lab = 1.8)
lines(xx, predict(fecquasi,
                  data.frame(Ht = xx),
                  type = 'response'),
      col = 'red',
      lty = 2)

plot(Repro ~ Plant_Height15,data = Regression.Data,
     xlab = 'Size (t+1)',
     ylab = 'Pr( Reproductive) (t+1)',
     cex.lab = 1.8)
xx<-seq(0, max(AllPlants1$Plant_Height15, na.rm=TRUE), 0.1)
lines(xx, predict(Repro.Glm, 
                  data.frame(Plant_Height15 = xx),
                  type = 'response'),
      lty = 2, 
      col = 'red')

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

# Build F Matrix-----------------------------------
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



# make a matrix of transitions for each growth distribution-------------

# I'm correcting for eviction using a 
# method that Maria Paniw recommended. This re-scales all of the columns
# based on the column sums
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

par(mfrow = c(1,2))

plot((1-colSums(Gmat_Cont_2)), type='l',ylim=c(-1e-12, 1e-12))
abline(h = 0, col = 'red', lty = 2)

plot((1-colSums(Gmat_CR_2)), type='l',ylim=c(-1e-12, 1e-12))
abline(h = 0, col = 'red', lty = 2)

# Maria's method handles eviction a bit better,
# so I'm going to use those kernels for the 
# rest of the analysis

# 
# for(i in 1:(S/5)) Gmat[1,i]<-Gmat[1,i]+1-sum(Gmat[,i])
# for(i in (S*4/5):S) Gmat[S,i]<-Gmat[S,i]+1-sum(Gmat[,i])
#  
source('C:/Users/sl13sise/Dropbox/ATSC 2018 participant folder/23.1.18/Rees/R code and Data/MatrixImage.R')

par(mfrow=c(1,1))
matrix.image(Gmat_Cont_2, xlab = 'size (t)', 
             ylab = 'size (t + 1)', 
             main = "Growth Kernel, Control")
matrix.image(Gmat_CR_2, xlab = 'size (t)',
             ylab = 'size (t + 1)',
             main = 'Growth Kernel, CR')

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
par(mfrow = c(1,2))
plot(colSums(Pmat_Cont), type = 'l')
plot(colSums(Pmat_CR), type = 'l')


# inspect survival + growth kernel
par(mfrow=c(1,1))
matrix.image(Pmat_Cont, main = "Survival+Growth Kernel, Control")
matrix.image(Pmat_CR, main = 'Survival+Growth Kernel, CR')

# build full K kernel-----------------------------------
Kstart <- rbind(Discrete.Corner, cbind(SB3.Emerge, SB2.Emerge))
K_Cont <- cbind(Kstart, rbind(SB3.Go, SB2.Go, Pmat_Cont))
K_CR <- cbind(Kstart, rbind(SB3.Go, SB2.Go, Pmat_CR))


# inspect kernels. This won't look like much for this IPM though because transitions in seed bank stages
# have probability of 1, while most other cells have probabilities between .01 and 1e-20.
matrix.image(K_Cont ^ 0.05, main = "Full Kernel, Control")
matrix.image(K_CR ^ 0.05, main = "Full Kernel, CR")


# Calculate lambdas and elasticities--------------------------------
eigen_Cont <- eigen(K_Cont)
lambda_Cont_obs <- max(Re(eigen_Cont$values))

eigen_CR <- eigen(K_CR)
lambda_CR_obs <- max(Re(eigen_CR$values))

lambda_Cont_obs
lambda_CR_obs

# Now, we need to bootstrap everything (hooray!!!).
# I doubt there will be significant differences, but you never
# know. 

# I'm going to set up some data frames outside of the 
# loop for my computer's sake, but this will still probably take a
# a very long time to run. For survival parameters, rather than
# refitting the models, I'm going to resample the output from
# the stan models.



