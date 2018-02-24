# Ligustrum IPM w/ bootstrapping code. This is a pared down version
# of the full script and does not include the models created during
# the model selection process, only the final ones!

rm(list = ls()) 
graphics.off()

library(dplyr)
library(magrittr)
library(stringr)
library(mgcv)

# setwd('C:/Users/sl13sise/Dropbox/invaders demography/Ligustrum/IPM/LigObt_RS_Server') # local
setwd("~/LigObt_RS_Server/") # server

source('IPM_Functions_Ligustrum.R') 

# load the data and do some basic restructuring
RAs <- read.csv("RA4R.csv")
# str(RAs)
AllPlants <- read.csv("Ligobt4R-8.8.15.csv",
                      stringsAsFactors = FALSE) %>%
  filter(Clone == 0 | is.na(Clone)) %>%
  filter(Trt != 'Herb')

AllPlants$Plant <- as.numeric(AllPlants$Plant)
AllPlants$growth <- AllPlants$Ht15 - AllPlants$Ht14
AllPlants <- AllPlants %>% arrange(Quadrat,Plant)

# Growth-------------------------------------------------------
# Linear models. 
AllBig <- filter(AllPlants, Trt == 'All')

growth.lms <- AllPlants %>% 
  filter(Trt != "All") %>% 
  group_by(Trt) %>%
  do(LM = lm(Ht15 ~ Ht14, data = rbind(., AllBig)))

# Survival across treatments----------------------------------------------------------------------
# Logistic regression first with all plants that are classed by treatment.
AllControl <- filter(AllPlants, Trt == 'Control')
AllCR <- filter(AllPlants, Trt == 'Comp')

ContLinGlm <- glm(Alive2015 ~ Ht14,
                  data = rbind(AllControl, AllBig),
                  family = binomial())
CRLinGlm <- glm(Alive2015 ~ Ht14,
                data = rbind(AllCR, AllBig),
                family = binomial())
ContQuadGlm <- glm(Alive2015 ~ Ht14 + I(Ht14^2),
                   data = rbind(AllControl, AllBig),
                   family = binomial())
CRQuadGlm <- glm(Alive2015 ~ Ht14 + I(Ht14^2),
                 data = rbind(AllCR, AllBig),
                 family = binomial())

# Fecundity-------------------------------------------------------------------------------

# first, we will estimate fecundity from our sample of 10 RAs
fec <- glm(Fruit ~ Ht,
           data = RAs,
           family = quasipoisson())
summary(fec)

#Pr(Reproductive) - Logistic regression ----------------------------------
AllPlants$Repro <- ifelse(AllPlants$Stg15 == "RA", 1, 0)

Regression.Data <- filter(AllPlants, Alive2015 != "NA")

Repro.Glm <- glm(Repro ~ Ht14,
                 data = Regression.Data,
                 family = binomial())


# Recruit size distribution----------------------------------------------------------------------
sdls <- filter(AllPlants, Stg14 == "SDL")
Sdl.mean <- mean(sdls$Ht14, na.rm = TRUE)
Sdl.SD <- sd(sdls$Ht14, na.rm = TRUE)

# Establishment Pr()
germ.prob <- 0.5067 

# Using this as baseline estimate, but will substitute
# from 0.01 - 1 in bootstrapping loop to estimate sensitivity
est.prob <- 0.15 

# Create fecundity parameters data frame-------------------------------
f.params <- data.frame(prob.repro.int = as.numeric(coefficients(Repro.Glm)[1]),
                       prob.repro.slope = as.numeric(coefficients(Repro.Glm)[2]),
                       recruit.size.mean = Sdl.mean,
                       recruit.size.sd = Sdl.SD,
                       germ = germ.prob,
                       est.prob = est.prob) 


# the size range must extend beyond the limits of the data
min.size <- min(AllPlants$Ht14, na.rm = TRUE) * .6
max.size <- max(AllPlants$Ht14, na.rm = TRUE) * 1.2
S <- 500 # Number of cells in matrix  

# matrix variables 
b <-  min.size + c(0:S)*(max.size - min.size)/S  # boundary points of mesh cells
Y <- 0.5* (b[1:S] + b[2:(S+1)])  # mid points of mesh cells 
h <- Y[2]-Y[1]  # cell widths

# make a matrix of transitions for each growth distribution-------------
Gmat_comp <- h * (outer(Y, Y,
                        FUN = grow.prob,
                        model = growth.lms$LM[[1]]))

Gmat_cont <- h * (outer(Y, Y, 
                        FUN = grow.prob,
                        model = growth.lms$LM[[2]]))

Gmat_cont <- Gmat_cont/matrix(as.vector(apply(Gmat_cont,
                                              2,
                                              sum)),
                              nrow = S,
                              ncol = S,
                              byrow = TRUE)

Gmat_comp <- Gmat_comp/matrix(as.vector(apply(Gmat_comp,
                                              2,
                                              sum)),
                              nrow = S,
                              ncol = S,
                              byrow = TRUE)

s.cont_Lin <- predict.surv(Y, ContLinGlm)
s.cont_Quad <- predict.surv(Y, ContQuadGlm)
s.comp_Lin <- predict.surv(Y, CRLinGlm)
s.comp_Quad <- predict.surv(Y, CRQuadGlm)

P_Cont_Lin <- P_Cont_Quad <- P_Comp_Lin <- P_Comp_Quad <- matrix(0, S, S)

# make a P matrix (growth and survival)
for(i in seq_len(S)){
  P_Cont_Lin[ ,i] <- Gmat_cont[ ,i] * s.cont_Lin[i] 
  P_Cont_Quad[ ,i] <- Gmat_cont[ ,i] * s.cont_Quad[i]
  P_Comp_Lin[ ,i] <- Gmat_comp[ ,i] * s.comp_Lin[i]
  P_Comp_Quad[ ,i] <- Gmat_comp[ ,i] * s.comp_Quad[i]
}
# Build F Matrix-----------------------------------

FRow <- c(0, SB.Go(Y, params =  f.params, 
                   repro.model = Repro.Glm,
                   fec.model = fec))
FCol <- h * SB.Emerge(Y, f.params)

# build full K kernel -----------------------------------
K_cont_Lin <- rbind(FRow, cbind(FCol, P_Cont_Lin))
K_cont_Quad <- rbind(FRow, cbind(FCol, P_Cont_Quad))
K_comp_Lin <- rbind(FRow, cbind(FCol, P_Comp_Lin))
K_comp_Quad <- rbind(FRow, cbind(FCol, P_Comp_Quad))

# Calculate lambdas--------------------------------
eigen_cont_Lin <- eigen(K_cont_Lin)
lambda_cont_Lin <- max(Re(eigen_cont_Lin$values))

eigen_cont_Quad <- eigen(K_cont_Quad)
lambda_cont_Quad <- max(Re(eigen_cont_Quad$values))

eigen_comp_Lin <- eigen(K_comp_Lin)
lambda_comp_Lin <- max(Re(eigen_comp_Lin$values))

eigen_comp_Quad <- eigen(K_comp_Quad)
lambda_comp_Quad <- max(Re(eigen_comp_Quad$values))

# boot strapping and sensitivity to establishment probability-------------------------
# Additionally, at for one of these iterations, we need to store
# the boot strapped values for parameters so we can generate the
# vital rate figure panels
AllLambdas<-list(est.prob = seq(.01, 1, .01),
                 lambda_comp_Lin = rep(0,100),
                 lambda_comp_Lin_up = rep(0, 100),
                 lambda_comp_Lin_lo = rep(0, 100),
                 lambda_cont_Lin = rep(0,100),
                 lambda_cont_Lin_up = rep(0, 100),
                 lambda_cont_Lin_lo = rep(0, 100),
                 lambda_comp_Quad = rep(0,100),
                 lambda_comp_Quad_up = rep(0, 100),
                 lambda_comp_Quad_lo = rep(0, 100),
                 lambda_cont_Quad = rep(0,100),
                 lambda_cont_Quad_up = rep(0, 100),
                 lambda_cont_Quad_lo = rep(0, 100))

OutputValues <- list(GrowSlope_CR = c(coef(growth.lms$LM[[1]])[2],
                                      rep(NA, 1000)),
                     GrowInt_CR = c(coef(growth.lms$LM[[1]])[1],
                                    rep(NA, 1000)),
                     GrowSlope_Cont = c(coef(growth.lms$LM[[2]])[2],
                                        rep(NA, 1000)),
                     GrowInt_Cont = c(coef(growth.lms$LM[[2]])[1],
                                      rep(NA, 1000)),
                     RepSlope = c(coef(Repro.Glm)[2],
                                  rep(NA, 1000)),
                     RepInt = c(coef(Repro.Glm)[1],
                                rep(NA, 1000)),
                     RecMean = c(Sdl.mean,
                                 rep(NA, 1000)),
                     RecSD = c(Sdl.SD,
                               rep(NA, 1000)),
                     SurvInt_Cont_Lin = c(coef(ContLinGlm)[1],
                                          rep(NA, 1000)),
                     SurvSlope_Cont_Lin = c(coef(ContLinGlm)[2],
                                            rep(NA, 1000)),
                     SurvInt_Cont_Quad = c(coef(ContQuadGlm)[1],
                                           rep(NA, 1000)),
                     SurvSlope_Cont_Quad = c(coef(ContQuadGlm)[2],
                                             rep(NA, 1000)),
                     SurvSlope2_Cont_Quad = c(coef(ContQuadGlm)[3],
                                              rep(NA, 1000)),
                     SurvInt_CR_Lin = c(coef(CRLinGlm)[1],
                                        rep(NA, 1000)),
                     SurvSlope_CR_Lin = c(coef(CRLinGlm)[2],
                                          rep(NA, 1000)),
                     SurvInt_CR_Quad = c(coef(CRQuadGlm)[1],
                                         rep(NA, 1000)),
                     SurvSlope_CR_Quad = c(coef(CRQuadGlm)[2],
                                           rep(NA, 1000)),
                     SurvSelope2_CR_Quad = c(coef(CRQuadGlm)[3],
                                             rep(NA, 1000)),
                     Obs_Boot = c('Observed', rep('Boot', 1000)))

AllCR <- filter(AllPlants, Trt == 'Comp')
AllCont <- filter(AllPlants, Trt == 'Control')
AllBig <- filter(AllPlants, Trt == 'All')

for(i in unique(AllLambdas$est.prob)){
  it <- i * 100
  boot_lambdas_cont_Lin <- rep(NA, 1000)
  boot_lambdas_comp_Lin <- rep(NA, 1000)
  boot_lambdas_cont_Quad <- rep(NA, 1000)
  boot_lambdas_comp_Quad <- rep(NA, 1000)
  
  # substitute in new establishment probability
  f.params$est.prob <- i
  
  # Reconstruct the emergence column. This is the only parameter
  # affected by establishment probability, so first we test out
  # sensitivity to that. Then, bootstrap the data set with the
  # new establishment probability 
  FCol <- h * SB.Emerge(Y, f.params)
  
  # # build full K kernel -----------------------------------
  K_cont_Lin <- rbind(FRow, cbind(FCol, P_Cont_Lin))
  K_cont_Quad <- rbind(FRow, cbind(FCol, P_Cont_Quad))
  K_comp_Lin <- rbind(FRow, cbind(FCol, P_Comp_Lin))
  K_comp_Quad <- rbind(FRow, cbind(FCol, P_Comp_Quad))
  
  eigen_cont_Lin <- eigen(K_cont_Lin)
  lambda_cont_Lin <- max(Re(eigen_cont_Lin$values))
  AllLambdas$lambda_cont_Lin[it] <- lambda_cont_Lin
  
  eigen_comp_Lin <- eigen(K_comp_Lin)
  lambda_comp_Lin <- max(Re(eigen_comp_Lin$values))
  AllLambdas$lambda_comp_Lin[it] <- lambda_comp_Lin
  
  eigen_cont_Quad <- eigen(K_cont_Quad)
  lambda_cont_Quad <- max(Re(eigen_cont_Quad$values))
  AllLambdas$lambda_cont_Quad[it] <- lambda_cont_Quad
  
  eigen_comp_Quad <- eigen(K_comp_Quad)
  lambda_comp_Quad <- max(Re(eigen_comp_Quad$values))
  AllLambdas$lambda_comp_Quad[it] <- lambda_comp_Quad
  
  for(j in 1:1000){
    
    # Resample raw data sets, re-fit survival, growth and pr(Repro values)
    x1 <- sample(1:dim(AllControl)[1], dim(AllControl)[1], replace = TRUE)
    x2 <- sample(1:dim(AllCR)[1], dim(AllCR)[1], replace = TRUE)
    x3 <- sample(1:dim(AllBig)[1], dim(AllBig)[1], replace = TRUE)
    
    BootCont <- AllCont[x1, ]
    BootCR <- AllCR[x2, ]
    BootBig <- AllBig[x3, ]
    
    # Growth
    BootCRGrow <- lm(Ht15 ~ Ht14, data = rbind(BootCR, BootBig))
    BootContGrow <- lm(Ht15 ~ Ht14, data = rbind(BootCont, BootBig))
    
    # Survival
    BootCrSurvLin <- glm(Alive2015 ~ Ht14, 
                         data = rbind(BootCR, BootBig),
                         family = binomial())
    BootContSurvLin <- glm(Alive2015 ~ Ht14, 
                           data = rbind(BootCont, BootBig),
                           family = binomial())
    BootCrSurvQuad <- glm(Alive2015 ~ Ht14 + I(Ht14^2),
                          data = rbind(BootCR, BootBig),
                          family = binomial())
    BootContSurvQuad <- glm(Alive2015 ~ Ht14 + I(Ht14^2),
                            data = rbind(BootCont, BootBig),
                            family = binomial())
    
    # pr(Repro) and seedling size distributions
    AllData <- rbind(BootCR, BootCont, BootBig)
    
    BootReproData <- filter(AllData, Alive2015 != 'NA')
    BootReproGLM <- glm(Repro ~ Ht14, 
                        data = BootReproData,
                        family = binomial())
    
    BootSdls <- filter(AllData, Stg14 == 'SDL')
    BootSdlMean <- mean(BootSdls$Ht14, na.rm = TRUE)
    BootSdlSD <- sd(BootSdls$Ht14, na.rm = TRUE)
    
    # Combine all the bootstrapped fecundity parameters, but use
    # the same establishment probability as the outer loop to generate
    # a confidence interval for each level of it.
    BootFecParams <- data.frame(prob.repro.slope = coef(BootReproGLM)[2],
                                prob.repro.int = coef(BootReproGLM)[1],
                                recruit.size.mean = BootSdlMean,
                                recruit.size.sd = BootSdlSD,
                                germ = germ.prob,
                                est.prob = i)
    
    
    boot_Gmat_comp <- h * (outer(Y, Y,
                                 FUN = grow.prob,
                                 model = BootCRGrow))
    
    boot_Gmat_cont <- h * (outer(Y, Y, 
                                 FUN = grow.prob,
                                 model = BootContGrow))
    
    boot_Gmat_cont <- boot_Gmat_cont/matrix(as.vector(apply(boot_Gmat_cont,
                                                            2,
                                                            sum)),
                                            nrow = S,
                                            ncol = S,
                                            byrow = TRUE)
    
    boot_Gmat_comp <- boot_Gmat_comp/matrix(as.vector(apply(boot_Gmat_comp,
                                                            2,
                                                            sum)),
                                            nrow = S,
                                            ncol = S,
                                            byrow = TRUE)
    
    boot_s.cont_Lin <- predict.surv(Y, BootContSurvLin)
    boot_s.cont_Quad <- predict.surv(Y, BootContSurvQuad)
    boot_s.comp_Lin <- predict.surv(Y, BootCrSurvLin)
    boot_s.comp_Quad <- predict.surv(Y, BootCrSurvQuad)
    
    boot_P_Cont_Lin <- boot_P_Cont_Quad <- boot_P_Comp_Lin <- boot_P_Comp_Quad <- matrix(0, S, S)
    
    # make a P matrix (growth and survival)
    for(x in seq_len(S)){
      boot_P_Cont_Lin[ ,x] <- boot_Gmat_cont[ ,x] * boot_s.cont_Lin[x] 
      boot_P_Cont_Quad[ ,x] <- boot_Gmat_cont[ ,x] * boot_s.cont_Quad[x]
      boot_P_Comp_Lin[ ,x] <- boot_Gmat_comp[ ,x] * boot_s.comp_Lin[x]
      boot_P_Comp_Quad[ ,x] <- boot_Gmat_comp[ ,x] * boot_s.comp_Quad[x]
    }
    # Build fecundity vectors
    boot_FRow <- c(0, SB.Go(Y, BootFecParams,
                            repro.model = BootReproGLM,
                            fec.model = fec))
    boot_FCol <- h * SB.Emerge(Y, BootFecParams)
    
    # Build full kernel
    boot_K_cont_Lin <- rbind(boot_FRow, cbind(boot_FCol, boot_P_Cont_Lin))
    boot_K_cont_Quad <- rbind(boot_FRow, cbind(boot_FCol, boot_P_Cont_Quad))
    boot_K_comp_Lin <- rbind(boot_FRow, cbind(boot_FCol, boot_P_Comp_Lin))
    boot_K_comp_Quad <- rbind(boot_FRow, cbind(boot_FCol, boot_P_Comp_Quad))
    
    # Calculate lambdas--------------------------------
    boot_eigen_cont_Lin <- eigen(boot_K_cont_Lin)
    boot_lambdas_cont_Lin[j] <- max(Re(boot_eigen_cont_Lin$values))
    
    boot_eigen_cont_Quad <- eigen(boot_K_cont_Quad)
    boot_lambdas_cont_Quad[j] <- max(Re(boot_eigen_cont_Quad$values))
    
    boot_eigen_comp_Lin <- eigen(boot_K_comp_Lin)
    boot_lambdas_comp_Lin[j] <- max(Re(boot_eigen_comp_Lin$values))
    
    boot_eigen_comp_Quad <- eigen(boot_K_comp_Quad)
    boot_lambdas_comp_Quad[j] <- max(Re(boot_eigen_comp_Quad$values))
    
    # Store boot strap values for vital rates that don't depend
    # on EstProb. We only do this once because doing so 100
    # is wildly unnecessary
    if(it == 1) {
      OutputValues$GrowSlope_CR[j + 1] <- coef(BootCRGrow)[2]
      OutputValues$GrowInt_CR[j + 1] <- coef(BootCRGrow)[1]
      OutputValues$GrowSlope_Cont[j + 1] <- coef(BootContGrow)[2]
      OutputValues$GrowInt_Cont[j + 1] <- coef(BootContGrow)[1]
      OutputValues$RepSlope[j + 1] <- coef(BootReproGLM)[2]
      OutputValues$RepInt[j + 1] <- coef(BootReproGLM)[1]
      OutputValues$RecMean[j + 1] <- BootSdlMean
      OutputValues$RecSD[j + 1] <- BootSdlSD
      OutputValues$SurvInt_Cont_Lin[j + 1] <- coef(BootContSurvLin)[1]
      OutputValues$SurvSlope_Cont_Lin[j + 1] <- coef(BootContSurvLin)[2]
      OutputValues$SurvInt_Cont_Quad[j + 1] <- coef(BootContSurvQuad)[1]
      OutputValues$SurvSlope_Cont_Quad[j + 1] <- coef(BootContSurvQuad)[2]
      OutputValues$SurvSlope2_Cont_Quad[j + 1] <- coef(BootContSurvQuad)[3]
      OutputValues$SurvInt_CR_Lin[j + 1] <- coef(BootCrSurvLin)[1]
      OutputValues$SurvSlope_CR_Lin[j + 1] <- coef(BootCrSurvLin)[2]
      OutputValues$SurvInt_CR_Quad[j + 1] <- coef(BootCrSurvQuad)[1]
      OutputValues$SurvSlope_CR_Quad[j + 1] <- coef(BootCrSurvQuad)[2]
      OutputValues$SurvSelope2_CR_Quad[j + 1] <- coef(BootCrSurvQuad)[3]
      OutputValues$Obs_Boot[j + 1] <- 'Boot'
    }
    
  } # inner bootstrapping loop
  
  # Store confidence intervals from inner loop
  AllLambdas$lambda_comp_Lin_up[it] <- boot_lambdas_comp_Lin %>%
    sort %>%
    .[975]
  AllLambdas$lambda_comp_Lin_lo[it] <- boot_lambdas_comp_Lin %>%
    sort %>%
    .[25]
  AllLambdas$lambda_comp_Quad_up[it] <- boot_lambdas_comp_Quad %>%
    sort %>%
    .[975]
  AllLambdas$lambda_comp_Quad_lo[it] <- boot_lambdas_comp_Quad %>%
    sort %>%
    .[25]
  
  AllLambdas$lambda_cont_Lin_up[it] <- boot_lambdas_cont_Lin %>%
    sort %>%
    .[975]
  AllLambdas$lambda_cont_Lin_lo[it] <- boot_lambdas_cont_Lin %>%
    sort %>%
    .[25]
  AllLambdas$lambda_cont_Quad_up[it] <- boot_lambdas_cont_Quad %>%
    sort %>%
    .[975]
  AllLambdas$lambda_cont_Quad_lo[it] <- boot_lambdas_cont_Quad %>%
    sort %>%
    .[25]
  
  if(it %% 5 == 0){
    base::message(paste0(it, "% of data processed"))
  }
  
  
  
} # Estprob perturbation loop

OutputValues %>%
  as_tibble() %>%
  write.csv(., file = 'Ligustrum_Bootstrapped_Vital_Rates.csv',
            row.names = FALSE)

write.csv(AllLambdas, 'Bootstrapping_Output_Ligustrum.csv',
          row.names = FALSE)

library(gmailr)
sig <- source("gmailr_signature.R")


draft <- mime() %>%
  to(c("levisc8@gmail.com")) %>%
  from("samlevin.rstudio@gmail.com") %>%
  subject("Ligustrum Bootstrapping complete") %>%
  text_body(paste("Script complete, results on RStudio Server",
                  sig$value, sep = '\n\n')) %>%
  send_message()
