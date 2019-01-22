# Ligustrum IPM w/ bootstrapping code. This is a pared down version
# of the full script and does not include the models created during
# the model selection process, only the final ones! 

# It is intended to run
# on iDiv's RStudio Server or on UFZ's EVE Cluster, as the bootstrapping takes
# quite a bit of time. Computing the transfer function of the 
# establishment probability is also quite time consuming. Thus, that portion
# is commented out so that you can review how the calculations were done, but the 
# code will use the stored outputs from those computations rather re-performing them
# 
# rm(list = ls()) 
# graphics.off()

library(dplyr)
library(magrittr)
library(stringr)
library(mgcv)
library(ggplot2)
library(tidyr)

source('Ligustrum_IPM/IPM_Functions_Ligustrum.R') 

# load the data and do some basic restructuring

RAs <- read.csv("Ligustrum_IPM/Ligustrum_RA_Clean.csv")
AllPlants <- read.csv("Ligustrum_IPM/Ligustrum_Clean.csv",
                      stringsAsFactors = FALSE)

AllPlants$Plant <- as.numeric(AllPlants$Plant)
AllPlants$growth <- AllPlants$HeightNext - AllPlants$Height
AllPlants <- AllPlants %>% arrange(Plot,Plant)

# Growth-------------------------------------------------------
# Linear models. 
AllBig <- filter(AllPlants, Treatment == 'All')

growth.lms <- AllPlants %>% 
  filter(Treatment != "All") %>% 
  group_by(Treatment) %>%
  do(LM = lm(HeightNext ~ Height, data = rbind(., AllBig)))

# Survival across treatments----------------------------------------------------------------------
# Logistic regression first with all plants that are classed by treatment.
AllControl <- filter(AllPlants, Treatment == 'Control')
AllCR <- filter(AllPlants, Treatment == 'Comp')

ContLinGlm <- glm(Survival ~ Height,
                  data = rbind(AllControl, AllBig),
                  family = binomial())
CRLinGlm <- glm(Survival ~ Height,
                data = rbind(AllCR, AllBig),
                family = binomial())
ContQuadGlm <- glm(Survival ~ Height + I(Height^2),
                   data = rbind(AllControl, AllBig),
                   family = binomial())
CRQuadGlm <- glm(Survival ~ Height + I(Height^2),
                 data = rbind(AllCR, AllBig),
                 family = binomial())

# Fecundity-------------------------------------------------------------------------------

# first, we will estimate fecundity from our sample of 10 RAs
fec <- glm(Seeds ~ Height,
           data = RAs,
           family = quasipoisson())
summary(fec)

#Pr(Reproductive) - Logistic regression ----------------------------------
AllPlants$Repro <- ifelse(AllPlants$Stage15 == "RA", 1, 0)

Regression.Data <- filter(AllPlants, Survival != "NA")

Repro.Glm <- glm(Repro ~ Height,
                 data = Regression.Data,
                 family = binomial())


# Recruit size distribution----------------------------------------------------------------------
sdls <- filter(AllPlants, Stage14 == "SDL")
Sdl.mean <- mean(sdls$Height, na.rm = TRUE)
Sdl.SD <- sd(sdls$Height, na.rm = TRUE)

# Establishment Pr.  See manuscript for details

germ.ligustrum <- read.csv('Germination/Clean_Germ.csv', 
                           stringsAsFactors = FALSE) %>%
  filter(Species == 'Ligustrum')

germ.prob <- mean(c(germ.ligustrum$Surface_Germ_Prop, 
                    germ.ligustrum$Buried_Germ_Prop))

# Using this as baseline estimate, but will substitute
# from 0.01 - 1 in bootstrapping loop to estimate lambda for all possible values
est.prob <- 0.15 

# Create fecundity parameters data frame-------------------------------
f.params <- data.frame(prob.repro.int = as.numeric(coefficients(Repro.Glm)[1]),
                       prob.repro.slope = as.numeric(coefficients(Repro.Glm)[2]),
                       recruit.size.mean = Sdl.mean,
                       recruit.size.sd = Sdl.SD,
                       germ = germ.prob,
                       est.prob = est.prob) 


# the size range must extend beyond the limits of the data
min.size <- min(AllPlants$Height, na.rm = TRUE) * .6
max.size <- max(AllPlants$Height, na.rm = TRUE) * 1.2
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

# This is commented out in lieu of using the outputs that were already computed
AllLambdas <- list(est.prob = seq(.01, 1, .01),
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

AllCR <- filter(AllPlants, Treatment == 'Comp')
AllCont <- filter(AllPlants, Treatment == 'Control')
AllBig <- filter(AllPlants, Treatment == 'All')

# If you do have the necessary computational resources, then uncomment the loop
# below to replicate the results. Otherwise, skip down to reading in the 
# computed values

# for(i in unique(AllLambdas$est.prob)){
#   it <- i * 100
#   boot_lambdas_cont_Lin <- rep(NA, 1000)
#   boot_lambdas_comp_Lin <- rep(NA, 1000)
#   boot_lambdas_cont_Quad <- rep(NA, 1000)
#   boot_lambdas_comp_Quad <- rep(NA, 1000)
#   
#   # substitute in new establishment probability
#   f.params$est.prob <- i
#   
#   # Reconstruct the emergence column. This is the only parameter
#   # affected by establishment probability, so first we test out
#   # sensitivity to that. Then, bootstrap the data set with the
#   # new establishment probability 
#   FCol <- h * SB.Emerge(Y, f.params)
#   
#   # # build full K kernel -----------------------------------
#   K_cont_Lin <- rbind(FRow, cbind(FCol, P_Cont_Lin))
#   K_cont_Quad <- rbind(FRow, cbind(FCol, P_Cont_Quad))
#   K_comp_Lin <- rbind(FRow, cbind(FCol, P_Comp_Lin))
#   K_comp_Quad <- rbind(FRow, cbind(FCol, P_Comp_Quad))
#   
#   eigen_cont_Lin <- eigen(K_cont_Lin)
#   lambda_cont_Lin <- max(Re(eigen_cont_Lin$values))
#   AllLambdas$lambda_cont_Lin[it] <- lambda_cont_Lin
#   
#   eigen_comp_Lin <- eigen(K_comp_Lin)
#   lambda_comp_Lin <- max(Re(eigen_comp_Lin$values))
#   AllLambdas$lambda_comp_Lin[it] <- lambda_comp_Lin
#   
#   eigen_cont_Quad <- eigen(K_cont_Quad)
#   lambda_cont_Quad <- max(Re(eigen_cont_Quad$values))
#   AllLambdas$lambda_cont_Quad[it] <- lambda_cont_Quad
#   
#   eigen_comp_Quad <- eigen(K_comp_Quad)
#   lambda_comp_Quad <- max(Re(eigen_comp_Quad$values))
#   AllLambdas$lambda_comp_Quad[it] <- lambda_comp_Quad
#   
#   # Once we have an observed lambda for each value est.prob, we need to bootstrap
#   # the dataset 1000 times to generate confidence intervals. We create new regressions
#   # from each new data set and then construct P and F kernels for them
#   for(j in 1:1000){
#     
#     # Resample raw data sets, re-fit survival, growth and pr(Repro values)
#     x1 <- sample(1:dim(AllControl)[1], dim(AllControl)[1], replace = TRUE)
#     x2 <- sample(1:dim(AllCR)[1], dim(AllCR)[1], replace = TRUE)
#     x3 <- sample(1:dim(AllBig)[1], dim(AllBig)[1], replace = TRUE)
#     
#     BootCont <- AllCont[x1, ]
#     BootCR <- AllCR[x2, ]
#     BootBig <- AllBig[x3, ]
#     
#     # Growth
#     BootCRGrow <- lm(HeightNext ~ Height, data = rbind(BootCR, BootBig))
#     BootContGrow <- lm(HeightNext ~ Height, data = rbind(BootCont, BootBig))
#     
#     # Survival
#     BootCrSurvLin <- glm(Survival ~ Height, 
#                          data = rbind(BootCR, BootBig),
#                          family = binomial())
#     BootContSurvLin <- glm(Survival ~ Height, 
#                            data = rbind(BootCont, BootBig),
#                            family = binomial())
#     BootCrSurvQuad <- glm(Survival ~ Height + I(Height^2),
#                           data = rbind(BootCR, BootBig),
#                           family = binomial())
#     BootContSurvQuad <- glm(Survival ~ Height + I(Height^2),
#                             data = rbind(BootCont, BootBig),
#                             family = binomial())
#     
#     # pr(Repro) and seedling size distributions
#     AllData <- rbind(BootCR, BootCont, BootBig)
#     
#     BootReproData <- filter(AllData, Survival != 'NA')
#     BootReproGLM <- glm(Repro ~ Height, 
#                         data = BootReproData,
#                         family = binomial())
#     
#     BootSdls <- filter(AllData, Stage14 == 'SDL')
#     BootSdlMean <- mean(BootSdls$Height, na.rm = TRUE)
#     BootSdlSD <- sd(BootSdls$Height, na.rm = TRUE)
#     
#     # Combine all the bootstrapped fecundity parameters, but use
#     # the same establishment probability as the outer loop to generate
#     # a confidence interval for each level of it.
#     BootFecParams <- data.frame(prob.repro.slope = coef(BootReproGLM)[2],
#                                 prob.repro.int = coef(BootReproGLM)[1],
#                                 recruit.size.mean = BootSdlMean,
#                                 recruit.size.sd = BootSdlSD,
#                                 germ = germ.prob,
#                                 est.prob = i)
#     
#     
#     boot_Gmat_comp <- h * (outer(Y, Y,
#                                  FUN = grow.prob,
#                                  model = BootCRGrow))
#     
#     boot_Gmat_cont <- h * (outer(Y, Y, 
#                                  FUN = grow.prob,
#                                  model = BootContGrow))
#     
#     boot_Gmat_cont <- boot_Gmat_cont/matrix(as.vector(apply(boot_Gmat_cont,
#                                                             2,
#                                                             sum)),
#                                             nrow = S,
#                                             ncol = S,
#                                             byrow = TRUE)
#     
#     boot_Gmat_comp <- boot_Gmat_comp/matrix(as.vector(apply(boot_Gmat_comp,
#                                                             2,
#                                                             sum)),
#                                             nrow = S,
#                                             ncol = S,
#                                             byrow = TRUE)
#     
#     boot_s.cont_Lin <- predict.surv(Y, BootContSurvLin)
#     boot_s.cont_Quad <- predict.surv(Y, BootContSurvQuad)
#     boot_s.comp_Lin <- predict.surv(Y, BootCrSurvLin)
#     boot_s.comp_Quad <- predict.surv(Y, BootCrSurvQuad)
#     
#     boot_P_Cont_Lin <- boot_P_Cont_Quad <- boot_P_Comp_Lin <- boot_P_Comp_Quad <- matrix(0, S, S)
#     
#     # make a P matrix (growth and survival)
#     for(x in seq_len(S)){
#       boot_P_Cont_Lin[ ,x] <- boot_Gmat_cont[ ,x] * boot_s.cont_Lin[x] 
#       boot_P_Cont_Quad[ ,x] <- boot_Gmat_cont[ ,x] * boot_s.cont_Quad[x]
#       boot_P_Comp_Lin[ ,x] <- boot_Gmat_comp[ ,x] * boot_s.comp_Lin[x]
#       boot_P_Comp_Quad[ ,x] <- boot_Gmat_comp[ ,x] * boot_s.comp_Quad[x]
#     }
#     # Build fecundity vectors
#     boot_FRow <- c(0, SB.Go(Y, BootFecParams,
#                             repro.model = BootReproGLM,
#                             fec.model = fec))
#     boot_FCol <- h * SB.Emerge(Y, BootFecParams)
#     
#     # Build full kernel
#     boot_K_cont_Lin <- rbind(boot_FRow, cbind(boot_FCol, boot_P_Cont_Lin))
#     boot_K_cont_Quad <- rbind(boot_FRow, cbind(boot_FCol, boot_P_Cont_Quad))
#     boot_K_comp_Lin <- rbind(boot_FRow, cbind(boot_FCol, boot_P_Comp_Lin))
#     boot_K_comp_Quad <- rbind(boot_FRow, cbind(boot_FCol, boot_P_Comp_Quad))
#     
#     # Calculate lambdas--------------------------------
#     boot_eigen_cont_Lin <- eigen(boot_K_cont_Lin)
#     boot_lambdas_cont_Lin[j] <- max(Re(boot_eigen_cont_Lin$values))
#     
#     boot_eigen_cont_Quad <- eigen(boot_K_cont_Quad)
#     boot_lambdas_cont_Quad[j] <- max(Re(boot_eigen_cont_Quad$values))
#     
#     boot_eigen_comp_Lin <- eigen(boot_K_comp_Lin)
#     boot_lambdas_comp_Lin[j] <- max(Re(boot_eigen_comp_Lin$values))
#     
#     boot_eigen_comp_Quad <- eigen(boot_K_comp_Quad)
#     boot_lambdas_comp_Quad[j] <- max(Re(boot_eigen_comp_Quad$values))
#     
#     # Store boot strap values for vital rates that don't depend
#     # on EstProb. We only do this once, as we do not need to know
#     # the values of all 100,000 bootstrap samples to get a confidence interval 
#     # estimate 
#     
#     if(it == 1) {
#       OutputValues$GrowSlope_CR[j + 1] <- coef(BootCRGrow)[2]
#       OutputValues$GrowInt_CR[j + 1] <- coef(BootCRGrow)[1]
#       OutputValues$GrowSlope_Cont[j + 1] <- coef(BootContGrow)[2]
#       OutputValues$GrowInt_Cont[j + 1] <- coef(BootContGrow)[1]
#       OutputValues$RepSlope[j + 1] <- coef(BootReproGLM)[2]
#       OutputValues$RepInt[j + 1] <- coef(BootReproGLM)[1]
#       OutputValues$RecMean[j + 1] <- BootSdlMean
#       OutputValues$RecSD[j + 1] <- BootSdlSD
#       OutputValues$SurvInt_Cont_Lin[j + 1] <- coef(BootContSurvLin)[1]
#       OutputValues$SurvSlope_Cont_Lin[j + 1] <- coef(BootContSurvLin)[2]
#       OutputValues$SurvInt_Cont_Quad[j + 1] <- coef(BootContSurvQuad)[1]
#       OutputValues$SurvSlope_Cont_Quad[j + 1] <- coef(BootContSurvQuad)[2]
#       OutputValues$SurvSlope2_Cont_Quad[j + 1] <- coef(BootContSurvQuad)[3]
#       OutputValues$SurvInt_CR_Lin[j + 1] <- coef(BootCrSurvLin)[1]
#       OutputValues$SurvSlope_CR_Lin[j + 1] <- coef(BootCrSurvLin)[2]
#       OutputValues$SurvInt_CR_Quad[j + 1] <- coef(BootCrSurvQuad)[1]
#       OutputValues$SurvSlope_CR_Quad[j + 1] <- coef(BootCrSurvQuad)[2]
#       OutputValues$SurvSelope2_CR_Quad[j + 1] <- coef(BootCrSurvQuad)[3]
#       OutputValues$Obs_Boot[j + 1] <- 'Boot'
#     }
#     
#   } # inner bootstrapping loop
#   
#   # Store confidence intervals from inner loop
#   AllLambdas$lambda_comp_Lin_up[it] <- boot_lambdas_comp_Lin %>%
#     sort %>%
#     .[975]
#   AllLambdas$lambda_comp_Lin_lo[it] <- boot_lambdas_comp_Lin %>%
#     sort %>%
#     .[25]
#   AllLambdas$lambda_comp_Quad_up[it] <- boot_lambdas_comp_Quad %>%
#     sort %>%
#     .[975]
#   AllLambdas$lambda_comp_Quad_lo[it] <- boot_lambdas_comp_Quad %>%
#     sort %>%
#     .[25]
#   
#   AllLambdas$lambda_cont_Lin_up[it] <- boot_lambdas_cont_Lin %>%
#     sort %>%
#     .[975]
#   AllLambdas$lambda_cont_Lin_lo[it] <- boot_lambdas_cont_Lin %>%
#     sort %>%
#     .[25]
#   AllLambdas$lambda_cont_Quad_up[it] <- boot_lambdas_cont_Quad %>%
#     sort %>%
#     .[975]
#   AllLambdas$lambda_cont_Quad_lo[it] <- boot_lambdas_cont_Quad %>%
#     sort %>%
#     .[25]
#   
#   if(it %% 5 == 0){
#     base::message(paste0(it, "% of data processed"))
#   }
#   
#   
#   
# } # Estprob perturbation loop
#
# Fix a typo in the name of a column
# names(OutputValues)[18] <- 'SurvSlope2_CR_Quad'

# OutputValues %>%
#   as_tibble() %>%
#   write.csv(., file = 'Ligustrum_IPM/Ligustrum_Bootstrap_Output.csv',
#             row.names = FALSE)
# 
# write.csv(AllLambdas, 'Ligustrum_IPM/Ligustrum_All_Lambdas.csv',
#           row.names = FALSE)
# 

# Read in the stored outputs from the loop above
AllLambdas <- read.csv('Ligustrum_IPM/Ligustrum_All_Lambdas.csv',
                       stringsAsFactors = FALSE) %>%
  .[-c(7), ]

# The vital rate data need some re-shaping to work in ggplot2
OutputValues <- read.csv('Ligustrum_IPM/Ligustrum_Bootstrap_Output.csv',
                         stringsAsFactors = FALSE) %>%
  gather(key = 'Variable', value = 'Value', -Obs_Boot) %>%
  mutate(Trt = vapply(.$Variable, 
                      FUN = function(x) str_split(x, "_")[[1]][2],
                      FUN.VALUE = ''),
         SurvModel = vapply(.$Variable,
                            FUN = function(x) str_split(x, '_')[[1]][3],
                            FUN.VALUE = '')) %>%
  group_by(Variable, Trt) %>%
  arrange(desc(Value)) %>%
  summarise(Obs = Value[Obs_Boot == 'Observed'],
            UpCI = Value[26],
            LoCI = Value[976]) %>%
  ungroup() %>%
  mutate(Variable = vapply(.$Variable,
                           FUN = function(x) str_split(x, '_')[[1]][1],
                           FUN.VALUE = ''))

write.csv(OutputValues, file = 'Ligustrum_IPM/Ligustrum_Summarized_Output.csv', 
          row.names = FALSE)


ggplot(AllLambdas, aes(x = est.prob)) +
  theme(panel.background = element_blank(), # Set up plot thematic elements
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.15, 0.9), # legend in top left corner
        legend.title = element_text(size = 22), # make legend text bigger
        legend.text = element_text(size = 18),
        axis.title.y = element_text(size = 22,
                                    margin = margin(t = 0, # pad the axis title with a bit of whitespace
                                                    r = 20,
                                                    b = 0,
                                                    l = 0)),
        axis.title.x = element_text(size = 22,
                                    margin = margin(t = 20,
                                                    r = 0,
                                                    b = 0,
                                                    l = 0)),
        axis.line = element_line(colour = 'black', # make axis lines a bit bolder
                                 size = 2),
        axis.ticks = element_line(size = 1.2),
        axis.text = element_text(size = 16)) + 
  geom_line(aes(y = lambda_comp_Quad,
                color = 'CR'),
            size = 1.5,
            alpha = 1,
            linetype = 'dashed') + 
  geom_hline(yintercept = 1,
             color = 'grey50',
             alpha = 0.7,
             size = 2,
             linetype = 'dashed') + 
  geom_ribbon(aes(ymin = lambda_comp_Quad_lo, # add shaded confidence intervals
                  ymax = lambda_comp_Quad_up),
              fill = 'green',
              alpha = 0.2) + 
  geom_line(aes(y = lambda_cont_Quad,
                color = 'Control'),
            alpha = 1,
            size = 1.5,
            linetype = 'dashed') + 
  geom_ribbon(aes(ymin = lambda_cont_Quad_lo,
                  ymax = lambda_cont_Quad_up),
              fill = 'black',
              alpha = 0.2) + 
  scale_y_continuous(expression(paste('Per-capita Growth Rate (', 
                                      lambda, ')')),
                     breaks = seq(1, 1.8, 0.2)) + 
  scale_x_continuous(expression(paste('Establishment Probability (',
                                      italic(E[p]), ')'))) + 
  scale_color_manual('Treatment',
                     breaks = c('Control', 'CR'),
                     values = c('black', 'green')) + 
  geom_text(aes(x = 1, y = 1.8, label = 'B'), size = 16) 



ggsave(filename = 'Lambda_EstP.png',
       path = 'Ligustrum_IPM/',
       height = 8,
       width = 8,
       units = 'in')

# brief sanity check
all(OutputValues$Obs > OutputValues$LoCI &
      OutputValues$Obs < OutputValues$UpCI)

# Fix a dumb typo from my bootstrapping code.
OutputValues$Variable[OutputValues$Variable == 'SurvSelope2'] <- 'SurvSlope2'


# Start renaming parameters so figure titles look pretty
OutputValues$Trt[is.na(OutputValues$Trt)] <- 'Pooled'
OutputValues$Trt[OutputValues$Trt == 'Cont'] <- 'Control'
OutputValues$Variable[OutputValues$Variable == 'RecMean'] <- "paste(mu, ' Recruit Size Mean' )"
OutputValues$Variable[OutputValues$Variable == 'RecSD'] <- "paste(sigma, ' Recruit Size SD')"
OutputValues$Variable[OutputValues$Variable == 'RepInt'] <- "paste( italic(f[p](x)), ' Intercept')"
OutputValues$Variable[OutputValues$Variable == 'RepSlope'] <- "paste(italic(f[p](x)), ' Slope')"
OutputValues$Variable[OutputValues$Variable == 'SurvInt'] <- "paste(italic(s(x)),' Intercept')"
OutputValues$Variable[OutputValues$Variable == 'SurvSlope'] <- "paste(italic(s(x)),' Linear Term')"
OutputValues$Variable[OutputValues$Variable == 'SurvSlope2'] <- "paste(italic(s(x)),' Quadratic Term')"
OutputValues$Variable[OutputValues$Variable == 'GrowInt'] <- "paste(italic(g(y,x)),' Intercept')"
OutputValues$Variable[OutputValues$Variable == 'GrowSlope'] <- "paste(italic(g(y,x)),' Slope')"



ggplot(data = OutputValues,
       aes(x = Trt)) + 
  geom_point(aes(y = Obs, 
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
        # legend.position = c(0.85, 0.1), # legend in top bottom right corner
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
  scale_color_manual('A\nTreatment',
                     breaks = c('Control', 'CR', "Pooled"),
                     values = c('black', 'green', 'blue'))

ggsave(filename = 'Ligustrum_Vital_Rate_Coefficients.png',
       path = 'Ligustrum_IPM/',
       height = 9,
       width = 11,
       unit = 'in')




