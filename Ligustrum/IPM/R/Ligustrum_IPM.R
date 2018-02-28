#####################################################################################
### Ligustrum Exploration and IPM ######################################
#####################################################################################

# All plants are modeled as ramets. will switch to genets if needed
# 1. Growth: splines look good when separated by treatment and then add in big plants.
#    sticking with those for now. Tried a linear model but found that slopes were all < 1 which implied
#    plants shrink on average. That didn't sound reasonable too me.
# Rae: Why use a 200 cm cutoff for adding all large? It appears that 165-175 might be a slightly better place to start.
# Then you could say all RAs were used in estimating growth for each treatment. This is something we might want to discuss more.
# The larger plants seem to drive the slope where as the smaller plants drive the intercept. Are we biasing the slope by lumping?
#   
# 2. Survival: Doesn't look like much of a treatment effect, but will model them separately for now
#    to see if it affects lambda. I have also written code to combine all treatments so combining them later
#    would not take very long. 
#     -For now, keeping plants <200 separate, combining above that.
# 3. Fecundity: We have fruit production, probability of reproduction
#    recruit size distribution (assumed no treatment effect since we have no data on that from 2015). 
#     -Use Spline for fruit production, binom for Pr(repro), normal for recruit size dist.
#       All are pooled across trts
#
# 4. Clones: Competitor removal seems to produce smaller clones than the others. Will investigate other parameters
#     -No relationship between size of parents and whether they clone. will be a fixed mean probability
#     -Marginally significant relationship between size and number of clones in quasipoisson. Normal poisson
#      model appears to be underdispersed, and is non-significant. Visual examination appears to indicate
#      there is a relationship, but there are large outliers too. Next, tried Splines and Zero-Truncated Poisson.
#      Zero-Truncated Poisson works, and actually makes sense, since a plant that is considered clonal cannot
#      produce less than 1 clone. Using 0-truncated Pois for now.
#     -Clone size distributions look different, modeling separately
#   OMITTING CLONES FOR NOW DUE TO POSSIBLE ISSUES WITH UNDERLYING DATA COLLECTION METHODS

rm(list = ls()) 
graphics.off()

### LIBRARIES ###
# library(rstan)
library(dplyr)
library(magrittr)
library(stringr)
library(mgcv)


source('Ligustrum/IPM/R/IPM_Functions_Ligustrum.r') 

# load the data and do some basic restructuring
RAs<-read.csv("Ligustrum/IPM/Data/LO_RA_Clean.csv")
# str(RAs)
AllPlants<-read.csv("Ligustrum/IPM/Data/LO_Clean.csv",
                    stringsAsFactors = FALSE) 

AllPlants$Plant_Number <- as.numeric(AllPlants$Plant_Number)
# str(AllPlants)
AllPlants$growth <- AllPlants$Plant_Height15 - AllPlants$Plant_Height14
AllPlants <- AllPlants %>% arrange(Plot,Plant_Number)




# Growth---------------------------------------------------------------------------------------
# Linear models. 
Big.Plants <- filter(AllPlants, Treatment == "All")

growth.lms <- AllPlants %>% 
  filter(Treatment!="All") %>% 
  group_by(Treatment) %>%
  do(LM = lm(Plant_Height15 ~ Plant_Height14, data = rbind(., Big.Plants)),
     SumLM = summary(lm(Plant_Height15 ~ Plant_Height14, data = rbind(. ,Big.Plants))),
     GAM = gam(Plant_Height15 ~ s(Plant_Height14), data = rbind(., Big.Plants)))


par(mfrow = c(1, 1))
xx <- seq(0, max(AllPlants$Plant_Height14, na.rm = TRUE), .1)
plot(Plant_Height15 ~ Plant_Height14, data = AllPlants)
abline(growth.lms$LM[[2]], col = 'black')
abline(growth.lms$LM[[1]], col = 'green')
lines(xx, predict(growth.lms$GAM[[1]], 
                  data.frame(Plant_Height14 = xx),
                  type = 'response'),
      col = 'green', lty = 2)
lines(xx, predict(growth.lms$GAM[[2]], 
                  data.frame(Plant_Height14 = xx),
                  type = 'response'),
      col = 'black', lty = 2)

abline(1:500,1:500,col='red')
legend('topleft',c("Control", "CR",
                   "Control GAM", "CR GAM",
                   "1:1 Growth Line"),
       col = c('black', 'green', 
               'black', 'green',
               'red'),
       lty = c(1, 1, 2, 2, 1))

growth.lms$SumLM[[1]]
growth.lms$SumLM[[2]]

AIC(growth.lms$LM[[1]], growth.lms$GAM[[1]])
AIC(growth.lms$LM[[2]], growth.lms$GAM[[2]])

AllGrow <- lm(Plant_Height15 ~ Plant_Height14, data = AllPlants)
summary(AllGrow)

# Survival across treatments----------------------------------------------------------------------
# Logistic regression first with all plants that are classed by treatment.
AllSmall <- filter(AllPlants, Treatment != 'All')
AllBig <- filter(AllPlants, Treatment == 'All')
AllControl <- filter(AllPlants, Treatment == 'Control')
AllCR <- filter(AllPlants, Treatment == 'Comp')


# There is no way to model all of them 
# together, but we'll see what their regressions look like separately.
ContLinearGlm <- glm(Survival ~ Plant_Height14,
                     data = rbind(AllControl, AllBig),
                     family = binomial())
CRLinearGlm <- glm(Survival ~ Plant_Height14,
                   data = rbind(AllCR, AllBig),
                   family = binomial())

# Both return warnings of probabilities = 0|1. trying Quadratic fits now
ContQuadGlm <- glm(Survival ~ Plant_Height14 + I(Plant_Height14^2),
                     data = rbind(AllControl, AllBig),
                     family = binomial())
CRQuadGlm <- glm(Survival ~ Plant_Height14 + I(Plant_Height14^2),
                   data = rbind(AllCR, AllBig),
                   family = binomial())
summary(ContLinearGlm)
summary(CRLinearGlm)
summary(ContQuadGlm)
summary(CRQuadGlm)

AIC(ContLinearGlm, ContQuadGlm)
AIC(CRLinearGlm,CRQuadGlm)
# Plot fits to see how they look
xx <- seq(0, max(AllPlants$Plant_Height15, na.rm = TRUE) + 50, 1)
plot(Survival ~ Plant_Height14, data = AllPlants,
     xlim = c(0, 600))
lines(xx, predict(ContLinearGlm,
                  data.frame(Plant_Height14 = xx),
                  type = 'response'),
      col = 'black', 
      lty = 1)
lines(xx, predict(ContQuadGlm, 
                  data.frame(Plant_Height14 = xx),
                  type = 'response'),
      col = 'black',
      lty = 2)
lines(xx, predict(CRLinearGlm,
                  data.frame(Plant_Height14 = xx),
                  type = 'response'),
      col = 'green',
      lty = 1)
lines(xx, predict(CRQuadGlm,
                  data.frame(Plant_Height14 = xx),
                  type = 'response'),
      col = 'green',
      lty = 2)
legend('bottomright',
       c('Control', 'Control Quad', 'CR', 'CR Quad'),
       col = c('black', 'black', 'green', 'green'),
       lty = c(1, 2, 1, 2))


# Fecundity-------------------------------------------------------------------------------
# we are missing a recruit size distribution from year two. For now, I am using the initial
# seedling size distribution and seeing how that goes. We can figure out something else
# later. Modeled without a seedbank based on work on congeners by Shelton and Cain 2002

# first, we will estimate fecundity from our sample of 10 RAs
fec <- glm(Seeds ~ Height,
           data = RAs,
           family = quasipoisson())
summary(fec)
# normal poisson is very overdispersed, trying lm and splines next
feclm <- lm(Seeds ~ Height,
            data = RAs)
summary(feclm)

fec.gam <- gam(Seeds ~ s(Height),
               data = RAs)
summary(fec.gam)

par(mfrow = c(1,1))
xx <- seq(0, max(RAs$Height), 1)
plot(Seeds ~ Height, data = RAs)
lines(xx, predict(fec, data.frame(Height = xx), type = 'response'),
      col = 1)
abline(feclm, col = 2)
lines(xx, predict(fec.gam, data.frame(Height = xx), type = 'response'),
      col = 3)
legend('topleft', c('Poisson','Linear','Spline'), 
       col = c(1, 2, 3), lty = 1)

# Rae: The spline is returning a negative intercept.

#Pr(Reproductive) - Logistic regression ----------------------------------
AllPlants$Repro <- ifelse(AllPlants$Stage15 == "RA", 1, 0)

Regression.Data <- filter(AllPlants, Survival != "NA")

Repro.Glm <- glm(Repro ~ Plant_Height14,
                 data = Regression.Data,
                 family = binomial())
summary(Repro.Glm)
plot(Repro ~ Plant_Height14,
     data = Regression.Data)
xx<-seq(0, max(AllPlants$Plant_Height14, na.rm = TRUE), 0.1)
lines(xx, predict(Repro.Glm, data.frame(Plant_Height14 = xx), type = 'response'),
      lty = 2, col = 'red')




# Recruit size distribution----------------------------------------------------------------------
# Not really sure how to do this with no data from 2015. Using 2014 size distribution for now,
# but there's really no way to separate out the treatments, and I think this is the only place
# we could potentially see treatment effect. However, we may be able to justify this by pointing out
# that there was no difference in growth or survival for small plants between the two censuses,
# so maybe we're ok
sdls <- filter(AllPlants, Stage14 == "SDL")
Sdl.mean <- mean(sdls$Plant_Height14, na.rm = TRUE)
Sdl.SD <- sd(sdls$Plant_Height14, na.rm = TRUE)

xx <- seq(0, max(sdls$Plant_Height14, na.rm = TRUE), 0.5)
hist(sdls$Plant_Height14,
     breaks = 1,
     border = 'white',
     xlab = 'Height 2014',
     freq = FALSE,
     ylim = c(0, max(dnorm(xx, Sdl.mean, Sdl.SD)) + 0.05))
lines(xx, dnorm(xx, Sdl.mean, Sdl.SD), 
      col = 'red', lty = 2)

# no treatment effect on seedling growth or survival, so I think we can justify saying no treatment
# effect on seedling size distribution
sdl.tab <- table(sdls$Treatment, sdls$Survival)
sdl.chisq <- chisq.test(sdl.tab, correct = FALSE)
sdl.chisq
sdl.tab
sdl.chisq$exp

sdl.grow.aov <- aov(growth ~ Treatment,
                    data = AllPlants[AllPlants$Stage14 == "SDL", ])
summary(sdl.grow.aov)
boxplot(growth ~ Treatment,
        data = AllPlants[AllPlants$Stage14 == "SDL", ])

# Figures for Appendix
xx <- seq(0, max(AllPlants$Plant_Height14, na.rm = TRUE), 0.1)
par(mfrow = c(2,3), xpd = TRUE)
plot(Plant_Height15 ~ Plant_Height14, data = AllPlants,
     xlab = 'Size (t)',
     ylab = 'Size (t+1)')
lines(xx, predict(growth.lms$LM[[2]],
                  data.frame(Plant_Height14 = xx), 
                  type = 'response'),
      col = 'black',
      lty = 2)
lines(xx, predict(growth.lms$LM[[1]], 
                  data.frame(Plant_Height14 = xx),
                  type = 'response'),
      col = 'green',
      lty = 2)

xx <- seq(0, max(AllPlants$Plant_Height14, na.rm = TRUE) + 30, 0.1)
plot(Survival ~ Plant_Height14, data = AllPlants,
     xlab = 'Size (t)',
     ylab = 'Survival (t+1)',
     xlim = c(0, max(AllPlants$Plant_Height14, na.rm = TRUE) + 30))
lines(xx, predict.surv(xx, CRLinearGlm),
      col = 'green', lty = 1)
lines(xx, predict.surv(xx, CRQuadGlm),
      col = 'green', lty = 2)
lines(xx, predict.surv(xx, ContLinearGlm),
      col = 'black', lty = 1)
lines(xx, predict.surv(xx, ContQuadGlm),
      col = 'black', lty = 2)
legend('bottomright', c('CR Linear',
                        'CR Quad',
                        'Cont Linear',
                        'Cont Quad'),
       lty = c(1, 2, 1, 2),
       col = c('green', 'green',
               'black', 'black'))

xx <- seq(min(RAs$Height), max(RAs$Height) - 5, 1)
plot(Seeds ~ Height, data = RAs,
     xlab = 'Size (t+1)',
     ylab = 'Seeds (t+1)')
lines(xx, predict(fec, 
                  data.frame(Height = xx),
                  type = 'response'),
      col = 'red',
      lty = 2)

xx <- seq(0, max(AllPlants$Plant_Height14, na.rm = TRUE), .1)

plot(Repro ~ Plant_Height14, data = AllPlants,
     xlab = 'Size (t)',
     ylab = 'Reproductive (t+1)')
lines(xx, predict(Repro.Glm, data.frame(Plant_Height14 = xx),
                  type = 'response'),
      col = 'red',
      lty = 2)

xx <- seq(0, max(sdls$Plant_Height14, na.rm = TRUE), 0.5)
hist(sdls$Plant_Height14,
     breaks = 1,
     border = 'white',
     xlab = 'Recruit size (t)',
     freq = FALSE,
     ylim = c(0, max(dnorm(xx, Sdl.mean, Sdl.SD)) + 0.05),
     main = '')
lines(xx, dnorm(xx, Sdl.mean, Sdl.SD),
      col = 'red',
      lty = 2)

legend(x = 33.25, y = 0.23,
       'Panel A',
       bty = 'n',
       xpd = NA,
       cex = 1.3)
legend(x = 30, y = 0.2, c('Treatment','Control', 'Competitor removal',
                           'Pooled values'),
       lty = c(NA, 2, 2, 2), col = c(NA,'black', 'green', 'red'),
       xpd = NA, bty = 'n',
       cex = 1.3)

# Clonal reproduction-----------------------------------

# For.Cloning<-filter(AllPlants,!is.na(Clone))
# # add column indicating parent plant
# for(i in 1:length(For.Cloning$Quadrat)){
#   if(For.Cloning$Clone[i] !=0){
#     For.Cloning$Parent[i] <- For.Cloning$Plant_Number[i]%>%str_sub(1,-2)
#   }else{
#     For.Cloning$Parent[i]<-NA
#   }
# }
# 
# For.Cloning$Parent<-as.numeric(For.Cloning$Parent)
# For.Cloning$QuadPlant<-paste(For.Cloning$Quadrat,For.Cloning$Plant_Number,sep="-")
# 
# # calculate number of clones produced by each parent and add that to For.Cloning data frame
# Clone.Count<-For.Cloning%>%filter(Parent!="NA")%>%group_by(Quadrat,Parent)%>%summarise(Clone.N=n())
# Clone.Count$QuadPlant<-paste(Clone.Count$Quadrat,Clone.Count$Parent,sep="-")
# 
# For.Cloning<-left_join(For.Cloning,Clone.Count,by=c("QuadPlant"))%>%select(-Parent.x,
#                                                                            -Quadrat.y,-Parent.y,-QuadPlant)
# 
# For.Cloning$Clonal<-ifelse(!is.na(For.Cloning$Clone.N),1,0)
# For.Cloning$Clonal[For.Cloning$Clone==1]<-NA #the clones themselves wouldn't be expected to be clonal...
# 
# # Pr(Clonal) by treatment-------------------------------
# tab<-table(For.Cloning$Treatment,For.Cloning$Clonal)
# chisq.test(tab)
# tab
# chisq.test(tab)$exp
# # Pr(clonal) does not differ by treatment, so we will lump them all together.
# Pr.Clonal.Reg<-glm(Clonal~Plant_Height14,data=For.Cloning, family=binomial())
# summary(Pr.Clonal.Reg)
# 
# xx<-seq(0,max(For.Cloning$Plant_Height14,na.rm=T),1)
# plot(Clonal~Plant_Height14,data=For.Cloning)
# lines(xx,predict(Pr.Clonal.Reg,data.frame(Plant_Height14=xx),type='response'),col='red')
# 
# # it would appear that initial size is not a useful predictor of whether or not a plant clones, so this will be
# # a discrete parameter
# Pr.Clonal<-mean(For.Cloning$Clonal,na.rm=T)
# 
# # number of clones by size and treatment
# Clone.Treatment.N<-aov(Clone.N~Treatment,data=For.Cloning)
# summary(Clone.Treatment.N)
# 
# # Treatment has no strong effect on number of clones produced. Looking to see if size is a reasonable
# # predictor when treatments are lumped together
# 
# par(mfrow=c(1,1))
# plot(Clone.N~Plant_Height14,data=For.Cloning,xlim=c(0,200))
# Clone.N.Pois<-glm(Clone.N~Plant_Height14,data=For.Cloning,family=quasipoisson())
# Clone.N.Spline<-gam(Clone.N~s(Plant_Height14),data=For.Cloning)
# Clone.N.ZTPois<-vglm(Clone.N~Plant_Height14,data=For.Cloning,family=pospoisson)
# summary(Clone.N.Pois)
# summary(Clone.N.Spline)
# summary(Clone.N.ZTPois)
# lines(xx,predict(Clone.N.Spline,data.frame(Plant_Height14=xx),type='response'),col='green')
# lines(xx,predict(Clone.N.Pois,data.frame(Plant_Height14=xx),type='response'),col='red')
# lines(xx,predict(Clone.N.ZTPois,data.frame(Plant_Height14=xx),type='response'),col='blue')
# legend('topleft',c('Spline','Quasi-Poisson','Zero-truncated Poisson'),col=c('green','red','blue'),lty=1)
# 
# # Clone size distribution-----------------
# Clones<-filter(For.Cloning,Clone==1)
# 
# cont.clone.mean<-mean(Clones$Clone.Height[Clones$Treatment=="Control"])
# cont.clone.sd<-sd(Clones$Clone.Height[Clones$Treatment=="Control"])
# comp.clone.mean<-mean(Clones$Clone.Height[Clones$Treatment=="Comp"])
# comp.clone.sd<-sd(Clones$Clone.Height[Clones$Treatment=="Comp"])
# herb.clone.mean<-mean(Clones$Clone.Height[Clones$Treatment=="Herb"])
# herb.clone.sd<-mean(Clones$Clone.Height[Clones$Treatment=="Herb"])
# 
# par(mfrow=c(1,1))
# xx<-seq(0,max(Clones$Clone.Height),1)
# hist(Clones$Clone.Height,main="",breaks=1,border='white',freq=F,ylim=c(0,.03))
# lines(xx,dnorm(xx,cont.clone.mean,cont.clone.sd),col=1)
# lines(xx,dnorm(xx,comp.clone.mean,comp.clone.sd),col=2)
# lines(xx,dnorm(xx,herb.clone.mean,herb.clone.sd),col=3)
# legend('topright',c('Control','CR','HR'),col=c(1,2,3),lty=1)
# 
# clone.size.aov<-aov(Clone.Height~Treatment,data=Clones)
# summary(clone.size.aov)
# TukeyHSD(clone.size.aov)
# print(model.tables(clone.size.aov,'mean'))
# boxplot(Clone.Height~Treatment,data=Clones)


# Discrete parameters------------------------------------

# Establishment Pr()
# germination will come from our experiment and will be the average of all non-acid
# treated cells. Establishment will have to come from somewhere else and I am really 
# not sure where given the dearth of literature on the demography of this species.
# Ramula 2008 used simulated parameters in lieu of certain vital rates, so maybe 
# we should investigate that further. One idea I have is to try a variety of 
# values and see which ones make the most sense. However we'd still be guessing 
# at the most sensible value since we don't have much by way of recruitment info in 
# 2015

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
min.size <- min(AllPlants$Plant_Height14, na.rm = TRUE) * .6
max.size <- max(AllPlants$Plant_Height14, na.rm = TRUE) * 1.2
S <- 500 # Number of cells in matrix  

# matrix variables 
b <-  min.size + c(0:S)*(max.size - min.size)/S  # boundary points of mesh cells
Y <- 0.5* (b[1:S] + b[2:(S+1)])  # mid points of mesh cells 
h <- Y[2]-Y[1]  # cell widths

# make a matrix of transitions for each growth distribution-------------
source('C:/Users/sl13sise/Dropbox/ATSC 2018 participant folder/23.1.18/Rees/R code and Data/MatrixImage.R')

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

matrix.image(Gmat_cont, main = 'Control')
matrix.image(Gmat_comp, main = 'CR')

par(mfrow=c(1,2))
plot((1-colSums(Gmat_comp)),type='l',main="CR")
abline(h=0,col='red',lty=2)
plot((1-colSums(Gmat_cont)),type='l', main="Cont")
abline(h=0,col='red',lty=2)

s.cont_Lin <- diag(predict.surv(Y, ContLinearGlm))
s.cont_Quad <- diag(predict.surv(Y, ContQuadGlm))
s.comp_Lin <- diag(predict.surv(Y, CRLinearGlm))
s.comp_Quad <- diag(predict.surv(Y, CRQuadGlm))

# make a P matrix (growth and survival)
P_Cont_Lin <- Gmat_cont %*% s.cont_Lin 
P_Cont_Quad <- Gmat_cont %*% s.cont_Quad
P_Comp_Lin <- Gmat_comp %*% s.comp_Lin
P_Comp_Quad <- Gmat_comp %*% s.comp_Quad

# inspect survival + growth kernel
matrix.image(P_Cont_Lin, main = 'Control, Linear SurvModel')
matrix.image(P_Cont_Quad, main = 'Control, Quadratic SurvModel')
matrix.image(P_Comp_Lin, main = 'CR, Linear SurvModel')
matrix.image(P_Comp_Quad, main = 'CR, Quadratic SurvModel')

# Build F Matrix-----------------------------------

par(mfrow = c(1,1))
FRow <- c(0, SB.Go(Y,
                   f.params,
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

lambda_cont_Lin
lambda_cont_Quad
lambda_comp_Lin
lambda_comp_Quad

