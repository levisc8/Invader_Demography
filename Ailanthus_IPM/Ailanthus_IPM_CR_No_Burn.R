### Infile growth, survival, and reproduction data ###

### Infile growth, survival, and reproduction data ###
Ailanthus4R <- read.csv("Ailanthus_IPM/Ailanthus_Clean.csv", header = TRUE)


### Infile germination data ###
Ailanthus4R_Germ <- read.csv("Ailanthus_IPM/Ailanthus_Germ_Clean.csv")
head(Ailanthus4R_Germ)

## Seedling data
AASeedlings <- read.csv("Ailanthus_IPM/Ailanthus_Seedlings_Clean.csv",
                        stringsAsFactors = FALSE)

### Subset data ###
# PopData is used for growth and survival analyses
# Good_PopData is used to model size distribution of measured clones
# PopDataSDL is used to model size distribution of seedlings
# PopDataRA is used to model the relationship between size and probability of reproduction
# PopDataSeeds is used to model relationship between fecundity and size; it only includes reproducing individuals

###########################################################################
# Competitor removal without fire
PopData1 <- Ailanthus4R[Ailanthus4R$Treatment == 'Comp' | Ailanthus4R$Treatment == 'All',]
PopData <- PopData1[PopData1$Burn != 'Y',]

PopDataSDL <- Ailanthus4R[Ailanthus4R$Treatment == 'Comp' & Ailanthus4R$Burn != 'Y' & Ailanthus4R$StageNext == 'SDL',]

##########################################################################
# These subsets are used for all treatments
Good_PopData <- PopData[PopData$Quality =='good',]
Nonclones <- Good_PopData[Good_PopData$Stage !='gone' & Good_PopData$StageNext != 'SDL' & Good_PopData$Stage != 'SDL',]
Clones1 <- PopData[PopData$StageNext != 'gone',]
Clones2 <- Clones1[Clones1$StageNext != 'SDL' & Clones1$Stage != 'SDL',]
Clones <- Clones2[Clones2$Stage == 'gone' & Clones2$Quality == 'good',]

PopDataRA1 = Ailanthus4R[Ailanthus4R$StageNext != 'RA',]  
PopDataRA <- PopDataRA1[PopDataRA1$StageNext != 'Male',]
PopDataSeeds <- Ailanthus4R[Ailanthus4R$StageNext == 'RA' | Ailanthus4R$StageNext == 'Female',]

Germ <- Ailanthus4R_Germ[Ailanthus4R_Germ$Study == 'Germ', ]
Viable <- Ailanthus4R_Germ[Ailanthus4R_Germ$Study == 'Viable',]



p.vec=data.frame(
  surv.int=NA,           # Intercept from logistic regression of survival
  surv.slope=NA,         # Slope from logistic regression of survival
  growth.int=NA,         # Intercept from linear regression of growth
  growth.slope=NA,       # Slope from linear regression of growth
  growth.sd=NA,          # Residual sd from linear regression of growth
  seed.int=NA,           # Intercept from Poisson regression of seed number and size
  seed.slope=NA,         # Slope from Poisson regression of seed number and size
  prob.repro.int=NA,     # Intercept from logistic regression of reproduction
  prob.repro.slope=NA,   # Slope from logistic regression of reproduction
  recruit.size.mean=NA,  # Mean recruit size
  recruit.size.sd=NA,    # Standard deviation of recruit size
  establishment.prob=NA, # Probability of establishment
  clonal.size.mean=NA,   # Mean clone size
  clonal.size.sd=NA,     # Standard deviation of clone size
  clonal.prob=NA,        # Probability of producting a clone
  stay.seedbank.prob=NA, # Probability of staying in the seedbank
  go.seedbank.prob=NA,   # Probablity of entering the seedbank
  germ.seedbank.prob=NA  # Probablity of germinating from the seedbank
)

# 1. Survival: Logistic regression
surv.reg <- glm(Survival ~ Size, data = PopData, family = binomial())
summary(surv.reg)

p.vec$surv.int <- coefficients(surv.reg)[1]
p.vec$surv.slope <- coefficients(surv.reg)[2]

# 2. Growth: Linear regression with x-intercept set to zero
growth.reg <- lm(SizeNext ~ 0 + Size, data = PopData)
summary(growth.reg)
p.vec$growth.int <- 0  
p.vec$growth.slope <- coefficients(growth.reg)[1] 
p.vec$growth.sd <- sd(resid(growth.reg))

# 3. Seeds: Divided into two parameters
# Logistic regression for probability of being reproductive
prob.repro.reg <- glm(Reproductive ~ SizeNext, data = PopDataRA, family = binomial())
summary(prob.repro.reg)
p.vec$prob.repro.int <- coefficients(prob.repro.reg)[1]
p.vec$prob.repro.slope <- coefficients(prob.repro.reg)[2]

#  Poisson regression for relationship between SizeNext and number of seeds
seed.reg=glm(Seeds ~ SizeNext, data = PopDataSeeds, family = poisson())
summary(seed.reg)
p.vec$seed.int <- coefficients(seed.reg)[1]
p.vec$seed.slope <- coefficients(seed.reg)[2]

# 4. Size distribution of seedlings
p.vec$recruit.size.mean <- mean(PopDataSDL$SizeNext)
p.vec$recruit.size.sd <- sd(PopDataSDL$SizeNext)

# 5. Staying in continuous stage (i.e., does not enter seedbank)
k <- 1

v <- mean(Germ$v)
 
g <- mean(Germ$g)
# 
# ### Infile seedling count data ###
# AASeedlings <- read.csv("AA SDL Counts by plot.csv", header = TRUE)
# 
# # Subset data
# SDL_NoBurn <- AASeedlings[AASeedlings$Burn == 'N',] 

# From Crandall and Knight (2018): "q is the total effect of fire (direct and indirect) on seed and early seedling survivorship.
# Because the difference in seedling number was highest for the burned plots, we set q=1 for these 
# plots (the highest possible value). For unburned plots, q=0.160585 (because difference in seedling 
# number for unburned plots was 0.160585*difference in seedling number for the burned plots), which 
# takes into account that unburned plants have lower recruitment compared to burned plants" 
q <- 0.160585

p.vec$establishment.prob <- v * g * q

# 6.  Size distribution of clones
p.vec$clonal.size.mean <- mean(Clones$SizeNext, na.rm = TRUE)
p.vec$clonal.size.sd <- sd(Clones$SizeNext, na.rm = TRUE) * 1.05

# 7. Probability of clonal growth
# Estimated by dividing the number of new recruits (i.e., Clones) in 2013 
# by the number of NRA and RA (i.e., Nonclones) in 2012
p.vec$clonal.prob <- length(Clones$SizeNext)/length(Nonclones$Size)
CompN_clonal.prob <- length(Clones$SizeNext)/length(Nonclones$Size)

# 8.  Probabilities associated with discrete stage (i.e., seedbank)
# Probability of staying in the seedbank
p.vec$stay.seedbank.prob <- v * (1 - g) * k

# Probablity of entering the seedbank from continuous stage (more below)
p.vec$go.seedbank.prob <- v * (1 - g) * k

# Probablity of germinating from the seedbank (more below)
p.vec$germ.seedbank.prob <- v * g * q

### Define functions to describe life history ###

# 1. Survival probability function
s.x.CompN <- function(x,params) {
  u <- exp(params$surv.int+params$surv.slope * x)
  return(u / (1 + u))
}

# 2. Growth function
g.yx.CompN <- function(xp,x,params) {
  dnorm(xp, mean = params$growth.int + params$growth.slope * x, sd = params$growth.sd) 
}

# 3. Reproduction function
p.repro.x.CompN <- function(x,params) {
  u <- exp(params$prob.repro.int + params$prob.repro.slope * x)
  return(u / (1 + u))
}

# Fecundity model while taking into account some seeds move into 
# the discrete portion of the model (i.e., seedbank)
f.yx.CompN <- function(xp,x,params) {
  p.repro.x.CompN(x,params)*
    params$establishment.prob *
    dnorm(xp, mean = params$recruit.size.mean, sd = params$recruit.size.sd) *
    exp(params$seed.int + params$seed.slope * x)
}

# 4. Seedling size distribution
d.x.CompN <- function(xp,x,params) {
  dnorm(xp, mean = p.vec$recruit.size.mean, sd = p.vec$recruit.size.sd) * p.vec$germ.seedbank.prob
}

# 5. Seeds entering seedbank
e.x.CompN <- function(xp,x,params) {
  p.repro.x.CompN(x,params) * exp(params$seed.int + params$seed.slope * x) * params$go.seedbank.prob
}

# 6. Clonal function
c.yx.CompN <- function(xp,x,params) {
    params$clonal.prob *
    dnorm(xp, mean = params$clonal.size.mean, sd = params$clonal.size.sd)
}

### Combine vital rate functions to build the IPM matrix
# Set number of classes (n.size) or points for midpoint rule approximation
min.size <- 0.9 * min(c(PopData$Size, PopData$SizeNext), na.rm = T)
max.size_CompN <- 1.1 * max(c(PopData$Size, PopData$SizeNext), na.rm = T)
n.size <- 100 # number of cells in the matrix 
b <- min.size + c(0:n.size) * (max.size_CompN - min.size) / n.size # boundary points
y <- 0.5 * (b[1:n.size] + b[2:(n.size + 1)]) # mesh points
h <- y[2] - y[1] # step size

### Make IPM matrices ###

# Growth matrix
G_CompN <- h * outer(y, y, g.yx.CompN, params = p.vec)
# Larger individuals are evicted (see Williams et al. 2012), so return the evicted individuals to the 
# cells at the boundaries where they were evicted (i.e., rerout growth to sizes outside the allowed range 
# to the extreme sizes avoiding eviction).
for(i in 1:(n.size/2)) G_CompN[1,i] <- G_CompN[1,i] + 1 - sum(G_CompN[,i])
for(i in (n.size/2+1):n.size) G_CompN[n.size,i] <- G_CompN[n.size,i] + 1 - sum(G_CompN[,i])

# Survival vector
S_CompN <- s.x.CompN(y, params = p.vec) 
SMatrix <- array(0, dim = c(100,99))
S_CompN1 <- cbind(S_CompN, SMatrix)          

# Fecundity martix
F_CompN <- array(0, dim=c(1 + n.size, 1 + n.size))
F_CompN[2:(1 + n.size), 2:(1 + n.size)] <- h * outer(y, y, f.yx.CompN, params = p.vec)
# e.x.CompN tells seeds to go into the seedbank
F_CompN[1, 2:(1 + n.size)] <- e.x.CompN(y, y, params = p.vec)

# Survival/growth matrix
# The first row and column are seedbank transitions
P_CompN <- array(0, dim = c(1 + n.size, 1 + n.size))

# Add probability of staying in seedbank 
P_CompN[1,1] <- p.vec$stay.seedbank.prob 
P_CompN1 <- P_CompN

# Add probability of germinating and moving to continuous
P_CompN[2:(1 + n.size), 1] <- d.x.CompN(y, params = p.vec)

# Build growth/survival matrix including discrete seedbank stage
for(i in 2:(1 + n.size))
  P_CompN[i, 2:(1 + n.size)] <- S_CompN * G_CompN[(i - 1),]

C_CompN <- array(0,dim=c(1 + n.size, 1 + n.size))
C_CompN[2:(1 + n.size), 2:(1 + n.size)] <- h * outer(y, y, c.yx.CompN, params = p.vec) # clonal growth matrix

# Build matrix
K_CompN <- P_CompN + C_CompN + F_CompN

# Calculate lambda
CompN_lambda <- Re(eigen(K_CompN)$values[1])
CompN_lambda

## Bootstrapping code
# Compute annual changes in seedlings numbers. This is used to calculate 
# the establishment probability. The loop below calculates the absolute difference
# in each plot and stores it in a new data frame.
AASdls <- data.frame(Burn = rep(NA, length(unique(AASeedlings$Plot))),
                     Treatment = rep(NA, length(unique(AASeedlings$Plot))),
                     Treatment_X_Burn = rep(NA, length(unique(AASeedlings$Plot))),
                     Plot = rep(NA, length(unique(AASeedlings$Plot))),
                     Change = rep(NA, length(unique(AASeedlings$Plot))))

it <- 1
for(i in unique(AASeedlings$Plot)) {
  
  temp <- subset(AASeedlings, Plot == i)
  
  sdl_12 <- temp$SDL_Count[temp$Year == 2012]
  sdl_13 <- temp$SDL_Count[temp$Year == 2013]
  
  AASdls[it, 1:4] <- temp[1, 1:4]
  AASdls[it, 'Change'] <- abs(sdl_13 - sdl_12)
  it <- it + 1
}


AASeedlingsN = AASdls[AASdls$Burn == 'N',]
AASeedlingsY = AASdls[AASdls$Burn == 'Y',]

### Subset data ###
# PopData is file for growth and survival analyses; uncomment the treatment of interest
# Good.PopData is used to model size distribution of only the good clones
# PopDataSDL is used to model size distribution and establishment probability
# PopDataRA is used to model the relationship between fecundity and size of RAs only
###########################################################################

### Start bootstrapping
CompN_bootlambda = rep(NA, 1000)
p.vec.all <- lapply(p.vec,
                    function(x) x[1])



for(j in 1:1000) 
{
  
  ###########################################################################
  # Competitor removal without fire
  PopData1 = Ailanthus4R[Ailanthus4R$Treatment == 'Comp' | Ailanthus4R$Treatment == 'All',]
  PopDataRAW = PopData1[PopData1$Burn != 'Y' & PopData1$StageNext != 'SDL',]
  
  Clones1 = subset(PopData1, Burn != 'Y' & StageNext != 'gone')
  Clones2 = subset(Clones1, StageNext != 'SDL' & Stage != 'SDL')
  ALLClones = subset(Clones2, Stage == 'gone')
  
  PopDataSDL = subset(Ailanthus4R, (Treatment != 'Herb' & Treatment != 'Control' & Burn != 'Y' & StageNext == 'SDL'))
  
  ###########################################################################
  
  # Use these subsets for all treatments
  PopDataSeeds = Ailanthus4R[Ailanthus4R$StageNext == 'RA' | Ailanthus4R$StageNext == 'Female',]
  
  Germ = subset(Ailanthus4R_Germ, Study == 'Germ')
  

  ###########################################################################
  
  ### start the loop for seedling data
  nseedlings = length(PopDataSDL[ ,1]) 
  x1 = sample(1:nseedlings, nseedlings, replace = TRUE)
  
  # bootstrap dataset
  PopDataSDLBoot = PopDataSDL[x1[1], ]
  for(i in 2:nseedlings){ 
    PopDataSDLBoot = rbind(PopDataSDLBoot, PopDataSDL[x1[i], ])
  }
  PopDataSDL=PopDataSDLBoot
  
  ### start the loop for seeds data
  nSeeds = length(PopDataSeeds[,1])
  x2 = sample(1:nSeeds, nSeeds, replace = TRUE)
  
  # bootstrap dataset
  PopDataSeedsBoot = PopDataSeeds[x2[1], ]
  for(i in 2:nSeeds){ 
    PopDataSeedsBoot = rbind(PopDataSeedsBoot,PopDataSeeds[x2[i], ])
  }
  PopDataSeeds=PopDataSeedsBoot
  
  ### start the loop for NRA and RA data
  nNRA = length(PopDataRAW[,1])
  x3 = sample(1:nNRA, nNRA, replace = TRUE)
  
  # bootstrap dataset
  PopDataBootNRA = PopDataRAW[x3[1], ]
  for(i in 2:nNRA){ 
    PopDataBootNRA = rbind(PopDataBootNRA,PopDataRAW[x3[i], ])
  } 
  
  ### start  loop for clones data
  nClones = length(ALLClones[,1])
  x4 = sample(1:nClones, nClones, replace = TRUE)  
  
  # bootstrap dataset
  ClonesBoot = ALLClones[x4[1], ]
  for(i in 2:nClones){ 
    ClonesBoot = rbind(ClonesBoot,ALLClones[x4[i], ])
  }    
  ALLClones=ClonesBoot
  
  ### start  loop for seedling change data
  nAASeedlings = length(AASeedlingsN[,1])
  x5 = sample(1:nAASeedlings, nAASeedlings, replace = TRUE)  
  
  # bootstrap each dataset (for q- direct and indirect)
  BootAASeedlingsN = AASeedlingsN[x5[1], ]
  for(i in 2:nAASeedlings){ 
    BootAASeedlingsN = rbind(BootAASeedlingsN,AASeedlingsN[x5[i], ])
  } 
  BootAASeedlingsY = AASeedlingsY[x5[1], ]
  for(i in 2:nAASeedlings){ 
    BootAASeedlingsY = rbind(BootAASeedlingsY,AASeedlingsY[x5[i], ])
  }  
  
  # Combine Clones, NRA, and RA (including those outside of plots) for new PopData
  PopData = rbind(PopDataBootNRA, PopDataSDLBoot)
  
  ### These subsets are used for all treatments
  Good.PopData = subset(PopData, Quality =='good') 
  Nonclones = subset(PopData, Stage !='gone' & Quality == 'good')
  Clones = subset(ALLClones, Quality == 'good')
  PopDataRA = PopData[PopData$StageNext != 'RA' & PopData$StageNext != 'Male',]
  
  ### start  loop for germination data
  nGerm = length(Germ[,1])
  x4 = sample(1:nGerm, nGerm, replace = TRUE)  
  
  # bootstrap dataset
  GermBoot = Germ[x4[1], ]
  for(i in 2:nGerm){ 
    GermBoot = rbind(GermBoot,Germ[x4[i], ])
  }    
  Germ=GermBoot
  
  
  p.vec=data.frame(
    surv.int=NA,           # Intercept from logistic regression of survival
    surv.slope=NA,         # Slope from logistic regression of survival
    growth.int=NA,         # Intercept from linear regression of growth
    growth.slope=NA,       # Slope from linear regression of growth
    growth.sd=NA,          # Residual sd from linear regression of growth
    seed.int=NA,           # Intercept from Poisson regression of seed number and size
    seed.slope=NA,         # Slope from Poisson regression of seed number and size
    prob.repro.int=NA,     # Intercept from logistic regression of reproduction
    prob.repro.slope=NA,   # Slope from logistic regression of reproduction
    recruit.size.mean=NA,  # Mean recruit size
    recruit.size.sd=NA,    # Standard deviation of recruit size
    establishment.prob=NA, # Probability of establishment
    clonal.size.mean=NA,   # Mean clone size
    clonal.size.sd=NA,     # Standard deviation of clone size
    clonal.prob=NA,        # Probability of producting a clone
    stay.seedbank.prob=NA, # Probability of staying in the seedbank
    go.seedbank.prob=NA,   # Probablity of entering the seedbank
    germ.seedbank.prob=NA  # Probablity of germinating from the seedbank
  )
  
  # 1. Survival: logistic regression
  surv.reg=glm(Survival~Size,data=PopData,family=binomial())
  
  p.vec$surv.int=coefficients(surv.reg)[1]
  p.vec$surv.slope=coefficients(surv.reg)[2]
  
  # 2. Growth: linear regression with x-intercept set to zero
  growth.reg=lm(SizeNext~0+Size,data=PopData)
  p.vec$growth.int=0  # If intercept is above one naturally, change to: coefficients(growth.reg)[1]
  p.vec$growth.slope=coefficients(growth.reg)[1]  # And change this to: coefficients(growth.reg)[2] 
  p.vec$growth.sd=sd(resid(growth.reg))
  
  # 3. Seeds: Divided into two parameters
  # Logistic regression for probability of being reproductive
  prob.repro.reg=glm(Reproductive~SizeNext,data=PopDataRA,family=binomial())
  
  ### The warnings in the output result from fitting a binomial to the survival and probability of reproduction data. 
  ### "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  ### This error has bypassed using an if/else statement; when the probabilities were numerically
  ### 0 or 1, I used the original values instead of the resulting inflated values for slope and intercept.
  ### If the original model is changed, these hard-coded values will also have to be changed.
  mu=exp(PopData$Reproductive)
  eps=10 * .Machine$double.eps
  
  if (any(mu > 1 - eps) || any(mu < eps)) {
    p.vec$prob.repro.int=-27.8
    p.vec$prob.repro.slope=0.247 
  } else {
    p.vec$prob.repro.int=coefficients(prob.repro.reg)[1]
    p.vec$prob.repro.slope=coefficients(prob.repro.reg)[2]
  }
  
  #  Poisson regression for relationship between SizeNext and number of seeds
  seed.reg=glm(Seeds~SizeNext,data=PopDataSeeds,family=poisson())
  p.vec$seed.int=coefficients(seed.reg)[1]
  p.vec$seed.slope=coefficients(seed.reg)[2]
  
  # Mean and standard deviation of seedling size for use in a normal distribution. 
  # Assume that offspring size is independent of maternal size,
  # so we only need to describe the distribution of offspring sizes in terms of its mean and variance.
  
  # 4. Size distribution of seedlings
  p.vec$recruit.size.mean=mean(PopDataSDL$SizeNext)
  p.vec$recruit.size.sd=sd(PopDataSDL$SizeNext)
  
  # 5. Staying in continuous stage (i.e., does not enter seedbank)
  k=1
  
  v=mean(Germ$v)
  
  g=mean(Germ$g)
  
  # Subset data
  SDL_NoBurn = AASeedlings[AASeedlings$Burn == 'N',] 
  
  # For unburned plots
  e = (mean(BootAASeedlingsN$Change) + 10) / (mean(BootAASeedlingsY$Change) + 10)
  
  p.vec$establishment.prob = v*g*e
  
  # 6.  Size distribution of clones
  p.vec$clonal.size.mean=mean(Clones$SizeNext, na.rm = TRUE)
  p.vec$clonal.size.sd=sd(Clones$SizeNext, na.rm = TRUE)
  
  # 7. Probability of clonal growth
  # Estimated by dividing the number of new recruits (i.e., Clones) in 2013 
  # by the number of NRA and RA (i.e., Nonclones) in 2012
  p.vec$clonal.prob=length(Clones$SizeNext)/length(Nonclones$Size)
  CompN_clonal.prob=length(Clones$SizeNext)/length(Nonclones$Size)
  
  # 8.  Probabilities associated with discrete stage (i.e., seedbank)
  # Probability of staying in the seedbank
  p.vec$stay.seedbank.prob=v*(1-g)*k
  
  # Probablity of entering the seedbank from continuous stage (more below)
  p.vec$go.seedbank.prob=v*(1-g)*k
  
  # Probablity of germinating from the seedbank (more below)
  p.vec$germ.seedbank.prob=v*g*e
  
  ### Define functions to describe life history ###
  
  # 1. Survival probability function
  s.x.CompN=function(x,params) {
    u=exp(params$surv.int+params$surv.slope*x)
    return(u/(1+u))
  }
  
  # 2. Growth function
  g.yx.CompN=function(xp,x,params) {
    dnorm(xp,mean=params$growth.int+params$growth.slope*x,sd=params$growth.sd) 
  }
  
  # 3. Reproduction function
  p.repro.x.CompN=function(x,params) {
    u=exp(params$prob.repro.int+params$prob.repro.slope*x)
    return(u/(1+u))
  }
  
  # Fecundity model while taking into account some seeds move into 
  # the discrete portion of the model (i.e., seedbank)
  f.yx.CompN=function(xp,x,params) {
    p.repro.x.CompN(x,params)*
      params$establishment.prob*
      dnorm(xp,mean=params$recruit.size.mean,sd=params$recruit.size.sd)*
      exp(params$seed.int+params$seed.slope*x)
  }
  
  # 4. Seedling size distribution
  d.x.CompN=function(xp,x,params) {
    dnorm(xp, mean=p.vec$recruit.size.mean,sd=p.vec$recruit.size.sd)*p.vec$germ.seedbank.prob
  }
  
  # 5. Seeds entering seedbank
  e.x.CompN=function(xp,x,params) {
    p.repro.x.CompN(x,params)*exp(params$seed.int+params$seed.slope*x)*params$go.seedbank.prob
  }
  
  # 6. Clonal function
  c.yx.CompN=function(xp,x,params) {
    params$clonal.prob*
      dnorm(xp,mean=params$clonal.size.mean,sd=params$clonal.size.sd)
  }
  
  ### Combine vital rate functions to build the discretized IPM kernal (i.e., IPM matrix)
  # Set number of classes (n.size), or points for midpoint rule approximation
  # Define the boundary points (b; the edges of the cells defining the matrix),
  # mesh points (y; the centers of the cells defining the matrix and the points at which
  # the matrix is evaluated for the midpoint rule of numerical integration, and
  # step size (h; the widths of the cells)
  # The integration limits (min.size and max.size_CompN) span the range of sizes observed in the data set.
  
  min.size=0.9*min(c(PopData$Size,PopData$SizeNext),na.rm=T)  # Use values slightly above and below limits
  max.size_CompN=1.1*max(c(PopData$Size,PopData$SizeNext),na.rm=T)
  n.size=100 # number of cells in the matrix 
  b=min.size+c(0:n.size)*(max.size_CompN-min.size)/n.size # boundary points
  y=0.5*(b[1:n.size]+b[2:(n.size+1)]) # mesh points
  h=y[2]-y[1] # step size
  
  ### Make IPM matrices ###
  # The function outer() evaluates the matrix at all pairwise combinations of the two 
  # vectors y and y and returns matrices representing the kernel components for growth 
  # and fecundity, respectively. For the numerical integration, we’re using the midpoint 
  # rule estimate the area under a curve. The midpoint rule assumes a rectangular 
  # approximation. The heights of the rectangles are given by the outer function and 
  # the width of the rectangles is h.
  # The result is n.size x n.size cell discretization of the kernel, K.
  
  # Growth matrix
  G_CompN=h*outer(y,y,g.yx.CompN,params=p.vec)
  # Larger individuals are evicted (see Williams et al. 2012), so return the evicted individuals to the 
  # cells at the boundaries where they were evicted (i.e., rerout growth to sizes outside the allowed range 
  # to the extreme sizes avoiding eviction).
  for(i in 1:(n.size/2)) G_CompN[1,i]=G_CompN[1,i]+1-sum(G_CompN[,i])
  for(i in (n.size/2+1):n.size) G_CompN[n.size,i]=G_CompN[n.size,i]+1-sum(G_CompN[,i])
  
  # Survival vector
  S_CompN=s.x.CompN(y,params=p.vec) 
  SMatrix = array(0,dim=c(100,99))
  S_CompN1= cbind(S_CompN, SMatrix)          
  
  # Fecundity martix
  F_CompN = array(0,dim=c(1+n.size,1+n.size))
  F_CompN[2:(1+n.size),2:(1+n.size)] = h*outer(y,y,f.yx.CompN,params=p.vec)
  # e.x.CompN tells seeds to go into the seedbank
  F_CompN[1,2:(1+n.size)] = e.x.CompN(y,y,params=p.vec)
  
  # Survival/growth matrix
  # The first row and column are seedbank transitions
  P_CompN = array(0,dim=c(1+n.size,1+n.size))
  
  # Add probability of staying in seedbank 
  P_CompN[1,1] = p.vec$stay.seedbank.prob 
  P_CompN1 = P_CompN
  
  # Add probability of germinating and moving to continuous
  P_CompN[2:(1+n.size),1] = d.x.CompN(y,params=p.vec)
  
  # Build growth/survival matrix including discrete seedbank stage
  for(i in 2:(1+n.size))
    P_CompN[i,2:(1+n.size)] = S_CompN*G_CompN[(i-1),]
  
  C_CompN = array(0,dim=c(1+n.size,1+n.size))
  C_CompN[2:(1+n.size),2:(1+n.size)] = h*outer(y,y,c.yx.CompN,params=p.vec) # clonal growth matrix
  
  # Build complete matrix
  K_CompN=P_CompN+C_CompN+F_CompN
  
  ### Calculate lambda
  CompN_lambda_boot =Re(eigen(K_CompN)$values[1])
  
  CompN_bootlambda[j] = CompN_lambda_boot
  
  p.vec.all <- lapply(seq_along(p.vec.all), function(x) p.vec.all[[x]] <- c(p.vec.all[[x]], p.vec[1, x]))
  
}
# name everything and now do a little data manipulation
names(p.vec.all) <- names(p.vec)
p.vec.all$lambda <- c(CompN_lambda, CompN_bootlambda)


saveRDS(p.vec.all, file = 'Ailanthus_IPM/CR_Bootstrap_Output.rds')
