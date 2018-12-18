### Infile growth, survival, and reproduction data ###
Ailanthus4R <- read.csv("Ailanthus4R.csv", header = TRUE)
head(Ailanthus4R)

### Infile seedling count data ###
AASeedlings <- read.csv("AA SDL Counts by plot.csv", header = TRUE)

### Infile germination data ###
Ailanthus4R_Germ <- read.csv("Ailanthus4R_Germ.csv")
head(Ailanthus4R_Germ)

### Subset data ###
# PopData is used for growth and survival analyses
# Good_PopData is used to model size distribution of measured clones
# PopDataSDL is used to model size distribution of seedlings
# PopDataRA is used to model the relationship between size and probability of reproduction
# PopDataSeeds is used to model relationship between fecundity and size; it only includes reproducing individuals

###########################################################################
# Control without fire
PopData1 <- Ailanthus4R[Ailanthus4R$Trt == 'Control' | Ailanthus4R$Trt == 'All',]
PopData <- PopData1[PopData1$Burn != 'Y',]

PopDataSDL <- Ailanthus4R[Ailanthus4R$Trt == 'Control' & Ailanthus4R$Burn != 'Y' & Ailanthus4R$Stage2013 == 'SDL',]

##########################################################################
# These subsets are used for all treatments
Good_PopData <- PopData[PopData$Quality2013 =='good',]
Nonclones <- Good_PopData[Good_PopData$Stage2012 !='gone' & Good_PopData$Stage2013 != 'SDL' & Good_PopData$Stage2012 != 'SDL',]
Clones1 <- PopData[PopData$Stage2013 != 'gone',]
Clones2 <- Clones1[Clones1$Stage2013 != 'SDL' & Clones$Stage2012 != 'SDL',]
Clones <- Clones2[Clones2$Stage2012 == 'gone' & Clones2$Quality2013 == 'good',]

PopDataRA1 <- PopData[PopData$Stage2013 != 'RA',] 
PopDataRA <- PopDataRA1[PopDataRA1$Stage2013 != 'Boy',]
PopDataSeeds <- Ailanthus4R[Ailanthus4R$Stage2013 == 'RA' | Ailanthus4R$Stage2013 == 'girl',]

Germ <- Ailanthus4R_Germ[Ailanthus4R_Germ$Study == 'Germ',]
Viable <- Ailanthus4R_Germ[Ailanthus4R_Germ$Study == 'Viable',]
FireY <- Ailanthus4R_Germ[Ailanthus4R_Germ$Study == 'Fire' & Ailanthus4R_Germ$Location == '0' & Ailanthus4R_Germ$Burn == 'Yes',]
FireN <- Ailanthus4R_Germ[Ailanthus4R_Germ$Study == 'Fire' & Ailanthus4R_Germ$Location == '0' & Ailanthus4R_Germ$Burn == 'No',]

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
surv.reg <- glm(Survival ~ Size, data = PopData, family = quasibinomial())
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
prob.repro.reg <- glm(Fecundity ~ SizeNext, data = PopDataRA, family = binomial())
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


# Subset data
SDL_NoBurn <- AASeedlings[AASeedlings$Burn == 'N',] 

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
ContN_clonal.prob <- length(Clones$SizeNext)/length(Nonclones$Size)

# 8.  Probabilities associated with discrete stage (i.e., seedbank)
# Probability of staying in the seedbank
p.vec$stay.seedbank.prob <- v * (1 - g) * k

# Probablity of entering the seedbank from continuous stage (more below)
p.vec$go.seedbank.prob <- v * (1 - g) * k

# Probablity of germinating from the seedbank (more below)
p.vec$germ.seedbank.prob <- v * g * q

### Define functions to describe life history ###

# 1. Survival probability function
s.x.CONTN <- function(x,params) {
  u <- exp(params$surv.int+params$surv.slope * x)
  return(u / (1 + u))
}

# 2. Growth function
g.yx.CONTN <- function(xp,x,params) {
  dnorm(xp, mean = params$growth.int + params$growth.slope * x, sd = params$growth.sd) 
}

# 3. Reproduction function
p.repro.x.CONTN <- function(x,params) {
  u <- exp(params$prob.repro.int + params$prob.repro.slope * x)
  return(u / (1 + u))
}

# Fecundity model while taking into account some seeds move into 
# the discrete portion of the model (i.e., seedbank)
f.yx.CONTN <- function(xp,x,params) {
  p.repro.x.CONTN(x,params)*
    params$establishment.prob *
    dnorm(xp, mean = params$recruit.size.mean, sd = params$recruit.size.sd) *
    exp(params$seed.int + params$seed.slope * x)
}

# 4. Seedling size distribution
d.x.CONTN <- function(xp,x,params) {
  dnorm(xp, mean = p.vec$recruit.size.mean, sd = p.vec$recruit.size.sd) * p.vec$germ.seedbank.prob
}

# 5. Seeds entering seedbank
e.x.CONTN <- function(xp,x,params) {
  p.repro.x.CONTN(x,params) * exp(params$seed.int + params$seed.slope * x) * params$go.seedbank.prob
}

# 6. Clonal function
c.yx.CONTN <- function(xp,x,params) {
    params$clonal.prob *
    dnorm(xp, mean = params$clonal.size.mean, sd = params$clonal.size.sd)
}

### Combine vital rate functions to build the IPM matrix
# Set number of classes (n.size) or points for midpoint rule approximation
min.size <- 0.9 * min(c(PopData$Size, PopData$SizeNext), na.rm = T)
max.size_ContN <- 1.1 * max(c(PopData$Size, PopData$SizeNext), na.rm = T)
n.size <- 100 # number of cells in the matrix 
b <- min.size + c(0:n.size) * (max.size_ContN - min.size) / n.size # boundary points
y <- 0.5 * (b[1:n.size] + b[2:(n.size + 1)]) # mesh points
h <- y[2] - y[1] # step size

### Make IPM matrices ###

# Growth matrix
G_ContN <- h * outer(y, y, g.yx.CONTN, params = p.vec)
# Larger individuals are evicted (see Williams et al. 2012), so return the evicted individuals to the 
# cells at the boundaries where they were evicted (i.e., rerout growth to sizes outside the allowed range 
# to the extreme sizes avoiding eviction).
for(i in 1:(n.size/2)) G_ContN[1,i] <- G_ContN[1,i] + 1 - sum(G_ContN[,i])
for(i in (n.size/2+1):n.size) G_ContN[n.size,i] <- G_ContN[n.size,i] + 1 - sum(G_ContN[,i])

# Survival vector
S_ContN <- s.x.CONTN(y, params = p.vec) 
SMatrix <- array(0, dim = c(100,99))
S_ContN1 <- cbind(S_ContN, SMatrix)          

# Fecundity martix
F_ContN <- array(0, dim=c(1 + n.size, 1 + n.size))
F_ContN[2:(1 + n.size), 2:(1 + n.size)] <- h * outer(y, y, f.yx.CONTN, params = p.vec)
# e.x.CONTN tells seeds to go into the seedbank
F_ContN[1, 2:(1 + n.size)] <- e.x.CONTN(y, y, params = p.vec)

# Survival/growth matrix
# The first row and column are seedbank transitions
P_ContN <- array(0, dim = c(1 + n.size, 1 + n.size))

# Add probability of staying in seedbank 
P_ContN[1,1] <- p.vec$stay.seedbank.prob 
P_ContN1 <- P_ContN

# Add probability of germinating and moving to continuous
P_ContN[2:(1 + n.size), 1] <- d.x.CONTN(y, params = p.vec)

# Build growth/survival matrix including discrete seedbank stage
for(i in 2:(1 + n.size))
  P_ContN[i, 2:(1 + n.size)] <- S_ContN * G_ContN[(i - 1),]

C_ContN <- array(0,dim=c(1 + n.size, 1 + n.size))
C_ContN[2:(1 + n.size), 2:(1 + n.size)] <- h * outer(y, y, c.yx.CONTN, params = p.vec) # clonal growth matrix

# Build matrix
K_ContN <- P_ContN + C_ContN + F_ContN

# Calculate lambda
ContN_lambda <- Re(eigen(K_ContN)$values[1])
ContN_lambda