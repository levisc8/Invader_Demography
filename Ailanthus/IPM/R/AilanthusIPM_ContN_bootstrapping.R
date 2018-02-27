### Infile growth, survival, and reproduction data ###
Ailanthus4R = read.csv("Ailanthus/IPM/Data/Ailanthus4R.csv", header = TRUE)
head(Ailanthus4R)

### Infile germination data ###
Ailanthus4R_Germ = read.csv("Ailanthus/IPM/Data/Ailanthus4R_Germ.csv")

### Infile seedling count data ###
AASeedlings = read.csv("Ailanthus/IPM/Data/AA SDL Changes by plot.csv", header = TRUE)
AASeedlingsN = AASeedlings[AASeedlings$Burn == 'N',]
AASeedlingsY = AASeedlings[AASeedlings$Burn == 'Y',]

### Subset data ###
# PopData is file for growth and survival analyses; uncomment the treatment of interest
# Good.PopData is used to model size distribution of only the good clones
# PopDataSDL is used to model size distribution and establishment probability
# PopDataRA is used to model the relationship between fecundity and size of RAs only
# UNCOMMENT THE TREATMENT OF INTEREST BEFORE RUNNING CODE.
# POPDATA AND POPDATASDL ARE UNIQUE FOR EACH TREATMENT.
# WHEN SELECTING A TREATMENT FOR ANALYSIS, ALSO SEE LINES 170-177.
###########################################################################

### Start bootstrapping
ContN_bootstrap_clones = rep(NA,1000)
ContN_bootlambda = rep(NA, 1000) # This will save all the lambdas from each bootstrap. For now I fill in a column with 1000 rows with the value 'NA'.  Later they are replacee with data.
for(j in 1:1000) 
{
  
  ###########################################################################
  # Competitor removal without fire
  PopData1 = Ailanthus4R[Ailanthus4R$Trt == 'Control' | Ailanthus4R$Trt == 'All',]
  PopDataRAW = PopData1[PopData1$Burn != 'Y' & PopData1$Stage2013 != 'SDL',]
  
  Clones1 = subset(PopData1, Burn != 'Y' & Stage2013 != 'gone')
  Clones2 = subset(Clones1, Stage2013 != 'SDL' & Stage2012 != 'SDL')
  ALLClones = subset(Clones2, Stage2012 == 'gone')
  
  PopDataSDL = subset(Ailanthus4R, (Trt != 'Herb' & Trt != 'Comp' & Burn != 'Y' & Stage2013 == 'SDL'))
  
  ###########################################################################
  
  # Use these subsets for all treatments
  PopDataSeeds = Ailanthus4R[Ailanthus4R$Stage2013 == 'RA' | Ailanthus4R$Stage2013 == 'girl',]
  
  Germ = subset(Ailanthus4R_Germ, Study == 'Germ')
  
  FireY = subset(Ailanthus4R_Germ, Study == 'Fire' & Location == '0' & Burn == 'Yes')
  FireN = subset(Ailanthus4R_Germ, Study == 'Fire' & Location == '0' & Burn == 'No')
  
  ###########################################################################
  
  ### start the loop for seedling data
  nseedlings = length(PopDataSDL[,1]) 
  x1 = sample(1:nseedlings, nseedlings, replace = TRUE)
  
  # bootstrap dataset
  PopDataSDLBoot = PopDataSDL[x1[1],c('Plot','Site','Trt', 'Burn', 'Trt2', 'Stage2012', 'Number2012', 'MaxD2012', 'Size', 'Number2013', 'Stage2013', 'MaxD2013', 'Quality2013', 'SizeNext', 'Survival', 'Seeds_All', 'Seeds', 'Fecundity')]
  for(i in 2:nseedlings){ 
    PopDataSDLBoot = rbind(PopDataSDLBoot,PopDataSDL[x1[i],c('Plot','Site','Trt', 'Burn', 'Trt2', 'Stage2012', 'Number2012', 'MaxD2012', 'Size', 'Number2013', 'Stage2013', 'MaxD2013', 'Quality2013', 'SizeNext', 'Survival', 'Seeds_All', 'Seeds', 'Fecundity')])
  }
  PopDataSDL=PopDataSDLBoot
  
  ### start the loop for seeds data
  nSeeds = length(PopDataSeeds[,1])
  x2 = sample(1:nSeeds, nSeeds, replace = TRUE)
  
  # bootstrap dataset
  PopDataSeedsBoot = PopDataSeeds[x2[1],c('Plot','Site','Trt', 'Burn', 'Trt2', 'Stage2012', 'Number2012', 'MaxD2012', 'Size', 'Number2013', 'Stage2013', 'MaxD2013', 'Quality2013', 'SizeNext', 'Survival', 'Seeds_All', 'Seeds', 'Fecundity')]
  for(i in 2:nSeeds){ 
    PopDataSeedsBoot = rbind(PopDataSeedsBoot,PopDataSeeds[x2[i],c('Plot','Site','Trt', 'Burn', 'Trt2', 'Stage2012', 'Number2012', 'MaxD2012', 'Size', 'Number2013', 'Stage2013', 'MaxD2013', 'Quality2013', 'SizeNext', 'Survival', 'Seeds_All', 'Seeds', 'Fecundity')])
  }
  PopDataSeeds=PopDataSeedsBoot
  
  ### start the loop for NRA and RA data
  nNRA = length(PopDataRAW[,1])
  x3 = sample(1:nNRA, nNRA, replace = TRUE)
  
  # bootstrap dataset
  PopDataBootNRA = PopDataRAW[x3[1],c('Plot','Site','Trt', 'Burn', 'Trt2', 'Stage2012', 'Number2012', 'MaxD2012', 'Size', 'Number2013', 'Stage2013', 'MaxD2013', 'Quality2013', 'SizeNext', 'Survival', 'Seeds_All', 'Seeds', 'Fecundity')]
  for(i in 2:nNRA){ 
    PopDataBootNRA = rbind(PopDataBootNRA,PopDataRAW[x3[i],c('Plot','Site','Trt', 'Burn', 'Trt2', 'Stage2012', 'Number2012', 'MaxD2012', 'Size', 'Number2013', 'Stage2013', 'MaxD2013', 'Quality2013', 'SizeNext', 'Survival', 'Seeds_All', 'Seeds', 'Fecundity')])
  } 
  
  ### start  loop for clones data
  nClones = length(ALLClones[,1])
  x4 = sample(1:nClones, nClones, replace = TRUE)  
  
  # bootstrap dataset
  ClonesBoot = ALLClones[x4[1],c('Plot','Site','Trt', 'Burn', 'Trt2', 'Stage2012', 'Number2012', 'MaxD2012', 'Size', 'Number2013', 'Stage2013', 'MaxD2013', 'Quality2013', 'SizeNext', 'Survival', 'Seeds_All', 'Seeds', 'Fecundity')]
  for(i in 2:nClones){ 
    ClonesBoot = rbind(ClonesBoot,ALLClones[x4[i],c('Plot','Site','Trt', 'Burn', 'Trt2', 'Stage2012', 'Number2012', 'MaxD2012', 'Size', 'Number2013', 'Stage2013', 'MaxD2013', 'Quality2013', 'SizeNext', 'Survival', 'Seeds_All', 'Seeds', 'Fecundity')])
  }    
  ALLClones=ClonesBoot
  
  ### start  loop for seedling change data
  nAASeedlings = length(AASeedlingsN[,1])
  x5 = sample(1:nAASeedlings, nAASeedlings, replace = TRUE)  
  
  # bootstrap each dataset (for q- direct and indirect)
  BootAASeedlingsN = AASeedlingsN[x5[1],c('Burn',  'Trt',  'Trt_Burn',  'Plot_Value',	'SDL_no_2012', 'SDL_no_2013', 'Change')]
  for(i in 2:nAASeedlings){ 
    BootAASeedlingsN = rbind(BootAASeedlingsN,AASeedlingsN[x5[i],c('Burn',  'Trt',  'Trt_Burn',	'Plot_Value',	'SDL_no_2012', 'SDL_no_2013', 'Change')])
  } 
  BootAASeedlingsY = AASeedlingsY[x5[1],c('Burn',  'Trt',  'Trt_Burn',	'Plot_Value',	'SDL_no_2012', 'SDL_no_2013', 'Change')]
  for(i in 2:nAASeedlings){ 
    BootAASeedlingsY = rbind(BootAASeedlingsY,AASeedlingsY[x5[i],c('Burn',  'Trt',  'Trt_Burn',	'Plot_Value',	'SDL_no_2012', 'SDL_no_2013', 'Change')])
  } 
  
  # Combine Clones, NRA, and RA (including those outside of plots) for new PopData
  PopData = rbind(PopDataBootNRA, PopDataSDLBoot)
  
  ### These subsets are used for all treatments
  Good.PopData = subset(PopData, Quality2013 =='good') 
  Nonclones = subset(PopData, Stage2012 !='gone' & Quality2013 == 'good')
  Clones = subset(ALLClones, Quality2013 == 'good')
  PopDataRA = PopData[PopData$Stage2013 != 'RA' & PopData$Stage2013 != 'Boy',]
  
  ### start  loop for germination data
  nGerm = length(Germ[,1])
  x4 = sample(1:nGerm, nGerm, replace = TRUE)  
  
  # bootstrap dataset
  GermBoot = Germ[x4[1],c('Study',  'ID',  'Burn',	'Protected',	'Location',	'Total',	'Germinated',	'Proportion',	'Viable',	'v',	'g')]
  for(i in 2:nGerm){ 
    GermBoot = rbind(GermBoot,Germ[x4[i],c('Study',  'ID',  'Burn',	'Protected',	'Location',	'Total',	'Germinated',	'Proportion',	'Viable',	'v',	'g')])
  }    
  Germ=GermBoot
  
  ### start  loop for burn germination data
  nFireY = length(FireY[,1])
  x5 = sample(1:nFireY, nFireY, replace = TRUE)  
  
  # bootstrap dataset
  FireYBoot = FireY[x5[1],c('Study',  'ID',  'Burn',	'Protected',	'Location',	'Total',	'Germinated',	'Proportion',	'Viable',	'v',	'g')]
  for(i in 2:nFireY){ 
    FireYBoot = rbind(FireYBoot,FireY[x5[i],c('Study',  'ID',  'Burn',	'Protected',	'Location',	'Total',	'Germinated',	'Proportion',	'Viable',	'v',	'g')])
  }    
  FireY=FireYBoot
  
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
  surv.reg=glm(Survival~Size,data=PopData,family=quasibinomial())
  
  mu=exp(PopData$Survival)
  eps=10 * .Machine$double.eps
  
  if (any(mu > 1 - eps) || any(mu < eps)) {
    p.vec$surv.int=0.944
    p.vec$surv.slope=0.0629 
  } else {
    p.vec$surv.int=coefficients(prob.repro.reg)[1]
    p.vec$surv.slope=coefficients(prob.repro.reg)[2]
  }
  
  # 2. Growth: linear regression with x-intercept set to zero
  growth.reg=lm(SizeNext~0+Size,data=PopData)
  p.vec$growth.int=0  # If intercept is above one naturally, change to: coefficients(growth.reg)[1]
  p.vec$growth.slope=coefficients(growth.reg)[1]  # And change this to: coefficients(growth.reg)[2] 
  p.vec$growth.sd=sd(resid(growth.reg))
  
  # 3. Seeds: Divided into two parameters
  # Logistic regression for probability of being reproductive
  prob.repro.reg=glm(Fecundity~SizeNext,data=PopDataRA,family=binomial())
  
  ### The warnings in the output result from fitting a binomial to the survival and probability of reproduction data. 
  ### "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  ### This error has bypassed using an if/else statement; when the probabilities were numerically
  ### 0 or 1, I used the original values instead of the resulting inflated values for slope and intercept.
  ### If the original model is changed, these hard-coded values will also have to be changed.
  mu=exp(PopData$Fecundity)
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
  
  ### Infile seedling count data ###
  AASeedlings = read.csv("Ailanthus/IPM/Data/Copy of AA SDL Counts by plot.csv", header = TRUE)
  
  # Subset data
  SDL_NoBurn = AASeedlings[AASeedlings$Burn == 'N',] 
  
  # Plot data
  SDL_no.reg=lm(SDL_no~Year,data=SDL_NoBurn)
  
  # For unburned plots
  e = (mean(BootAASeedlingsN$Change)+10)/(mean(BootAASeedlingsY$Change)+10)
  
  p.vec$establishment.prob=v*g*e
  
  # 6.  Size distribution of clones
  p.vec$clonal.size.mean=mean(Clones$SizeNext, na.rm = TRUE)
  p.vec$clonal.size.sd=sd(Clones$SizeNext, na.rm = TRUE)
  
  # 7. Probability of clonal growth
  # Estimated by dividing the number of new recruits (i.e., Clones) in 2013 
  # by the number of NRA and RA (i.e., Nonclones) in 2012
  p.vec$clonal.prob=length(Clones$SizeNext)/length(Nonclones$Size)
  ContN_clonal.prob=length(Clones$SizeNext)/length(Nonclones$Size)
  
  # 8.  Probabilities associated with discrete stage (i.e., seedbank)
  # Probability of staying in the seedbank
  p.vec$stay.seedbank.prob=v*(1-g)*k
  
  # Probablity of entering the seedbank from continuous stage (more below)
  p.vec$go.seedbank.prob=v*(1-g)*k
  
  # Probablity of germinating from the seedbank (more below)
  p.vec$germ.seedbank.prob=v*g*e
  
  ### Define functions to describe life history ###
  
  # 1. Survival probability function
  s.x.ContN=function(x,params) {
    u=exp(params$surv.int+params$surv.slope*x)
    return(u/(1+u))
  }
  
  # 2. Growth function
  g.yx.ContN=function(xp,x,params) {
    dnorm(xp,mean=params$growth.int+params$growth.slope*x,sd=params$growth.sd) 
  }
  
  # 3. Reproduction function
  p.repro.x.ContN=function(x,params) {
    u=exp(params$prob.repro.int+params$prob.repro.slope*x)
    return(u/(1+u))
  }
  
  # Fecundity model while taking into account some seeds move into 
  # the discrete portion of the model (i.e., seedbank)
  f.yx.ContN=function(xp,x,params) {
    p.repro.x.ContN(x,params)*
      params$establishment.prob*
      dnorm(xp,mean=params$recruit.size.mean,sd=params$recruit.size.sd)*
      exp(params$seed.int+params$seed.slope*x)
  }
  
  # 4. Seedling size distribution
  d.x.ContN=function(xp,x,params) {
    dnorm(xp, mean=p.vec$recruit.size.mean,sd=p.vec$recruit.size.sd)*p.vec$germ.seedbank.prob
  }
  
  # 5. Seeds entering seedbank
  e.x.ContN=function(xp,x,params) {
    p.repro.x.ContN(x,params)*exp(params$seed.int+params$seed.slope*x)*params$go.seedbank.prob
  }
  
  # 6. Clonal function
  c.yx.ContN=function(xp,x,params) {
    params$clonal.prob*
      dnorm(xp,mean=params$clonal.size.mean,sd=params$clonal.size.sd)
  }
  
  ### Combine vital rate functions to build the discretized IPM kernal (i.e., IPM matrix)
  # Set number of classes (n.size), or points for midpoint rule approximation
  # Define the boundary points (b; the edges of the cells defining the matrix),
  # mesh points (y; the centers of the cells defining the matrix and the points at which
  # the matrix is evaluated for the midpoint rule of numerical integration, and
  # step size (h; the widths of the cells)
  # The integration limits (min.size and max.size_ContN) span the range of sizes observed in the data set.
  
  min.size=0.9*min(c(PopData$Size,PopData$SizeNext),na.rm=T)  # Use values slightly above and below limits
  max.size_ContN=1.1*max(c(PopData$Size,PopData$SizeNext),na.rm=T)
  n.size=100 # number of cells in the matrix 
  b=min.size+c(0:n.size)*(max.size_ContN-min.size)/n.size # boundary points
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
  G_ContN=h*outer(y,y,g.yx.ContN,params=p.vec)
  # Larger individuals are evicted (see Williams et al. 2012), so return the evicted individuals to the 
  # cells at the boundaries where they were evicted (i.e., rerout growth to sizes outside the allowed range 
  # to the extreme sizes avoiding eviction).
  for(i in 1:(n.size/2)) G_ContN[1,i]=G_ContN[1,i]+1-sum(G_ContN[,i])
  for(i in (n.size/2+1):n.size) G_ContN[n.size,i]=G_ContN[n.size,i]+1-sum(G_ContN[,i])
  
  # Survival vector
  S_ContN=s.x.ContN(y,params=p.vec) 
  SMatrix = array(0,dim=c(100,99))
  S_ContN1= cbind(S_ContN, SMatrix)          
  
  # Fecundity martix
  F_ContN = array(0,dim=c(1+n.size,1+n.size))
  F_ContN[2:(1+n.size),2:(1+n.size)] = h*outer(y,y,f.yx.ContN,params=p.vec)
  # e.x.ContN tells seeds to go into the seedbank
  F_ContN[1,2:(1+n.size)] = e.x.ContN(y,y,params=p.vec)
  
  # Survival/growth matrix
  # The first row and column are seedbank transitions
  P_ContN = array(0,dim=c(1+n.size,1+n.size))
  
  # Add probability of staying in seedbank 
  P_ContN[1,1] = p.vec$stay.seedbank.prob 
  P_ContN1 = P_ContN
  
  # Add probability of germinating and moving to continuous
  P_ContN[2:(1+n.size),1] = d.x.ContN(y,params=p.vec)
  
  # Build growth/survival matrix including discrete seedbank stage
  for(i in 2:(1+n.size))
    P_ContN[i,2:(1+n.size)] = S_ContN*G_ContN[(i-1),]
  
  C_ContN = array(0,dim=c(1+n.size,1+n.size))
  C_ContN[2:(1+n.size),2:(1+n.size)] = h*outer(y,y,c.yx.ContN,params=p.vec) # clonal growth matrix
  
  # Build complete matrix
  K_ContN=P_ContN+C_ContN+F_ContN
  
  ### Calculate lambda
  ContN_lambda_boot =Re(eigen(K_ContN)$values[1])
  
  ContN_bootstrap_clones[j] = p.vec$clonal.prob
  ContN_bootlambda[j] = ContN_lambda_boot
}

ContN_bootstrap_clones = sort(ContN_bootstrap_clones)
ContN_clones_CI = c(ContN_bootstrap_clones[25], ContN_bootstrap_clones[975])

ContN_bootlambda = sort(ContN_bootlambda)

ContN_lambda_CI = c(ContN_bootlambda[25], ContN_bootlambda[975])

### End bootstrapping
###########################################################################
###########################################################################
###########################################################################