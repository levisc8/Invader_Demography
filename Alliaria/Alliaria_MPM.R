
##Clear everything
rm(list=ls(all=TRUE)) 
cat("\014") 

#import data into R studio
library(dplyr)
ks <- read.csv("C:/Users/sl13sise/Dropbox/invaders demography/Alliaria/GM4R.csv",
               stringsAsFactors = F) %>%
  filter(Treatment != 'Herb')

# Calculate number of plots per treatment
Plot.N <- ks %>%
  group_by(Treatment,Density,Burn) %>% 
  summarise(plot.n = length(unique(Plot)))

# calculate seedlings per plot to test for density dependence
sdlPerPlot <- ks %>% 
  group_by(Treatment, Density, Plot) %>%
  summarise(N.SDL = n(),
            N.RA = sum(RA, na.rm = T),
            s2 = N.RA/N.SDL,
            Fec = mean(Seeds,na.rm = T))

# visually examine relationship
plot(s2 ~ N.SDL, data = sdlPerPlot)
plot(Fec ~ N.RA, data = sdlPerPlot)

# cutting off plots above 200 seedlings/plot
idx <- sdlPerPlot$Plot[sdlPerPlot$N.SDL < 200]

paramsD <- ks[ks$Plot %in% idx, ] %>%
  group_by(Treatment) %>% 
  summarise(plot.n = length(unique(Plot)),
            N.SDL = n(),
            N.RA = sum(RA, na.rm = T),
            s2 = N.RA/N.SDL) 

Fec <- ks %>% 
  group_by(Treatment) %>% 
  summarise(Fec = mean(Seeds, na.rm = T))



paramsD[paramsD$Treatment == 'Comp', 'Fec'] <- Fec[Fec$Treatment == 'Comp',
                                                   'Fec']
paramsD[paramsD$Treatment == 'Control', 'Fec'] <- Fec[Fec$Treatment == 'Control',
                                                      'Fec']


##Control
s1_c <- 1
s2_c <- paramsD[paramsD$Treatment == "Control", "s2"] %>% as.numeric
f_c <- paramsD[paramsD$Treatment == "Control", "Fec"] %>% as.numeric

##Comp Rem
s1_cr <- 1
s2_cr <- paramsD[paramsD$Treatment=="Comp","s2"] %>% as.numeric
f_cr <- paramsD[paramsD$Treatment=="Comp","Fec"] %>% as.numeric


##Germination and seed viability taken from Pardini et al 2009
G1 <- 0.5503
G2 <- 0.3171
v <- 0.8228

A_cont <- matrix(c(1 - G2, 0, f_c * v * (1 - G1),
                   G2 * s1_c, 0, f_c * G1 * s1_c,
                   0, s2_c, 0),
                 nrow = 3, byrow = TRUE, 
                 dimnames = list(c("SB", "Rosette", "RA"),
                                 c("SB", "Rosette", "RA")))

A_cr <- matrix(c(1 - G2, 0, f_cr * v * (1 - G1),
                 G2 * s1_cr, 0, f_cr * G1 * s1_cr,
                 0, s2_cr, 0),
               nrow = 3, byrow = TRUE, 
               dimnames = list(c("SB", "Rosette", "RA"),
                               c("SB", "Rosette", "RA")))
 

ev <- eigen(A_cont)

lmax <- which(Re(ev$values) == max(Re(ev$values)))
lmax <- lmax[1]  
lambda_cont <- Re(ev$values[lmax])
lambda_cont

ev <- eigen(A_cr)

lmax <- which(Re(ev$values) == max(Re(ev$values)))
lmax <- lmax[1]  
lambda_cr <- Re(ev$values[lmax])
lambda_cr

values=c(s2_c,s2_cr,
         f_c, f_cr,
         lambda_cont,lambda_cr)


####################################################
#Bootstrapping code

nreps <- 1000

boot_s2_c <- rep(NA, nreps)
boot_s2_cr <- rep(NA, nreps)
boot_f_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)
boot_l_c <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)

C.RAs <- filter(ks, Treatment == 'Control' & RA == 1)
CR.RAs <- filter(ks, Treatment == 'Comp' & RA == 1)

CPlants <- filter(ks, Treatment == 'Control' & Plot %in% idx)
CRPlants <- filter(ks, Treatment == 'Comp' & Plot %in% idx)

n.c.ra <- dim(C.RAs)[1]
n.cr.ra <- dim(CR.RAs)[1]
n.c.plants <- dim(CPlants)[1]
n.cr.plants <- dim(CRPlants)[1]


#start the loop
for(j in 1:nreps) {
  
  #x3 is a column that picks n integers between the numbers 1 and n (because we have n plants in the data set).  
  #replace=TRUE means that we are sampling with replacement (once we choose a number, we can choose it again).
  xCP <- sample(1:n.c.plants, n.c.plants,
                replace = T)
  xCRP <- sample(1:n.cr.plants, n.cr.plants,
                 replace = T)
  xCRA <- sample(1:n.c.ra, n.c.ra,
                 replace = TRUE)
  xCRRA <- sample(1:n.cr.ra, n.cr.ra,
                  replace = TRUE)
  
  bootSurvData <- rbind(CPlants[xCP, ],
                        CRPlants[xCRP, ])
  bootFecData <- rbind(C.RAs[xCRA, ],
                       CR.RAs[xCRRA, ])
  paramsD1 <- bootSurvData %>% group_by(Treatment) %>% summarise(
    N.SDL = n(),
    N.RA = sum(RA, na.rm = T),
    s2 = N.RA/N.SDL
  )
  
  Fec1 <- bootFecData %>% group_by(Treatment) %>% summarise(
    Fec = mean(Seeds, na.rm = T)
  )
  
  paramsD1[paramsD1$Treatment == 'Comp',
           'Fec'] <- Fec1[Fec1$Treatment == 'Comp',
                         'Fec']
  
  paramsD1[paramsD1$Treatment == 'Control',
           'Fec'] <- Fec1[Fec1$Treatment == 'Control',
                         'Fec']
  
  # replace with observed means if no vital rate cannot be calculated from bootstrap
  # sample
  for(i in 1:dim(paramsD1)[1]){
    for(x in 1:dim(paramsD1)[2]){
      if(is.na(paramsD1[i, x])){
        paramsD1[i, x] <- paramsD[i, x]
      }
    }
  }
  

  ##Control
  s1_c=1
  s2_c=paramsD1[paramsD1$Treatment=="Control","s2"] %>% as.numeric
  f_c=paramsD1[paramsD1$Treatment=="Control","Fec"] %>% as.numeric
  boot_s2_c[j]=s2_c
  boot_f_c[j]=f_c
  
  ##Comp Rem
  s1_cr=1
  s2_cr=paramsD1[paramsD1$Treatment=="Comp","s2"] %>% as.numeric
  f_cr=paramsD1[paramsD1$Treatment=="Comp","Fec"] %>% as.numeric
  boot_s2_cr[j]=s2_cr
  boot_f_cr[j]=f_cr
  
  
  ##Bootstrapped matrices!
  G1=.5503
  G2=.3171
  v=.8228
  
  A_cont=matrix(c(1-G2,0,f_c*v*(1-G1),
                  G2*s1_c,0,f_c*G1*s1_c,
                  0,s2_c,0),
                nrow=3, byrow=TRUE, 
                dimnames=list(c("SB","Rosette","RA"),c("SB","Rosette","RA")))
  
  A_cr=matrix(c(1-G2,0,f_cr*v*(1-G1),
                G2*s1_cr,0,f_cr*G1*s1_cr,
                0,s2_cr,0),
              nrow=3, byrow=TRUE, 
              dimnames=list(c("SB","Rosette","RA"),c("SB","Rosette","RA")))

  ev <- eigen(A_cont)
  
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lmax <- lmax[1]  
  lambda_cont <- Re(ev$values[lmax])
   boot_l_c[j]=lambda_cont
  
  ev <- eigen(A_cr)
  
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lmax <- lmax[1]  
  lambda_cr <- Re(ev$values[lmax])
  boot_l_cr[j]=lambda_cr
  
}


boot_s2_c=sort(boot_s2_c)
boot_s2_cr=sort(boot_s2_cr)

boot_f_c=sort(boot_f_c)
boot_f_cr=sort(boot_f_cr)

boot_l_c=sort(boot_l_c)
boot_l_cr=sort(boot_l_cr)

##Extract Confidence Intervals
lower=c(boot_s2_c[25], boot_s2_cr[25],
        boot_f_c[25], boot_f_cr[25], 
        boot_l_c[25], boot_l_cr[25])

upper=c( boot_s2_c[975], boot_s2_cr[975],
        boot_f_c[975], boot_f_cr[975], 
        boot_l_c[975], boot_l_cr[975])

results=data.frame(values, lower, upper) 

results$Trtvalue=c("Control","CR")
results$Trtvalue=factor(results$Trtvalue,
                        levels=c("Control","CR"))
S2=results[1:2,]
fec=results[3:4,]
lambda=results[5:6,]


par(mfrow = c(1, 3), mar = c(5,6,4,2) + 0.15)

plot(as.integer(S2$Trtvalue), S2$values,
     pch = 1,
     ylim = c(0, 1),
     axes = FALSE,
     main = "",
     xlab = "", 
     ylab = "Mean Survival (S)",
     cex.lab = 1.8)
arrows(as.integer(S2$Trtvalue), S2$lower, 
       as.integer(S2$Trtvalue),S2$upper, 
       length = 0.05,
       angle = 90,
       code = 3)
mtext(c("Control","CR"),
      side = 1,
      line = 1, at = c(1, 2),
      cex = 1.0)
axis(2, cex.axis = 1)
box(lwd = 2)

plot(as.integer(fec$Trtvalue), fec$values,
     pch = 1,
     ylim = c(0, 600),
     axes = FALSE,
     main = "", 
     xlab = "Treatment",
     ylab = "Mean Fecundity (F)",
     cex.lab = 1.8)
arrows(as.integer(fec$Trtvalue), fec$lower, 
       as.integer(fec$Trtvalue),fec$upper, 
       length = 0.05,
       angle = 90,
       code = 3)
mtext(c("Control","CR"),
      side = 1,
      line = 1,
      at = c(1, 2),
      cex = 1.0)
axis(2, cex.axis = 1)
box(lwd = 2)

plot(as.integer(lambda$Trtvalue), lambda$values,
     pch = 1,
     ylim = c(0, 3),
     axes = FALSE,
     main = "", 
     xlab = "",
     ylab = expression(paste('Lambda (', lambda, ')')),
     cex.lab = 1.8)
arrows(as.integer(lambda$Trtvalue), lambda$lower, 
       as.integer(lambda$Trtvalue),lambda$upper,
       length = 0.05,
       angle = 90,
       code = 3)
mtext(c("Control", "CR"),
      side = 1,
      line = 1,
      at = c(1, 2),
      cex = 1.0)
axis(2, cex.axis = 1)
box(lwd = 2)
