##Clear environment and console
rm(list=ls(all=TRUE)) 
cat("\014")

library(dplyr)

#import data into R studio
setwd('C:/Users/sl13sise/Dropbox/invaders demography')
ks=read.csv("Carduus/tenative.carduus.model.file.csv") %>%
   filter(Treatment != 'Herbivore')

ks$Seeds=ks$Inflor15*88.5

# check for density dependence
ddTest <- group_by(ks,Treatment, Quad) %>% summarise(N=n(),
                                               Sb=mean(Sbetween,na.rm=TRUE),
                                               Seeds=mean(Seeds,na.rm=TRUE))

params=group_by(ks, Treatment) %>% summarise(
             N=n(),
             Sb=mean(Sbetween,na.rm=TRUE),
             Seeds=mean(Seeds,na.rm=TRUE))

##Control
s_c=params[params$Treatment=="Control","Sb"] %>% as.numeric
f_c=params[params$Treatment=="Control","Seeds"] %>% as.numeric

##Comp Rem
s_cr=params[params$Treatment=="Competitor","Sb"] %>% as.numeric 
f_cr=params[params$Treatment=="Competitor","Seeds"] %>% as.numeric


#seed viability and germination rate, Bolting probability set to 1 for the
##Time being, can be more precisely calculated later....
Gsb=.03
Ssb=.2597
B=1
e=.2333
Gi=.03

A_cont=matrix(c(Ssb, 0, f_c * e,
                Gsb, 0, f_c * Gi,
                0, s_c, 0),
              nrow = 3, byrow = TRUE, 
                  dimnames=list(c("SB", "Rosette", "RA"),
                                c("SB", "Rosette", "RA")))

A_cr=matrix(c(Ssb, 0, f_cr * e,
                  Gsb, 0, f_cr * Gi,
                  0, s_cr, 0),
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

values=c(s_c, s_cr,
         f_c, f_cr,
         lambda_cont, lambda_cr)


####################################################
#Bootstrapping code
#I do 1000 bootstrap replicates.  This will save all the data. 
##For now I fill in a column with 1000 rows with the value 'NA'.
##Later, I will replace the NAs with data.

nreps=1000
boot_s_c <- rep(NA, nreps)
boot_s_cr <- rep(NA, nreps)
boot_f_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)
boot_l_cont <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)

cplants=subset(ks,Treatment=="Control")
crplants=subset(ks,Treatment=="Competitor")


nc=length(cplants[,1])
ncr=length(crplants[,1])
#start the loop
for(j in 1:nreps) {
  
  #x3 is a column that picks n integers between the numbers 1 and n (because we have n plants in the data set).  
  #replace=TRUE means that we are sampling with replacement (once we choose a number, we can choose it again).
  
  x2=sample(1:nc, nc, replace = TRUE)
  x3=sample(1:ncr, ncr, replace = TRUE)
  
  bootc <- cplants[x2, ]
  bootcr <- crplants[x3, ]
  
  bootdata <- rbind(bootc, bootcr)
  
  params1 <- group_by(bootdata,Treatment) %>% summarise(
               N = length(Treatment),
               Sb = mean(Sbetween, na.rm = TRUE),
               Seeds = mean(Seeds, na.rm = TRUE))
  
  ##Replaces all NaNs from summary function with observed parameter means
  ## otherwise, bootstrap code will not work.

  for(k in 1:dim(params1)[1]){
    for(l in 1:dim(params1)[2]){
      if(is.na(params1[k, l])){
        params1[k, l] <- params[k, l]
      }
    }
  }
 
  ##Control
  s_c <- params1[params1$Treatment == "Control","Sb"] %>% as.numeric
  f_c <- params1[params1$Treatment == "Control","Seeds"] %>% as.numeric
  boot_s_c[j] <- s_c
  boot_f_c[j] <- f_c
  
  ##Comp Rem
  s_cr <- params1[params1$Treatment == "Competitor","Sb"] %>% as.numeric
  f_cr <- params1[params1$Treatment == "Competitor","Seeds"] %>% as.numeric
  boot_s_cr[j] <- s_cr
  boot_f_cr[j] <- f_cr
  
  
  ##Bootstrapped matrices!
  Gsb=.03
  Ssb=.2597
  B=1
  e=.2333
  Gi=.03
  
  A_cont=matrix(c(Ssb,0,f_c*e,
                  Gsb,0,f_c*Gi,
                  0,s_c,0),
                nrow=3, byrow=TRUE, 
                dimnames=list(c("SB","Rosette","RA"),c("SB","Rosette","RA")))
  
  A_cr=matrix(c(Ssb,0,f_cr*e,
                Gsb,0,f_cr*Gi,
                0,s_cr,0),
              nrow=3, byrow=TRUE, 
              dimnames=list(c("SB","Rosette","RA"),c("SB","Rosette","RA")))
             
                 
             ev <- eigen(A_cont)
             
             lmax <- which(Re(ev$values) == max(Re(ev$values)))
             lmax <- lmax[1]  
             lambda_cont <- Re(ev$values[lmax])
             boot_l_cont[j]=lambda_cont
             
             ev <- eigen(A_cr)
             
             lmax <- which(Re(ev$values) == max(Re(ev$values)))
             lmax <- lmax[1]  
             lambda_cr <- Re(ev$values[lmax])
             boot_l_cr[j]=lambda_cr
             
}

boot_s_c=sort(boot_s_c)
boot_s_cr=sort(boot_s_cr)
boot_f_c=sort(boot_f_c)
boot_f_cr=sort(boot_f_cr)
boot_l_cont=sort(boot_l_cont)
boot_l_cr=sort(boot_l_cr)

##Extract Confidence Intervals
lower=c(boot_s_c[25], boot_s_cr[25],
        boot_f_c[25], boot_f_cr[25], 
        boot_l_cont[25], boot_l_cr[25])

upper=c(boot_s_c[975], boot_s_cr[975],
        boot_f_c[975], boot_f_cr[975], 
        boot_l_cont[975], boot_l_cr[975])

results=data.frame(values, lower, upper) 

results$Trtvalue=c("Control","CR")
results$Trtvalue=factor(results$Trtvalue,
                        levels=c("Control","CR"))
S1=results[1:2,]
fec=results[3:4,]
lambda=results[5:6,]


par(mfrow = c(1, 3), mar = c(5,6,4,2) + 0.1)
plot(as.integer(S1$Trtvalue), S1$values,
     pch = 1, 
     ylim = c(0,1), 
     axes = FALSE,
     xlab = "",
     ylab = expression(paste("Mean Survival (S "['j']*')')),
     cex.lab=1.8)

arrows(as.integer(S1$Trtvalue), S1$lower, 
       as.integer(S1$Trtvalue),S1$upper, 
       length=0.05, angle=90, code=3)
mtext(c("Control","CR"),side=1,line=1,at=c(1,2), cex=1.0)
axis(2, cex.axis=1)
box(lwd=2)


plot(as.integer(fec$Trtvalue), fec$values,pch = 1,
     ylim = c(0,1500), axes = FALSE, xlab = "Treatment",
     ylab = "Mean Fecundity (F)", cex.lab = 1.8)
arrows(as.integer(fec$Trtvalue), fec$lower, 
       as.integer(fec$Trtvalue),fec$upper, 
       length = 0.05, angle = 90, code = 3)
mtext(c("Control", "CR"),
      side = 1, 
      line = 1,
      at = c(1, 2),
      cex = 1.0)
axis(2, cex.axis = 1)
box(lwd = 2)

plot(as.integer(lambda$Trtvalue), lambda$values,
     pch = 1,
     ylim = c(0,6), 
     axes = FALSE,
     xlab = "", 
     ylab = expression(paste('Lambda (', lambda, ')', sep = "")),
     cex.lab = 1.8)

arrows(as.integer(lambda$Trtvalue), lambda$lower, 
       as.integer(lambda$Trtvalue),lambda$upper, 
       length = 0.05, angle = 90, code = 3)
mtext(c("Control","CR"),
      side = 1,
      line = 1,
      at = c(1, 2),
      cex = 1.0)
axis(2, cex.axis = 1)
box(lwd = 2)

