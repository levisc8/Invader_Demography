##A quick note: I've tried to add as detailed annotations as possible, but
##the following code for Potentilla is a bit confusing,
##so if you have any questions (or find any typos), let me know!

### clear everything ###
rm(list=ls(all=TRUE)) 
cat("\014") 
library(dplyr)

#import data into R studio
ks <- read.csv("C:/Users/sl13sise/Dropbox/invaders demography/Potentilla/potrec20154R AG.csv") %>%
  filter(Treatment != 'Herb')

##We found an average of 87.3 seeds/fruit in the field
ks$Seeds <- ks$Fruits15 * 87.3

##We also have a separate seedling survival data file, which makes this code a bit more complicated
##if you have it all in one file, then you won't need a couple steps outlined here.
sdlsurv <- read.csv("C:/Users/sl13sise/Dropbox/invaders demography/Potentilla/Potentilla Seedlings 2015 AG.csv")%>%
  filter(Treatment != 'Herb')


##Subset seedlings into treatments for use later
csdl <- subset(sdlsurv, Treatment == "Control")
crsdl <- subset(sdlsurv, Treatment == "Comp")

params <- ks %>%
  group_by(Treatment, Stg14) %>%
  summarise(N = length(Treatment),
            Survival = mean(Alive15, na.rm = TRUE),
            NRA.15 = length(Treatment[Stg15 == "NRA"]),
            RA.15 = length(Treatment[Stg15 == "RA"]),
            N15 = length(Treatment[Alive15 != 0]),
            pf = RA.15 / N15)

params <- as.data.frame(params)

##extract model paramater values from "params" dataframe we created above
##A quick reference guide for vital rates
##S1_T=Sdl survival, see notes in Herbivore removal for additional details
##S2_T=NRA survival
##S3_T=RA Survival
##Pfsdl_T=Probability that a seedling in 2014 flowers in 2015
##Pfnra_T=probability that an NRA in 2014 flowers in 2015
##Pfra_T=probablity that an RA in 2014 flowers in 2015
##F_T=seed production in 2015

##Control
s1a_c=mean(csdl$Alive.815,na.rm=T)
s1b_c=params[params$Stg14=="SDL"&params$Treatment=="Control","Survival"]
s1_c=s1a_c*s1b_c
s2_c=params[params$Stg14=="NRA"&params$Treatment=="Control","Survival"]
s3_c=params[params$Stg14=="RA"&params$Treatment=="Control","Survival"]
pfsdl_c=params[params$Stg14=="SDL"&params$Treatment=="Control","pf"]
pfnra_c=params[params$Stg14=="NRA"&params$Treatment=="Control","pf"]
pfra_c=params[params$Stg14=="RA"&params$Treatment=="Control","pf"]
f_c=mean(ks[ks$Treatment=="Control","Seeds"],na.rm=T)

##Comp Rem
s1a_cr=mean(crsdl$Alive.815,na.rm=T)
s1b_cr=params[params$Stg14=="SDL"&params$Treatment=="Comp","Survival"]
s1_cr=s1b_cr*s1a_cr
s2_cr=params[params$Stg14=="NRA"&params$Treatment=="Comp","Survival"]
s3_cr=params[params$Stg14=="RA"&params$Treatment=="Comp","Survival"]
pfsdl_cr=params[params$Stg14=="SDL"&params$Treatment=="Comp","pf"]
pfnra_cr=params[params$Stg14=="NRA"&params$Treatment=="Comp","pf"]
pfra_cr=params[params$Stg14=="RA"&params$Treatment=="Comp","pf"]
f_cr=mean(ks[ks$Treatment=="Comp","Seeds"],na.rm=T)

#seed viability and germination rate
##All rates are based on Kiemnec and McInnis 2009 study
Gi=.177 #Immediate germination
v=.218 #Viability
Ssb=.18 #Survival of a seed in the seedbank
Gsb=.025 #Germination rate of seeds in seedbank

A_cont=matrix(c(Ssb,0,0,f_c*v*(1-Gi),
                Gsb,0,0,f_c*v*Gi,
                0,s1_c*(1-pfsdl_c),s2_c*(1-pfnra_c),s3_c*(1-pfra_c),
                0,s1_c*pfsdl_c,s2_c*pfnra_c,s3_c*pfra_c),
              nrow=4, byrow=TRUE, dimnames=list(c("SB","SDL","NRA","RA"),
                                                c("SB","SDL","NRA","RA")))
A_cr=matrix(c(Ssb,0,0,f_cr*v*(1-Gi),
              Gsb,0,0,f_cr*v*Gi,
              0,s1_cr*(1-pfsdl_cr),s2_cr*(1-pfnra_cr),s3_cr*(1-pfra_cr),
              0,s1_cr*pfsdl_cr,s2_cr*pfnra_cr,s3_cr*pfra_cr),
            nrow=4, byrow=TRUE, dimnames=list(c("SB","SDL","NRA","RA"),
                                              c("SB","SDL","NRA","RA")))

##Obtain lambda values from matrices

ev <- eigen(A_cont)

lmax <- which(Re(ev$values) == max(Re(ev$values)))
lmax <- lmax[1]  
lambda_c <- Re(ev$values[lmax])
lambda_c

ev <- eigen(A_cr)

lmax <- which(Re(ev$values) == max(Re(ev$values)))
lmax <- lmax[1]  
lambda_cr <- Re(ev$values[lmax])
lambda_cr


##Store vital rates observed in field data for later
values=c(s1_c, s1_cr,
         s2_c, s2_cr,
         s3_c, s3_cr,
         pfsdl_c, pfsdl_cr,
         pfnra_c, pfnra_cr,
         pfra_c, pfra_cr,
         f_c, f_cr,
         lambda_c, lambda_cr)


####################################################
#Bootstrapping code
##We bootstrap the data set 1000 times, running a new matrix for each one, then store
##the bootstrapped vital rates in the vectors created below. 
nreps=1000
boot_s1_c <- rep(NA, nreps)
boot_s1_cr <- rep(NA, nreps)
boot_s2_c=rep(NA,nreps)
boot_s2_cr=rep(NA,nreps)
boot_s3_c=rep(NA,nreps)
boot_s3_cr=rep(NA,nreps)

boot_pfsdl_c=rep(NA,nreps)
boot_pfsdl_cr=rep(NA,nreps)
boot_pfnra_c=rep(NA,nreps)
boot_pfnra_cr=rep(NA,nreps)
boot_pfra_c=rep(NA,nreps)
boot_pfra_cr=rep(NA,nreps)

boot_f_c=rep(NA,nreps)
boot_f_cr=rep(NA,nreps)

boot_l_c <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)

##Subset by treatment and stage class so that each one can be bootstrapped separately
##This ensures that each bootstrapped data frame contains the same number of Sdls,NRAs, and RAs
##from each treatment as our field data does. Each one is sampled with replacement though,
##so the dataset is different in each of the 1000 runs.
ccsdl=subset(ks,Treatment=="Control"&Stg14=="SDL")
crcrsdl=subset(ks,Treatment=="Comp"&Stg14=="SDL")
cnra=subset(ks,Treatment=="Control"&Stg14=="NRA")
crnra=subset(ks,Treatment=="Comp"&Stg14=="NRA")
cra=subset(ks,Treatment=="Control"&Stg14=="RA")
crra=subset(ks,Treatment=="Comp"&Stg14=="RA")

##calculate number of plants in each stage class-treatment combination, so R knows
##how any to put into each one when bootstrapping

#NOTE: if you aren't using s1a and s1b, then you do not need the first three lines here (highlighted lines).
##In other words, if you don't have the equivalent of the "Sdlsurv"data frame, then DO NOT
##include these lines
ncsdl=length(csdl[,c(1)])
ncrsdl=length(crsdl[,1])
ncnra=length(cnra[,1])
ncrnra=length(crnra[,1])
ncra=length(cra[,1])
ncrra=length(crra[,1])
nncsdl=length(ccsdl[,1])
nncrsdl=length(crcrsdl[,1])

#start the loop
for(j in 1:nreps) 
{
  
  #xN is a column that picks n integers between the numbers 1 and n (because we have n plants in the data set). 
  ##NOTE: as noted above, depending on you have your data formatted, X1-3 may not be needed. If not, remove them 
  
  x1=sample(1:ncsdl,ncsdl,replace=T)
  x2=sample(1:ncrsdl,ncrsdl,replace=T)
  x4=sample(1:ncnra,ncnra,replace=T)
  x5=sample(1:ncrnra,ncrnra,replace=T)
  x7=sample(1:ncra,ncra,replace=T)
  x8=sample(1:ncrra,ncrra,replace=T)
  x10=sample(1:nncsdl,nncsdl,replace=T)
  x11=sample(1:nncrsdl,nncrsdl,replace=T)

  
  ##bootstrap each one separately, then recombine them for use later
  ##NOTE: highlighted lines. If you have all of the seedling data you need in one file, remove the highlighted lines
  bootcsdl <- csdl[x1,]
  bootcrsdl=crsdl[x2,]
  bootcnra=cnra[x4,]
  bootcrnra=crnra[x5,]
  bootcra=cra[x7,]
  bootcrra=crra[x8,]
  dbootcsdl=ccsdl[x10,]
  dbootcrsdl=crcrsdl[x11,]

  bootdata=rbind(bootcnra,
                 bootcrnra,
                 bootcra,
                 bootcrra,
                 dbootcsdl,
                 dbootcrsdl)
  ##Delete the bootsdl data frame if you do not have an equivalent to the Sdlsurv data file
  bootsdl=rbind(bootcsdl,bootcrsdl)

  params1=summarise(group_by(bootdata,Treatment,Stg14),
                   N=length(Treatment),
                   Survival=mean(Alive15,na.rm=T),
                   NRA.15=length(Treatment[Stg15=="NRA"]),
                   RA.15=length(Treatment[Stg15=="RA"]),
                   N15=length(Treatment[Alive15!=0]),
                   pf=RA.15/N15)
  
  params1=as.data.frame(params1)
  ##These next few lines go through the bootstrapped vital rate data frame
  ##and remove any NaNs by replacing them with the value we observed in our field data.
  ##This won't alter the confidence intervals, and keeps the loop running smoothly.
  ##Occasionally when R is resampling, it will only choose dead plants, which 
  ##causes the loop to stop and return an error. This prevents that from happening.

  for(k in 1:length(params$Treatment)){
    for(l in 1:length(params)){
      if(is.na(params1[k,l])==T){
        params1[k,l]=params[k,l]
      }
    }
  }
  
  #extract model paramater values from "params" dataframe we created above!
  ##Control
  s1a_c=mean(bootcsdl$Alive.815,na.rm=T)
  s1b_c=params1[params1$Stg14=="SDL"&params1$Treatment=="Control","Survival"]
  s1_c=s1a_c*s1b_c
  s2_c=params1[params1$Stg14=="NRA"&params1$Treatment=="Control","Survival"]
  s3_c=params1[params1$Stg14=="RA"&params1$Treatment=="Control","Survival"]
  pfsdl_c=params1[params1$Stg14=="SDL"&params1$Treatment=="Control","pf"]
  pfnra_c=params1[params1$Stg14=="NRA"&params1$Treatment=="Control","pf"]
  pfra_c=params1[params1$Stg14=="RA"&params1$Treatment=="Control","pf"]
  f_c=mean(bootdata[bootdata$Treatment=="Control","Seeds"],na.rm=T)
  
  boot_s1_c[j]=s1_c
  boot_s2_c[j]=s2_c
  boot_s3_c[j]=s3_c
  boot_pfsdl_c[j]=pfsdl_c
  boot_pfnra_c[j]=pfnra_c
  boot_pfra_c[j]=pfra_c
  boot_f_c[j]=f_c
  
  
  ##Comp Rem
  s1a_cr=mean(bootcrsdl$Alive.815,na.rm=T)
  s1b_cr=params1[params1$Stg14=="SDL"&params1$Treatment=="Comp","Survival"]
  s1_cr=s1a_cr*s1b_cr
  s2_cr=params1[params1$Stg14=="NRA"&params1$Treatment=="Comp","Survival"]
  s3_cr=params1[params1$Stg14=="RA"&params1$Treatment=="Comp","Survival"]
  pfsdl_cr=params1[params1$Stg14=="SDL"&params1$Treatment=="Comp","pf"]
  pfnra_cr=params1[params1$Stg14=="NRA"&params1$Treatment=="Comp","pf"]
  pfra_cr=params1[params1$Stg14=="RA"&params1$Treatment=="Comp","pf"]
  f_cr=mean(bootdata[bootdata$Treatment=="Comp","Seeds"],na.rm=T)
  
  boot_s1_cr[j]=s1_cr
  boot_s2_cr[j]=s2_cr
  boot_s3_cr[j]=s3_cr
  boot_pfsdl_cr[j]=pfsdl_cr
  boot_pfnra_cr[j]=pfnra_cr
  boot_pfra_cr[j]=pfra_cr
  boot_f_cr[j]=f_cr
  
  ##These matrices are the same as above, but use the boot strapped vital rates instead
  ##instead of values from field observations.
  Gi=.177
  v=.218
  Ssb=.18
  Gsb=.025
  
  A_cont=matrix(c(Ssb,0,0,f_c*v*(1-Gi),
                  Gsb,0,0,f_c*v*Gi,
                  0,s1_c*(1-pfsdl_c),s2_c*(1-pfnra_c),s3_c*(1-pfra_c),
                  0,s1_c*pfsdl_c,s2_c*pfnra_c,s3_c*pfra_c),
                nrow=4, byrow=TRUE, dimnames=list(c("SB","SDL","NRA","RA"),
                                                  c("SB","SDL","NRA","RA")))
  A_cr=matrix(c(Ssb,0,0,f_cr*v*(1-Gi),
                Gsb,0,0,f_cr*v*Gi,
                0,s1_cr*(1-pfsdl_cr),s2_cr*(1-pfnra_cr),s3_cr*(1-pfra_cr),
                0,s1_cr*pfsdl_cr,s2_cr*pfnra_cr,s3_cr*pfra_cr),
              nrow=4, byrow=TRUE, dimnames=list(c("SB","SDL","NRA","RA"),
                                                c("SB","SDL","NRA","RA")))
  
 
  ev <- eigen(A_cont)
  
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lmax <- lmax[1]  
  lambda_c <- Re(ev$values[lmax])
  boot_l_c[j]=lambda_c
  
  ev <- eigen(A_cr)
  
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lmax <- lmax[1]  
  lambda_cr <- Re(ev$values[lmax])
  boot_l_cr[j]=lambda_cr
  
}

##Sort vectors with bootstrapped vital rates
boot_s1_c=sort(boot_s1_c)
boot_s1_cr=sort(boot_s1_cr)
boot_s2_c=sort(boot_s2_c)
boot_s2_cr=sort(boot_s2_cr)
boot_s3_c=sort(boot_s3_c)
boot_s3_cr=sort(boot_s3_cr)

boot_pfsdl_c=sort(boot_pfsdl_c)
boot_pfsdl_cr=sort(boot_pfsdl_cr)
boot_pfnra_c=sort(boot_pfnra_c)
boot_pfnra_cr=sort(boot_pfra_cr)
boot_pfra_c=sort(boot_pfra_c)
boot_pfra_cr=sort(boot_pfra_cr)

boot_f_c=sort(boot_f_c)
boot_f_cr=sort(boot_f_cr)
boot_l_c=sort(boot_l_c)
boot_l_cr=sort(boot_l_cr)

##Extract confidence intervals
lower=c(boot_s1_c[25], boot_s1_cr[25],
        boot_s2_c[25], boot_s2_cr[25],
        boot_s3_c[25], boot_s3_cr[25],
        boot_pfsdl_c[25], boot_pfsdl_cr[25],
        boot_pfnra_c[25], boot_pfnra_cr[25],
        boot_pfra_c[25], boot_pfra_cr[25],
        boot_f_c[25], boot_f_cr[25], 
        boot_l_c[25], boot_l_cr[25])


upper=c(boot_s1_c[975], boot_s1_cr[975],
        boot_s2_c[975], boot_s2_cr[975],
        boot_s3_c[975], boot_s3_cr[975],
        boot_pfsdl_c[975], boot_pfsdl_cr[975],
        boot_pfnra_c[975], boot_pfnra_cr[975],
        boot_pfra_c[975], boot_pfra_cr[975],
        boot_f_c[975], boot_f_cr[975], 
        boot_l_c[975], boot_l_cr[975])

##Create data frame of observed values and bootstrapped CIs
results=data.frame(values, lower, upper) 

##Add treatment column to results DF
results$Trtvalue=c("Control","CR")
results$Trtvalue=factor(results$Trtvalue,
                        levels=c("Control","CR"))

##Subset by vital rate to make graphing easier
s1=results[1:2, ]
s2=results[3:4, ]
s3=results[5:6, ]
PFSdl=results[7:8, ]
PFNra=results[9:10, ]
PFRa=results[11:12, ]
fec=results[13:14, ]
lambda=results[15:16, ]

##Graph results
par(mfrow=c(2,4))
plot(as.integer(s1$Trtvalue), 
     s1$values,pch=1,
     ylim=c(0,1), 
     axes=FALSE,
     main="", 
     xlab="",
     ylab="SDL Survival",
     cex.lab=1.5)
arrows(as.integer(s1$Trtvalue), s1$lower,
       as.integer(s1$Trtvalue),s1$upper, 
       length=0.05, 
       angle=90, 
       code=3)
mtext(c("Control","CR"),
      side=1,
      line=1,
      at=c(1,2),
      cex=.8)
axis(2, cex.axis=1)
box(lwd=2)

plot(as.integer(s2$Trtvalue), s2$values,
     pch=1, 
     ylim=c(0,1), 
     axes=FALSE,
     main="",
     xlab="",
     ylab="NRA Survival",
     cex.lab=1.5)
arrows(as.integer(s2$Trtvalue), s2$lower,
       as.integer(s2$Trtvalue),s2$upper,
       length=0.05, 
       angle=90, 
       code=3)
mtext(c("Control","CR"),
      side=1,
      line=1,
      at=c(1,2), 
      cex=.8)
axis(2, cex.axis=1)
box(lwd=2)

plot(as.integer(s3$Trtvalue), s3$values,
     pch=1,
     ylim=c(0,1), 
     axes=FALSE,
     main="",
     xlab="",
     ylab="RA Survival",
     cex.lab=1.5)
arrows(as.integer(s3$Trtvalue), s3$lower,
       as.integer(s3$Trtvalue),s3$upper,
       length=0.05,
       angle=90, 
       code=3)
mtext(c("Control","CR"),
      side=1,
      line=1,
      at=c(1,2),
      cex=.8)
axis(2, cex.axis=1)
box(lwd=2)

plot(as.integer(PFSdl$Trtvalue), PFSdl$values,
     pch=1, 
     ylim=c(0,1),
     axes=FALSE,
     main="",
     xlab="", 
     ylab="P(Flowering - Seedling)",
     cex.lab=1.5)
arrows(as.integer(PFSdl$Trtvalue), PFSdl$lower, 
       as.integer(PFSdl$Trtvalue),PFSdl$upper,
       length=0.05, 
       angle=90, 
       code=3)
mtext(c("Control","CR"),
      side=1,
      line=1,
      at=c(1,2),
      cex=.8)
axis(2, cex.axis=1)
box(lwd=2)

plot(as.integer(PFNra$Trtvalue), PFNra$values,
     pch=1, 
     ylim=c(0,1), 
     axes=FALSE,
     main="", 
     xlab="", 
     ylab="P(Flowering - NRAs)",
     cex.lab=1.5)
arrows(as.integer(PFNra$Trtvalue), PFNra$lower,
       as.integer(PFNra$Trtvalue),PFNra$upper,
       length=0.05,
       angle=90, 
       code=3)
mtext(c("Control","CR"),
      side=1,
      line=1,
      at=c(1,2),
      cex=.8)
axis(2, cex.axis=1)
box(lwd=2)


plot(as.integer(PFRa$Trtvalue), PFRa$values,
     pch=1,
     ylim=c(0,1),
     axes=FALSE,
     main="", 
     xlab="Treatment", 
     ylab="P(Flowering - RAs)",
     cex.lab=1.5)
arrows(as.integer(PFRa$Trtvalue), PFRa$lower, 
       as.integer(PFRa$Trtvalue),PFRa$upper, 
       length=0.05,
       angle=90, 
       code=3)
mtext(c("Control","CR"),
      side=1,
      line=1,
      at=c(1,2), 
      cex=.8)
axis(2, cex.axis=1)
box(lwd=2)


plot(as.integer(fec$Trtvalue), fec$values,
     pch=1,
     ylim=c(0,3000), 
     axes=FALSE,
     main="", 
     xlab="Treatment", 
     ylab="Mean Fecundity",
     cex.lab=1.5)
arrows(as.integer(fec$Trtvalue), fec$lower, 
       as.integer(fec$Trtvalue),fec$upper,
       length=0.05, 
       angle=90, 
       code=3)
mtext(c("Control","CR"),
      side=1,
      line=1,
      at=c(1,2),
      cex=.8)
axis(2, cex.axis=1)
box(lwd=2)

plot(as.integer(lambda$Trtvalue), lambda$values,
     pch=1,
     ylim=c(0,max(lambda$upper)+1),
     axes=FALSE,
     main="", 
     xlab="", 
     ylab="Lambda",
     cex.lab=1.5)
arrows(as.integer(lambda$Trtvalue), lambda$lower, 
       as.integer(lambda$Trtvalue),lambda$upper,
       length=0.05, 
       angle=90, 
       code=3)
mtext(c("Control", "CR"),
      side=1,
      line=1,
      at=c(1,2),
      cex=.8)
axis(2, cex.axis=1)
box(lwd=2)

