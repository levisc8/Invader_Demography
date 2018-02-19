##K striata model with 2015 survival and 2014 fecundity!
#import data into R studio
rm(list=ls(all=TRUE)) 
cat("\014")

library(dplyr)

#import data into R studio
ks15=read.csv("C:/Users/sl13sise/Dropbox/invaders demography/Kummerowia/kstra2015.2 AG.csv") %>%
  filter(Treatment != 'Herb')
ks14=read.csv("C:/Users/sl13sise/Dropbox/invaders demography/Kummerowia/kummerowia AG.csv") %>%
  filter(Treatment != 'Herb')

survival=summarise(group_by(ks15,Treatment),
                 Survival=mean(Alive,na.rm=T))
fec=summarise(group_by(ks14,Treatment),
                 Seeds=mean(f,na.rm=T))

params=full_join(survival,fec,by="Treatment")

params=as.data.frame(params)

s_c=params[params$Treatment=="Control","Survival"]
f_c=params[params$Treatment=="Control","Seeds"]

s_cr=params[params$Treatment=="Comp","Survival"]
f_cr=params[params$Treatment=="Comp","Seeds"]

V=0.6622*0.92 #seed viability
G0=0.92 #proportion of viable seeds that germinate immediately (mean of above and buried seeds, 2014/15 germ exp)
G1=0.04 #proportion of viable seeds that germinate after one year in Seed bank
G2=0.04 #proportion of viable seeds that germinate after two years in Seed bank

A_cont=matrix(c(0, 0, (f_c*V*G2),
                1, 0, (f_c*V*G1),
                0, 1, (s_c*f_c*V*G0)),
              nrow=3, byrow=TRUE, 
              dimnames=list(c("seedbank2", "seedbank1", "plant"),
                            c("seedbank2", "seedbank1", "plant")))
A_cr=matrix(c(0, 0, (f_cr*V*G2),
              1, 0, (f_cr*V*G1),
              0, 1, (s_cr*f_cr*V*G0)),
            nrow=3, byrow=TRUE, 
            dimnames=list(c("seedbank2", "seedbank1", "plant"),
                          c("seedbank2", "seedbank1", "plant")))

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

values=c(s_c, s_cr, 
         f_c, f_cr,
         lambda_c, lambda_cr)


####################################################
#Bootstrapping code
nreps=1000
boot_s_c <- rep(NA, nreps)
boot_s_cr <- rep(NA, nreps)
boot_f_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)
boot_l_c <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)
es_crcontrol <- rep(NA, nreps)

cr15=subset(ks15,Treatment=="Comp")
cr14=subset(ks14,Treatment=="Comp")
cont15=subset(ks15,Treatment=="Control")
cont14=subset(ks14,Treatment=="Control")
n2=length(cr15[,1])
n3=length(cont15[,1])
n5=length(cr14[,1])
n6=length(cont14[,1])

for(j in 1:nreps) 
{
  
  #x3 is a column that picks n integers between the numbers 1 and n (because we have n plants in the data set).  
  #replace=TRUE means that we are sampling with replacement (once we choose a number, we can choose it again).
  
  x2=sample(1:n2,n2,replace=T)
  x3=sample(1:n3,n3,replace=T)
  x5=sample(1:n5,n5,replace=T)
  x6=sample(1:n6,n6,replace=T)
  
  bootcomp15=cr15[x2,]
  bootc15=cont15[x3,]
  bootcomp14=cr14[x5,]
  bootcont14=cont14[x6,]
  bootdata15=rbind(bootcomp15,bootc15)
  bootdata14=rbind(bootcomp14,bootcont14)
  
  Survival1=summarise(group_by(bootdata15,Treatment),
                    Survival=mean(Alive,na.rm=T))
  Fec1=summarise(group_by(bootdata14,Treatment),
                 Seeds=mean(f,na.rm=T))
  params1=full_join(Survival1,Fec1,by="Treatment")
  
  params1=as.data.frame(params1)
  
  for(k in 1:length(params$Treatment)){
    for(l in 1:length(params)){
      if(is.na(params1[k,l])==T){
        params1[k,l]=params[k,l]
      }
    }
  }

  s_c=params1[params1$Treatment=="Control","Survival"]
  f_c=params1[params1$Treatment=="Control","Seeds"]
  boot_s_c[j]=s_c
  boot_f_c[j]=f_c
  
  s_cr=params1[params1$Treatment=="Comp","Survival"]
  f_cr=params1[params1$Treatment=="Comp","Seeds"]
  boot_s_cr[j]=s_cr
  boot_f_cr[j]=f_cr
  
  V=0.6622*0.92 #seed viability
  G0=0.92 #proportion of viable seeds that germinate immediately (mean of above and buried seeds, 2014/15 germ exp)
  G1=0.04 #proportion of viable seeds that germinate after one year in Seed bank
  G2=0.04 #proportion of viable seeds that germinate after two years in Seed bank
  
   A_cont=matrix(c(0, 0, (f_c*V*G2),
                  1, 0, (f_c*V*G1),
                  0, 1, (s_c*f_c*V*G0)),
                 nrow=3,
                 byrow=TRUE,
                 dimnames=list(c("seedbank2", "seedbank1", "plant"),
                               c("seedbank2", "seedbank1", "plant")))
  A_cr=matrix(c(0, 0, (f_cr*V*G2),
                1, 0, (f_cr*V*G1),
                0, 1, (s_cr*f_cr*V*G0)),
              nrow=3,
              byrow=TRUE,
              dimnames=list(c("seedbank2", "seedbank1", "plant"),
                            c("seedbank2", "seedbank1", "plant")))

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
es_crcontrol <- log(boot_l_cr+1)-log(boot_l_c+1)
escr=sort(es_crcontrol)

es_crcontrol2=log(boot_l_cr+.5)-log(boot_l_c+.5)
escr2=sort(es_crcontrol2)


boot_s_c=sort(boot_s_c)
boot_s_cr=sort(boot_s_cr)
boot_f_c=sort(boot_f_c)
boot_f_cr=sort(boot_f_cr)
boot_l_c=sort(boot_l_c)
boot_l_cr=sort(boot_l_cr)


lower=c(boot_s_c[25], boot_s_cr[25],
        boot_f_c[25], boot_f_cr[25], 
        boot_l_c[25], boot_l_cr[25])
upper=c(boot_s_c[975], boot_s_cr[975],
        boot_f_c[975], boot_f_cr[975], 
        boot_l_c[975], boot_l_cr[975])

results=data.frame(values, lower, upper) 

results$Trtvalue=c("Control","CR")
results$Trtvalue=factor(results$Trtvalue,
                        levels=c("Control","CR"))
surv=results[1:2,]
fec=results[3:4,]
lambda=results[5:6,]


par(mfrow = c(1, 3), mar = c(5,6,4,2) + 0.2)
plot(as.integer(surv$Trtvalue), surv$values,
     pch = 1, 
     ylim = c(0,1), 
     axes = FALSE,
     main = "",
     xlab = "",
     ylab = "Mean Survival",
     cex.lab = 1.8)
arrows(as.integer(surv$Trtvalue), surv$lower,
       as.integer(surv$Trtvalue),surv$upper,
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

plot(as.integer(fec$Trtvalue), fec$values,
     pch = 1, 
     ylim = c(0, max(fec$upper) + 20), 
     axes = FALSE,
     main = "", 
     xlab = "Treatment",
     ylab = "Mean Fecundity",
     cex.lab = 1.8)
arrows(as.integer(fec$Trtvalue), fec$lower, 
       as.integer(fec$Trtvalue), fec$upper,
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
     ylim = c(0, max(lambda$upper) + 3), 
     axes = FALSE,
     main = "",
     xlab = "", 
     ylab = expression(paste('Lambda (', lambda, ')')),
     cex.lab = 1.8)
arrows(as.integer(lambda$Trtvalue), lambda$lower, 
       as.integer(lambda$Trtvalue), lambda$upper, 
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

effectsizes <- c(mean(escr), escr[25], escr[975])
effectsizes2 <- c(mean(escr2), escr2[25], escr2[975])
