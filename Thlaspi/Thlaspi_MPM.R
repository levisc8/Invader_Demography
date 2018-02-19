##Clear environment and console
rm(list=ls(all=TRUE)) 
cat("\014") 
library(dplyr)

##Import data
ks=read.csv("C:/Users/sl13sise/Dropbox/invaders demography/Thlaspi/ThlaspiD.csv",
            stringsAsFactors = FALSE) %>%
  filter(Treatment == 'Control' |
           Treatment == 'Comp')

params <- ks %>%
  group_by(Treatment) %>%
  summarise(Survival = mean(Alive, na.rm = TRUE),
            Fruit = mean(Fruit, na.rm = TRUE),
            nAlive = sum(Alive, na.rm = TRUE))

params <- as.data.frame(params)

F_C <- params$Fruit[params$Treatment == 'Control'] * 3.76 # mean seeds/fruit
F_CR <- params$Fruit[params$Treatment == 'Comp'] * 3.76 # mean seeds/fruit

S_C <- params$Survival[params$Treatment == 'Control']
S_CR <- params$Survival[params$Treatment == 'Comp']


# All rates below taken from Burns et al 2012/2013
r = 0.67 # remain in seedbank
e = 0.12 # emerge from seedbank
g = 11 / (400 * .95)
v = 0.95


A_c <- matrix(c(r, S_C * F_C * (1 - g) * v,
                e, S_C * F_C * g * v),
              byrow = TRUE, nrow = 2)
 
A_cr <- matrix(c(r, S_CR * F_CR * (1 - g) * v,
                 e, S_CR * F_CR * g * v),
               nrow = 2, byrow = TRUE)


ev <- eigen(A_c)

lmax <- which(Re(ev$values) == max(Re(ev$values)))
lmax <- lmax[1]  
lambda_c <- Re(ev$values[lmax])
lambda_c

ev <- eigen(A_cr)

lmax <- which(Re(ev$values) == max(Re(ev$values)))
lmax <- lmax[1]  
lambda_cr <- Re(ev$values[lmax])
lambda_cr



values <- c(S_C, S_CR,
            F_C, F_CR,
            lambda_c, lambda_cr)

####################################################
#Bootstrapping code
nreps <- 1000
boot_s_c <- rep(NA, nreps)
boot_s_cr <- rep(NA, nreps)

boot_f_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)

boot_l_c <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)

#subset plants into treatments
c <- subset(ks, Treatment == "Control")
cr <- subset(ks, Treatment == "Comp")

nc = length(c[ ,1])
ncr = length(cr[ ,1])

#start the loop
for(j in 1:nreps) {
  
  x1 <- sample(1:nc, nc, replace = TRUE)
  x3 <- sample(1:ncr, ncr, replace = TRUE)

  bootc <- c[x1, ]
  bootcr <- cr[x3, ]
  bootdata <- rbind(bootc, bootcr)
  
  params1 <- bootdata %>% 
    group_by(Treatment) %>%
    summarise(Survival = mean(Alive, na.rm = TRUE),
             Fruit = mean(Fruit, na.rm = TRUE),
             nAlive = sum(Alive, na.rm = TRUE))
  
  params1 <- as.data.frame(params1)
  
  for(k in 1:length(params1$Treatment)){
    for(l in 1:length(params1)){
      if(is.na(params1[k, l])){
        params1[k, l] <- params[k, l]
      }
    }
  }
  
  
  S_C <- params1[params1$Treatment == "Control", "Survival"]
  S_CR <- params1[params1$Treatment == "Comp", "Survival"]
  
  F_C <- params1[params1$Treatment == "Control", "Fruit"] * 3.76
  F_CR <- params1[params1$Treatment == "Comp", "Fruit"] * 3.76
  
  boot_s_c[j] = S_C
  boot_s_cr[j] = S_CR
  
  boot_f_c[j] = F_C
  boot_f_cr[j] = F_CR
  
  r = 0.67 # remain in seedbank
  e = 0.12 # emerge from seedbank
  g = 11 / (400 * .95)
  v = 0.95
  
  
  A_c <- matrix(c(r, S_C * F_C * (1 - g) * v,
                  e, S_C * F_C * g * v),
                byrow = TRUE, nrow = 2)
  
  A_cr <- matrix(c(r, S_CR * F_CR * (1 - g) * v,
                   e, S_CR * F_CR * g * v),
                 nrow = 2, byrow = TRUE)
  
  ev <- eigen(A_c)
  
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

boot_s_c <- sort(boot_s_c)
boot_s_cr <- sort(boot_s_cr)

boot_f_c <- sort(boot_f_c)
boot_f_cr <- sort(boot_f_cr)

boot_l_c <- sort(boot_l_c)
boot_l_cr <- sort(boot_l_cr)

lower <- c(boot_s_c[25], boot_s_cr[25],
           boot_f_c[25], boot_f_cr[25],
           boot_l_c[25], boot_l_cr[25])
upper <- c(boot_s_c[975], boot_s_cr[975],
           boot_f_c[975], boot_f_cr[975],
           boot_l_c[975], boot_l_cr[975])

results <- data.frame(values, lower, upper) 

results$Trtvalue <- c("Control", "CR")
results$Trtvalue <- factor(results$Trtvalue,
                           levels = c("Control", "CR"))

surv <- results[1:2, ]
fec <- results[3:4, ]
lambda <- results[5:6, ]

par(mfrow = c(1, 3), mar = c(5,6,4,2) + 0.2)
plot(as.integer(surv$Trtvalue),surv$values,
     pch=1, 
     ylim=c(0,1), 
     axes=FALSE,
     main="", 
     xlab="", 
     ylab="Survival",
     cex.lab=1.8)
arrows(as.integer(surv$Trtvalue), surv$lower, 
       as.integer(surv$Trtvalue),surv$upper,
       length=0.05, 
       angle=90, 
       code=3)
mtext(c("Control","CR"),
      side=1,
      line=1,
      at=c(1,2), 
      cex=1.0)
axis(2, cex.axis=1)
box(lwd=2)

plot(as.integer(fec$Trtvalue), fec$values,
     pch=1, 
     ylim=c(0, max(fec$upper) + 3), 
     axes=FALSE,
     main="", 
     xlab = 'Treatment',
     ylab="Mean Fecundity",
     cex.lab=1.8)
arrows(as.integer(fec$Trtvalue), fec$lower, 
       as.integer(fec$Trtvalue),fec$upper, 
       length=0.05,
       angle=90,
       code=3)
mtext(c("Control","CR"),
      side=1,
      line=1,
      at=c(1,2), 
      cex=1.0)
axis(2, cex.axis=1)
box(lwd=2)

plot(as.integer(lambda$Trtvalue), lambda$values,
     pch=1,
     ylim=c(0,max(lambda$upper)+1),
     axes=FALSE,
     main="", 
     xlab="", 
     ylab=expression(paste('Lambda (', lambda, ')')),
     cex.lab=1.8)
arrows(as.integer(lambda$Trtvalue), lambda$lower, 
       as.integer(lambda$Trtvalue),lambda$upper, 
       length=0.05,
       angle=90, 
       code=3)
mtext(c("Control","CR"),
      side=1,
      line=1,
      at=c(1,2), 
      cex=1.0)
axis(2, cex.axis=1)
box(lwd=2)

