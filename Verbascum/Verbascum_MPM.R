##Clear console and environment
rm(list = ls()) 
cat("\014") 

#import data into R studio
ks <- read.csv("C:/Users/sl13sise/Dropbox/invaders demography/Verbascum/Herbdemo4R.csv",
            stringsAsFactors = FALSE)

ks <- ks[ ,-c(1, 3, 9:16)]
names(ks) <- c("Treatment","Plant","Diameter","Notes","Therb","Alive1014",
            "Seeds","Alive15","Notes2")

library(dplyr)
ks <- filter(ks, Treatment != 'Herb')
params <- ks %>%
  group_by(Treatment) %>%
  summarise(N = length(Treatment),
            Sb = mean(Alive15, na.rm = TRUE),
            Seeds = mean(Seeds, na.rm = TRUE))

params <- as.data.frame(params)

#extract model paramater values from "params" dataframe we created above!


##Control
s_c <- params[params$Treatment == "Control", "Sb"]
f_c <- params[params$Treatment == "Control", "Seeds"]

##Comp Rem
s_cr <- params[params$Treatment == "Comp","Sb"]
f_cr <- params[params$Treatment == "Comp","Seeds"]


#Gi=proportion of seeds germinating immediately
#Gsb=proportion of seeds germinating out of seed bank
#v=seed viability
#e=proportion of seeds that don't germinate immediately and enter seedbank
#Ssb=probability that a seed that enters the seed bank will survive to t+1
Gsb <- 0.0025
v <- 0.51111
Ssb <- 0.91
Gi <- 0.0083
e <- 1 - Gi



A_cont <- matrix(c(Ssb * (1 - Gsb), 0, f_c * v * e,
                   Ssb * Gsb, 0, f_c * v * Gi,
                   0, s_c, 0),
                  nrow = 3, byrow = TRUE,
                 dimnames = list(c("SB", "Rosette", "RA"),
                                 c("SB", "Rosette", "RA")))

A_cr <- matrix(c(Ssb * (1 - Gsb), 0, f_cr * v * e,
                 Ssb * Gsb, 0, f_cr * v * Gi,
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

values <- c(s_c, s_cr, 
            f_c, f_cr,
            lambda_cont, lambda_cr)


####################################################
#Bootstrapping code
#I do 1000 bootstrap replicates.  This will save all the data. 
##For now I fill in a column with 1000 rows with the value 'NA'.
##Later, I will replace the NAs with data.

nreps <- 1000
boot_s_c <- rep(NA,nreps)
boot_s_cr <- rep(NA,nreps)
boot_f_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)
boot_l_c <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)

cplants <- filter(ks, Treatment == "Control")
crplants <- filter(ks, Treatment == "Comp")

nc <- length(cplants[ ,1])
ncr <- length(crplants[ ,1])

n <- length(ks[ ,c(1)])
#start the loop
for(j in seq_len(nreps)) {
  
  #x3 is a column that picks n integers between the numbers 1 and n 
  #(because we have n plants in the data set).  
  #replace=TRUE means that we are sampling with replacement 
  #(once we choose a number, we can choose it again).
  
  x1 <- sample(1:nc, nc, replace = TRUE)
  x2 <- sample(1:ncr, ncr, replace = TRUE)
  
  bootc <- cplants[x1, ]
  bootcr <- crplants[x2, ]
  bootdata <- rbind(bootc, bootcr)
  
  params1 <- bootdata %>%
    group_by(Treatment) %>%
    summarise(N = length(Treatment),
              Sb = mean(Alive15, na.rm = TRUE),
              Seeds = mean(Seeds, na.rm = TRUE))
  
  params1 <- as.data.frame(params1)
  
  ##These lines changes all NaN's to observed parameter means in the params data frame
  ##I'm not sure if this is the best way to actually deal with the problem, but it 
  ##Keeps the loop from shutting down if resampling is fucky. 
  
  for(k in 1:length(params1$Treatment)){
    for(l in 1:length(params1)){
      if(is.na(params1[k, l])){
        params1[k, l] <- params[k, l]
      }
    }
  }
  #extract model paramater values from "params1" dataframe we created above!

  ##Control
  s_c <- params1[params1$Treatment == "Control", "Sb"]
  f_c <- params1[params1$Treatment == "Control", "Seeds"]
  boot_s_c[j] <- s_c
  boot_f_c[j] <- f_c
  
  ##Comp Rem
  s_cr <- params1[params1$Treatment == "Comp", "Sb"]
  f_cr <- params1[params1$Treatment == "Comp", "Seeds"]
  boot_s_cr[j] <- s_cr
  boot_f_cr[j] <- f_cr
  
  
  Gi <- 0.0083
  Gsb <- 0.0025
  v <- 0.51111
  Ssb <- 0.91 # Lincoln nebraska site Bernside 1996
  e <- 1-Gi
  
  A_cont <- matrix(c(Ssb * (1 - Gsb), 0, f_c * v * e,
                     Ssb * Gsb, 0, f_c * v * Gi,
                     0, s_c, 0),
                   nrow = 3, byrow = TRUE,
                   dimnames = list(c("SB", "Rosette", "RA"),
                                   c("SB", "Rosette", "RA")))
  
  A_cr <- matrix(c(Ssb * (1 - Gsb), 0, f_cr * v * e,
                   Ssb * Gsb, 0, f_cr * v * Gi,
                   0, s_cr, 0),
                 nrow = 3, byrow = TRUE,
                 dimnames = list(c("SB", "Rosette", "RA"),
                                 c("SB", "Rosette", "RA")))             
             
             ev <- eigen(A_cont)
             
             lmax <- which(Re(ev$values) == max(Re(ev$values)))
             lmax <- lmax[1]  
             lambda_cont <- Re(ev$values[lmax])
             boot_l_c[j] <- lambda_cont
             
             ev <- eigen(A_cr)
             
             lmax <- which(Re(ev$values) == max(Re(ev$values)))
             lmax <- lmax[1]  
             lambda_cr <- Re(ev$values[lmax])
             boot_l_cr[j] <- lambda_cr
             
}


boot_s_c <- sort(boot_s_c)
boot_s_cr <- sort(boot_s_cr)
boot_f_c <- sort(boot_f_c)
boot_f_cr <- sort(boot_f_cr)
boot_l_c <- sort(boot_l_c)
boot_l_cr <- sort(boot_l_cr)

##Extract Confidence Intervals
lower <- c(boot_s_c[25], boot_s_cr[25],
           boot_f_c[25], boot_f_cr[25], 
           boot_l_c[25], boot_l_cr[25])

upper <- c(boot_s_c[975], boot_s_cr[975],
           boot_f_c[975], boot_f_cr[975], 
           boot_l_c[975], boot_l_cr[975])

results <- data.frame(values, lower, upper) 

results$Trtvalue <- c("Control", "CR")
results$Trtvalue <- factor(results$Trtvalue,
                           levels = c("Control","CR"))
S1 <- results[1:2,]
fec <- results[3:4,]
lambda <- results[5:6,]


par(mfrow = c(1, 3), 
    mar = c(5,6,4,2) + 0.2)
plot(as.integer(S1$Trtvalue), S1$values,
     pch = 1, 
     ylim = c(0,1), 
     axes = FALSE,
     main = "",
     xlab = "", 
     ylab = "Survival",
     cex.lab = 1.8)
arrows(as.integer(S1$Trtvalue), S1$lower, 
       as.integer(S1$Trtvalue), S1$upper, 
       length = 0.05, 
       angle = 90,
       code = 3)
mtext(c("Control","CR"),
      side = 1,
      line = 1,
      at = c(1,2), 
      cex = 1.0)
axis(2, cex.axis = 1)
box(lwd = 2)

plot(as.integer(fec$Trtvalue), fec$values,
     pch = 1, 
     ylim = c(0, 85000),
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
      at = c(1,2),
      cex = 1.0)
axis(2, cex.axis = 1)
box(lwd = 2)

plot(as.integer(lambda$Trtvalue), lambda$values,
     pch = 1, 
     ylim = c(0,15), 
     axes = FALSE,
     main = "",
     xlab = "Treatment", 
     ylab = expression(paste('Lambda (', lambda, ')')),
     cex.lab = 1.8)
arrows(as.integer(lambda$Trtvalue), lambda$lower, 
       as.integer(lambda$Trtvalue), lambda$upper, 
       length = 0.05, 
       angle = 90, 
       code = 3)
mtext(c("Control","CR"),
      side = 1,
      line = 1,
      at = c(1,2), 
      cex = 1.0)
axis(2, cex.axis = 1)
box(lwd = 2)
