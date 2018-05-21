##Clear environment and console
rm(list = ls()) 

library(dplyr)
library(ggplot2)

#import data into R studio
ks <- read.csv("Carduus/CN_Clean.csv") 
ks$Seeds <- ks$Flower_Number * 88.5

# check for density dependence
ddTest <- ks %>% 
  group_by(Treatment, Plot) %>% 
  summarise(N = n(),
            Sb = mean(Survival, na.rm = TRUE),
            Seeds = mean(Seeds, na.rm = TRUE))

params <- ks %>% 
  group_by(Treatment) %>%
  summarise(N = n(),
            Sb = mean(Survival, na.rm = TRUE),
            Seeds = mean(Seeds, na.rm = TRUE))

##Control
s_c <- params[params$Treatment == "Control","Sb"] %>% as.numeric
f_c <- params[params$Treatment == "Control","Seeds"] %>% as.numeric

##Comp Rem
s_cr <- params[params$Treatment == "Competitor","Sb"] %>% as.numeric 
f_cr <- params[params$Treatment == "Competitor","Seeds"] %>% as.numeric


#seed viability and germination rate,
Gsb <- 0.03
Ssb <- 0.2597
e <- 0.2333
Gi <- 0.03

A_cont <- matrix(c(Ssb, 0, f_c * e,
                   Gsb, 0, f_c * Gi,
                   0, s_c, 0),
                 nrow = 3, byrow = TRUE, 
                 dimnames=list(c("SB", "Rosette", "RA"),
                               c("SB", "Rosette", "RA")))

A_cr <- matrix(c(Ssb, 0, f_cr * e,
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

values <- c(s_c, s_cr,
            f_c, f_cr,
            lambda_cont, lambda_cr)


####################################################
#Bootstrapping code
#I do 1000 bootstrap replicates.  This will save all the data. 
##For now I fill in a column with 1000 rows with the value 'NA'.
##Later, I will replace the NAs with data.

nreps <- 1000
boot_s_c <- rep(NA, nreps)
boot_s_cr <- rep(NA, nreps)
boot_f_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)
boot_l_cont <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)

cplants <- filter(ks, Treatment == "Control")
crplants <- filter(ks, Treatment == "Competitor")


nc <- length(cplants[ ,1])
ncr <- length(crplants[ ,1])
#start the loop
for(j in seq_len(nreps)) {
  
  #x3 is a column that picks n integers between the numbers 1 and n (because we have n plants in the data set).  
  #replace=TRUE means that we are sampling with replacement (once we choose a number, we can choose it again).
  
  x2 <- sample(1:nc, nc, replace = TRUE)
  x3 <- sample(1:ncr, ncr, replace = TRUE)
  
  bootc <- cplants[x2, ]
  bootcr <- crplants[x3, ]
  
  bootdata <- rbind(bootc, bootcr)
  
  params1 <- bootdata %>% 
    group_by(Treatment) %>%
    summarise(N = length(Treatment),
              Sb = mean(Survival, na.rm = TRUE),
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
  Gsb <- 0.03
  Ssb <- 0.2597
  e <- 0.2333
  Gi <- 0.03
  
  A_cont <- matrix(c(Ssb, 0, f_c*e,
                     Gsb, 0, f_c*Gi,
                     0, s_c, 0),
                   nrow = 3,
                   byrow = TRUE, 
                   dimnames = list(c("SB","Rosette","RA"),
                                   c("SB","Rosette","RA")))
  
  A_cr <- matrix(c(Ssb, 0, f_cr*e,
                   Gsb, 0, f_cr*Gi,
                   0, s_cr, 0),
                 nrow = 3,
                 byrow = TRUE, 
                 dimnames = list(c("SB","Rosette","RA"),
                                 c("SB","Rosette","RA")))
             
                 
             ev <- eigen(A_cont)
             
             lmax <- which(Re(ev$values) == max(Re(ev$values)))
             lmax <- lmax[1]  
             lambda_cont <- Re(ev$values[lmax])
             boot_l_cont[j] <- lambda_cont
             
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
boot_l_cont <- sort(boot_l_cont)
boot_l_cr <- sort(boot_l_cr)

##Extract Confidence Intervals
lower <- c(boot_s_c[25], boot_s_cr[25],
           boot_f_c[25], boot_f_cr[25], 
           boot_l_cont[25], boot_l_cr[25])

upper <- c(boot_s_c[975], boot_s_cr[975],
           boot_f_c[975], boot_f_cr[975], 
           boot_l_cont[975], boot_l_cr[975])

results <- tibble(values, lower, upper) 


results$text <- NA
results$y <- NA

if(results$upper[1] < results$lower[2] |
   results$upper[2] < results$lower[1]) {
  results$text[1:2] <- '***'
} else {
  results$text[1:2] <- 'NS'
}

results$y[1:2] <- .99 * max(results$upper[1:2])

if(results$upper[3] < results$lower[4] |
   results$upper[4] < results$lower[3]) {
  results$text[3:4] <- '***'
} else {
  results$text[3:4] <- 'NS'
}

results$y[3:4] <- .99 * max(results$upper[3:4])


if(results$upper[5] < results$lower[6] |
   results$upper[6] < results$lower[5]) {
  results$text[5:6] <- '***'
} else {
  results$text[5:6] <- 'NS'
}
results$y[5:6] <- .99 * max(results$upper[5:6])



results$Trt <- c("Control","CR")
results$Var <- factor(c('paste(italic(s[j]))', 
                        'paste(italic(s[j]))',
                        'paste(italic(f))',
                        'paste(italic(f))',
                        'lambda',
                        'lambda'),
                      levels = c('paste(italic(s[j]))',
                                 'paste(italic(f))',
                                 'lambda'),
                      ordered = TRUE)

ggplot(data = results,
       aes(x = Trt)) +
  geom_point(aes(y = values,
                 color = Trt),
             size = 4.5) + 
  facet_wrap(~Var,
             scales = 'free',
             labeller = label_parsed) + 
  geom_linerange(aes(ymin = lower,
                     ymax = upper,
                     color = Trt),
                 size = 1.25) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = 'black',
                                    fill = NA),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 18,
                                    margin = margin(t = 20,
                                                    r = 0,
                                                    l = 0,
                                                    b = 0)),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  scale_x_discrete('Treatment') + 
  scale_y_continuous('') + 
  scale_color_manual('Treatment',
                     breaks = c('Control', 'CR'),
                     values = c('black', 'green')) + 
  geom_text(x = 1.5,
            y = results$y,
            label = results$text)


ggsave(filename = 'Carduus_VR_Panel.png',
       path = 'Carduus',
       height = 5,
       width = 8,
       unit = 'in')
