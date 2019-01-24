
# rm(list = ls()) 

#import data into R studio
library(dplyr)
library(ggplot2)

ks <- read.csv('Alliaria_MPM/Alliaria_Census_Clean.csv',
                stringsAsFactors = FALSE)

# Calculate vital rates by treatment
paramsD <- ks %>%
  group_by(Treatment) %>%
  summarise(plot.n = length(unique(Plot)),
            N.SDL = n(),
            N.RA = sum(RA, na.rm = TRUE),
            s2 = N.RA/N.SDL)

Fec <- ks %>% 
  group_by(Treatment) %>% 
  summarise(Fec = mean(Seeds, na.rm = TRUE))


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

# The matrix  and eigenvalues
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

values <- c(s2_c,s2_cr,
            f_c, f_cr,
            lambda_cont,lambda_cr)


####################################################
#Bootstrapping code

# Set up vectors to hold boot strapped vital rates

nreps <- 1000

boot_s2_c <- rep(NA, nreps)
boot_s2_cr <- rep(NA, nreps)
boot_f_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)
boot_l_c <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)

# Split data into distinct data frames
C.RAs <- filter(ks, Treatment == 'Control' & RA == 1)
CR.RAs <- filter(ks, Treatment == 'Comp' & RA == 1)

CPlants <- filter(ks, Treatment == 'Control')
CRPlants <- filter(ks, Treatment == 'Comp')

# Compute number of plants for resampling (see next comment)
n.c.ra <- dim(C.RAs)[1]
n.cr.ra <- dim(CR.RAs)[1]
n.c.plants <- dim(CPlants)[1]
n.cr.plants <- dim(CRPlants)[1]


#start the loop
for(j in 1:nreps) {
  
  # xTrt is a vector that picks n integers between the numbers 1 and n (because we 
  # have n plants in the data set).  
  # replace=TRUE means that we are sampling with replacement (once we choose a 
  # number, we can choose it again).
  xCP <- sample(1:n.c.plants, 
                n.c.plants,
                replace = TRUE)
  xCRP <- sample(1:n.cr.plants, 
                 n.cr.plants,
                 replace = TRUE)
  xCRA <- sample(1:n.c.ra, 
                 n.c.ra,
                 replace = TRUE)
  xCRRA <- sample(1:n.cr.ra, 
                  n.cr.ra,
                  replace = TRUE)
  
  bootSurvData <- rbind(CPlants[xCP, ],
                        CRPlants[xCRP, ])
  bootFecData <- rbind(C.RAs[xCRA, ],
                       CR.RAs[xCRRA, ])
  
  # boot strapped vital rates by treatment
  paramsD1 <- bootSurvData %>% 
    group_by(Treatment) %>% 
    summarise(N.SDL = n(),
              N.RA = sum(RA, na.rm = TRUE),
              s2 = N.RA/N.SDL)
  
  Fec1 <- bootFecData %>% 
    group_by(Treatment) %>% 
    summarise(Fec = mean(Seeds,
                         na.rm = TRUE))
  
  paramsD1[paramsD1$Treatment == 'Comp','Fec'] <- 
    Fec1[Fec1$Treatment == 'Comp','Fec']
  
  paramsD1[paramsD1$Treatment == 'Control', 'Fec'] <- 
    Fec1[Fec1$Treatment == 'Control', 'Fec']
  
  # replace with observed means if no vital rate can be calculated from bootstrap
  # sample. This occurs if the boot strapping sampler happens to, for example, 
  # only pick dead RAs leading to NaNs in the fecundity values
  
  for(i in 1:dim(paramsD1)[1]){
    for(x in 1:dim(paramsD1)[2]){
      if(is.na(paramsD1[i, x])){
        paramsD1[i, x] <- paramsD[i, x]
      }
    }
  }
  

  ##Control
  s1_c <- 1
  s2_c <- paramsD1[paramsD1$Treatment == "Control", "s2"] %>% as.numeric
  f_c <- paramsD1[paramsD1$Treatment == "Control", "Fec"] %>% as.numeric
  boot_s2_c[j] <- s2_c
  boot_f_c[j] <- f_c
  
  ##Comp Rem
  s1_cr <- 1
  s2_cr <- paramsD1[paramsD1$Treatment == "Comp", "s2"] %>% as.numeric
  f_cr <- paramsD1[paramsD1$Treatment == "Comp", "Fec"] %>% as.numeric
  boot_s2_cr[j] <- s2_cr
  boot_f_cr[j] <- f_cr
  
  
  ##Bootstrapped matrices!
  G1 <- 0.5503
  G2 <- 0.3171
  v <- 0.8228
  
  A_cont <- matrix(c(1-G2, 0, f_c*v*(1-G1),
                     G2*s1_c, 0, f_c*G1*s1_c,
                     0, s2_c, 0),
                nrow = 3,
                byrow = TRUE, 
                dimnames = list(c("SB","Rosette","RA"),
                                c("SB","Rosette","RA")))
  
  A_cr <- matrix(c(1-G2, 0, f_cr*v*(1-G1),
                   G2*s1_cr, 0, f_cr*G1*s1_cr,
                   0, s2_cr, 0),
              nrow = 3,
              byrow = TRUE, 
              dimnames = list(c("SB","Rosette","RA"),
                              c("SB","Rosette","RA")))

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


boot_s2_c <- sort(boot_s2_c)
boot_s2_cr <- sort(boot_s2_cr)

boot_f_c <- sort(boot_f_c)
boot_f_cr <- sort(boot_f_cr)

boot_l_c <- sort(boot_l_c)
boot_l_cr <- sort(boot_l_cr)

##Extract Confidence Intervals
lower <- c(boot_s2_c[25], boot_s2_cr[25],
           boot_f_c[25], boot_f_cr[25], 
           boot_l_c[25], boot_l_cr[25])

upper <- c(boot_s2_c[975], boot_s2_cr[975],
           boot_f_c[975], boot_f_cr[975], 
           boot_l_c[975], boot_l_cr[975])

results <- tibble(values, lower, upper) 

results$text <- NA
results$y <- NA

# Insert significance codes. This is a bit convoluted, but it works
if(results$upper[1] < results$lower[2] |
   results$upper[2] < results$lower[1]) {
  results$text[1:2] <- '***'
} else {
  results$text[1:2] <- 'NS'
}

# This looks weird, but it just calculates where to put the significance code
# on the plot's y-axis. In this case, I put them along the top of each plot 
# just inside the plot margin.

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

# Create pretty plot label expressions
results$Trt <- c('Control', 'CR')
results$Var <- factor(c('paste(italic(s))', 
                        'paste(italic(s))',
                        'paste(italic(f))',
                        'paste(italic(f))',
                        'lambda',
                        'lambda'),
                      levels = c('paste(italic(s))',
                                 'paste(italic(f))',
                                 'lambda'),
                      ordered = TRUE)

# horizontal line at lambda = 1
results$yints <- NA_integer_
results$yints[results$Var == 'lambda'] <- 1

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
  geom_hline(data = results,
             aes(yintercept = yints),
             color = 'grey50',
             alpha = 0.7,
             linetype = 'dashed',
             size = 2) + 
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

ggsave(filename = 'Alliaria_VR_Panel.png',
       path = 'Alliaria_MPM/',
       height = 5,
       width = 8,
       unit = 'in')
