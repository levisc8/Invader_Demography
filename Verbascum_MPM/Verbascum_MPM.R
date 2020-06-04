##Clear console and environment
# rm(list = ls()) 

library(dplyr)
library(ggplot2)

# import data into R studio
ks <- read.csv("Verbascum_MPM/Verbascum_Census_Clean.csv",
            stringsAsFactors = FALSE) 

# Extract vital rate parameters
params <- ks %>%
  group_by(Treatment) %>%
  summarise(N = length(Treatment),
            Sb = mean(Survival, na.rm = TRUE),
            Seeds = mean(Seeds, na.rm = TRUE)) %>%
  as.data.frame()

# extract model paramater values from "params" dataframe we created above

# Control
s_c <- params[params$Treatment == "Control", "Sb"]
f_c <- params[params$Treatment == "Control", "Seeds"]

# Comp Rem
s_cr <- params[params$Treatment == "Comp","Sb"]
f_cr <- params[params$Treatment == "Comp","Seeds"]


# Gi = proportion of seeds germinating immediately
# Gsb = proportion of seeds germinating out of seed bank
# v = seed viability
# e = proportion of seeds that don't germinate immediately 
# and enter seedbank
# Ssb = probability that a seed that enters the seed bank 
# will survive to t+1

GermData <- read.csv('Germination/Clean_Germ.csv',
                 stringsAsFactors = FALSE) %>%
  filter(Species == 'Verbascum')


Gsb <- 0.0025
v <- mean(GermData$Surface_Germ_Prop)
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
# Bootstrapping code
# This will save all the data. 

nreps <- 1000
boot_s_c <- rep(NA,nreps)
boot_s_cr <- rep(NA,nreps)
boot_f_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)
boot_l_c <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)

# Subset plants into treatments and calculate numbers/treatment
cplants <- filter(ks, Treatment == "Control")
crplants <- filter(ks, Treatment == "Comp")

nc <- length(cplants[ ,1])
ncr <- length(crplants[ ,1])

# start the loop
for(j in seq_len(nreps)) {
  
  # x3 is a column that picks n integers between the numbers 1 and n 
  # (because we have n plants in the data set).  
  # replace = TRUE means that we are sampling with replacement 
  # (once we choose a number, we can choose it again).
  
  x1 <- sample(1:nc, nc, replace = TRUE)
  x2 <- sample(1:ncr, ncr, replace = TRUE)
  
  bootc <- cplants[x1, ]
  bootcr <- crplants[x2, ]
  bootdata <- rbind(bootc, bootcr)
  
  # recalculate parameter values w/ bootstrapped data
  params1 <- bootdata %>%
    group_by(Treatment) %>%
    summarise(N = length(Treatment),
              Sb = mean(Survival, na.rm = TRUE),
              Seeds = mean(Seeds, na.rm = TRUE))
  
  params1 <- as.data.frame(params1)
  
  # These lines changes all NaN's to observed parameter means 
  # in the params data frame
  
  for(k in 1:length(params1$Treatment)){
    for(l in 1:length(params1)){
      if(is.na(params1[k, l])){
        params1[k, l] <- params[k, l]
      }
    }
  }
  
  # Control
  s_c <- params1[params1$Treatment == "Control", "Sb"]
  f_c <- params1[params1$Treatment == "Control", "Seeds"]
  boot_s_c[j] <- s_c
  boot_f_c[j] <- f_c
  
  # Comp Rem
  s_cr <- params1[params1$Treatment == "Comp", "Sb"]
  f_cr <- params1[params1$Treatment == "Comp", "Seeds"]
  boot_s_cr[j] <- s_cr
  boot_f_cr[j] <- f_cr
  
  
  Gi <- 0.0083
  Gsb <- 0.0025
  v <- mean(GermData$Surface_Germ_Prop)
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

var_escr <- var(log(boot_l_cr + 0.5) - log(boot_l_c + 0.5))

lambdas <- data.frame(lambda_c = c(values[5], boot_l_c),
                      lambda_cr = c(values[6], boot_l_cr),
                      boot_obs = c("observed", rep("bootstrap", 1000)))

saveRDS(lambdas, file = '../Data/bootstrap_lambdas/Verbascum_lambdas.rds')


boot_l_c <- sort(boot_l_c)
boot_l_cr <- sort(boot_l_cr)

# Extract Confidence Intervals
lower <- c(boot_s_c[25], boot_s_cr[25],
           boot_f_c[25], boot_f_cr[25], 
           boot_l_c[25], boot_l_cr[25])

upper <- c(boot_s_c[975], boot_s_cr[975],
           boot_f_c[975], boot_f_cr[975], 
           boot_l_c[975], boot_l_cr[975])

results <- tibble(values, lower, upper) 

results$text <- NA
results$y <- NA

# Add in significance codes

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

# Make pretty labels for facet_wrap()
results$Trt <- rep(c('Control', 'CR'), 3)
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

# Add in horizontal line at lambda = 1
results$yints <- NA_integer_
results$yints[results$Var == 'lambda'] <- 1

# plot the vital rates
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



ggsave(filename = 'Verbascum_VR_Panel.png',
       path = 'Verbascum_MPM/',
       height = 5,
       width = 8,
       unit = 'in')
