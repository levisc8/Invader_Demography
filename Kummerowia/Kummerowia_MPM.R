##K striata model with 2015 survival and 2014 fecundity
#import data into R studio
rm(list = ls()) 
cat("\014")

library(dplyr)
library(ggplot2)

#import data into R studio
ks15 <- read.csv("Kummerowia/kstra2015.2 AG.csv") %>%
  filter(Treatment != 'Herb')
ks14 <- read.csv("Kummerowia/kummerowia AG.csv") %>%
  filter(Treatment != 'Herb')

survival <- ks15 %>%
  group_by(Treatment) %>%
  summarise(Survival = mean(Alive, na.rm = TRUE))

params <- ks14 %>%
  group_by(Treatment) %>%
  summarise(Seeds = mean(f, na.rm=TRUE)) %>%
  full_join(survival, ., by = 'Treatment') %>%
  as.data.frame()

s_c <- params[params$Treatment=="Control", "Survival"]
f_c <- params[params$Treatment=="Control", "Seeds"]

s_cr <- params[params$Treatment == "Comp", "Survival"]
f_cr <- params[params$Treatment == "Comp", "Seeds"]

V <- 0.6622*0.92 #seed viability
G0 <- 0.92 #proportion of viable seeds that germinate immediately (mean of above and buried seeds, 2014/15 germ exp)
G1 <- 0.04 #proportion of viable seeds that germinate after one year in Seed bank
G2 <- 0.04 #proportion of viable seeds that germinate after two years in Seed bank

A_cont <- matrix(c(0, 0, (f_c*V*G2),
                   1, 0, (f_c*V*G1),
                   0, 1, (s_c*f_c*V*G0)),
                 nrow = 3,
                 byrow = TRUE, 
                 dimnames = list(c("seedbank2", "seedbank1", "plant"),
                                 c("seedbank2", "seedbank1", "plant")))
A_cr <- matrix(c(0, 0, (f_cr*V*G2),
              1, 0, (f_cr*V*G1),
              0, 1, (s_cr*f_cr*V*G0)),
            nrow = 3,
            byrow = TRUE, 
            dimnames = list(c("seedbank2", "seedbank1", "plant"),
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
nreps <- 1000
boot_s_c <- rep(NA, nreps)
boot_s_cr <- rep(NA, nreps)
boot_f_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)
boot_l_c <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)
es_crcontrol <- rep(NA, nreps)

cr15 <- filter(ks15, Treatment == "Comp")
cr14 <- filter(ks14, Treatment == "Comp")
cont15 <- filter(ks15, Treatment == "Control")
cont14 <- filter(ks14, Treatment == "Control")
n2 <- dim(cr15)[1]
n3 <- dim(cont15)[1]
n5 <- dim(cr14)[1]
n6 <- dim(cont14)[1]

for(j in 1:nreps) 
{
  
  #x3 is a column that picks n integers between the numbers 1 and n (because we have n plants in the data set).  
  #replace=TRUE means that we are sampling with replacement (once we choose a number, we can choose it again).
  
  x2 <- sample(1:n2,  n2, replace = TRUE)
  x3 <- sample(1:n3, n3, replace = TRUE)
  x5 <- sample(1:n5, n5, replace = TRUE)
  x6 <- sample(1:n6, n6, replace = TRUE)
  
  bootcomp15 <- cr15[x2, ]
  bootc15 <- cont15[x3, ]
  bootcomp14 <- cr14[x5, ]
  bootcont14 <- cont14[x6, ]
  bootdata15 <- rbind(bootcomp15, bootc15)
  bootdata14 <- rbind(bootcomp14, bootcont14)
  
  Survival1 <- bootdata15 %>%
    group_by(Treatment) %>%
    summarise(Survival = mean(Alive, na.rm = TRUE))

  params1 <- bootdata14 %>%
    group_by(Treatment) %>%
    summarise(Seeds = mean(f, na.rm = TRUE)) %>%
    full_join(Survival1, ., by = 'Treatment') %>%
    as.data.frame()
  
  for(k in 1:length(params$Treatment)){
    for(l in 1:length(params)){
      if(is.na(params1[k, l]) == TRUE){
        params1[k,l] <- params[k, l]
      }
    }
  }

  s_c <- params1[params1$Treatment == "Control","Survival"]
  f_c <- params1[params1$Treatment == "Control","Seeds"]
  boot_s_c[j] <- s_c
  boot_f_c[j] <- f_c
  
  s_cr <- params1[params1$Treatment == "Comp","Survival"]
  f_cr <- params1[params1$Treatment == "Comp","Seeds"]
  boot_s_cr[j] <- s_cr
  boot_f_cr[j] <- f_cr
  
  V <- 0.6622*0.92 #seed viability
  G0 <- 0.92 #proportion of viable seeds that germinate immediately (mean of above and buried seeds, 2014/15 germ exp)
  G1 <- 0.04 #proportion of viable seeds that germinate after one year in Seed bank
  G2 <- 0.04 #proportion of viable seeds that germinate after two years in Seed bank
  
  A_cont <- matrix(c(0, 0, (f_c*V*G2),
                     1, 0, (f_c*V*G1),
                     0, 1, (s_c*f_c*V*G0)),
                   nrow = 3,
                   byrow = TRUE,
                   dimnames = list(c("seedbank2", "seedbank1", "plant"),
                                   c("seedbank2", "seedbank1", "plant")))
  A_cr <- matrix(c(0, 0, (f_cr*V*G2),
                   1, 0, (f_cr*V*G1),
                   0, 1, (s_cr*f_cr*V*G0)),
                 nrow = 3,
                 byrow = TRUE,
                 dimnames = list(c("seedbank2", "seedbank1", "plant"),
                                 c("seedbank2", "seedbank1", "plant")))

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

es_crcontrol2 <- log(boot_l_cr+.5)-log(boot_l_c+.5)
escr2 <- sort(es_crcontrol2)


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

results <- tibble(values, lower, upper) 

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
                     values = c('black', 'green'))

ggsave(filename = 'Kummerowia_VR_Panel.png',
       path = 'Kummerowia',
       height = 5,
       width = 8,
       unit = 'in')

