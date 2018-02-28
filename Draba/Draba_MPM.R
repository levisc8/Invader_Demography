##Clear environment
rm(list = ls()) 
##Import data
ks <- read.csv("Draba/DV_Clean.csv")

ks$Seeds <- ks$Fruits * 17.2 

#subset plants into treatments

c <- subset(ks, Treatment == "Control")
cr <- subset(ks, Treatment == "Comp")


library(dplyr)
library(ggplot2)
params <-  ks %>%
  group_by(Treatment) %>%
  summarise(Survival = mean(Survival, na.rm = TRUE),
            Seeds = mean(Seeds, na.rm = TRUE))

params <- as.data.frame(params)

names(params) <- c("Treatment","Survival","Seeds")

s_c <- params[params$Treatment == "Control", "Survival"]
s_cr <- params[params$Treatment == "Comp", "Survival"]

f_c <- params[params$Treatment == "Control","Seeds"]
f_cr <- params[params$Treatment == "Comp","Seeds"]


v0 <- 0.428 # max of symonides 84b
v1 <- 0.18
v2 <- 0.05
v3 <- 0.02
v4 <- 0.005
G <- 0.057/v0 #proportion of viable seeds that germinate immediately

A_cont <- matrix(c(0, v3 * (1-G), 0, 0, 0,
                   0, 0, v2 * (1-G), 0, 0,
                   0, 0, 0, v1 * (1-G), 0,
                   0, 0, 0, 0, s_c * f_c * v0 * (1-G),  
                   v4 * G, v3 * G, v2 * G, v1 * G, s_c * f_c * v0 * G),
                 nrow = 5, 
                 ncol = 5, 
                 byrow = TRUE,
                 dimnames = list(c('SB4', 'SB3', 'SB2', 'SB1', "P"),
                                 c('SB4', 'SB3', 'SB2', 'SB1', "P")))
A_cr <- matrix(c(0, v3 * (1-G), 0, 0, 0,
                 0, 0, v2 * (1-G), 0, 0,
                 0, 0, 0, v1 * (1-G), 0,
                 0, 0, 0, 0, s_cr * f_cr * v0 * (1-G),  
                 v4 * G, v3 * G, v2 * G, v1 * G, s_cr * f_cr * v0 * G),
               nrow = 5,
               ncol = 5,
               byrow = TRUE,
               dimnames = list(c('SB4', 'SB3', 'SB2', 'SB1', "P"),
                               c('SB4', 'SB3', 'SB2', 'SB1', "P")))

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

values <- c(s_c, s_cr,
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

nc <- length(c[ ,1])
ncr <- length(cr[ ,1])


#start the loop
for(i in 1:nreps) {
  
  x1 <- sample(1:nc, nc, replace = TRUE)
  x2 <- sample(1:ncr, ncr, replace = TRUE)
 
  bootc <- c[x1, ]
  bootcr <- cr[x2, ]
  bootdata <- rbind(bootc,bootcr)
  
  params1 <- bootdata %>%
    group_by(Treatment) %>%
    summarise(Survival = mean(Survival, na.rm = TRUE),
              Seeds = mean(Seeds, na.rm = TRUE))
  
  params1 <- as.data.frame(params1)
  names(params1) <- c("Treatment","Survival","Seeds")

  for(k in 1:length(params1$Treatment)){
    for(j in 1:length(params1)){
      if(is.na(params1[k, j])){
        params1[k, j] <- params[k, j]
       }
    }
  }  

                
  s_c <- params1[params1$Treatment == "Control", "Survival"]
  s_cr <- params1[params1$Treatment == "Comp", "Survival"]
  
  f_c <- params1[params1$Treatment == "Control","Seeds"]
  f_cr <- params1[params1$Treatment == "Comp","Seeds"]
  
  boot_s_c[i] <- s_c
  boot_s_cr[i] <- s_cr
 
  boot_f_c[i] <- f_c
  boot_f_cr[i] <- f_cr
  
  v <- 0.613 #seed viability
  G0 <- 0.057 #proportion of viable seeds that germinate immediately
  
  A_cont <- matrix(c(0, v3 * (1-G), 0, 0, 0,
                     0, 0, v2 * (1-G), 0, 0,
                     0, 0, 0, v1 * (1-G), 0,
                     0, 0, 0, 0, s_c * f_c * v0 * (1-G),  
                     v4 * G, v3 * G, v2 * G, v1 * G, s_c * f_c * v0 * G),
                   nrow = 5, 
                   ncol = 5,
                   byrow = TRUE,
                   dimnames = list(c('SB4', 'SB3', 'SB2', 'SB1', "P"),
                                   c('SB4', 'SB3', 'SB2', 'SB1', "P")))
  A_cr <- matrix(c(0, v3 * (1-G), 0, 0, 0,
                   0, 0, v2 * (1-G), 0, 0,
                   0, 0, 0, v1 * (1-G), 0,
                   0, 0, 0, 0, s_cr * f_cr * v0 * (1-G),  
                   v4 * G, v3 * G, v2 * G, v1 * G, s_cr * f_cr * v0 * G),
                 nrow = 5, 
                 ncol = 5, 
                 byrow = TRUE,
                 dimnames = list(c('SB4', 'SB3', 'SB2', 'SB1', "P"),
                                 c('SB4', 'SB3', 'SB2', 'SB1', "P")))

  
  ev <- eigen(A_cont)
  
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lmax <- lmax[1]  
  lambda_c <- Re(ev$values[lmax])
  boot_l_c[i] <-lambda_c
  
  ev <- eigen(A_cr)
  
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lmax <- lmax[1]  
  lambda_cr <- Re(ev$values[lmax])
  boot_l_cr[i] <- lambda_cr
 
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

results <- tibble(values, lower, upper) 

##Use the code below to graph results for surival, fecundity, and lambda
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

ggsave(filename = 'Draba_VR_Panel.png',
       path = 'Draba',
       height = 5,
       width = 8,
       unit = 'in')
