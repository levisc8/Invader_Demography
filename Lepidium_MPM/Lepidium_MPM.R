# Clear environment and console
# rm(list = ls()) 

# Import data and add seed counts. Lepidium always makes 2 seeds/fruit
library(dplyr)
library(ggplot2)
ks <- read.csv("Lepidium_MPM/Lepidium_Census_Clean.csv") %>%
  mutate(Seeds = Fruit * 2) # 2 seeds/fruit 

# Compute other parameter values by treatment
params <- ks %>%
  group_by(Treatment) %>%
  summarise(Survival = mean(Survival, na.rm = TRUE),
            Seeds = mean(Seeds, na.rm = TRUE)) %>%
  as.data.frame()

s_c <- params[params$Treatment == "Control","Survival"]
s_cr <- params[params$Treatment == "Comp","Survival"]

f_c <- params[params$Treatment == "Control","Seeds"]
f_cr <- params[params$Treatment == "Comp","Seeds"]

G0 <- (mean(c(.017, .049, .068)) + 0.02)/2 #proportion of viable seeds that germinate immediately

# The "matrix" for this model isn't really a matrix, it's just an equation. 
# I've used the same set up as the other models for consistency's sake.
# The dominant eigenvalue extracted below should be identical to the value
# of each "matrix"

A_c <- s_c*f_c*G0

A_cr <- s_cr*f_cr*G0


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

values <- c(s_c,s_cr,
            f_c, f_cr,
            lambda_c, lambda_cr)


####################################################
# Bootstrapping code
# Set up vectors to store bootstrapped values
nreps <- 1000
boot_s_c <- rep(NA, nreps)
boot_s_cr <- rep(NA, nreps)
boot_f_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)
boot_l_c <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)

#subset plants into treatments and compute how many are in each one
cont <- subset(ks, Treatment == "Control")
cr <- subset(ks, Treatment == "Comp")

nc <- dim(cont)[1]
ncr <- dim(cr)[1]

#start the loop
for(j in seq_len(nreps)) {
  
  # x1 and x2 are integer vectors used to sample with replacement.
  x1 <- sample(1:nc, nc, replace = TRUE)
  x2 <- sample(1:ncr, ncr, replace = TRUE)

  bootc <- cont[x1, ]
  bootcr <- cr[x2, ]
  bootdata <- rbind(bootc,bootcr)
  
  # Compute vital rates based on bootstrapped data. If a value is missing,
  # then set it to the observed value of the vital rate in question.
  params1 <- bootdata %>%
    group_by(Treatment) %>%
    summarise(Survival = mean(Survival, na.rm = TRUE),
              Seeds = mean(Seeds, na.rm = TRUE)) %>%
    as.data.frame()

  for(k in 1:length(params1$Treatment)){
    for(l in 1:length(params1)){
      if(is.na(params1[k, l]) == TRUE){
        params1[k, l] <- params[k, l]
      }
    }
  }
  s_c <- params1[params1$Treatment == "Control","Survival"]
  s_cr <- params1[params1$Treatment == "Comp","Survival"]
  
  f_c <- params1[params1$Treatment == "Control","Seeds"]
  f_cr <- params1[params1$Treatment == "Comp","Seeds"]
  
  boot_s_c[j] <- s_c
  boot_s_cr[j] <- s_cr
  boot_f_c[j] <- f_c
  boot_f_cr[j] <- f_cr
  
  G0 <- (mean(c(.017, .049, .068)) + 0.02)/2 #proportion of viable seeds that germinate immediately
  
  

  A_c <- s_c*f_c*G0
  
  A_cr <- s_cr*f_cr*G0
  
  ev <- eigen(A_c)
  
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lmax <- lmax[1]  
  lambda_c <- Re(ev$values[lmax])
  boot_l_c[j] <- lambda_c
  
  ev <- eigen(A_cr)
  
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lmax <- lmax[1]  
  lambda_cr <- Re(ev$values[lmax])
  boot_l_cr[j] <- lambda_cr
  
}

# Extract the confidence intervals for each parameter.
boot_s_c <- sort(boot_s_c)
boot_s_cr <- sort(boot_s_cr)
boot_f_c <- sort(boot_f_c)
boot_f_cr <- sort(boot_f_cr)
boot_l_c <- sort(boot_l_c)
boot_l_cr <- sort(boot_l_cr)

var_escr <- var(log(boot_l_cr + 0.5) - log(boot_l_c + 0.5))


lower <- c(boot_s_c[25],boot_s_cr[25],
           boot_f_c[25],boot_f_cr[25], 
           boot_l_c[25],boot_l_cr[25])
upper <- c(boot_s_c[975],boot_s_cr[975],
           boot_f_c[975],boot_f_cr[975], 
           boot_l_c[975],boot_l_cr[975])

results <- tibble(values, lower, upper) 

results$text <- NA
results$y <- NA

# Determine significance code and where to place it on each figure
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


# Ugly expressions make pretty labels! For use in labeller = label_parsed in 
# facet_wrap(). Ordered = TRUE ensures that the plots are drawn in the correct
# order (e.g. survival, fecundity, lambda from left to right)
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


ggsave(filename = 'Lepidium_VR_Panel.png',
       path = 'Lepidium_MPM/',
       height = 5,
       width = 8,
       unit = 'in')
