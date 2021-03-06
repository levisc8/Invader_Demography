# Clear environment and console
# rm(list = ls())

library(dplyr)
library(ggplot2)

# Import data
ks <- read.csv("Thlaspi_MPM/Thlaspi_Census_Clean.csv",
               stringsAsFactors = FALSE)

# Calculate parameter values by treatment
params <- ks %>%
  group_by(Treatment) %>%
  summarise(Survival = mean(Survival, na.rm = TRUE),
            Fruit = mean(Fruit, na.rm = TRUE),
            nAlive = sum(Survival, na.rm = TRUE)) %>%
  as.data.frame()

# mean seeds/fruit
F_C <- params$Fruit[params$Treatment == 'Control'] * 3.76 
F_CR <- params$Fruit[params$Treatment == 'Comp'] * 3.76 

S_C <- params$Survival[params$Treatment == 'Control']
S_CR <- params$Survival[params$Treatment == 'Comp']


# All rates below taken from Burns et al 2012/2013
r  <-  0.67 # remain in seedbank
e  <-  0.12 # emerge from seedbank
g  <-  11 / (400 * .95)
v  <-  0.95


# Create matrices and extract dominant eigenvalues

A_c <- matrix(c(r, S_C * F_C * (1 - g) * v,
                e, S_C * F_C * g * v),
              byrow = TRUE, 
              nrow = 2)
 
A_cr <- matrix(c(r, S_CR * F_CR * (1 - g) * v,
                 e, S_CR * F_CR * g * v),
               nrow = 2,
               byrow = TRUE)


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
# Bootstrapping code
nreps <- 1000
boot_s_c <- rep(NA, nreps)
boot_s_cr <- rep(NA, nreps)

boot_f_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)

boot_l_c <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)

# subset plants into treatments
cont <- filter(ks, Treatment == "Control")
cr <- filter(ks, Treatment == "Comp")

# Calculate number of plants per treatment
nc <- dim(cont)[1]
ncr <- dim(cr)[1]

# start the loop
for(j in 1:nreps) {
  
  # Create bootstrapping vectors and subset the data using them
  x1 <- sample(1:nc, nc, replace = TRUE)
  x3 <- sample(1:ncr, ncr, replace = TRUE)

  bootc <- cont[x1, ]
  bootcr <- cr[x3, ]
  bootdata <- rbind(bootc, bootcr)
  
  # re-calculate parameter values with bootstrapped data
  params1 <- bootdata %>% 
    group_by(Treatment) %>%
    summarise(Survival = mean(Survival, na.rm = TRUE),
              Fruit = mean(Fruit, na.rm = TRUE),
              nAlive = sum(Survival, na.rm = TRUE)) %>%
    as.data.frame()
  

  # If any parameters are given NAs, then substitute in the values
  # we observed in the field. 
  
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
  
  boot_s_c[j] <- S_C
  boot_s_cr[j] <- S_CR
  
  boot_f_c[j] <- F_C
  boot_f_cr[j] <- F_CR
  
  r <- 0.67 # remain in seedbank
  e <- 0.12 # emerge from seedbank
  g <- 11 / (400 * .95)
  v <- 0.95
  
  # Bootstrapped matrices
  
  A_c <- matrix(c(r, S_C * F_C * (1 - g) * v,
                  e, S_C * F_C * g * v),
                byrow = TRUE,
                nrow = 2)
  
  A_cr <- matrix(c(r, S_CR * F_CR * (1 - g) * v,
                   e, S_CR * F_CR * g * v),
                 nrow = 2,
                 byrow = TRUE)
  
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

# Extract confidence intervals for each parameter
boot_s_c <- sort(boot_s_c)
boot_s_cr <- sort(boot_s_cr)

boot_f_c <- sort(boot_f_c)
boot_f_cr <- sort(boot_f_cr)

var_escr <- var(log(boot_l_cr + 0.5) - log(boot_l_c + 0.5))

lambdas <- data.frame(lambda_c = c(values[5], boot_l_c),
                      lambda_cr = c(values[6], boot_l_cr),
                      boot_obs = c("observed", rep("bootstrap", 1000)))

saveRDS(lambdas, file = '../Data/bootstrap_lambdas/Thlaspi_lambdas.rds')

boot_l_c <- sort(boot_l_c)
boot_l_cr <- sort(boot_l_cr)

lower <- c(boot_s_c[25], boot_s_cr[25],
           boot_f_c[25], boot_f_cr[25],
           boot_l_c[25], boot_l_cr[25])
upper <- c(boot_s_c[975], boot_s_cr[975],
           boot_f_c[975], boot_f_cr[975],
           boot_l_c[975], boot_l_cr[975])

results <- tibble(values, lower, upper) 

results$text <- NA
results$y <- NA

# Determine significance and place code on plot
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

# make pretty labels for facet_wrap()
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

# insert horizontal line at lambda = 1
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

ggsave(filename = 'Thlaspi_VR_Panel.png',
       path = 'Thlaspi_MPM/',
       height = 5,
       width = 8,
       unit = 'in')
