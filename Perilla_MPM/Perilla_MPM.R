# rm(list = ls()) 

library(dplyr)
library(ggplot2)

# read in data
ks <- read.csv("Perilla_MPM/Perilla_Census_Clean.csv",
               stringsAsFactors = FALSE)

cont <- filter(ks, Treatment == "Control")
cr <- filter(ks, Treatment == "Comp")

#calculate vital rates for each treatment

params <- ks %>%
  group_by(Treatment) %>%
  summarise(Survival = mean(Survival, na.rm = TRUE),
            Seeds = mean(Seeds, na.rm = TRUE)) %>%
  as.data.frame()


##Extract vital rates from dataframe

s_cr <- params[params$Treatment == "Comp", "Survival"]
f_cr <- params[params$Treatment == "Comp", "Seeds"]

s_c <- params[params$Treatment == "Control", "Survival"]
f_c <- params[params$Treatment == "Control", "Seeds"]

# calculate germination rate. See manuscript for details.
GermData <- read.csv('Germination/Clean_Germ.csv', 
                     stringsAsFactors = FALSE) %>%
  filter(Species == 'Perilla')

G0 <- mean(c(0.22, mean(GermData$Surface_Germ_Prop)))

# This species is being modeled as an annual with no seed bank
# so there isn't really a matrix for it. 

A_cr <- f_cr * s_cr * G0
A_c <- f_c * s_c * G0

# The code below is a bit redundant, but i retained the structure from
# other models so it won't look confusing when you see it for other species

eig_c <- eigen(A_c)
eig_cr <- eigen(A_cr)

l_cr <- max(Re(eig_cr$values))

l_c <- max(Re(eig_c$values))

l_cr
l_c


##Store observed values for vital rates for use later on
values <- c(s_cr,s_c,
            f_cr,f_c,
            l_cr,l_c)


# We bootstrap the data set 1000 times, creating a new matrix for each one, 
# then store the bootstrapped vital rates in the vectors created below. 
nreps <- 1000
boot_s_cr <- rep(NA, nreps) 
boot_s_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)
boot_f_c <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)
boot_l_c <- rep(NA, nreps)


# These numbers are used to tell R how many times to sample with replacement for
# each bootstrap simulation we run
nc <- dim(cont)[1]
ncr <- dim(cr)[1]

# start the loop
for(j in seq_len(nreps)) {
  
  # xN is a column that picks n integers between the numbers 1 and n 
  # (because we have n plants in the data set).  
  x1 <- sample(1:nc, nc, replace = TRUE)
  x3 <- sample(1:ncr, ncr, replace = TRUE)
  
  
  # The next four lines bootstrap each treatment separately, then recombines them
  # into a single dataframe
  boot_c <- cont[x1, ]
  boot_cr <- cr[x3, ]
  bootdata <- rbind(boot_c,boot_cr)
  
  # Vital rates out of the bootstrapped data frame, except with
  # a slightly different name

  params1 <- bootdata %>%
    group_by(Treatment) %>%
    summarise(Survival = mean(Survival, na.rm = TRUE),
              Seeds = mean(Seeds, na.rm = TRUE)) %>%
    as.data.frame()
  
  # These next few lines go through the bootstrapped vital rate data frame
  # and remove any NaNs by replacing them with the value we observed in our 
  # field data. Occasionally when R is resampling, it will only choose dead plants. 
  # This prevents that from happening.

  for(k in 1:length(params1$Treatment)){
    for(l in 1:length(params1)){
      if(is.na(params1[k, l])){
        params1[k, l] <- params[k, l]
      }
    }
  }
  
  # New vital rates for the equations below, as well as a few lines to store
  # the bootstrapped vital rate values from each run in the appropriate vector.
  # The vectors are then sorted to obtain confidence intervals.
  
  s_cr <- params1[params1$Treatment == "Comp", "Survival"]
  f_cr <- params1[params1$Treatment == "Comp", "Seeds"]
  boot_s_cr[j] <- s_cr
  boot_f_cr[j] <- f_cr
  
  s_c <- params1[params1$Treatment == "Control", "Survival"]
  f_c <- params1[params1$Treatment == "Control", "Seeds"]
  boot_s_c[j] <- s_c
  boot_f_c[j] <- f_c

  G0 <- mean(c(0.22, mean(GermData$Surface_Germ_Prop)))
  
  # This species is being modeled as an annual with no seed bank
  # so there isn't really a matrix for it. 
  
   A_cr <- f_cr * s_cr * G0
  A_c <- f_c * s_c * G0
  
  
  boot_l_cr[j] <- A_cr

  boot_l_c[j] <- A_c
}

# sorting vectors containing bootstraped vital rates
boot_s_cr <- sort(boot_s_cr)
boot_s_c <- sort(boot_s_c)
boot_f_cr <- sort(boot_f_cr)
boot_f_c <- sort(boot_f_c)
boot_l_cr <- sort(boot_l_cr)
boot_l_c <- sort(boot_l_c)

var_escr <- var(log(boot_l_cr + 0.5) - log(boot_l_c + 0.5))

# creating vector of upper and lower confidence intervals for each vital rate
lower <- c(boot_s_cr[25], boot_s_c[25],
           boot_f_cr[25], boot_f_c[25], 
           boot_l_cr[25], boot_l_c[25])
upper <- c(boot_s_cr[975], boot_s_c[975],
           boot_f_cr[975], boot_f_c[975],
           boot_l_cr[975], boot_l_c[975])


results <- tibble(values, lower, upper) 


results$text <- NA
results$y <- NA

# This determines whether or not there are significant differences in the vital
# rate and places the significance code on the plot
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

# Makes for prettier labels when used w/ facet_wrap()
results$Treatment <- c('Control', 'CR')
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
       aes(x = Treatment)) +
  geom_point(aes(y = values,
                 color = Treatment),
             size = 4.5) + 
  facet_wrap(~Var,
             scales = 'free',
             labeller = label_parsed) + 
  geom_linerange(aes(ymin = lower,
                     ymax = upper,
                     color = Treatment),
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
                     values = c('black', 'green'))  + 
  geom_text(x = 1.5,
            y = results$y,
            label = results$text)

ggsave(filename = 'Perilla_VR_Panel.png',
       path = 'Perilla_MPM/',
       height = 5,
       width = 8,
       unit = 'in')

