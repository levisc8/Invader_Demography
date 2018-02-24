rm(list = ls()) 
cat("\014") 
library(dplyr)
library(ggplot2)

ks <- read.csv("Perilla/Perilla4R2013 AG.csv",
               stringsAsFactors = FALSE) %>%
  filter(Treatment != 'Herb')

cont <- filter(ks, Treatment == "Control")
cr <- filter(ks, Treatment == "Comp")

#calculate vital rates for each treatment

params <- ks %>%
  group_by(Treatment) %>%
  summarise(Survival = mean(Alive, na.rm = TRUE),
            Seeds = mean(Seeds, na.rm = TRUE)) %>%
  as.data.frame()


##Extract vital rates from dataframe

s_cr <- params[params$Treatment == "Comp", "Survival"]
f_cr <- params[params$Treatment == "Comp", "Seeds"]

s_c <- params[params$Treatment == "Control", "Survival"]
f_c <- params[params$Treatment == "Control", "Seeds"]


G0 <- mean(c(0.22, 0.489))

##This species is being modeled as an annual with no seed bank
##so there isn't really a matrix for it. 

A_cr <- f_cr * s_cr * G0
A_c <- f_c * s_c * G0

##The code below is a bit redundant, but i retained the structure from
##other models so it won't look confusing when you see it for other species
l_cr <- A_cr

l_c <- A_c

##These lambda values are pretty high,
##but when we thinned the plots down to 15 plants/plot, we removed
##all of the density dependent mortality (which our model can't account
##for), which might mean our estimates of survival are a bit higher than they normally would be.
##It's probably ok, because we are only looking at effect sizes anyway.
l_cr
l_c


##Store observed values for vital rates for use later on
values <- c(s_cr,s_c,
            f_cr,f_c,
            l_cr,l_c)


##We bootstrap the data set 1000 times, running a new matrix for each one, then store
##the bootstrapped vital rates in the vectors created below. 
nreps <- 1000
boot_s_cr <- rep(NA, nreps) 
boot_s_c <- rep(NA, nreps)
boot_f_cr <- rep(NA, nreps)
boot_f_c <- rep(NA, nreps)
boot_l_cr <- rep(NA, nreps)
boot_l_c <- rep(NA, nreps)


##These numbers are used to tell R how many times to sample with replacement for
##each bootstrap simulation we run
nc <- dim(cont)[1]
ncr <- dim(cr)[1]

#start the loop
for(j in seq_len(nreps)) {
  
  #xN is a column that picks n integers between the numbers 1 and n (because we have n plants in the data set).  
  x1 <- sample(1:nc, nc, replace = TRUE)
  x3 <- sample(1:ncr, ncr, replace = TRUE)
  
  
  ##The next four lines bootstrap each treatment separately, then recombines them
  ##into a larger data frame that goes into ddply()
  boot_c <- cont[x1, ]
  boot_cr <- cr[x3, ]
  bootdata <- rbind(boot_c,boot_cr)
  
  # Vital rates out of the bootstrapped data frame, except with
  # a slightly different name
  # You'll see that the new name comes in handy below
  
  params1 <- bootdata %>%
    group_by(Treatment) %>%
    summarise(Survival = mean(Alive, na.rm = TRUE),
              Seeds = mean(Seeds, na.rm = TRUE)) %>%
    as.data.frame()
  
  ##These next few lines go through the bootstrapped vital rate data frame
  ##and remove any NaNs by replacing them with the value we observed in our field data.
  ##This won't alter the confidence intervals, and keeps the loop running smoothly.
  ##Occasionally when R is resampling, it will only choose dead plants, which screws up the loop
  ##causing it to stop and return an error. This prevents that from happening.
  ##I doubt it would be a problem for this species, but it has given me fits with Teucrium
  ##and a few others, so I've retained it here just in case, and to limit any alterations 
  ##to methodology.
  
  for(k in 1:length(params1$Treatment)){
    for(l in 1:length(params1)){
      if(is.na(params1[k, l])){
        params1[k, l] <- params[k, l]
      }
    }
  }
  
  ##New vital rates for the equations below, as well as a few lines to store
  ##the bootstrapped vital rate values from each run in the appropriate vector.
  ##The vectors are then sorted to obtain confidence intervals.
  
  s_cr <- params1[params1$Treatment == "Comp", "Survival"]
  f_cr <- params1[params1$Treatment == "Comp", "Seeds"]
  boot_s_cr[j] <- s_cr
  boot_f_cr[j] <- f_cr
  
  s_c <- params1[params1$Treatment == "Control", "Survival"]
  f_c <- params1[params1$Treatment == "Control", "Seeds"]
  boot_s_c[j] <- s_c
  boot_f_c[j] <- f_c

  G0 <- mean(c(0.22, 0.489))
  
  ##This species is being modeled as an annual with no seed bank
  ##so there isn't really a matrix for it. 
  
   A_cr <- f_cr * s_cr * G0
  A_c <- f_c * s_c * G0
  
  
  boot_l_cr[j] <- A_cr

  boot_l_c[j] <- A_c
}

##sorting vectors containing bootstraped vital rates
boot_s_cr <- sort(boot_s_cr)
boot_s_c <- sort(boot_s_c)
boot_f_cr <- sort(boot_f_cr)
boot_f_c <- sort(boot_f_c)
boot_l_cr <- sort(boot_l_cr)
boot_l_c <- sort(boot_l_c)

##creating vector of upper and lower confidence intervals for each vital rate
lower <- c(boot_s_cr[25], boot_s_c[25],
        boot_f_cr[25], boot_f_c[25], 
        boot_l_cr[25], boot_l_c[25])
upper <- c(boot_s_cr[975], boot_s_c[975],
        boot_f_cr[975], boot_f_c[975],
        boot_l_cr[975], boot_l_c[975])


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

ggsave(filename = 'Perilla_VR_Panel.png',
       path = 'Perilla',
       height = 5,
       width = 8,
       unit = 'in')

