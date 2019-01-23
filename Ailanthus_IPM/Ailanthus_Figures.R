library(dplyr)
library(ggplot2)
library(tidyr)

# source('Ailanthus_IPM/Ailanthus_IPM_CR_No_Burn.R')
# source('Ailanthus_IPM/Ailanthus_IPM_Control_No_Burn.R')

CROutput <- readRDS('Ailanthus_IPM/CR_Bootstrap_Output.rds') %>%
  as_tibble() %>%
  select(-growth.int, -prob.repro.slope, -prob.repro.int) # remove parameters which did not vary

CRCleaned <- CROutput %>%
  .[-1, ] %>% # remove top row to get CI's out
  gather(key = 'Variable', value = 'Value') %>%
  group_by(Variable) %>%
  arrange(desc(Value)) %>%
  summarise(Obs = NA_real_,
            UpCI = Value[25],
            LoCI = Value[975]) %>%
  mutate(Treatment = 'CR')

for(i in seq_len(dim(CROutput)[2])) {
  var <- names(CROutput)[i]
  
  CRCleaned[CRCleaned$Variable == var, 'Obs'] <- CROutput[1, i]
}



ContOutput <- readRDS('Ailanthus_IPM/Cont_Bootstrap_Output.rds')%>%
  as_tibble() %>%
  select(-growth.int, -prob.repro.slope, -prob.repro.int) # remove parameters which did not vary


ContCleaned <- ContOutput %>%
  .[-1, ] %>% # remove top row to get CI's out
  gather(key = 'Variable', value = 'Value') %>%
  group_by(Variable) %>%
  arrange(desc(Value)) %>%
  summarise(Obs = NA_real_,
            UpCI = Value[25],
            LoCI = Value[975]) %>% 
  mutate(Treatment = 'Control')

for(i in seq_len(dim(ContOutput)[2])) {
  var <- names(ContOutput)[i]
  ContCleaned[ContCleaned$Variable == var, 'Obs'] <- ContOutput[1, i]
}

AllParams <- rbind(ContCleaned, CRCleaned)

write.csv(AllParams, file = 'Ailanthus_IPM/Ailanthus_Summarized_Output.csv',
          row.names = FALSE)

# Now, some ugly brute force renaming of variables so plot labels look pretty

AllParams$Variable[AllParams$Variable == 'clonal.prob'] <- "paste(h[d])"
AllParams$Variable[AllParams$Variable == 'clonal.size.mean'] <- "paste(mu[c], ' Clone Size Mean')"
AllParams$Variable[AllParams$Variable == 'clonal.size.sd'] <- "paste(sigma[c], ' Clone Size SD')"
AllParams$Variable[AllParams$Variable == 'establishment.prob'] <- "paste('Establishment Probability')"
AllParams$Variable[AllParams$Variable == 'germ.seedbank.prob'] <- "paste(italic(g[s]))"
AllParams$Variable[AllParams$Variable == 'go.seedbank.prob'] <- "paste(italic(g[b]))"
AllParams$Variable[AllParams$Variable == 'growth.sd'] <- "paste(sigma[g], ' Growth SD')"
AllParams$Variable[AllParams$Variable == 'growth.slope'] <- "paste(italic(g(y,x)), ' Slope')"
AllParams$Variable[AllParams$Variable == 'recruit.size.mean'] <- "paste(mu[r], ' Recruit Size Mean')"
AllParams$Variable[AllParams$Variable == 'recruit.size.sd'] <- "paste(sigma[r], ' Recruit Size SD')"
AllParams$Variable[AllParams$Variable == 'seed.int'] <- "paste(italic(f[s](x)), ' Intercept')"
AllParams$Variable[AllParams$Variable == 'seed.slope'] <- "paste(italic(f[s](x)), ' Slope')"
AllParams$Variable[AllParams$Variable == 'stay.seedbank.prob'] <- "paste(s[b])"
AllParams$Variable[AllParams$Variable == 'surv.int'] <- "paste(italic(s(x)), ' Intercept')"
AllParams$Variable[AllParams$Variable == 'surv.slope'] <- "paste(italic(s(x)), ' Slope')"

VR_Plot <- ggplot(AllParams, aes(x = Treatment)) + 
  geom_point(aes(y = Obs,
                 color = Treatment),
             size = 4.5) +
  geom_linerange(aes(ymax = UpCI,
                     ymin = LoCI,
                     color = Treatment),
                 size = 1.25) + 
  facet_wrap(~Variable,
             scales = 'free',
             labeller = label_parsed) + 
  theme(panel.background = element_blank(), # Set up plot thematic elements
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # legend.position = c(0.85, 0.1), # legend in bottom right corner
        legend.title = element_text(size = 22), # make legend text bigger
        legend.text = element_text(size = 18),
        legend.background = element_blank(),
        axis.title.y = element_text(size = 22,
                                    margin = margin(t = 0, # pad the axis title with a bit of whitespace
                                                    r = 20,
                                                    b = 0,
                                                    l = 0)),
        axis.title.x = element_text(size = 22,
                                    margin = margin(t = 20,
                                                    r = 0,
                                                    b = 10,
                                                    l = 0)),
        axis.line = element_line(colour = 'black', # make axis lines a bit bolder
                                 size = 2),
        axis.ticks = element_line(size = 1.2),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = 'white')) +
  scale_y_continuous('Observed Value + Confidence Intervals') + 
  scale_x_discrete('Treatment ') +
  scale_color_manual('Treatment',
                     breaks = c('Control', 'CR'),
                     values = c('black', 'green'))

VR_Plot

ggsave('Ailanthus_VR_Plot.png',
       device = 'png',
       path = 'Ailanthus_IPM/',
       height = 11,
       width = 17,
       units = 'in',
       dpi = 400)
