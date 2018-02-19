rm(list = ls())


library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# Read in the output from the boot strapping run on iDiv's RStudio server.
LigObtData <- read.csv('Ligustrum/IPM/Data/Bootstrapping_Output_Ligustrum.csv',
                       stringsAsFactors = FALSE) %>%
  .[-c(7), ]

LigObtPlot <- ggplot(LigObtData, aes(x = est.prob)) +
  theme(panel.background = element_blank(), # Set up plot thematic elements
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.15, 0.9), # legend in top left corner
        legend.title = element_text(size = 22), # make legend text bigger
        legend.text = element_text(size = 18),
        axis.title.y = element_text(size = 22,
                                    margin = margin(t = 0, # pad the axis title with a bit of whitespace
                                                    r = 20,
                                                    b = 0,
                                                    l = 0)),
        axis.title.x = element_text(size = 22,
                                    margin = margin(t = 20,
                                                    r = 0,
                                                    b = 0,
                                                    l = 0)),
        axis.line = element_line(colour = 'black', # make axis lines a bit bolder
                                 size = 2),
        axis.ticks = element_line(size = 1.2),
        axis.text = element_text(size = 16)) + 
  geom_line(aes(y = lambda_comp_Quad,
                color = 'CR'),
            size = 1.5,
            alpha = 1,
            linetype = 'dashed') + 
  geom_ribbon(aes(ymin = lambda_comp_Quad_lo, # add shaded confidence intervals
                  ymax = lambda_comp_Quad_up),
              fill = 'green',
              alpha = 0.2) + 
  geom_line(aes(y = lambda_cont_Quad,
                color = 'Control'),
            alpha = 1,
            size = 1.5,
            linetype = 'dashed') + 
  geom_ribbon(aes(ymin = lambda_cont_Quad_lo,
                  ymax = lambda_cont_Quad_up),
              fill = 'black',
              alpha = 0.2) + 
  scale_y_continuous(expression(paste('Per-capita Growth Rate (', lambda, ')')),
                     breaks = seq(1, 1.8, 0.2)) + 
  scale_x_continuous(expression(paste('Establishment Probability (', italic(E[p]), ')'))) + 
  scale_color_manual('Treatment',
                     breaks = c('Control', 'CR'),
                     values = c('black', 'green')) + 
  geom_text(aes(x = 1, y = 1.8, label = 'B'), size = 16) 


LigObtPlot

# ggsave(filename = 'Lambda_EstP.png', 
#        path = 'Ligustrum/IPM/Figures',
#        height = 8, 
#        width = 8, 
#        units = 'in')

# Ligustrum done! Now on to Euonymus. This one is a bit different, since
# we don't have some continuous variable that we are interested in. Doing
# a normal dot and whisker plot
# Euonymus needs a bit of tidying. Extracting the confidence intervals as well.
# We'll add in the observed values next.


# I've saved the output as an RData object. I plan to use that, the IPM Script
# and the raw data itself as the submission to BDJ. That way, one can 
# use the script and raw data to re analyze as they see fit, but also just
# skip to the summarized data if that is all they want.

# EuoAlaData <- read.csv('Euonymus/IPM/EuoAla_RS_Folder/BootStrap_Output_Euonymus.csv',
#                        stringsAsFactors = FALSE) %>%
#   gather(key = 'Variable', value = 'Value') %>%
#   mutate(Trt = vapply(.$Variable, 
#                       FUN = function(x) str_split(x, '_')[[1]][2],
#                       FUN.VALUE = '')) %>%
#   group_by(Variable, Trt) %>%
#   arrange(desc(Value)) %>%
#   summarise(obs = NA,
#             UpCI = Value[25],
#             LoCI = Value[975])


# # Insert observed values into bootstrap output
# source('Euonymus/IPM/R/Euonymus_IPM.R')
# 
# EuoAlaData$obs[1] <- lambda_Cont_obs
# EuoAlaData$obs[2] <- lambda_CR_obs
# EuoAlaData$obs[3] <- Sdl.mean
# EuoAlaData$obs[4] <- Sdl.SD
# EuoAlaData$obs[5] <- f.params$prob.repro.int
# EuoAlaData$obs[6] <- f.params$prob.repro.slope
# 
# SurvCIsCR <- fixef(BrmSurvQuad_CR)[ ,-2]
# SurvCIsCont <- fixef(BrmSurvQuad_Cont)[ ,-2]
# 
# EuoAlaData$Variable <- c('Lambda', 'Lambda',
#                          'Recruit Size Mean',
#                          'Recruit Size SD',
#                          'Repro Intercept',
#                          'Repro Slope')
# 
# 
# SurvSummary <- tibble(Variable = rep(c('Surv Intercept',
#                                        'Surv Linear Term',
#                                        'Surv Quadratic Term'), 2),
#                       Trt = c(rep('Cont', 3),
#                               rep('CR', 3)),
#                       obs = c(SurvCIsCont[ ,1],
#                               SurvCIsCR[ ,1]),
#                       UpCI = c(SurvCIsCont[ ,3],
#                                SurvCIsCR[ ,3]),
#                       LoCI = c(SurvCIsCont[ ,2],
#                                SurvCIsCR[ ,2]))
# 
# 
# PlotData <- rbind.data.frame(SurvSummary,
#                              EuoAlaData, 
#                              stringsAsFactors = FALSE)

# save(PlotData, file = 'Euonymus/IPM/Data/BootStrap_Summary_Data.RData')

# Load in saved parameter estimates and CIs
load('Euonymus/IPM/Data/BootStrap_Summary_Data.RData')


VR_Plot <- ggplot(data = PlotData,
                     aes(x = Trt)) + 
  geom_point(aes(y = obs, 
                 color = Trt),
             size = 2.5) + 
  geom_linerange(aes(ymin = LoCI,
                    ymax = UpCI,
                    color = Trt)) + 
  facet_wrap(~Variable,
             scales = 'free') +  
  theme(panel.background = element_blank(), # Set up plot thematic elements
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9, 0.1), # legend in top left corner
        legend.title = element_text(size = 22), # make legend text bigger
        legend.text = element_text(size = 18),
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
        axis.text = element_text(size = 16)) +
  scale_y_continuous('Observed Value + Confidence Intervals') + 
  scale_x_discrete('Treatment (if applicable)') +
  scale_color_discrete('Treatment')
        
VR_Plot    

ggsave(filename = 'Euonymus_Vital_Rate_Coefficients.png',
       path = 'Euonymus/IPM/Figures',
       height = 8,
       width = 8,
       unit = 'in')
