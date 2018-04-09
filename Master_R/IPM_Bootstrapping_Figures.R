rm(list = ls())


library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# Read in the output from the boot strapping run on iDiv's RStudio server.
LigObtData <- read.csv('Ligustrum/IPM/Data/Bootstrapping_Output_Ligustrum.csv',
                       stringsAsFactors = FALSE) %>%
  .[-c(7), ]

ggplot(LigObtData, aes(x = est.prob)) +
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



ggsave(filename = 'Lambda_EstP.png',
       path = 'Ligustrum/IPM/Figures',
       height = 8,
       width = 8,
       units = 'in')

# Now make the boot strapped vital rate figures
LigObtVRData <- read.csv('Ligustrum/IPM/Data/Ligustrum_Bootstrapped_Vital_Rates.csv',
                         stringsAsFactors = FALSE) %>%
  gather(key = 'Variable', value = 'Value', -Obs_Boot) %>%
  mutate(Trt = vapply(.$Variable, 
                      FUN = function(x) str_split(x, "_")[[1]][2],
                      FUN.VALUE = ''),
         SurvModel = vapply(.$Variable,
                            FUN = function(x) str_split(x, '_')[[1]][3],
                            FUN.VALUE = '')) %>%
  group_by(Variable, Trt) %>%
  arrange(desc(Value)) %>%
  summarise(Obs = Value[Obs_Boot == 'Observed'],
            UpCI = Value[26],
            LoCI = Value[976]) %>%
  ungroup() %>%
  mutate(Variable = vapply(.$Variable,
                           FUN = function(x) str_split(x, '_')[[1]][1],
                           FUN.VALUE = ''))

# brief sanity check
all(LigObtVRData$Obs > LigObtVRData$LoCI &
      LigObtVRData$Obs < LigObtVRData$UpCI)

# Fix a dumb typo from my bootstrapping code.
LigObtVRData$Variable[LigObtVRData$Variable == 'SurvSelope2'] <- 'SurvSlope2'


# Start renaming parameters so figure titles look pretty
LigObtVRData$Trt[is.na(LigObtVRData$Trt)] <- 'Pooled'
LigObtVRData$Trt[LigObtVRData$Trt == 'Cont'] <- 'Control'
LigObtVRData$Variable[LigObtVRData$Variable == 'RecMean'] <- "paste(mu, ' Recruit Size Mean' )"
LigObtVRData$Variable[LigObtVRData$Variable == 'RecSD'] <- "paste(sigma, ' Recruit Size SD')"
LigObtVRData$Variable[LigObtVRData$Variable == 'RepInt'] <- "paste( italic(f[p](x)), ' Intercept')"
LigObtVRData$Variable[LigObtVRData$Variable == 'RepSlope'] <- "paste(italic(f[p](x)), ' Slope')"
LigObtVRData$Variable[LigObtVRData$Variable == 'SurvInt'] <- "paste(italic(s(x)),' Intercept')"
LigObtVRData$Variable[LigObtVRData$Variable == 'SurvSlope'] <- "paste(italic(s(x)),' Linear Term')"
LigObtVRData$Variable[LigObtVRData$Variable == 'SurvSlope2'] <- "paste(italic(s(x)),' Quadratic Term')"
LigObtVRData$Variable[LigObtVRData$Variable == 'GrowInt'] <- "paste(italic(g(y,x)),' Intercept')"
LigObtVRData$Variable[LigObtVRData$Variable == 'GrowSlope'] <- "paste(italic(g(y,x)),' Slope')"



ggplot(data = LigObtVRData,
       aes(x = Trt)) + 
  geom_point(aes(y = Obs, 
                 color = Trt),
             size = 4.5) + 
  geom_linerange(aes(ymin = LoCI,
                     ymax = UpCI,
                     color = Trt),
                 size = 1.25) + 
  facet_wrap(~Variable,
             scales = 'free',
             labeller = label_parsed) +  
  theme(panel.background = element_blank(), # Set up plot thematic elements
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # legend.position = c(0.85, 0.1), # legend in top bottom right corner
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
  scale_x_discrete('Treatment (if applicable)') +
  scale_color_manual('A\nTreatment',
                     breaks = c('Control', 'CR', "Pooled"),
                     values = c('black', 'green', 'blue'))

ggsave(filename = 'Ligustrum_Vital_Rate_Coefficients.png',
       path = 'Ligustrum/IPM/Figures',
       height = 9,
       width = 11,
       unit = 'in')


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

PlotData$Trt[is.na(PlotData$Trt)] <- 'Pooled'
PlotData$Trt[PlotData$Trt == 'Cont'] <- 'Control'
PlotData$Variable[PlotData$Variable == 'Recruit Size Mean'] <- "paste(mu, ' Recruit Size Mean' )"
PlotData$Variable[PlotData$Variable == 'Recruit Size SD'] <- "paste(sigma, ' Recruit Size SD')"
PlotData$Variable[PlotData$Variable == 'Repro Intercept'] <- "paste( italic(f[p](x)), ' Intercept')"
PlotData$Variable[PlotData$Variable == 'Repro Slope'] <- "paste(italic(f[p](x)), ' Slope')"
PlotData$Variable[PlotData$Variable == 'Surv Intercept'] <- "paste(italic(s(x)),' Intercept')"
PlotData$Variable[PlotData$Variable == 'Surv Linear Term'] <- "paste(italic(s(x)),' Linear Term')"
PlotData$Variable[PlotData$Variable == 'Surv Quadratic Term'] <- "paste(italic(s(x)),' Quadratic Term')"
PlotData$Variable[PlotData$Variable == 'Lambda'] <- 'lambda'


VR_Plot <- ggplot(data = PlotData,
                     aes(x = Trt)) + 
  geom_point(aes(y = obs, 
                 color = Trt),
             size = 4.5) + 
  geom_linerange(aes(ymin = LoCI,
                    ymax = UpCI,
                    color = Trt),
                 size = 1.25) + 
  facet_wrap(~Variable,
             scales = 'free',
             labeller = label_parsed) +  
  theme(panel.background = element_blank(), # Set up plot thematic elements
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.1), # legend in top bottom right corner
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
  scale_x_discrete('Treatment (if applicable)') +
  scale_color_manual('Treatment',
                     breaks = c('Control', 'CR', "Pooled"),
                     values = c('black', 'green', 'blue'))
        
VR_Plot    

ggsave(filename = 'Euonymus_Vital_Rate_Coefficients.png',
       path = 'Euonymus/IPM/Figures',
       height = 8,
       width = 10,
       unit = 'in')


# LonMaaData <- read.csv('Lonicera/IPM/Lonicera_RS_Folder/BootStrap_Output_Lonicera.csv',
#                        stringsAsFactors = FALSE) %>%
#   gather(key = 'Variable', value = 'Value', -Boot_Obs) %>%
#   mutate(Trt = vapply(.$Variable,
#                       FUN = function(x) str_split(x, '_')[[1]][2],
#                       FUN.VALUE = ''),
#          SurvModel = vapply(.$Variable,
#                             FUN = function(x) str_split(x, '_')[[1]][3],
#                             FUN.VALUE = '')) %>%
#   group_by(Variable, Trt, SurvModel) %>%
#   arrange(desc(Value)) %>%
#   summarise(obs = Value[Boot_Obs == 'Observed'],
#             UpCI = Value[25],
#             LoCI = Value[975]) %>%
#   ungroup() %>%
#   mutate(Variable = vapply(.$Variable,
#                            FUN = function(x) str_split(x, '_')[[1]][1],
#                            FUN.VALUE = ''))

# Function to re-run just the part of the IPM that calculates survival
# https://stackoverflow.com/questions/12214963/source-only-part-of-a-file
# sourcePartial <- function(file, start, end, ...) {
#   file.lines <- scan(file, what = character(),
#                      skip = start - 1, nlines = end-start + 1,
#                      sep = '\n')
#   file.lines.collapsed <- paste(file.lines, collapse='\n')
#   source(textConnection(file.lines.collapsed), ...)
# }
# 
# 
# sourcePartial('Lonicera/IPM/R/Lonicera_IPM.R', start = 6, end = 136)
# 
# SurvData <- tibble(Variable = rep(c('SurvInt', 'SurvSlope', 
#                                       'SurvInt', 'SurvSlope', 'SurvSlope2'), 2),
#                        Trt = c(rep('Cont', 5),
#                                rep('CR', 5)),
#                        SurvModel = c(rep('Lin', 2),
#                                      rep('Quad', 3),
#                                      rep('Lin', 2),
#                                      rep('Quad', 3)))
# 
# 
# PlotData <- rbind(fixef(ContSurvBRM_Lin),
#                    fixef(ContSurvBRM_Quad),
#                    fixef(CompSurvBRM_Lin), 
#                    fixef(CompSurvBRM_Quad)) %>%
#   as_tibble() %>%
#   .[ ,-2] %>%
#   setNames(c('obs', 'LoCI', 'UpCI')) %>%
#   select(obs, UpCI, LoCI) %>%
#   cbind(SurvData, .) %>%
#   rbind(LonMaaData)
# save(PlotData, file = 'Lonicera/IPM/Data/Lonicera_Vital_Rate_Summary.RData')
load('Lonicera/IPM/Data/Lonicera_Vital_Rate_Summary.RData')

PlotData$Trt[PlotData$Trt == 'None'] <- 'Pooled'
PlotData$Trt[PlotData$Trt == 'Cont'] <- 'Control'
PlotData$Variable[PlotData$Variable == 'GrowInt'] <- "paste(italic(g(y,x)), ' Intercept')"
PlotData$Variable[PlotData$Variable == 'GrowSlope'] <- "paste(italic(g(y,x)), ' Slope')"
PlotData$Variable[PlotData$Variable == 'RecMean'] <- "paste(mu, ' Recruit Size Mean' )"
PlotData$Variable[PlotData$Variable == 'RecSD'] <- "paste(sigma, ' Recruit Size SD')"
PlotData$Variable[PlotData$Variable == 'RepInt'] <- "paste( italic(f[p](x)), ' Intercept')"
PlotData$Variable[PlotData$Variable == 'RepSlope'] <- "paste(italic(f[p](x)), ' Slope')"
PlotData$Variable[PlotData$Variable == 'SeedInt'] <- "paste(italic(f[s](x)), ' Intercept')"
PlotData$Variable[PlotData$Variable == 'SeedSlope'] <- "paste(italic(f[s](x)),' Slope')"
PlotData$Variable[PlotData$Variable == 'SurvInt'] <- "paste(italic(s(x)),' Intercept')"
PlotData$Variable[PlotData$Variable == 'SurvSlope'] <- "paste(italic(s(x)),' Linear Term')"
PlotData$Variable[PlotData$Variable == 'SurvSlope2'] <- "paste(italic(s(x)),' Quadratic Term')"
PlotData$Variable[PlotData$Variable == 'Lambda'] <- 'lambda'

# Create data set for surivival insensitive vital rates
PlotData2 <- filter(PlotData, is.na(SurvModel))

LM_SurvIns_Plot <- ggplot(data = PlotData2,
                          aes(x = Trt)) + 
  geom_point(aes(y = obs, 
                 color = Trt),
             size = 4.5) + 
  geom_linerange(aes(ymin = LoCI,
                     ymax = UpCI,
                     color = Trt),
                 size = 1.25) + 
  facet_wrap(~Variable,
             scales = 'free',
             labeller = label_parsed) +
  theme(panel.background = element_blank(), # Set up plot thematic elements
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.87, 0.1), # legend in top bottom right corner
        legend.title = element_text(size = 22), # make legend text bigger
        legend.text = element_text(size = 17),
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
    scale_color_manual('A\nTreatment',
                       breaks = c('Control', 'CR', 'Pooled'),
                       values = c('black', 'green', 'blue')) + 
  scale_y_continuous('Observed Value + Confidence Intervals') + 
  scale_x_discrete('Treatment (if applicable)') 

LM_SurvIns_Plot

ggsave(filename = 'Lonicera_Survival_Insensitive_Vital_Rates.png',
       device = 'png',
       width = 10,
       height = 8,
       units = 'in',
       path = 'Lonicera/IPM/Figures')

# Create two dummy rows to occupy the right side of the 
# x-axis for SurvSlope^2
DummyRows <- tibble(Variable = rep("paste(italic(s(x)),' Quadratic Term')", 2),
                    Trt = rep(NA, 2),
                    SurvModel = rep(NA, 2),
                    obs = rep(NA, 2),
                    UpCI = rep(NA_real_, 2),
                    LoCI = rep(NA_real_, 2))

PlotData1 <- filter(PlotData, !is.na(SurvModel)) %>%
  rbind(., DummyRows) %>%
  mutate(Dummy = paste(.$Trt, .$SurvModel, sep = '-')) 
PlotData1$Dummy <- gsub('NA-NA', 'Legend', PlotData1$Dummy) 
PlotData1$SurvModel[PlotData1$SurvModel == 'Lin'] <- 'Linear'
PlotData1$SurvModel[PlotData1$SurvModel == 'Quad'] <- 'Quadratic'

PlotData1$Facet <- relevel(as.factor(PlotData1$Variable), ref = 'lambda')

PlotData1 %>%
  as_tibble() %>%
  ggplot(data = .,
         aes(x = Dummy)) + 
  geom_point(aes(y = obs, 
                 shape = SurvModel,
                 color = Trt),
             size = 4.5) + 
  geom_linerange(aes(ymax = UpCI,
                     ymin = LoCI,
                     color = Trt),
                 size = 1.25) +
  facet_wrap(~Facet,
             scales = 'free',
             labeller = label_parsed) +  
  theme(panel.background = element_blank(), # Set up plot thematic elements
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.91, 0.2), # legend in top bottom right corner
        legend.title = element_text(size = 20), # make legend text bigger
        legend.text = element_text(size = 17),
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
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = 'white')) +
  scale_y_continuous('Observed Value + Confidence Intervals') + 
  scale_x_discrete('Treatment and Survival Model') +
  labs(color = 'B\nTreatment', shape = 'Survival Model') +
  scale_color_manual(breaks = c('Control', 'CR'),
                     values = c('black', 'green')) + 
  scale_shape_manual(breaks = c('Quadratic', 'Linear'),
                     values = c(17,19))

ggsave(filename = 'Lonicera_Survival_Sensitive_Vital_Rates.png',
       device = 'png',
       width = 12,
       height = 8,
       units = 'in',
       path = 'Lonicera/IPM/Figures')
