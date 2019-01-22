# Script to compile all demographic models into a single place.
# This will then be used to modify the tyson.Rdata object in FunPhylo

rm(list = ls(all.names = TRUE)) 
library(FunPhylo)
library(dplyr)

data(tyson)

# Function to restore a relatively vanilla environment each time a 
# model script is sourced. 
.detachAllPackages <- function(keep = NULL) {
  
  KeepPackages <- c("package:stats",
                    "package:graphics",
                    "package:grDevices",
                    "package:utils",
                    "package:datasets",
                    "package:methods",
                    "package:base", 
                    keep)
  
  AllPackages <- search()[ifelse(unlist(gregexpr("package:",
                                                 search())) == 1 ,
                                 TRUE,
                                 FALSE)]
  
  PackageList <- setdiff(AllPackages, KeepPackages)
  
  if(length(PackageList) > 0) {
    lapply(PackageList,  
           detach, 
           character.only = TRUE)
  }
  
}

# Update CR Biomass with fully reproducible script using raw, plot-level data.
# last iteration done by hand in excel was a mess and this will be much better
source('CR_Biomass/Biomass_Cleaning.R')


# Create hidden object to store data since each script will clear all non-hidden objects
# in global environment. This script will not work if any of the sourced files contain
# a call to rm(list = ls(all = TRUE)) (emphasis on the all = TRUE part)

.DemoTable <- tibble(Species = tyson$demo.data$Species,
                    ESCR = NA_real_,
                    ESCR_2 = NA_real_,
                    CRBM = Means$Log_Standardized_Biomass,
                    Habitat = tyson$demo.data$Habitat)


# Going alphabetically through species. First up is Ailanthus
# Ailanthus--------------------
source('Ailanthus_IPM/Ailanthus_Figures.R') # Contains bootstrapping for both!


# Test if confidence intervals overlap. If not, effect size is LLR,
# if so, effect size is 0
if(min(CompN_lambda_CI) > max(ContN_lambda_CI) |
   max(CompN_lambda_CI) < min(ContN_lambda_CI)){
  AA_LRR <- AA_LRR2 <-  log(CompN_lambda + 0.5) - log(ContN_lambda + 0.5)
} else {
  AA_LRR <- 0
  AA_LRR2 <-  log(CompN_lambda + 0.5) - log(ContN_lambda + 0.5)
}

.DemoTable$ESCR[1] <- AA_LRR
.DemoTable$ESCR_2[1] <- AA_LRR2
# Alliaria-----------------
.detachAllPackages(keep = c('package:dplyr',  
                            'package:FunPhylo',
                            'package:ggplot2'))

source('Alliaria_MPM/Alliaria_MPM.R')
if(results$lower[6] > results$upper[5] |
   results$upper[6] < results$lower[5]) {
  GM_LRR <- GM_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
} else {
  GM_LRR <- 0
  GM_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
}

.DemoTable$ESCR[2] <- GM_LRR
.DemoTable$ESCR_2[2] <- GM_LRR2

# Carduus------------------
.detachAllPackages(keep = c('package:dplyr', 
                            'package:FunPhylo', 
                            'package:ggplot2'))

source('Carduus_MPM/Carduus_MPM.R')
if(results$lower[6] > results$upper[5] |
   results$upper[6] < results$lower[5]) {
  CN_LRR <- CN_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
} else {
  CN_LRR <- 0
  CN_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
}

.DemoTable$ESCR[3] <- CN_LRR
.DemoTable$ESCR_2[3] <- CN_LRR2

# Draba-------------------
.detachAllPackages(keep = c('package:dplyr',  
                            'package:FunPhylo', 
                            'package:ggplot2'))

source('Draba_MPM/Draba_MPM.R')
if(results$lower[6] > results$upper[5] |
   results$upper[6] < results$lower[5]) {
  DV_LRR <- DV_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
} else {
  DV_LRR <- 0
  DV_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
}

.DemoTable$ESCR[4] <- DV_LRR
.DemoTable$ESCR_2[4] <- DV_LRR2
# Euonymus -----------------------------
# Ran this on iDiv's RStudio server and saved results. Takes a while,
# but one can run it locally too using the code below. For now, I skip that
# and use the .Rdata file with outputs for this

# source('Euonymus_IPM/Euonymus_BootStrap_IPM.R')
PlotData <- read.csv('Euonymus_IPM/Euonymus_Summarized_Output.csv',
                     stringsAsFactors = FALSE)
results <- filter(PlotData, Variable == 'Lambda')

if(results$LoCI[2] > results$UpCI[1] |
   results$UpCI[2] < results$LoCI[1]) {
  EA_LRR <- EA_LRR2 <- log(results$obs[2] + 0.5) - log(results$obs[1] + 0.5)
} else {
  EA_LRR <- 0
  EA_LRR2 <- log(results$obs[2] + 0.5) - log(results$obs[1] + 0.5)
}

.DemoTable$ESCR[5] <- EA_LRR
.DemoTable$ESCR_2[5] <- EA_LRR2

# Kummerowia ------------------
.detachAllPackages(keep = c('package:dplyr',    
                            'package:FunPhylo', 
                            'package:ggplot2'))

source('Kummerowia_MPM/Kummerowia_MPM.R')
if(results$lower[6] > results$upper[5] |
   results$upper[6] < results$lower[5]) {
  KS_LRR <- KS_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
} else {
  KS_LRR <- 0
  KS_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
}

.DemoTable$ESCR[6] <- KS_LRR
.DemoTable$ESCR_2[6] <- KS_LRR2

# Lepidium---------------------
.detachAllPackages(keep = c('package:dplyr',   
                            'package:FunPhylo', 
                            'package:ggplot2'))

source('Lepidium_MPM/Lepidium_MPM.R')
if(results$lower[6] > results$upper[5] |
   results$upper[6] < results$lower[5]) {
  LepCam_LRR <- LepCam_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
} else {
  LepCam_LRR <- 0
  LepCam_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
}

.DemoTable$ESCR[7] <- LepCam_LRR
.DemoTable$ESCR_2[7] <- LepCam_LRR2

# Lespedeza------------------
.detachAllPackages(keep = c('package:dplyr',  
                            'package:FunPhylo', 
                            'package:ggplot2'))

source('Lespedeza_MPM/Lespedeza_MPM.R')
if(results$LoCI[38] > results$UpCI[37] |
   results$UpCI[38] < results$LoCI[37]) {
  LesCun_LRR <- LesCun_LRR2 <- log(results$values[38] + 0.5) - log(results$values[37] + 0.5)
} else {
  LesCun_LRR <- 0
  LesCun_LRR2 <- log(results$values[38] + 0.5) - log(results$values[37] + 0.5)
}

.DemoTable$ESCR[8] <- LesCun_LRR
.DemoTable$ESCR_2[8] <- LesCun_LRR2

# Ligustrum------------------
# Same as for Euonymus
results <- read.csv('Ligustrum_IPM/Ligustrum_Bootstrap_Output.csv',
                    stringsAsFactors = FALSE)
if(any(results$lambda_comp_Quad_lo > results$lambda_cont_Quad_up) |
   any(results$lambda_comp_Quad_up < results$lambda_cont_Quad_lo)) {
  LO_LRR <- LO_LRR2 <- log(mean(results$lambda_comp_Quad) + 0.5) -
                       log(mean(results$lambda_cont_Quad) + 0.5)
} else {
  LO_LRR <- 0
  LO_LRR2 <- log(mean(results$lambda_comp_Quad) + 0.5) -
             log(mean(results$lambda_cont_Quad) + 0.5)
}

.DemoTable$ESCR[9] <- LO_LRR
.DemoTable$ESCR_2[9] <- LO_LRR2

# Lonicera-------------------
PlotData <- read.csv('Lonicera_IPM/Lonicera_Summarized_Output.csv',
                     stringsAsFactors = FALSE) 
results <- filter(PlotData, Variable == 'Lambda' &
                    SurvModel == 'Quad')

if(results$LoCI[2] > results$UpCI[1] |
   results$UpCI[2] < results$LoCI[1]) {
  LM_LRR <- LM_LRR2 <- log(results$obs[2] + 0.5) - log(results$obs[1] + 0.5)
} else {
  LM_LRR <- 0
  LM_LRR2 <- log(results$obs[2] + 0.5) - log(results$obs[1] + 0.5)
}

.DemoTable$ESCR[10] <- LM_LRR
.DemoTable$ESCR_2[10] <- LM_LRR2

# Perilla-------------------------
.detachAllPackages(keep = c('package:dplyr', 
                            'package:FunPhylo',
                            'package:ggplot2'))

source('Perilla_MPM/Perilla_MPM.R')
if(results$lower[6] > results$upper[5] |
   results$upper[6] < results$lower[5]) {
  PF_LRR <- PF_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
} else {
  PF_LRR <- 0
  PF_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
}

.DemoTable$ESCR[11] <- PF_LRR
.DemoTable$ESCR_2[11] <- PF_LRR2

# Potentilla------------------------
.detachAllPackages(keep = c('package:dplyr',
                            'package:FunPhylo',
                            'package:ggplot2'))

source('Potentilla_MPM/Potentilla_MPM.R')
if(results$lower[16] > results$upper[15] |
   results$upper[16] < results$lower[15]) {
  PR_LRR <- PR_LRR2 <- log(results$values[16] + 0.5) - log(results$values[15] + 0.5)
} else {
  PR_LRR <- 0
  PR_LRR2 <- log(results$values[16] + 0.5) - log(results$values[15] + 0.5)
}

.DemoTable$ESCR[12] <- PR_LRR
.DemoTable$ESCR_2[12] <- PR_LRR2

# Thlaspi ---------------------
.detachAllPackages(keep = c('package:dplyr', 
                            'package:FunPhylo',
                            'package:ggplot2'))

source('Thlaspi_MPM/Thlaspi_MPM.R')
if(results$lower[6] > results$upper[5] |
   results$upper[6] < results$lower[5]) {
  TP_LRR <- TP_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
} else {
  TP_LRR <- 0
  TP_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
}

.DemoTable$ESCR[13] <- TP_LRR
.DemoTable$ESCR_2[13] <- TP_LRR2

# Verbascum-------------------
.detachAllPackages(keep = c('package:dplyr',
                            'package:FunPhylo',  
                            'package:ggplot2'))

source('Verbascum_MPM/Verbascum_MPM.R')
if(results$lower[6] > results$upper[5] |
   results$upper[6] < results$lower[5]) {
  VT_LRR <- VT_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
} else {
  VT_LRR <- 0
  VT_LRR2 <- log(results$values[6] + 0.5) - log(results$values[5] + 0.5)
}

.DemoTable$ESCR[14] <- VT_LRR
.DemoTable$ESCR_2[14] <- VT_LRR2

.DemoTable <- setNames(.DemoTable, c('Species', 
                                     'ESCR2',
                                     'ESCR',
                                     'CRBM',
                                     'Habitat'))

save(.DemoTable, file = 'Master_R/Cleaned_Demographic_Model_Output.rda')

# Move the hidden object into the tyson data object and save it!
tyson$demo.data <- .DemoTable
devtools::use_data(tyson,
                   pkg = '../../Tyson/Writing/Tyson_Traits_Analysis/Fun_Phylo_Package/Fun_Phylo_Package',
                   overwrite = TRUE)

# Done! :)

