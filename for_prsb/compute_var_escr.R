# Generate variance in ESCR for weighting in regression step

library(FunPhylo)

data(tyson)

demo_data <- tyson$demo.data
demo_data$var_escr <- NA_real_


# Ailanthus - Outputs are stored in RDS files so we don't need to re-run these
ail_cr <- readRDS('Ailanthus_IPM/CR_Bootstrap_Output.rds')
ail_co <- readRDS('Ailanthus_IPM/Cont_Bootstrap_Output.rds')

ail_cr_l <- ail_cr$lambda[-1]
ail_co_l <- ail_co$lambda[-1]

demo_data$var_escr[demo_data$Species == 'Ailanthus_altissima'] <- 
  var(log(ail_cr_l + 0.5) - log(ail_co_l + 0.5))

# Alliaria - need to re-run this one. var_escr is created in the script so just insert it
source('Alliaria_MPM/Alliaria_MPM.R')
demo_data$var_escr[demo_data$Species == 'Alliaria_petiolata'] <- var_escr


# Carduus - ibid
source('Carduus_MPM/Carduus_MPM.R')
demo_data$var_escr[demo_data$Species == 'Carduus_nutans'] <- var_escr

# Draba - ibid
source('Draba_MPM/Draba_MPM.R')
demo_data$var_escr[demo_data$Species == 'Draba_verna'] <- var_escr


# Euonymus - lambdas stored in csv

euo_all <- read.csv('Euonymus_IPM/Euonymus_Bootstrap_Output.csv',
                    stringsAsFactors = FALSE)

demo_data$var_escr[demo_data$Species == 'Euonymus_alatus'] <-
  var(log(euo_all$Lambda_CR + 0.5) - log(euo_all$Lambda_Cont + 0.5))

# Kummerowia - need to re-run this one. var_escr is created in the script so just insert it
source('Kummerowia_MPM/Kummerowia_MPM.R')
demo_data$var_escr[demo_data$Species == 'Kummerowia_striata'] <- var_escr


# Lepidium - need to re-run this one. var_escr is created in the script so just insert it
source('Lepidium_MPM/Lepidium_MPM.R')
demo_data$var_escr[demo_data$Species == 'Lepidium_campestre'] <- var_escr

# Lespedeza - need to re-run this one. var_escr is created in the script so just insert it
source('Lespedeza_MPM/Lespedeza_MPM.R')
demo_data$var_escr[demo_data$Species == 'Lespedeza_cuneata'] <- var_escr


# Ligustrum - outputs stored in csv file. 
lig_all <- read.csv("Ligustrum_IPM/Ligustrum_All_Lambdas.csv",
                    stringsAsFactors = FALSE)

var_escr <- var(log(lig_all$lambda_comp_Quad + 0.5) - log(lig_all$lambda_cont_Quad + 0.5))
demo_data$var_escr[demo_data$Species == 'Ligustrum_obtusifolium'] <- var_escr

# Lonicera - Outputs stored in csv file

lon_all <- read.csv('Lonicera_IPM/Lonicera_Bootstrap_Output.csv',
                    stringsAsFactors = FALSE)[-1, ]

var_escr <- var(log(lon_all$Lambda_CR_Quad + 0.5) - log(lon_all$Lambda_Cont_Quad + 0.5))

demo_data$var_escr[demo_data$Species == 'Lonicera_maackii'] <- var_escr

# Perilla - var_escr is created by MPM script

source('Perilla_MPM/Perilla_MPM.R')
demo_data$var_escr[demo_data$Species == 'Perilla_frutescens'] <- var_escr

# Potentilla
source("Potentilla_MPM/Potentilla_MPM.R")
demo_data$var_escr[demo_data$Species == 'Potentilla_recta'] <- var_escr

# Thlaspi
source("Thlaspi_MPM/Thlaspi_MPM.R")
demo_data$var_escr[demo_data$Species == 'Thlaspi_perfoliatum'] <- var_escr

# Verbascum
source('Verbascum_MPM/Verbascum_MPM.R')
demo_data$var_escr[demo_data$Species == 'Verbascum_thapsus'] <- var_escr


tyson$demo.data <- demo_data

save(tyson, 
     file = '../../Tyson/Writing/Tyson_Traits_Analysis/Fun_Phylo_Package/Fun_Phylo_Package/data/tyson.rda')
