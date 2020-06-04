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

lambdas <- data.frame(lambda_c = ail_co$lambda,
                      lambda_cr = ail_cr$lambda,
                      boot_obs = c("observed", rep("bootstrap", 1000)))

saveRDS(lambdas, file = "../Data/bootstrap_lambdas/Ailanthus_lambdas.rds")

rm(lambdas)

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

euo_pt_ests <- read.csv("Euonymus_IPM/Euonymus_Summarized_Output.csv",
                        stringsAsFactors = FALSE)

lambdas <- data.frame(lambda_c = c(euo_pt_ests$obs[7], euo_all$Lambda_Cont),
                      lambda_cr = c(euo_pt_ests$obs[8], euo_all$Lambda_CR),
                      boot_obs = c("observed", rep("bootstrap", 1000)))

saveRDS(lambdas, file = "../Data/bootstrap_lambdas/Euonymus_lambdas.rds")
rm(lambdas)

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

sim_l_c <- rnorm(1000, mean = mean(lig_all$lambda_cont_Quad), 
                 sd = (mean(lig_all$lambda_cont_Quad - lig_all$lambda_cont_Quad_lo)) / 1.96)

sim_l_cr <- rnorm(1000, mean = mean(lig_all$lambda_comp_Quad), 
                 sd = (mean(lig_all$lambda_comp_Quad - lig_all$lambda_comp_Quad_lo)) / 1.96)

lambdas <- data.frame(lambda_c = c(mean(lig_all$lambda_cont_Quad), sim_l_c),
                      lambda_cr = c(mean(lig_all$lambda_comp_Quad), sim_l_cr),
                      boot_obs = c("observed", rep("bootstrap", 1000)))

saveRDS(lambdas, file = "../Data/bootstrap_lambdas/Ligustrum_lambdas.rds")

rm(lambdas)

var_escr <- var(log(lig_all$lambda_comp_Quad + 0.5) - log(lig_all$lambda_cont_Quad + 0.5))
demo_data$var_escr[demo_data$Species == 'Ligustrum_obtusifolium'] <- var_escr

# Lonicera - Outputs stored in csv file. First row is for values from observed
# data, next 1000 are bootstrapped values.

lon_all <- read.csv('Lonicera_IPM/Lonicera_Bootstrap_Output.csv',
                    stringsAsFactors = FALSE)

lambdas <- data.frame(lambda_c = lon_all$Lambda_Cont_Quad,
                      lambda_cr = lon_all$Lambda_CR_Quad,
                      boot_obs = c("observed", rep("bootstrap", 1000)))

saveRDS(lambdas, file = "../Data/bootstrap_lambdas/Lonicera_lambdas.rds")

rm(lambdas)

var_escr <- var(log(lon_all$Lambda_CR_Quad[2:1001] + 0.5) - log(lon_all$Lambda_Cont_Quad[2:1001] + 0.5))

demo_data$var_escr[demo_data$Species == 'Lonicera_maackii'] <- var_escr

# Perilla - var_escr is created by MPM script

source('Perilla_MPM/Perilla_MPM.R')
demo_data$var_escr[demo_data$Species == 'Perilla_frutescens'] <- var_escr
demo_data$ESCR[demo_data$Species == 'Perilla_frutescens'] <- log(results$values[5] + 0.5) - log(results$values[6] + 0.5)

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
