# Count total number of control plots

rm(list = ls())
library(fs)
library(dplyr)

# Read in all of the raw data files
speciesRegex <- c('Ailanthus4R|GM_Clean|CN_Clean|DV_Clean|EA_Clean|KS_14_Clean|LepCam_Clean|LesCun_Clean|LO_Clean|LM_Clean|PF_Clean|PR_Clean|TP_Clean|VT_Clean')

FilePaths <- dir_ls(recursive = TRUE,
                    regexp = speciesRegex) %>%
  .[grepl('.csv', .)] %>%
  .[!grepl('_RS_', .)] %>%
  .[!grepl('_Germ', .)]

RawFiles <- lapply(FilePaths, function(x) read.csv(x, stringsAsFactors = FALSE))

# Check names for consistency
Names <- lapply(RawFiles, names)

# Fix ailanthus
names(RawFiles[[1]])[3] <- 'Treatment'

N_Plots <- lapply(RawFiles, function(x) {
  x %>%
    filter(Treatment == 'Control' |
             Treatment == 'Cont' |
             Treatment == 'ContN') %>%
    summarise(N = length(unique(Plot))) %>%
    unlist()
  
})


unlist(N_Plots) %>% sum()


# Now all plots in CR or control
All_Plots <- lapply(RawFiles, function(x) {
  x %>%
    filter(Treatment == 'Control' |
             Treatment == 'Cont' |
             Treatment == 'ContN' |
             Treatment == 'Competitor' | 
             Treatment == 'Comp' | 
             Treatment == 'CompN') %>%
    summarise(N = length(unique(Plot))) %>%
    unlist()
  
})


unlist(All_Plots) %>% sum()