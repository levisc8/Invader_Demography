# Count total number of control plots

rm(list = ls())
library(fs)
library(dplyr)

# Read in all of the raw data files
speciesRegex <- c('Ailanthus_Clean|Alliaria_Clean|Carduus_Clean|Draba_Clean|Euonymus_Clean|Kummerowia_14_Clean|Lepidium_Clean|Lespedeza_Clean|Ligustrum_Clean|Lonicera_Clean|Perilla_Clean|Potentilla_Clean|Thlaspi_Clean|Verbascum_Clean')

FilePaths <- dir_ls(recursive = TRUE,
                    regexp = speciesRegex) %>%
  .[grepl('.csv', .)] %>%
  .[!grepl('_RS_', .)] %>%
  .[!grepl('_Germ', .)]

RawFiles <- lapply(FilePaths, function(x) read.csv(x, stringsAsFactors = FALSE))

# Check names for consistency
Names <- lapply(RawFiles, names)

# Stop if Plot isn't found in all of them
test <- lapply(Names, function(x) any(grepl('Plot', x))) %>%
  unlist() %>% 
  all() %>%
stopifnot()


N_Plots <- lapply(RawFiles, function(x) {
  x %>%
    filter(Treatment == 'Control' |
             Treatment == 'Cont' |
             Treatment == 'ContN' |
             Treatment == 'Cont') %>%
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
