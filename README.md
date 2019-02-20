[![DOI](https://zenodo.org/badge/122057810.svg)](https://zenodo.org/badge/latestdoi/122057810)

## Population projection models for a variety of invasive species in the central United States

This repository hosts R and Stan code for M/I-PMs and a raw demographic data set for 14 non-native species in the central US. Data includes information under ambient conditions and a competition removal experiment. The data paper has been accepted by _Ecology_ and contains all relevant metadata. The only condition for using this is that you cite that paper.

## How to work with the data

Clone the repo or download the archive from [Zenodo](https://doi.org/10.5281/zenodo.2573062). This contains genus names for all the species and the type of model (MPM/IPM) separated by an underscore. Within each folder is a `.R` script named by the genus name and the model type, again separated by an underscore. Data sets are labelled with the genus name followed by `Clean.csv`, also separated by underscores (e.g. `Alliaria_Clean.csv`). Additional data sets are supplied when necessary to complete the analysis (e.g. germination data denoted with a `_Germ.csv` suffix, data sets for regressions to compute fecundity denoted with a `_Fec.csv` suffix). 