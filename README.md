# Space-time interactions in Bayesian disease mapping with recent tools: making things easier for practitioners
This repository contains the R code to fit the models described in the paper entitled *"Space-time interactions in Bayesian disease mapping with recent tools: making things easier for practitioners"*.

## Table of contents
[Data](##Data)  
[R code](##R-code)<br/>
[References](##References)


# Data
Female breast cancer mortality data (ICD-10 code 50) in Spanish provinces during the period 1990-2010.

- [**BreastCancer_data.Rdata**](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/BreastCancer_data.Rdata)
  
  This .Rdata contains the following objects
  - **_Data_**: `data.frame` object with the number of observed and expected cases (_'Counts'_ and _'Expected'_ variables, respectively) for each province (_'Area'_) and time period (_'Year'_) for female breast cancer mortality data.
  - **_Carto_ESP_**: `sf` object containing the spatial polygons of the Spanish provinces. The data contains a `data.frame` with 50 rows and  _'Area'_ (character vector of geographic identifiers), _'Name'_ (character vector of province names), _'Longitude'_ (numeric vector of longitude values), _'Latitude'_ (numeric vector of latitude values) and _'geometry'_ (sfc_MULTIPOLYGON) variables.
  - **_Rs_**: adjacency matrix.
	

## R code
R code to fit the spatio-temporal models described in the paper has been included [here](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R).
Only models for the set of hyperprior distributions H1 are shown (to fit the models with H2 and H3 hyperprior distributions slight modifications are required in the code). 
- [icar_models](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/icar_models) and [bym_models](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/bym_models) folders contain the Rscripts with the spatio-temporal models fitted with ICAR and BYM spatial priors using R-INLA, Nimble 1 and Nimble 2. 
- [run](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/run) folder contains the Rscripts to run these models.
- [tables_figures_paper.R](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/tables_figures_paper.R) contains the necessary functions to reproduce all the figures and tables of Spanish breast cancer mortality data analysis.

## Acknowledgements
This work has been supported by Project PID2020-113125RB-I00 (AEI).

## References
	Urdangarin, A., Goicoa, T., and Ugarte, M.D. (2021). Space-time interactions in Bayesian disease mapping with recent tools: making things easier for practitioners. Submitted to _Statistical Methods in Medical Research_ accepted, January 2022.
