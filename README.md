# Space-time interactions in Bayesian disease mapping with recent tools: making things easier for practitioners
This repository contains the R code to fit the models described in the paper entitled *"Space-time interactions in Bayesian disease mapping with recent tools: making things easier for practitioners"*.

## Table of contents
[Data](##Data)  
[R code](##R-code)

[References](##References)



## Data
[BreastCancer_data.Rdata](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/BreastCancer_data.Rdata) containts breast cancer data in Spain during 1990-2010 used in this paper. 
This .Rdata file constains the following:

**Data:** a dataframe with the following variables
 - *Area*: provinces of Spain which are enumerated from 1 to 50
 - *Year*: year (from 1990 to 2010)
 - *Counts*: number of deaths of each province of Spain in each year due to breast cancer
 - *Expected*: expected cases of each province of Spain in each year

**Carto_ESP:** contains the cartography of each province of Spain

**Rs:** adjacency matrix


## R code
R code to fit the spatio-temporal models described in the paper has been included [here](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R).
Only models for the set of hyperprior distributions H1 are shown (to fit the models with H2 and H3 hyperprior distributions slight modifications are required in the code). 
- [icar_models](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/icar_models) and [bym_models](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/bym_models) folders contain the Rscripts with the spatio-temporal models fitted with ICAR and BYM spatial priors using R-INLA, Nimble 1 and Nimble 2. 
- [run](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/run) folder contains the Rscripts to run this models.
- [tables_figures_paper.R](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/tables_figures_paper.R) contains the necessary functions to reproduce all the figures and tables of Spanish breast cancer data analysis.

## Acknowledgements
This work has been supported by Project PID2020-113125RB-I00 (AEI).

## References
