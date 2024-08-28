# 2024-PSPM-TPC_TSR

# Differences in tri-trophic community responses to temperature-dependent vital rates, thermal niche mismatches and temperature-size rule

## Overview

This repository gathers the codes and generated data used in the study: 

Dijoux, S., Smalås, A., Primicerio, R., and Boukal, D. S. (2024). _ Differences in Tri-trophic community response to temperature-dependent vital rates, thermal niche mismatches and temperature-size rule_.
Manuscript submitted for publication.

Preprint version available on: [Authorea](https://www.authorea.com/users/381268/articles/1212671-tri-trophic-community-responses-to-temperature-dependent-vital-rates-thermal-niche-mismatches-and-temperature-size-rule)

Cite the code: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13382965.svg)](https://doi.org/10.5281/zenodo.13382965)

## Authors & Contact details:  
Samuel Dijoux (dijoux.samuel@gmail.com),  
Aslak Smalås (aslak.smalas@skandnat.no),  
Raul Primicerio (raul.primicerio@uit.no),  
David S. Boukal (dboukal@prf.jcu.cz).

Responsible of the repository: Samuel Dijoux.

## Cite the code: [![DOI](https://zenodo.org/badge/788494888.svg)](https://zenodo.org/doi/10.5281/zenodo.10993083)

## Brief summary of the study

This study investigates the direct influences of temperature
on the dynamic and structure of size-structured tri-trophic chain. We developed a physiologically-structured population model that implement temperature-dependencies in species vital rates (through temperature performance curves) and body sizes (temperature-size rule). Through multiple scenarios we investigate separate and joint effects of these actions on community structure along environmental gradients (temperature and habitat productivity) when applied in a single species (consumer or top predator) or multiple species at once. We then extend the model to investigate how thermal mismatches between consumer and predator further alter the community structure along these gradients.

### Layout
The present repository contains the models and scripts to run the simulations and to produce the figures in main text and supplementary information:

* **PSPM_TPC-TSR_PSPMequi_simulations.R**: Script to reproduce the simulations of the study.
* **PSPM_TPC-TSR_PSPMequi_data.RData**: R workspace containing the data generated by the simulations.
* **PSPM_TPC-TSR_FIG.R**: Script to reproduce the figures from the main study and supplementary information.
* **The models**
	* **PSPM_TPC-TSR_model.R** : tri-trophic chain model with physiologically-structured consumer population to implement separately and joint influences of temperacture performance curves and temperature-size rule in consumer traits, predator traits, or both consumer and predator traits.
	* **Demo_model.R** : model for demographic analysis of consumer life histories.
	* **PSPM_TPC-TSR_Thermal-Mismatch_model.R** : model extension of the PSPM implementing mismatch between consumer and predator thermal niches.

### Software requirements
* R version 4.3.2
* (optionnal) Rstudio
* Required R packages: PSPManalysis

