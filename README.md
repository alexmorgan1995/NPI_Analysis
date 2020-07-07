# Modelling Social Distancing Measures for COVID-19 Epidemic Control

Repository for R code used to model social distancing measures **WORK IN PROGRESS**.

## Running the Code
### R 
R code should run without further modification (after installing pre-requisite packages). 
The chosen directory where plots are saved can be altered by changing the `setwd()` function at the top of each script. Plots are automatically saved in the chosen working directory if each script is run. 

The analysis has been split into R files relating to the single internvetion, multi-intervention and figures for the supplementary material:

## Navigating the Repository 
### R
R code can be found in the `Enhanced Shielding` folder and is organised according to the figures found in the main text and the supplementary material **WORK IN PROGRESS**. Current R code is available for:

* **Single Intervention Analysis** - Modelling the effect of a single SDM intervention - `single_5_Scenario_Run.R`.
	* Baseline trajectory plots.
	* Single parameter sensitivity analyses.
 	* Multi-parameter sensitivity heatmap analyses. 
 
* **Multiple/Double Intervention Analysis** - Modelling the effect of a two SDM interventions - `multi_5_Scenario_Run.R`.
	* Multi-parameter sensitivity heatmap analyses (Trigger date/R0/Length). 
 
* **Supplementary Material for Single and Multiple Intervention Analysis** - `multi_5_Scenario_Run_supp.R`.
	* Multi-parameter sensitivity heatmap analyses - Exploring higher dimension parameter space.
	* Exploratory trajectory plots.

## Programs and Packages Used
COVID-19 modelling code was implemented using R (3.6.2) and R-Studio. ODEs were solved using the `desolve` (1.27.1) package in R and `odeint` in C++. Plotting in R was carried out using the `ggplot2` package (3.3.0). Dataframe manipulation was performed using `reshape2` (1.4.4). Finalised plot output was performed using `ggpubr` (0.2.4), `ggarrange` () and `Cairo` (1.5-10) packages. 

## Acknowledgements 
All analysis were performed by members and affiliates of Epigroup, University of Edinburgh: 
http://www.epigroup.biology.ed.ac.uk/
