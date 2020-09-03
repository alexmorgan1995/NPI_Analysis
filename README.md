# Modelling Social Distancing Measures for COVID-19 Epidemic Control

Repository for R code used to model social distancing measures **WORK IN PROGRESS**.

## Running the Code
### R 
R code should run without further modification (after installing pre-requisite packages). 
The chosen directory where plots are saved can be altered by changing the `setwd()` function at the top of each script. Plots are automatically saved in the chosen working directory if each script is run. 

The analysis has been split into R files relating to the single intervention, multi-intervention and figures for the supplementary material:

## Navigating the Repository 
### R
R code can be found in the `Enhanced Shielding` folder and is organised according to the figures found in the main text and the supplementary material **WORK IN PROGRESS**. Current R code is available for:

* **Single Intervention Analysis** - Modelling the effect of a single SDM intervention - `multi_5_Scenario_Run_FINAL.R`.
	* Baseline trajectory plots - Figure 1A.
	* Single parameter sensitivity analyses - Figure 1B.
 	* Multi-parameter sensitivity heatmap analyses - Figure 2. 
 
* **Multiple/Double Intervention Analysis** - Modelling the effect of a two SDM interventions - `multi_5_Scenario_Run_FINAL.R`.
	* Multi-parameter sensitivity heatmap analyses (Trigger date (t_p)/magnitude (c_min)/length (d_t)) - Figure 3. 
 
* **Supplementary Material for Single and Multiple Intervention Analysis** - `SUPPLEMENTARY_5_Scenario.R`.
	* Single-parameter sensitivity heatmap analyses - Exploring c_min.
	* Multi-parameter sensitivity heatmap analyses - Exploring d_t.

* **Supplementary Material for SEIR model example and optimising NPIs and sustainable control measures ** - `SUPPLEMENTARY_5_Scenario_conreduc_SEIR.R`.
	* Exploring the effect of an E compartment on model dynamics.
	* Exploring optimising intervention 1, followed by a constant reduction to beta(t) - representative of a later introduction of sustainable intervention measures.

## Programs and Packages Used
COVID-19 modelling code was implemented using R (3.6.2) and R-Studio. ODEs were solved using the `desolve` (1.27.1) package in R. Plotting in R was carried out using the `ggplot2` package (3.3.0). Dataframe manipulation was performed using `reshape2` (1.4.4). Finalised plot output was performed using `ggpubr` (0.2.4), `ggarrange` () and `Cairo` (1.5-10) packages. 

## Acknowledgements 
All analysis were performed by members and affiliates of Epigroup, University of Edinburgh: 
http://www.epigroup.biology.ed.ac.uk/
