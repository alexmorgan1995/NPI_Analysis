# Optimising time-limited non-pharmaceutical interventions for COVID-19 outbreak control.

Repository for R code used to model non-pharmaceutical interventions on a simulated UK COVID-19 outbreak. Modelled in the context of optimising NPIs.

## Running the Code
### R 
R code should run without further major modification (after installing pre-requisite packages). 

**Important** - Note that it is necessary to change `setwd` at the top of every script according to the user's preference.

The chosen directory where plots are saved can be altered by changing the `setwd()` function at the top of each script. Plots are automatically saved in the chosen working directory if each script is run.

**Important** - Note that supplementary code requires a folder labelled `supplementary` to be created in the working directory (defined by `setwd()`). Alternatively `setwd` can be altered for supplementary analysis and customised by the user. 

The analysis has been split into R files relating to the single intervention, multi-intervention and figures for the supplementary material:

## Navigating the Repository 
### R
R code is organised according to the figures found in the main text and the supplementary material. Current R code is available for:

* **Single Intervention Analysis** - Modelling the effect of a single SDM intervention - `multi_5_Scenario_Run_FINAL.R`.
	* Baseline trajectory plots - Figure 1A.
		*Same code was used to plot Figure S1 (keeping dt constant and altering cmin). 
	* Single parameter sensitivity analyses - Figure 1B.
 	* Multi-parameter sensitivity heatmap analysis - Figure 2. 
 
* **Multiple/Double Intervention Analysis** - Modelling the effect of a two SDM interventions - `multi_5_Scenario_Run_FINAL.R`.
	* Baseline double intervention trajectory plots - Figure S5.
		*Same code was used to plot Figure S2 (keeping dt1/dt2 constant and altering cmin1/cmin2). 
	* Multi-parameter sensitivity heatmap analyses (Trigger point 1(t_p1)/ Trigger point 2(t_p1) and Magnitude of intervention 1 (c_min1)/intervention 2 (c_min1)) - Figure 3. 
 
* **Supplementary Material for Single and Multiple Intervention Analysis** - `SUPPLEMENTARY_5_Scenario.R`.
	* Beta plots over time for the 5 intervention scenarios - Supplementary Table .
	* Single-parameter sensitivity heatmap analyses for all 5 scenarios - Exploring c_min - Figure S3-4
	* Multi-parameter sensitivity heatmap analyses tp1/tp2 and cmin1/cmin2 for all 5 scenarios - Exploring d_t - Figure S6-10 and Figure S11-15.

* **Supplementary Material for optimising NPIs + sustainable control measures** - `SUPPLEMENTARY_conreduc.R`.
	* Exploring optimisation of an intervention, in the context of a later, constant reduction to beta(t) - representative of an introduction of sustainable intervention measures indefinitely - Figure S16.

* **Supplementary Material for SEIR model** - `SUPPLEMENTARY_SEIR.R`.
	* Exploring the effect of an E compartment on model dynamics, both trajectory plots and the tp1/tp2 sensitivity analysis - Figure S17. 

* **Supplementary Material for SIRS model** - `SUPPLEMENTARY_SIRS.R`.
	* Exploring the effect of waning immunity on model dynamics, both the trajectory plot and the Scenario 1 sensitivity analysis (tp/dt) - Figure S18.  
	* We explored 3 average durations of immunity - 3, 6 and 12 months. 

## Programs and Packages Used
COVID-19 modelling code was implemented using R (3.6.2) and R-Studio. ODEs were solved using the `desolve` (1.27.1) package in R. Plotting in R was carried out using the `ggplot2` package (3.3.0). Dataframe manipulation was performed using `reshape2` (1.4.4). Finalised plot output was performed using `ggpubr` (0.2.4), `ggarrange` () and `Cairo` (1.5-10) packages. 

## Acknowledgements 
All analysis were performed by members and affiliates of Epigroup, University of Edinburgh: 
http://www.epigroup.biology.ed.ac.uk/
