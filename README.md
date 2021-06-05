# BSc thesis
This is the repository of the BSC thesis of Esmée van der Mark, Future Planet Studies, University of Amsterdam (2021). The project consists of extending a multitrophic, size-structured, bioenergetics model of the fish community in the Baltic Sea with temperature dependent vital rates and parameters in the model. 

This repository consists of multiple files. The most important one is the 'ExtendedTemperatureModel.py' file, which includes the Python code of my extended model. This file can be used to reproduce my figures which I used in my thesis as well. These figures are also available in this repository and they have the same numbering as in my thesis (Fig.3 - Fig.S3). 

My code can be used in several ways. At the moment, the code is set with the default values of ERmax and c (Ermax = -0.43 and c = 0.005). If you want to recreate the results of Ermax = 0.0 or c = 0.0, you can change these values in the 'parameter and constants' section at the top op the file. 

I have added a variable, called 'biomass_or_yield' , which you can also find in the 'parameter and constants' section, which you can change, ranging from numbers 0-3. This variable gives you the option which calculations you want to perform. These different calculations, for both low and high clupeid fishing mortalities (Fs) consists of:

0 : Calculating biomass levels, with temperature as bifurcation parameter (integrated from low to high temperature), with high Fc value and default parameter values

1 : Calculating biomass levels, with temperature as bifurcation parameter (integrated from high to low temperature), with high Fc value and default parameter values

2 : Calculating biomass levels, with cod fishing mortality (Fc) as bifurcation parameter, default parameter values as well

3 : Calculating mean annual yield (cod and clupeids), with Fc as bifurcation parameter, default parameter values

Setting this variable to one of these numbers will thus have as effect that only the calculations will be done which you want to calculate. These calculations can be found at the 'simulation' section of the file. For every simulation, I store the results in arrays, so I can at the end store them in CSV files, which can be seen in the 'storing the data in CSV files section' at the bottom of the file. 

These CSV files are then used in another Python file, called 'analysingResults.py', in which I have made the figures I used in my thesis. All the CSV files I have used for these figures are also available in this repository. There are 32 different CSV files and below I explain which I used for the different figures:

Fig3:

'results_low_clupeid_mortality_temp1'
'results_low_clupeid_mortality_temp2'
'results_low_clupeid_mortality_temp3'
'results_low_clupeid_mortality_temp4'
'results_high_clupeid_mortality_temp1'
'results_high_clupeid_mortality_temp2'
'results_high_clupeid_mortality_temp3'
'results_high_clupeid_mortality_temp4'

Fig4:

'results_low_clupeid_mortality2' 
'results_high_clupeid_mortality2'
'results_high_clupeid_mortality_reverse'
'results_low_clupeid_mortality_reverse'

Fig5:

'results yield low cod'
'results yield high cod'

Fig6:

'results_low_clupeid_mortality_temp1_Ermax0_2'
'results_low_clupeid_mortality_temp2_Ermax0_2'
'results_low_clupeid_mortality_temp3_Ermax0_2'
'results_low_clupeid_mortality_temp4_Ermax0_2'
'results_high_clupeid_mortality_temp1_Ermax0_2'
'results_high_clupeid_mortality_temp2_Ermax0_2'
'results_high_clupeid_mortality_temp3_Ermax0_2'
'results_high_clupeid_mortality_temp4_Ermax0_2'

FigS2:

'results yield low clup'
'results yield high clup'

FigS3:

'results_low_clupeid_mortality_temp1'
'results_low_clupeid_mortality_temp2'
'results_low_clupeid_mortality_temp3'
'results_low_clupeid_mortality_temp4'
'results_high_clupeid_mortality_temp1'
'results_high_clupeid_mortality_temp2'
'results_high_clupeid_mortality_temp3'
'results_high_clupeid_mortality_temp4'

'results_low_clupeid_mortality_temp1_c0'
'results_low_clupeid_mortality_temp2_c0'
'results_low_clupeid_mortality_temp3_c0'
'results_low_clupeid_mortality_temp4_c0'
'results_high_clupeid_mortality_temp1_c0'
'results_high_clupeid_mortality_temp2_c0'
'results_high_clupeid_mortality_temp3_c0'
'results_high_clupeid_mortality_temp4_c0'

