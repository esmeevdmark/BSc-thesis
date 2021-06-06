# Bachelor Thesis, Future Planet Studies
# Esmée van der Mark (12393894)
# Analysing the results of the extendedTemperatureModel.py
# 5/06/2021

# importing libraries 
import pandas as pd  
import matplotlib.pyplot as plt
import numpy as np

########################################################################################################################################################
###                                                  reading the CSV result files and saving them in this Python file                                ###
########################################################################################################################################################

# temperature bifurcation plots, low and high Fs, from low to high T
low = pd.read_csv('results_low_clupeid_mortality2', index_col = 0)
high = pd.read_csv('results_high_clupeid_mortality2', index_col = 0)

# temperature bifurcation plots, low and high Fs, from high to low T (reverse)
high_reverse = pd.read_csv('results_high_clupeid_mortality_reverse', index_col = 0)
low_reverse = pd.read_csv('results_low_clupeid_mortality_reverse', index_col = 0)

# Fc bifurcation plot (default, c = 0.005, Ermax = -0.43), both low and high Fs
low_clupeid_mortality_temperature1 = pd.read_csv('results_low_clupeid_mortality_temp1', index_col = 0)
low_clupeid_mortality_temperature2 = pd.read_csv('results_low_clupeid_mortality_temp2', index_col = 0)
low_clupeid_mortality_temperature4 = pd.read_csv('results_low_clupeid_mortality_temp3', index_col = 0)
low_clupeid_mortality_temperature3 = pd.read_csv('results_low_clupeid_mortality_temp4', index_col = 0)

high_clupeid_mortality_temperature1 = pd.read_csv('results_high_clupeid_mortality_temp1', index_col = 0)
high_clupeid_mortality_temperature2 = pd.read_csv('results_high_clupeid_mortality_temp2', index_col = 0)
high_clupeid_mortality_temperature4 = pd.read_csv('results_high_clupeid_mortality_temp3', index_col = 0)
high_clupeid_mortality_temperature3 = pd.read_csv('results_high_clupeid_mortality_temp4', index_col = 0)

# Fc bifurcation plot (c = 0, ERmax = -0.43)
lowMortality_temp1_c0 = pd.read_csv('results_low_clupeid_mortality_temp1_c0', index_col = 0)
lowMortality_temp2_c0 = pd.read_csv('results_low_clupeid_mortality_temp2_c0', index_col = 0)
lowMortality_temp3_c0 = pd.read_csv('results_low_clupeid_mortality_temp3_c0', index_col = 0)
lowMortality_temp4_c0 = pd.read_csv('results_low_clupeid_mortality_temp4_c0', index_col = 0)

highMortality_temp1_c0 = pd.read_csv('results_high_clupeid_mortality_temp1_c0', index_col = 0)
highMortality_temp2_c0 = pd.read_csv('results_high_clupeid_mortality_temp2_c0', index_col = 0)
highMortality_temp3_c0 = pd.read_csv('results_high_clupeid_mortality_temp3_c0', index_col = 0)
highMortality_temp4_c0 = pd.read_csv('results_high_clupeid_mortality_temp4_c0', index_col = 0)

# Fc bifurcation plot (Ermax = 0, c = 0.005)
lowMortality_temp1_Ermax0 = pd.read_csv('results_low_clupeid_mortality_temp1_Ermax0_2', index_col = 0)
lowMortality_temp2_Ermax0 = pd.read_csv('results_low_clupeid_mortality_temp2_Ermax0_2', index_col = 0)
lowMortality_temp3_Ermax0 = pd.read_csv('results_low_clupeid_mortality_temp3_Ermax0_2', index_col = 0)
lowMortality_temp4_Ermax0 = pd.read_csv('results_low_clupeid_mortality_temp4_Ermax0_2', index_col = 0)

highMortality_temp1_Ermax0 = pd.read_csv('results_high_clupeid_mortality_temp1_Ermax0_2', index_col = 0)
highMortality_temp2_Ermax0 = pd.read_csv('results_high_clupeid_mortality_temp2_Ermax0_2', index_col = 0)
highMortality_temp3_Ermax0 = pd.read_csv('results_high_clupeid_mortality_temp3_Ermax0_2', index_col = 0)
highMortality_temp4_Ermax0 = pd.read_csv('results_high_clupeid_mortality_temp4_Ermax0_2', index_col = 0)

# Mean annual yield (cod and clupeid) (Ermax = -0.43 , c = 0.005 )
yield_low_clup = pd.read_csv('results yield low clup')
yield_high_clup = pd.read_csv('results yield high clup')
yield_low_cod = pd.read_csv('results yield low cod')
yield_high_cod = pd.read_csv('results yield high cod')

###################################################################################################################################################
##                                                    Visualizing the results                                                                    ##
###################################################################################################################################################

# temperature arrays and Fc arrays to use for the plots for X-axis
temperature_range_reverse = np.linspace(20,6,160)
temperature_range = np.linspace(6,20,160)
codMortality = np.linspace(0, 1.4, 140)

############ Figure 3 in thesis: Subplot with Fc as bifurcation parameter (default parameter values) #################################

fig, ax  = plt.subplots(nrows = 4, ncols = 2, figsize = (6,4),
                                sharex='col', sharey= 'row')

ax[0,0].plot(codMortality, low_clupeid_mortality_temperature1['Cod adult'], label = '6.5 degrees C ' )
ax[0,0].plot(codMortality, low_clupeid_mortality_temperature2['Cod adult'],  label = '7.7 degrees C')
ax[0,0].plot(codMortality, low_clupeid_mortality_temperature3['Cod adult'],  label = '10 degrees C')
ax[0,0].plot(codMortality, low_clupeid_mortality_temperature4['Cod adult'],  label = '12 degrees C' )
ax[0,0].set_ylabel('Cod adult \n biomass ' r'$(kg/m^3)$', fontsize = 14 )
ax[0,0].set_title('Low-level clupeid fishing', fontsize = 16)

ax[0,1].plot(codMortality, high_clupeid_mortality_temperature1['Cod adult'],  label = '6.5 \u00b0C ')
ax[0,1].plot(codMortality, high_clupeid_mortality_temperature2['Cod adult'],  label = '7.7 \u00b0C' )
ax[0,1].plot(codMortality, high_clupeid_mortality_temperature3['Cod adult'],  label = '10 \u00b0C' )
ax[0,1].plot(codMortality, high_clupeid_mortality_temperature4['Cod adult'],  label = '12 \u00b0C' )
ax[0,1].set_title('High-level clupeid fishing', fontsize = 16)
ax[0,1].legend()

ax[1,0].plot(codMortality, low_clupeid_mortality_temperature1['Clupeid adult'],  label = '6.5 degrees C' )
ax[1,0].plot(codMortality, low_clupeid_mortality_temperature2['Clupeid adult'],  label = '7.7 degrees C' )
ax[1,0].plot(codMortality, low_clupeid_mortality_temperature3['Clupeid adult'],  label = '10 degrees C' )
ax[1,0].plot(codMortality, low_clupeid_mortality_temperature4['Clupeid adult'],  label = '12 degrees C' )
ax[1,0].set_ylabel('Clupeid adult \n biomass 'r'$(kg/m^3)$', fontsize = 14)

ax[1,1].plot(codMortality, high_clupeid_mortality_temperature1['Clupeid adult'],  label = '6.5 degrees C' )
ax[1,1].plot(codMortality, high_clupeid_mortality_temperature2['Clupeid adult'],  label = '7.7 degrees C' )
ax[1,1].plot(codMortality, high_clupeid_mortality_temperature3['Clupeid adult'],  label = '10 degrees C' )
ax[1,1].plot(codMortality, high_clupeid_mortality_temperature4['Clupeid adult'],  label = '12 degrees C' )

ax[2,0].plot(codMortality, low_clupeid_mortality_temperature1['Sj'],  label = '6.5 degrees C' )
ax[2,0].plot(codMortality, low_clupeid_mortality_temperature2['Sj'],  label = '7.7 degrees C' )
ax[2,0].plot(codMortality, low_clupeid_mortality_temperature3['Sj'],  label = '10 degrees C' )
ax[2,0].plot(codMortality, low_clupeid_mortality_temperature4['Sj'],  label = '12 degrees C' )
ax[2,0].set_ylabel('Clupeid juvenile \n biomass ' r'$(kg/m^3)$', fontsize = 14)

ax[2,1].plot(codMortality, high_clupeid_mortality_temperature1['Sj'],  label = '6.5 degrees C' )
ax[2,1].plot(codMortality, high_clupeid_mortality_temperature2['Sj'],  label = '7.7 degrees C' )
ax[2,1].plot(codMortality, high_clupeid_mortality_temperature3['Sj'],  label = '10 degrees C' )
ax[2,1].plot(codMortality, high_clupeid_mortality_temperature4['Sj'],  label = '12 degrees C' )

ax[3,0].plot(codMortality, low_clupeid_mortality_temperature1['Rs'],  label = '6.5 degrees C' )
ax[3,0].plot(codMortality, low_clupeid_mortality_temperature2['Rs'],  label = '7.7 degrees C' )
ax[3,0].plot(codMortality, low_clupeid_mortality_temperature3['Rs'],  label = '10 degrees C' )
ax[3,0].plot(codMortality, low_clupeid_mortality_temperature4['Rs'],  label = '12 degrees C' )
ax[3,0].set_ylabel('Resource for clupeids \n biomass ' r'$(kg/m^3)$', fontsize = 14)

ax[3,1].plot(codMortality, high_clupeid_mortality_temperature1['Rs'],  label = '6.5 degrees C' )
ax[3,1].plot(codMortality, high_clupeid_mortality_temperature2['Rs'],  label = '7.7 degrees C' )
ax[3,1].plot(codMortality, high_clupeid_mortality_temperature3['Rs'],  label = '10 degrees C' )
ax[3,1].plot(codMortality, high_clupeid_mortality_temperature4['Rs'],  label = '12 degrees C' )

fig.text(0.5, 0.04, 'Cod fishing mortality' r' $ F_{c}$' r' $(year^{-1})$', ha = 'center', fontsize = 16)

################## Figure 4 in thesis: biomass with T as bifurcation parameter (default parameter values) ##################################

fig, ax  = plt.subplots(nrows = 4, ncols = 2, figsize = (6,4),
                                sharex='col', sharey= 'row')

ax[0,0].plot(temperature_range, low['Cod adult'], label = 'low to high T' )
ax[0,0].plot(temperature_range_reverse, low_reverse['Cod adult'], '--' ,label = 'high to low T')
ax[0,0].set_ylabel('Cod adult \n biomass ' r'$(kg/m^3)$', fontsize = 14 )
ax[0,0].set_title('Low-level clupeid fishing', fontsize = 16)

ax[0,1].plot(temperature_range, high['Cod adult'],  label = 'low to high T')
ax[0,1].plot(temperature_range_reverse, high_reverse['Cod adult'], '--', label = 'high to low T' )
ax[0,1].set_title('High-level clupeid fishing', fontsize = 16)
ax[0,1].legend()

ax[1,0].plot(temperature_range, low['Clupeid adult'],  label = 'low to high T' )
ax[1,0].plot(temperature_range_reverse, low_reverse['Clupeid adult'], '--' , label = 'high to low T' )
ax[1,0].set_ylabel('Clupeid adult \n biomass 'r'$(kg/m^3)$', fontsize = 14)

ax[1,1].plot(temperature_range, high['Clupeid adult'],  label = 'low to high T' )
ax[1,1].plot(temperature_range_reverse, high_reverse['Clupeid adult'], '--', label = 'high to low T' )

ax[2,0].plot(temperature_range, low['Sj'],  label = 'low to high T' )
ax[2,0].plot(temperature_range_reverse, low_reverse['Sj'], '--' ,label = 'high to low T' )
ax[2,0].set_ylabel('Clupeid juvenile \n biomass ' r'$(kg/m^3)$', fontsize = 14)

ax[2,1].plot(temperature_range, high['Sj'],  label = '6.5 degrees C' )
ax[2,1].plot(temperature_range_reverse, high_reverse['Sj'], '--', label = '7.7 degrees C' )

ax[3,0].plot(temperature_range, low['Rs'],  label = 'low to high T' )
ax[3,0].plot(temperature_range_reverse, low_reverse['Rs'], '--' , label = 'high to low T' )
ax[3,0].set_ylabel('Resource for clupeids \n biomass ' r'$(kg/m^3)$', fontsize = 14)

ax[3,1].plot(temperature_range, high['Rs'],  label = 'low to high T' )
ax[3,1].plot(temperature_range_reverse, high_reverse['Rs'], '--', label = 'high to low T' )

fig.text(0.5, 0.04, 'Temperature (in \u00b0C)', ha = 'center', fontsize = 16)

################## Figure 5 in thesis : Cod mean annual yield (with default parameter values) ####################

fig, ax  = plt.subplots(nrows = 1, ncols = 2, figsize = (6,4),
                                sharex='col', sharey= 'row')

ax[0].plot(codMortality, yield_low_cod['temp1'],  label = '6.5 degrees C')
ax[0].plot(codMortality, yield_low_cod['temp2'],  label = '7.7 degrees C')
ax[0].plot(codMortality, yield_low_cod['temp3'],  label = '10 degrees C')
ax[0].plot(codMortality, yield_low_cod['temp4'],  label = '12 degrees C')
ax[0].set_title('Low-level clupeid \n fishing ', fontsize = 14)

ax[1].plot(codMortality, yield_high_cod['temp1'], label = '6.5 \u00b0C')
ax[1].plot(codMortality, yield_high_cod['temp2'],  label = '7.7 \u00b0C')
ax[1].plot(codMortality, yield_high_cod['temp3'],  label = '10 \u00b0C')
ax[1].plot(codMortality, yield_high_cod['temp4'],  label = '12 \u00b0C')
ax[1].set_title('High-level clupeid \n fishing ', fontsize = 14)
ax[1].legend()

fig.text(0.5, 0.04, 'Cod fishing mortality' r' $ F_{c}$' r' $(year^{-1})$', ha = 'center', fontsize = 16)
fig.text(0.04, 0.5, 'Cod mean annual yield ' r'$(kg/m^3)$', va = 'center', rotation = 'vertical', fontsize = 16)


################## Figure 6 in thesis: Biomass with Fc as bifurcation parameter (with ERmax = 0)

fig, ax  = plt.subplots(nrows = 4, ncols = 2, figsize = (6,4),
                                sharex='col', sharey= 'row')

ax[0,0].plot(codMortality, lowMortality_temp1_Ermax0['Cod adult'],  label = '6.5 degrees C' )
ax[0,0].plot(codMortality, lowMortality_temp2_Ermax0['Cod adult'],  label = '7.7 degrees C')
ax[0,0].plot(codMortality, lowMortality_temp3_Ermax0['Cod adult'],  label = '10 degrees C')
ax[0,0].plot(codMortality, lowMortality_temp4_Ermax0['Cod adult'],  label = '12 degrees C' )
ax[0,0].set_ylabel('Cod adult \n biomass ' r'$(kg/m^3)$', fontsize = 14 )
ax[0,0].set_title('Low-level clupeid fishing', fontsize = 16)

ax[0,1].plot(codMortality, highMortality_temp1_Ermax0['Cod adult'], label = '6.5 \u00b0C')
ax[0,1].plot(codMortality, highMortality_temp2_Ermax0['Cod adult'],  label = '7.7 \u00b0C' )
ax[0,1].plot(codMortality, highMortality_temp3_Ermax0['Cod adult'],  label = '10 \u00b0C' )
ax[0,1].plot(codMortality, highMortality_temp4_Ermax0['Cod adult'],  label = '12 \u00b0C' )
ax[0,1].set_title('High-level clupeid fishing', fontsize = 16)
ax[0,1].legend()

ax[1,0].plot(codMortality, lowMortality_temp1_Ermax0['Clupeid adult'],  label = '6.5 degrees C' )
ax[1,0].plot(codMortality, lowMortality_temp2_Ermax0['Clupeid adult'],  label = '7.7 degrees C' )
ax[1,0].plot(codMortality, lowMortality_temp3_Ermax0['Clupeid adult'],  label = '10 degrees C' )
ax[1,0].plot(codMortality, lowMortality_temp4_Ermax0['Clupeid adult'],  label = '12 degrees C' )
ax[1,0].set_ylabel('Clupeid adult \n biomass 'r'$(kg/m^3)$', fontsize = 14)

ax[1,1].plot(codMortality, highMortality_temp1_Ermax0['Clupeid adult'],  label = '6.5 degrees C' )
ax[1,1].plot(codMortality, highMortality_temp2_Ermax0['Clupeid adult'],  label = '7.7 degrees C' )
ax[1,1].plot(codMortality, highMortality_temp3_Ermax0['Clupeid adult'],  label = '10 degrees C' )
ax[1,1].plot(codMortality, highMortality_temp4_Ermax0['Clupeid adult'],  label = '12 degrees C' )

ax[2,0].plot(codMortality, lowMortality_temp1_Ermax0['Sj'],  label = '6.5 degrees C' )
ax[2,0].plot(codMortality, lowMortality_temp2_Ermax0['Sj'],  label = '7.7 degrees C' )
ax[2,0].plot(codMortality, lowMortality_temp3_Ermax0['Sj'],  label = '10 degrees C' )
ax[2,0].plot(codMortality, lowMortality_temp4_Ermax0['Sj'],  label = '12 degrees C' )
ax[2,0].set_ylabel('Clupeid juvenile \n biomass ' r'$(kg/m^3)$', fontsize = 14)

ax[2,1].plot(codMortality, highMortality_temp1_Ermax0['Sj'],  label = '6.5 degrees C' )
ax[2,1].plot(codMortality, highMortality_temp2_Ermax0['Sj'],  label = '7.7 degrees C' )
ax[2,1].plot(codMortality, highMortality_temp3_Ermax0['Sj'],  label = '10 degrees C' )
ax[2,1].plot(codMortality, highMortality_temp4_Ermax0['Sj'],  label = '12 degrees C' )

ax[3,0].plot(codMortality, lowMortality_temp1_Ermax0['Rs'],  label = '6.5 degrees C' )
ax[3,0].plot(codMortality, lowMortality_temp2_Ermax0['Rs'],  label = '7.7 degrees C' )
ax[3,0].plot(codMortality, lowMortality_temp3_Ermax0['Rs'],  label = '10 degrees C' )
ax[3,0].plot(codMortality, lowMortality_temp4_Ermax0['Rs'],  label = '12 degrees C' )
ax[3,0].set_ylabel('Resource for clupeids \n biomass ' r'$(kg/m^3)$', fontsize = 14)

ax[3,1].plot(codMortality, highMortality_temp1_Ermax0['Rs'],  label = '6.5 degrees C' )
ax[3,1].plot(codMortality, highMortality_temp2_Ermax0['Rs'],  label = '7.7 degrees C' )
ax[3,1].plot(codMortality, highMortality_temp3_Ermax0['Rs'],  label = '10 degrees C' )
ax[3,1].plot(codMortality, highMortality_temp4_Ermax0['Rs'],  label = '12 degrees C' )

fig.text(0.5, 0.04, 'Cod fishing mortality' r' $ F_{c}$' r' $(year^{-1})$', ha = 'center', fontsize = 16)

######################## Figure S2 in thesis : Clupeid mean annual yield (with default parameter values) ####################################

fig, ax  = plt.subplots(nrows = 1, ncols = 2, figsize = (6,4),
                                sharex='col', sharey= 'row')

ax[0].plot(codMortality, yield_low_clup['temp1'],  label = '6.5 degrees C')
ax[0].plot(codMortality, yield_low_clup['temp2'],  label = '7.7 degrees C')
ax[0].plot(codMortality, yield_low_clup['temp3'],  label = '10 degrees C')
ax[0].plot(codMortality, yield_low_clup['temp4'],  label = '12 degrees C')
ax[0].set_title('Low-level clupeid \n fishing ', fontsize = 14 )

ax[1].plot(codMortality, yield_high_clup['temp1'],  label = '6.5 \u00b0C')
ax[1].plot(codMortality, yield_high_clup['temp2'],  label = '7.7 \u00b0C')
ax[1].plot(codMortality, yield_high_clup['temp3'],  label = '10 \u00b0C')
ax[1].plot(codMortality, yield_high_clup['temp4'],  label = '12 \u00b0C')
ax[1].set_title('High-level clupeid \n fishing ', fontsize = 14) 
ax[1].legend()

fig.text(0.5, 0.04, 'Cod fishing mortality' r' $ F_{c}$' r' $(year^{-1})$', ha = 'center', fontsize = 16)
fig.text(0.04, 0.5, 'Clupeid mean annual yield ' r'$(kg/m^3)$', va = 'center', rotation = 'vertical', fontsize = 16)


###################### Figure S3 in thesis: Biomass with Fc as bifurcation parameter (c = 0.0) ################################################

fig, ax  = plt.subplots(nrows = 1, ncols = 2, figsize = (6,4),
                                sharex='col', sharey= 'row')

ax[0].plot(codMortality, low_clupeid_mortality_temperature1['Cod adult'], 'tab:blue',label = '6.5 degrees C ' )
ax[0].plot(codMortality, low_clupeid_mortality_temperature2['Cod adult'], 'tab:orange' ,label = '7.7 degrees C')
ax[0].plot(codMortality, low_clupeid_mortality_temperature3['Cod adult'], 'tab:green' ,label = '10 degrees C')
ax[0].plot(codMortality, low_clupeid_mortality_temperature4['Cod adult'], 'tab:red', label = '12 degrees C' )

ax[1].plot(codMortality, high_clupeid_mortality_temperature1['Cod adult'], 'tab:blue', label = '6.5 \u00b0C , c = 0.005')
ax[1].plot(codMortality, high_clupeid_mortality_temperature2['Cod adult'], 'tab:orange', label = '7.7 \u00b0C, c = 0.005' )
ax[1].plot(codMortality, high_clupeid_mortality_temperature3['Cod adult'], 'tab:green' ,  label = '10 \u00b0C, c = 0.005' )
ax[1].plot(codMortality, high_clupeid_mortality_temperature4['Cod adult'], 'tab:red', label = '12 \u00b0C, c = 0.005' )

ax[0].plot(codMortality, lowMortality_temp1_c0['Cod adult'], '--' ,color = 'tab:blue', label = '6.5 degrees C ' )
ax[0].plot(codMortality, lowMortality_temp2_c0['Cod adult'], '--', color ='tab:orange' ,label = '7.7 degrees C')
ax[0].plot(codMortality, lowMortality_temp3_c0['Cod adult'], '--', color ='tab:green' , label = '10 degrees C')
ax[0].plot(codMortality, lowMortality_temp4_c0['Cod adult'], '--', color ='tab:red', label = '12 degrees C' )
ax[0].set_ylabel('Cod adult \n biomass ' r'$(kg/m^3)$', fontsize = 14 )
ax[0].set_title('Low-level clupeid fishing', fontsize = 16)

ax[1].plot(codMortality, highMortality_temp1_c0['Cod adult'], '--', color ='tab:blue', label = '6.5 \u00b0C, c = 0.0 ')
ax[1].plot(codMortality, highMortality_temp2_c0['Cod adult'], '--', color ='tab:orange', label = '7.7 \u00b0C, c = 0.0' )
ax[1].plot(codMortality, highMortality_temp3_c0['Cod adult'], '--', color ='tab:green', label = '10 \u00b0C, c = 0.0' )
ax[1].plot(codMortality, highMortality_temp4_c0['Cod adult'], '--', color ='tab:red', label = '12 \u00b0C, c = 0.0' )
ax[1].set_title('High-level clupeid fishing', fontsize = 16)
ax[1].legend()

fig.text(0.5, 0.04, 'Cod fishing mortality' r' $ F_{c}$' r' $(year^{-1})$', ha = 'center', fontsize = 16)

plt.show()