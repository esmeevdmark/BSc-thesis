# Bachelor Thesis, Future Planet Studies
# Esmée van der Mark (12393894)
# extending the Soudijn model with the temperature dependancy of Lindmark et al. (2019)
# 5/06/2021

# libraries
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math 
from statistics import mean
import pandas as pd  

##########################################################################################################################################################
##                                                                  parameters and constants                                                           ###    
##########################################################################################################################################################                                                        

# reference temperature in Kelvin (so 10 degrees Celsius)
T0 = 283

# Boltzmann's constant 
k = 8.617332E-05 

# activation energies of metabolism (EM), maximum ingestion rate (EI), background mortality (Emu),
# maximum resource density (ERmax) and resource turnover rate (Edelta). 
# Note that ERmax can change between -0.43 and 0
EM = 0.594
EI = 0.594
Emu = 0.45
Edelta = 0.43
ERmax = -0.43

# Half-saturation density (g Vol-1)
H = 1

# Fisheries retention of clupeid juveniles (S) and cod juveniles
RhoS = 0.5
RhoC = 0.02

# Y = duration of growing season (250 days)
Y = 250

# parameter c can change between 0 and 0.005. The parameter c represents whether there is or is not 
# interactive effect between temperature and body size
c = 0.005

# Average size of clupeid juveniles (S_J), small adult clupeids (S_A), large adult clupeids (S_B),
S_J = 3.4
S_A = 12.7
S_B = 15.0

# Average size of cod juveniles (C_J), small adult cod (C_A), large adult cod (C_B),
C_J = 18.2
C_A = 350
C_B = 832

# if you want to calculate biomass, with temperature as bifurcation parameter (from low to high T) and high Fc , choose 0
# if you want to calculate biomass, with temperature aa bifurcation parameter (from high to low T) and high Fc, choose 1
# if you want to calculate biomass, with cod fishing mortality as bifucation parameter, choose 2
# if you want to calculate yield for both low and high Fs, choose 3
biomass_or_yield = 0

#########################################################################################################################################################
##                                                                    intitial values                                                                  ##
#########################################################################################################################################################

Rs1 = 32.7816744728
Rj1 = 0.4713946441
Ra1 =  0.7175924292
Sj1 = 10.4824915096
Sa1 =  69.8208217446
Sb1 =  89.7371125270
Sga1 = 0
Sgb1 = 0
Cj1 =  6.2772314150
Ca1 = 9.4197483660
Cb1 = 0.6284148535
Cga1 = 0
Cgb1 = 0

########################################################################################################################################################
##                                                                       functions                                                                    ##
########################################################################################################################################################

def Maturation(z, v, mu):
    """
    This function calculates the maturation rate.
    Input:
        z : float - ratio of initial to final body size
        v : float - food ingestion 
        mu : float - mortality rates

    Output:
        matrate: float - maturation rate

    """
    Max_Exp = 50.0
    logz = 0.0
    tmp = 0.0
    tres = 0.0
    matrate = 0.0 

    logz = math.log(z)

    if math.fabs(logz) < 1.0E-6:
        matrate = 0.0
    else:
        tres = mu / (1.0 - Max_Exp / logz)
        if v < tres:
            matrate = 0.0
        else:
            tmp = 1.0 - mu / v
            if math.fabs(tmp) < 1.0E-6:
                matrate = v * (tmp / 2 - 1 / logz)
            else:
                matrate = (v - mu) / (1.0 - math.exp(tmp * logz))

    return matrate


def function_growing_season(t, P):
    """
    This function calculates the biomasses in the growing season and consists of 13 ODEs for 13 different fish species / stages. 
    Input:
        t : array of floats - time
        P : array of floats - initial values 

    Output:
        dRs_dt : array of floats - resource for clupeids biomass
        dRj_dt : array of floats - resource for cod juveniles biomass
        dRa_dt : array of floats - resource for cod adults biomass
        dSj_dt : array of floats - juvenile clupeids biomass
        dSa_dt : array of floats - small adult clupeids biomass
        dSb_dt : array of floats - large adult clupeids biomass
        dSga_dt : array of floats - small adult reproductive storages biomass
        dSgb_dt : array of floats - large adult reproductive storages biomass
        dCj_dt : array of floats - juvenile cod biomass
        dCa_dt : array of floats - small adult cod biomass
        dCb_dt : array of floats - large adult cod biomass
        dCga_dt : array of floats - small adult reproductive storages biomass
        dCgb_dt : array of floats - large adult reproductive storages biomass

    """
    # parameter values for different life stages of clupeids and cods, stored in dictionaries
    Sj = {'I' : 0.23 * math.exp((EI * (T-T0))/(k*T*T0)) ,'Sigma':0.3 , 'T': 0.032 *  math.exp((EM * (T-T0))/(k*T*T0)) , 'kappa': 1, 'z': 0.05, 'mu':0.001 * math.exp((Emu * (T-T0))/(k*T*T0)), 'F': RhoS*Fs}
    Sa = {'I' : 0.078 * math.exp((EI * (T-T0))/(k*T*T0))  * (S_A/S_J) ** (c*(T-T0)),'Sigma':0.3 , 'T': 0.02 * math.exp((EM * (T-T0))/(k*T*T0)) * (S_A/S_J) ** (c*(T-T0)), 'kappa': 0.8, 'z': 0.7, 'mu':0.001 * math.exp((Emu * (T-T0))/(k*T*T0)), 'F': Fs}
    Sb = {'I' : 0.078 * math.exp((EI * (T-T0))/(k*T*T0)) * (S_B/S_J) ** (c*(T-T0)), 'Sigma':0.3 , 'T': 0.02 * math.exp((EM * (T-T0))/(k*T*T0)) * (S_B/S_J) ** (c*(T-T0)), 'kappa': 0.0, 'mu': 0.001 * math.exp((Emu * (T-T0))/(k*T*T0)) , 'F': Fs}
    Cj = {'I' : 0.08 * math.exp((EI * (T-T0))/(k*T*T0)), 'Sigma':0.3 , 'T': 0.015 * math.exp((EM * (T-T0))/(k*T*T0)), 'kappa': 1, 'z': 0.003, 'mu': 0.001 * math.exp((Emu * (T-T0))/(k*T*T0)), 'F': RhoC*Fc}
    Ca = {'I' : 0.022 * math.exp((EI * (T-T0))/(k*T*T0)) * (C_A/C_J) ** (c*(T-T0)), 'Sigma':0.4 , 'T': 0.006 * math.exp((EM * (T-T0))/(k*T*T0)) * (C_A/C_J) ** (c*(T-T0)) , 'kappa': 0.8, 'z': 0.125, 'mu': 0.001 * math.exp((Emu * (T-T0))/(k*T*T0)), 'F': Fc}
    Cb = {'I' : 0.022 * math.exp((EI * (T-T0))/(k*T*T0)) * (C_B/C_J) ** (c*(T-T0)), 'Sigma':0.4 , 'T': 0.006 * math.exp((EM * (T-T0))/(k*T*T0)) * (C_B/C_J) ** (c*(T-T0)), 'kappa': 0.0, 'mu': 0.001 * math.exp((Emu * (T-T0))/(k*T*T0)),'F': Fc}

    # cod population functions 
    Ecj = 0.8*P[1] + 0.2 * P[3]
    Eca = 0.5*P[2] + 0.3 * P[3] + 0.1*(P[4] + P[6]) + 0.1*(P[5] + P[7])
    Ecb = 0.2*P[2] + 0.25 * P[3] + 0.3*(P[4] + P[6]) + 0.25*(P[5] + P[7])

    Gcj = (Cj['I'] / (H + Ecj))
    Gca = (Ca['I'] / (H + Eca))
    Gcb = (Cb['I'] / (H + Ecb))

    # clupeid population functions
    Es = P[0]
    Psj = 0.2 * Gcj * P[8] + 0.3 * Gca * P[9] + 0.25 * Gcb * P[10]
    Psa = 0.1 * Gca * P[9] + 0.3 * Gcb * P[10]
    Psb = 0.1 * Gca * P[9] + 0.25 * Gcb * P[10]
    
    # resource functions 
    Gsj = (Sj['I'] / (H + Es))
    Gsa = (Sa['I'] / (H + Es))
    Gsb = (Sb['I'] / (H + Es))

    Gs = Gsj * P[3] + Gsa * P[4] + Gsb * P[5]
    Gj = 0.8 * Gcj * P[8] 
    Ga = 0.5 * Gca * P[9] + 0.2 * Gcb * P[10]

    # resource ODEs
    dRs_dt = Delta * (Rsmax - P[0]) - Gs * P[0]
    dRj_dt = Delta * (Rjmax - P[1]) - Gj * P[1]
    dRa_dt = Delta * (Ramax - P[2]) - Ga * P[2]
    
    # clupeid functions and ODEs 
    # Sj (juvenile clupeids)
    VSj = Sj['Sigma'] * Gsj * Es - Sj['T']
    dSj = Sj['mu'] + Sj['F'] + Psj
    GammaSj = Maturation(Sj['z'], VSj, dSj)
    dSj_dt = VSj * P[3] - GammaSj * P[3] - dSj * P[3]

    # Sa (small adult clupeids)
    VSa = Sa['Sigma'] * Gsa * Es - Sa['T']
    dSa = Sa['mu'] + Sa['F'] + Psa - min(VSa, 0)
    GammaSa = Maturation(Sa['z'], Sa['kappa'] * VSa, dSa)
    dSa_dt = GammaSj * P[3] + Sa['kappa'] * max(VSa, 0) * P[4] - GammaSa * P[4] - dSa * P[4]
    
    # Sb (big adult clupeids)
    VSb = Sb['Sigma'] * Gsb * Es - Sb['T']
    dSb = Sb['mu'] + Sb['F'] + Psb - min(VSb, 0)
    dSb_dt = GammaSa * P[4] - dSb * P[5]

    # Sga and Sgb 
    dSga_dt = (1 - Sa['kappa']) * max(VSa, 0) * P[4] - GammaSa * P[6] - dSa * P[6]
    dSgb_dt = GammaSa * P[6] + max(VSb, 0) * P[5] - dSb * P[7]

    # Cod functions and ODEs 
    # Cj (juvenile cod)
    VCj = Cj['Sigma'] * Gcj * Ecj - Cj['T']
    dCj = Cj['mu'] + Cj['F']
    GammaCj = Maturation(Cj['z'], VCj, dCj) 
    dCj_dt = VCj * P[8] - GammaCj * P[8] - dCj * P[8]

    # Ca (small adult cod)
    VCa = Ca['Sigma'] * Gca * Eca - Ca['T']
    dCa = Ca['mu'] + Ca['F'] - min(VCa, 0)
    GammaCa = Maturation(Ca['z'], Ca['kappa'] * VCa, dCa)
    dCa_dt = GammaCj * P[8] + Ca['kappa'] * max(VCa, 0) * P[9] - GammaCa * P[9] - dCa * P[9]

    # Cb (big adult cod)
    VCb = Cb['Sigma'] * Gcb * Ecb - Cb['T']
    dCb = Cb['mu'] + Cb['F'] - min(VCb, 0)
    dCb_dt = GammaCa * P[9] - dCb * P[10]

    # Cga and Cgb 
    dCga_dt = (1 - Ca['kappa']) * max(VCa, 0) * P[9] - GammaCa * P[11] - dCa * P[11]
    dCgb_dt = GammaCa * P[11] + max(VCb, 0) * P[10] - dCb * P[12]

    return dRs_dt, dRj_dt, dRa_dt, dSj_dt, dSa_dt, dSb_dt, dSga_dt, dSgb_dt, dCj_dt, dCa_dt, dCb_dt, dCga_dt, dCgb_dt

def function_reproduction_season(P):
    """
    Function for the reproduction season, to calculate how much reproduction will take place 
    and will thus be converted into juvenile biomass. At the same time, the biomass of the reproductive
    storages will be set to zero for the next growing season. 

    Input:
        P: array of floats - final biomass values of last growing season

    Output:
        Rs_tplus : float - Resource for clupeids biomass value after reproduction 
        Rj_tplus : float - Resource for juvenile cod biomass value after reproduction 
        Ra_tplus : float - Resource for adult cod biomass value after reproduction 
        Sj_tplus : float - Juvenile clupeids biomass value after reproduction 
        Sa_tplus : float - Small adult clupeids biomass value after reproduction 
        Sb_tplus : float - Large adult clupeids biomass value after reproduction 
        Sga_tplus : float - Small adult clupeids reproductive storages biomass value after reproduction 
        Sgb_tplus : float - Large adult clupeids reproductive storages biomass value after reproduction 
        Cj_tplus : float - Juvenile cod biomass value after reproduction 
        Ca_tplus : float - Small adult cod biomass value after reproduction 
        Cb_tplus : float - Large adult cod biomass value after reproduction 
        Cga_tplus : float - Small adult cod reproductive storages biomass value after reproduction 
        Cgb_tplus : float - Large adult cod reproductive storages biomass value after reproduction 

    """
    # resource reproduction funtions
    Rs_tplus = P[0]
    Rj_tplus = P[1]
    Ra_tplus = P[2]

    # clupeids reproduction functions
    Sj_tplus = P[3] + P[6] + P[7]
    Sa_tplus = P[4]
    Sb_tplus = P[5]
    Sga_tplus = 0
    Sgb_tplus = 0

    # cod reproduction functions 
    Cj_tplus = P[8] + P[11] + P[12]
    Ca_tplus = P[9]
    Cb_tplus = P[10]
    Cga_tplus = 0
    Cgb_tplus = 0 

    return Rs_tplus, Rj_tplus, Ra_tplus, Sj_tplus, Sa_tplus, Sb_tplus, Sga_tplus, Sgb_tplus, Cj_tplus, Ca_tplus, Cb_tplus, Cga_tplus, Cgb_tplus

######################################################################################################################################################
##                                                          simulation                                                                              ##
######################################################################################################################################################

temperature_range = np.linspace(6,20,160)
temperature_range_reverse = np.linspace(6,20,160)
temperatureKelvin = [279.5, 280.7, 283, 285]

# cod mortality (Fc)
codMortality = np.linspace(0, 1.4, 140)

# t_trans represents the first 50 years of an integration period
t_trans = 50

# t_measure represents the last year of an integration period 
t_measure = 99

### lists to store the data in, which eventually will be stored in CSV files 

# low clupeid mortality, with temperature as bifurcation parameter
# from low to high T
Rs_list_low = []
Rj_list_low  = []
Ra_list_low  = []
Sj_list_low  = []
Sa_list_low  = []
Sb_list_low  = []
Sga_list_low  = []
Sgb_list_low  = []
Cj_list_low  = []
Ca_list_low  = []
Cb_list_low  = []
Cga_list_low  = []
Cgb_list_low  = []
clupeid_adult_list_low  = []
cod_adult_list_low  = []

# high clupeid mortality, with temperature as bifurcation parameter
# from low to high T
Rs_list_high = []
Rj_list_high = []
Ra_list_high = []
Sj_list_high = []
Sa_list_high = []
Sb_list_high = []
Sga_list_high = []
Sgb_list_high = []
Cj_list_high = []
Ca_list_high = []
Cb_list_high = []
Cga_list_high = []
Cgb_list_high = []
clupeid_adult_list_high = []
cod_adult_list_high = []

# low clupeid mortality, with temperature as bifurcation parameter
# from high to low T
Rs_list_low_reverse = []
Rj_list_low_reverse  = []
Ra_list_low_reverse  = []
Sj_list_low_reverse  = []
Sa_list_low_reverse = []
Sb_list_low_reverse  = []
Sga_list_low_reverse  = []
Sgb_list_low_reverse  = []
Cj_list_low_reverse  = []
Ca_list_low_reverse  = []
Cb_list_low_reverse  = []
Cga_list_low_reverse  = []
Cgb_list_low_reverse  = []
clupeid_adult_list_low_reverse  = []
cod_adult_list_low_reverse  = []

# high clupeid mortality, with temperature as bifurcation parameter
# from high to low T 
Rs_list_high_reverse = []
Rj_list_high_reverse = []
Ra_list_high_reverse = []
Sj_list_high_reverse = []
Sa_list_high_reverse = []
Sb_list_high_reverse = []
Sga_list_high_reverse = []
Sgb_list_high_reverse = []
Cj_list_high_reverse = []
Ca_list_high_reverse = []
Cb_list_high_reverse = []
Cga_list_high_reverse = []
Cgb_list_high_reverse = []
clupeid_adult_list_high_reverse = []
cod_adult_list_high_reverse = []


# low clupeid mortality and 6.5 degrees Celcius (279.5 K)
Rs_list_mean279_5_low = []
Rj_list_mean279_5_low = []
Ra_list_mean279_5_low = []
Sj_list_mean279_5_low = []
Sa_list_mean279_5_low = []
Sb_list_mean279_5_low = []
Sga_list_mean279_5_low = []
Sgb_list_mean279_5_low = []
Cj_list_mean279_5_low = []
Ca_list_mean279_5_low = []
Cb_list_mean279_5_low = []
Cga_list_mean279_5_low = []
Cgb_list_mean279_5_low = []
clupeid_adult_list_mean279_5_low = []
cod_adult_list_mean279_5_low = []

# low clupeid mortality and 7.7 degrees Celcius (280.7 K)
Rs_list_mean280_7_low = []
Rj_list_mean280_7_low = []
Ra_list_mean280_7_low = []
Sj_list_mean280_7_low = []
Sa_list_mean280_7_low = []
Sb_list_mean280_7_low = []
Sga_list_mean280_7_low = []
Sgb_list_mean280_7_low = []
Cj_list_mean280_7_low = []
Ca_list_mean280_7_low = []
Cb_list_mean280_7_low = []
Cga_list_mean280_7_low = []
Cgb_list_mean280_7_low = []
clupeid_adult_list_mean280_7_low = []
cod_adult_list_mean280_7_low = []

# low clupeid mortality and 10 degrees Celcius (283 K)
Rs_list_mean283_low = []
Rj_list_mean283_low = []
Ra_list_mean283_low = []
Sj_list_mean283_low = []
Sa_list_mean283_low = []
Sb_list_mean283_low = []
Sga_list_mean283_low = []
Sgb_list_mean283_low = []
Cj_list_mean283_low = []
Ca_list_mean283_low = []
Cb_list_mean283_low = []
Cga_list_mean283_low = []
Cgb_list_mean283_low = []
clupeid_adult_list_mean283_low = []
cod_adult_list_mean283_low = []

# low clupeid mortality and 12 degrees Celcius (285 K)
Rs_list_mean285_low = []
Rj_list_mean285_low = []
Ra_list_mean285_low = []
Sj_list_mean285_low = []
Sa_list_mean285_low = []
Sb_list_mean285_low = []
Sga_list_mean285_low = []
Sgb_list_mean285_low = []
Cj_list_mean285_low = []
Ca_list_mean285_low = []
Cb_list_mean285_low = []
Cga_list_mean285_low = []
Cgb_list_mean285_low = []
clupeid_adult_list_mean285_low = []
cod_adult_list_mean285_low = []

# high clupeid mortality and 6.5 degrees Celcius (279.5 K)
Rs_list_mean279_5_high = []
Rj_list_mean279_5_high = []
Ra_list_mean279_5_high = []
Sj_list_mean279_5_high = []
Sa_list_mean279_5_high = []
Sb_list_mean279_5_high = []
Sga_list_mean279_5_high = []
Sgb_list_mean279_5_high = []
Cj_list_mean279_5_high = []
Ca_list_mean279_5_high = []
Cb_list_mean279_5_high = []
Cga_list_mean279_5_high = []
Cgb_list_mean279_5_high = []
clupeid_adult_list_mean279_5_high = []
cod_adult_list_mean279_5_high = []

# high clupeid mortality and 7.7 degrees Celcius (280.7 K)
Rs_list_mean280_7_high = []
Rj_list_mean280_7_high = []
Ra_list_mean280_7_high = []
Sj_list_mean280_7_high = []
Sa_list_mean280_7_high= []
Sb_list_mean280_7_high = []
Sga_list_mean280_7_high = []
Sgb_list_mean280_7_high = []
Cj_list_mean280_7_high = []
Ca_list_mean280_7_high = []
Cb_list_mean280_7_high = []
Cga_list_mean280_7_high = []
Cgb_list_mean280_7_high = []
clupeid_adult_list_mean280_7_high = []
cod_adult_list_mean280_7_high = []

# high clupeid mortality and 10 degrees Celcius (283 K)
Rs_list_mean283_high = []
Rj_list_mean283_high = []
Ra_list_mean283_high = []
Sj_list_mean283_high = []
Sa_list_mean283_high = []
Sb_list_mean283_high = []
Sga_list_mean283_high = []
Sgb_list_mean283_high = []
Cj_list_mean283_high = []
Ca_list_mean283_high = []
Cb_list_mean283_high = []
Cga_list_mean283_high = []
Cgb_list_mean283_high = []
clupeid_adult_list_mean283_high = []
cod_adult_list_mean283_high = []

# high clupeid mortality and 12 degrees Celcius (285 K)
Rs_list_mean285_high = []
Rj_list_mean285_high = []
Ra_list_mean285_high = []
Sj_list_mean285_high = []
Sa_list_mean285_high = []
Sb_list_mean285_high = []
Sga_list_mean285_high = []
Sgb_list_mean285_high = []
Cj_list_mean285_high = []
Ca_list_mean285_high = []
Cb_list_mean285_high = []
Cga_list_mean285_high = []
Cgb_list_mean285_high = []
clupeid_adult_list_mean285_high = []
cod_adult_list_mean285_high = []

# lists for mean annual total yield 
meanTotalYieldCod_temp1_low = []
meanTotalYieldCod_temp2_low = []
meanTotalYieldCod_temp3_low = []
meanTotalYieldCod_temp4_low = []

meanTotalYieldCod_temp1_high = []
meanTotalYieldCod_temp2_high = []
meanTotalYieldCod_temp3_high = []
meanTotalYieldCod_temp4_high = []

meanTotalYieldClup_temp1_low = []
meanTotalYieldClup_temp2_low = []
meanTotalYieldClup_temp3_low = []
meanTotalYieldClup_temp4_low = []

meanTotalYieldClup_temp1_high = []
meanTotalYieldClup_temp2_high = []
meanTotalYieldClup_temp3_high = []
meanTotalYieldClup_temp4_high = []

# calculating biomass with T as bifurcation parameter, high Fc, from low to high T
if biomass_or_yield == 0:
    # Fc is high 
    Fc = 0.002 

    # iterating over temperature, with low clupeid fishing mortality
    for temp in temperature_range:
        print(temp)
    
        # converting temperature from Celsius to Kelvin
        T = temp + 273 
        x = 0

        Delta = 0.1 * math.exp((Edelta * (T-T0))/(k*T*T0))
        Rsmax = 98 * math.exp((ERmax * (T-T0))/(k*T*T0))
        Rjmax = 1 * math.exp((ERmax * (T-T0))/(k*T*T0))
        Ramax = 0.75 * math.exp((ERmax * (T-T0))/(k*T*T0))

        # if it is the first year, the initial values will be reset 
        if x == 0:
            Rs1 = 32.7816744728
            Rj1 = 0.4713946441
            Ra1 =  0.7175924292
            Sj1 = 10.4824915096
            Sa1 =  69.8208217446
            Sb1 =  89.7371125270
            Sga1 = 0
            Sgb1 = 0
            Cj1 =  6.2772314150
            Ca1 = 9.4197483660
            Cb1 = 0.6284148535
            Cga1 = 0
            Cgb1 = 0

        # if it is not the first year, the initial values will be the biomass values of the last reproductive season
        if x > 0:
            Rs1, Rj1, Ra1, Sj1, Sa1, Sb1, Sga1, Sgb1, Cj1, Ca1, Cb1, Cga1, Cgb1 = P01

        # iterating over 100 years for each cod mortality value 
        for n in range(1, 100):

            # Fs = fishing mortality on clupeid (low = 0.2 per year = 0.0008 per day, high = 0.5 per year = 0.002 per day)
            Fs= 0.0008

            tn_minus= (n - 1) * Y
            tn_plus = n * Y 
            ts = np.linspace(tn_minus, tn_plus, Y)

            # ODE function , for 250 days, starting with the initial values 
            growing_season_1 = solve_ivp(function_growing_season, [tn_minus, tn_plus],  [Rs1, Rj1, Ra1, Sj1, Sa1, Sb1, Sga1, Sgb1, Cj1, Ca1, Cb1, Cga1, Cgb1], first_step = 0.01 ,max_step = 1.0, t_eval = ts, rtol = 1.0E-8, atol = 1.0E-10)
                
            # getting the useful data out of growing_season_1 and transposing it, so the rows and columns are placed in the right way
            data_1 = growing_season_1.y 
            data_1 = data_1.transpose()

            # the data of the first 50 years won't be stored 
            if n > t_trans: 
                if n == t_trans + 1:
                    all_growing_season_1 = data_1
                else:
                    all_growing_season_1 = np.vstack((all_growing_season_1, data_1))

            if n == t_measure: 
                break

            # using the reproduction function with the values of the last day of the growing season 
            P01 = function_reproduction_season(data_1[-1, :])

            # setting the new initial values for next year to be the reproduction values of this year 
            Rs1, Rj1, Ra1, Sj1, Sa1, Sb1, Sga1, Sgb1, Cj1, Ca1, Cb1, Cga1, Cgb1 = P01
                    
        # storing the data in lists 
        Rs_list_low.append(mean(all_growing_season_1[:, 0])*0.003)
        Rj_list_low.append(mean(all_growing_season_1[:, 1])*0.003)
        Ra_list_low.append(mean(all_growing_season_1[:, 2])*0.003)
        Sj_list_low.append(mean(all_growing_season_1[:, 3])*0.003)
        Sa_list_low.append(mean(all_growing_season_1[:, 4])*0.003)
        Sb_list_low.append(mean(all_growing_season_1[:, 5])*0.003)
        Sga_list_low.append(mean(all_growing_season_1[:, 6])*0.003)
        Sgb_list_low.append(mean(all_growing_season_1[:, 7])*0.003)
        Cj_list_low.append(mean(all_growing_season_1[:, 8])*0.003)
        Ca_list_low.append(mean(all_growing_season_1[:, 9])*0.003)
        Cb_list_low.append(mean(all_growing_season_1[:, 10])*0.003)
        Cga_list_low.append(mean(all_growing_season_1[:, 11])*0.003)
        Cgb_list_low.append(mean(all_growing_season_1[:, 12])*0.003)

        adult_clupeid_low = (mean(all_growing_season_1[:,4]) + mean(all_growing_season_1[:,5]) + mean(all_growing_season_1[:,6]) + mean(all_growing_season_1[:, 7])) * 0.003
        clupeid_adult_list_low.append(adult_clupeid_low)

        adult_cod_low = (mean(all_growing_season_1[:,9]) + mean(all_growing_season_1[:,10]) + mean(all_growing_season_1[:,11]) + mean(all_growing_season_1[:, 12])) * 0.003
        cod_adult_list_low.append(adult_cod_low)

        x = x + 1

    # iterating over temperature range, with high clupeid fishing mortality (Fs)
    for temp in temperature_range:

        # changing temperature from Celsius to Kelvin
        T = temp + 273
        x = 0
        
        Delta = 0.1 * math.exp((Edelta * (T-T0))/(k*T*T0))
        Rsmax = 98 * math.exp((ERmax * (T-T0))/(k*T*T0))
        Rjmax = 1 * math.exp((ERmax * (T-T0))/(k*T*T0))
        Ramax = 0.75 * math.exp((ERmax * (T-T0))/(k*T*T0))

        # if it is the first year, the initial values will be reset 
        if x == 0 :
            Rs2 = 32.7816744728
            Rj2 = 0.4713946441
            Ra2 =  0.7175924292
            Sj2 = 10.4824915096
            Sa2 =  69.8208217446
            Sb2 =  89.7371125270
            Sga2 = 0
            Sgb2 = 0
            Cj2 =  6.2772314150
            Ca2 = 9.4197483660
            Cb2 = 0.6284148535
            Cga2 = 0
            Cgb2 = 0

        # if it is not the first year, the initial values will be the biomass values of the last reproductive season
        if x > 0:
            Rs2, Rj2, Ra2, Sj2, Sa2, Sb2, Sga2, Sgb2, Cj2, Ca2, Cb2, Cga2, Cgb2 = P02

        # iterating over 100 years for each Fc value
        for n in range(1, 100):

            # Fs = fishing mortality on clupeid (low = 0.2 per year = 0.0008 per day, high = 0.5 per year = 0.002 per day)
            Fs = 0.002

            tn_minus= (n-1) * Y
            tn_plus = n * Y 
            ts = np.linspace(tn_minus,tn_plus,Y)

            growing_season_2 = solve_ivp(function_growing_season, [tn_minus, tn_plus],  [Rs2, Rj2, Ra2, Sj2, Sa2, Sb2, Sga2, Sgb2, Cj2, Ca2, Cb2, Cga2, Cgb2], first_step = 0.01, max_step = 1.0, t_eval = ts, rtol = 1.0E-8, atol = 1.0E-10)

            # getting the useful data out of growing_season_1 and transposing it, so the rows and columns are placed in the right way
            data_2 = growing_season_2.y 
            data_2 = data_2.transpose()

            # the data of the first 50 years won't be stored
            if n > t_trans:
                if n == t_trans + 1:
                    all_growing_season_2 = data_2
                else:
                    all_growing_season_2 = np.vstack((all_growing_season_2, data_2))
                    
            if n == t_measure: 
                break

            P02 = function_reproduction_season(data_2[-1, :])

            # setting the new initial values for next year to be the reproduction values of this year 
            Rs2, Rj2, Ra2, Sj2, Sa2, Sb2, Sga2, Sgb2, Cj2, Ca2, Cb2, Cga2, Cgb2 = P02
                    
        # storing the data in lists   
        Rs_list_high.append(mean(all_growing_season_2[:, 0])*0.003)
        Rj_list_high.append(mean(all_growing_season_2[:, 1])*0.003)
        Ra_list_high.append(mean(all_growing_season_2[:, 2])*0.003)
        Sj_list_high.append(mean(all_growing_season_2[:, 3])*0.003)
        Sa_list_high.append(mean(all_growing_season_2[:, 4])*0.003)
        Sb_list_high.append(mean(all_growing_season_2[:, 5])*0.003)
        Sga_list_high.append(mean(all_growing_season_2[:, 6])*0.003)
        Sgb_list_high.append(mean(all_growing_season_2[:, 7])*0.003)
        Cj_list_high.append(mean(all_growing_season_2[:, 8])*0.003)
        Ca_list_high.append(mean(all_growing_season_2[:, 9])*0.003)
        Cb_list_high.append(mean(all_growing_season_2[:, 10])*0.003)
        Cga_list_high.append(mean(all_growing_season_2[:, 11])*0.003)
        Cgb_list_high.append(mean(all_growing_season_2[:, 12])*0.003)

        adult_clupeid_high = (mean(all_growing_season_2[:,4]) + mean(all_growing_season_2[:,5]) + mean(all_growing_season_2[:,6]) + mean(all_growing_season_2[:, 7])) * 0.003
        clupeid_adult_list_high.append(adult_clupeid_high)

        adult_cod_high = (mean(all_growing_season_2[:,9]) + mean(all_growing_season_2[:,10]) + mean(all_growing_season_2[:,11]) + mean(all_growing_season_2[:, 12])) * 0.003
        cod_adult_list_high.append(adult_cod_high)

        x = x + 1

# calculating biomass with T as bifurcation parameter, high Fc, from high to low T
if biomass_or_yield == 1:
    # Fc is high 
    Fc = 0.003 

    # iterating over temperature, with low clupeid fishing mortality
    for temp in temperature_range_reverse:
    
        # converting temperature from Celsius to Kelvin
        T = temp + 273 
        x = 0

        Delta = 0.1 * math.exp((Edelta * (T-T0))/(k*T*T0))
        Rsmax = 98 * math.exp((ERmax * (T-T0))/(k*T*T0))
        Rjmax = 1 * math.exp((ERmax * (T-T0))/(k*T*T0))
        Ramax = 0.75 * math.exp((ERmax * (T-T0))/(k*T*T0))

        # if it is the first year, the initial values will be reset 
        if x == 0:
            Rs1 = 32.7816744728
            Rj1 = 0.4713946441
            Ra1 =  0.7175924292
            Sj1 = 10.4824915096
            Sa1 =  69.8208217446
            Sb1 =  89.7371125270
            Sga1 = 0
            Sgb1 = 0
            Cj1 =  6.2772314150
            Ca1 = 9.4197483660
            Cb1 = 0.6284148535
            Cga1 = 0
            Cgb1 = 0

        # if it is not the first year, the initial values will be the biomass values of the last reproductive season
        if x > 0:
            Rs1, Rj1, Ra1, Sj1, Sa1, Sb1, Sga1, Sgb1, Cj1, Ca1, Cb1, Cga1, Cgb1 = P01
            
            # if the adult cod (total) biomass is below 10**-5, the values will turn into 10**-5
            if (Ca1 + Cb1 + Cga1 + Cgb1) < 10**(-5):
                Ca1 = 10**(-5)
                Cb1 = 10**(-5)
                Cga1 = 10**(-5)
                Cgb1 = 10**(-5)

        # iterating over 100 years for each cod mortality value 
        for n in range(1, 100):

            # Fs = fishing mortality on clupeid (low = 0.2 per year = 0.0008 per day, high = 0.5 per year = 0.002 per day)
            Fs= 0.0008

            tn_minus= (n - 1) * Y
            tn_plus = n * Y 
            ts = np.linspace(tn_minus, tn_plus, Y)

            # ODE function , for 250 days, starting with the initial values 
            growing_season_1 = solve_ivp(function_growing_season, [tn_minus, tn_plus],  [Rs1, Rj1, Ra1, Sj1, Sa1, Sb1, Sga1, Sgb1, Cj1, Ca1, Cb1, Cga1, Cgb1], first_step = 0.01 ,max_step = 1.0, t_eval = ts, rtol = 1.0E-8, atol = 1.0E-10)
                
            # getting the useful data out of growing_season_1 and transposing it, so the rows and columns are placed in the right way
            data_1 = growing_season_1.y 
            data_1 = data_1.transpose()

            # the data of the first 50 years won't be stored 
            if n > t_trans: 
                if n == t_trans + 1:
                    all_growing_season_1 = data_1
                else:
                    all_growing_season_1 = np.vstack((all_growing_season_1, data_1))

            if n == t_measure: 
                break

            # using the reproduction function with the values of the last day of the growing season 
            P01 = function_reproduction_season(data_1[-1, :])

            # setting the new initial values for next year to be the reproduction values of this year 
            Rs1, Rj1, Ra1, Sj1, Sa1, Sb1, Sga1, Sgb1, Cj1, Ca1, Cb1, Cga1, Cgb1 = P01
                    
        # storing the data in lists 
        Rs_list_low_reverse.append(mean(all_growing_season_1[:, 0])*0.003)
        Rj_list_low_reverse.append(mean(all_growing_season_1[:, 1])*0.003)
        Ra_list_low_reverse.append(mean(all_growing_season_1[:, 2])*0.003)
        Sj_list_low_reverse.append(mean(all_growing_season_1[:, 3])*0.003)
        Sa_list_low_reverse.append(mean(all_growing_season_1[:, 4])*0.003)
        Sb_list_low_reverse.append(mean(all_growing_season_1[:, 5])*0.003)
        Sga_list_low_reverse.append(mean(all_growing_season_1[:, 6])*0.003)
        Sgb_list_low_reverse.append(mean(all_growing_season_1[:, 7])*0.003)
        Cj_list_low_reverse.append(mean(all_growing_season_1[:, 8])*0.003)
        Ca_list_low_reverse.append(mean(all_growing_season_1[:, 9])*0.003)
        Cb_list_low_reverse.append(mean(all_growing_season_1[:, 10])*0.003)
        Cga_list_low_reverse.append(mean(all_growing_season_1[:, 11])*0.003)
        Cgb_list_low_reverse.append(mean(all_growing_season_1[:, 12])*0.003)

        adult_clupeid_low_reverse = (mean(all_growing_season_1[:,4]) + mean(all_growing_season_1[:,5]) + mean(all_growing_season_1[:,6]) + mean(all_growing_season_1[:, 7])) * 0.003
        clupeid_adult_list_low_reverse.append(adult_clupeid_low_reverse)

        adult_cod_low_reverse = (mean(all_growing_season_1[:,9]) + mean(all_growing_season_1[:,10]) + mean(all_growing_season_1[:,11]) + mean(all_growing_season_1[:, 12])) * 0.003
        cod_adult_list_low_reverse.append(adult_cod_low_reverse)

        x = x + 1

    # iterating over temperature range, with high clupeid fishing mortality (Fs)
    for temp in temperature_range:
    
        # changing temperature from Celsius to Kelvin
        T = temp + 273
        x = 0
        
        Delta = 0.1 * math.exp((Edelta * (T-T0))/(k*T*T0))
        Rsmax = 98 * math.exp((ERmax * (T-T0))/(k*T*T0))
        Rjmax = 1 * math.exp((ERmax * (T-T0))/(k*T*T0))
        Ramax = 0.75 * math.exp((ERmax * (T-T0))/(k*T*T0))

        # if it is the first year, the initial values will be reset 
        if x == 0 :
            Rs2 = 32.7816744728
            Rj2 = 0.4713946441
            Ra2 =  0.7175924292
            Sj2 = 10.4824915096
            Sa2 =  69.8208217446
            Sb2 =  89.7371125270
            Sga2 = 0
            Sgb2 = 0
            Cj2 =  6.2772314150
            Ca2 = 9.4197483660
            Cb2 = 0.6284148535
            Cga2 = 0
            Cgb2 = 0

        # if it is not the first year, the initial values will be the biomass values of the last reproductive season
        if x > 0:
            Rs2, Rj2, Ra2, Sj2, Sa2, Sb2, Sga2, Sgb2, Cj2, Ca2, Cb2, Cga2, Cgb2 = P02
            
            # if the adult cod (total) biomass is below 10**-5, the values will turn into 10**-5
            if (Ca2 + Cb2 + Cga2 + Cgb2) < 10**(-5):
                Ca2 = 10**(-5)
                Cb2 = 10**(-5)
                Cga2 = 10**(-5)
                Cgb2 = 10**(-5)

        # iterating over 100 years for each Fc value
        for n in range(1, 100):

            # Fs = fishing mortality on clupeid (low = 0.2 per year = 0.0008 per day, high = 0.5 per year = 0.002 per day)
            Fs = 0.002

            tn_minus= (n-1) * Y
            tn_plus = n * Y 
            ts = np.linspace(tn_minus,tn_plus,Y)

            growing_season_2 = solve_ivp(function_growing_season, [tn_minus, tn_plus],  [Rs2, Rj2, Ra2, Sj2, Sa2, Sb2, Sga2, Sgb2, Cj2, Ca2, Cb2, Cga2, Cgb2], first_step = 0.01, max_step = 1.0, t_eval = ts, rtol = 1.0E-8, atol = 1.0E-10)

            # getting the useful data out of growing_season_1 and transposing it, so the rows and columns are placed in the right way
            data_2 = growing_season_2.y 
            data_2 = data_2.transpose()

            # the data of the first 50 years won't be stored
            if n > t_trans:
                if n == t_trans + 1:
                    all_growing_season_2 = data_2
                else:
                    all_growing_season_2 = np.vstack((all_growing_season_2, data_2))
                    
            if n == t_measure: 
                break

            P02 = function_reproduction_season(data_2[-1, :])

            # setting the new initial values for next year to be the reproduction values of this year 
            Rs2, Rj2, Ra2, Sj2, Sa2, Sb2, Sga2, Sgb2, Cj2, Ca2, Cb2, Cga2, Cgb2 = P02
                    
        # storing the data in lists   
        Rs_list_high_reverse.append(mean(all_growing_season_2[:, 0])*0.003)
        Rj_list_high_reverse.append(mean(all_growing_season_2[:, 1])*0.003)
        Ra_list_high_reverse.append(mean(all_growing_season_2[:, 2])*0.003)
        Sj_list_high_reverse.append(mean(all_growing_season_2[:, 3])*0.003)
        Sa_list_high_reverse.append(mean(all_growing_season_2[:, 4])*0.003)
        Sb_list_high_reverse.append(mean(all_growing_season_2[:, 5])*0.003)
        Sga_list_high_reverse.append(mean(all_growing_season_2[:, 6])*0.003)
        Sgb_list_high_reverse.append(mean(all_growing_season_2[:, 7])*0.003)
        Cj_list_high_reverse.append(mean(all_growing_season_2[:, 8])*0.003)
        Ca_list_high_reverse.append(mean(all_growing_season_2[:, 9])*0.003)
        Cb_list_high_reverse.append(mean(all_growing_season_2[:, 10])*0.003)
        Cga_list_high_reverse.append(mean(all_growing_season_2[:, 11])*0.003)
        Cgb_list_high_reverse.append(mean(all_growing_season_2[:, 12])*0.003)

        adult_clupeid_high_reverse = (mean(all_growing_season_2[:,4]) + mean(all_growing_season_2[:,5]) + mean(all_growing_season_2[:,6]) + mean(all_growing_season_2[:, 7])) * 0.003
        clupeid_adult_list_high_reverse.append(adult_clupeid_high_reverse)

        adult_cod_high_reverse = (mean(all_growing_season_2[:,9]) + mean(all_growing_season_2[:,10]) + mean(all_growing_season_2[:,11]) + mean(all_growing_season_2[:, 12])) * 0.003
        cod_adult_list_high_reverse.append(adult_cod_high_reverse)

        x = x + 1

# iterating over 4 different temperatures, for both low and high Fs 
# Fc as bifurcation parameter
if biomass_or_yield == 2:

    # iterating over 4 different T, for low Fs
    for temperature in temperatureKelvin:
        T = temperature
        x = 0

        for codMortalityPerYear in codMortality:
            # dividing the Fc value by 250, to get the Fc value per day
            Fc = codMortalityPerYear / 250 

            Delta = 0.1 * math.exp((Edelta * (T-T0))/(k*T*T0))
            Rsmax = 98 * math.exp((ERmax * (T-T0))/(k*T*T0))
            Rjmax = 1 * math.exp((ERmax * (T-T0))/(k*T*T0))
            Ramax = 0.75 * math.exp((ERmax * (T-T0))/(k*T*T0))

            # if it is the first year, the initial values will be reset 
            if x == 0:
                Rs1 = 32.7816744728
                Rj1 = 0.4713946441
                Ra1 =  0.7175924292
                Sj1 = 10.4824915096
                Sa1 =  69.8208217446
                Sb1 =  89.7371125270
                Sga1 = 0
                Sgb1 = 0
                Cj1 =  6.2772314150
                Ca1 = 9.4197483660
                Cb1 = 0.6284148535
                Cga1 = 0
                Cgb1 = 0

            # if it is not the first year, the initial values will be the biomass values of the last reproductive season
            if x > 0:
                Rs1, Rj1, Ra1, Sj1, Sa1, Sb1, Sga1, Sgb1, Cj1, Ca1, Cb1, Cga1, Cgb1 = P01

            # iterating over 100 years for each cod mortality value 
            for n in range(1, 100):

                # Fs = fishing mortality on clupeid (low = 0.2 per year = 0.0008 per day, high = 0.5 per year = 0.002 per day)
                Fs= 0.0008

                tn_minus= (n - 1) * Y
                tn_plus = n * Y 
                ts = np.linspace(tn_minus, tn_plus, Y)

                # ODE function , for 250 days, starting with the initial values 
                growing_season_1 = solve_ivp(function_growing_season, [tn_minus, tn_plus],  [Rs1, Rj1, Ra1, Sj1, Sa1, Sb1, Sga1, Sgb1, Cj1, Ca1, Cb1, Cga1, Cgb1], first_step = 0.01 ,max_step = 1.0, t_eval = ts, rtol = 1.0E-8, atol = 1.0E-10)
                
                # getting the useful data out of growing_season_1 and transposing it, so the rows and columns are placed in the right way
                data_1 = growing_season_1.y 
                data_1 = data_1.transpose()

                # the data of the first 50 years won't be stored 
                if n > t_trans: 
                    if n == t_trans + 1:
                        all_growing_season_1 = data_1
                    else:
                        all_growing_season_1 = np.vstack((all_growing_season_1, data_1))

                if n == t_measure: 
                    break

                # using the reproduction function with the values of the last day of the growing season 
                P01 = function_reproduction_season(data_1[-1, :])

                # setting the new initial values for next year to be the reproduction values of this year 
                Rs1, Rj1, Ra1, Sj1, Sa1, Sb1, Sga1, Sgb1, Cj1, Ca1, Cb1, Cga1, Cgb1 = P01
                    
            # storing the data in lists 
            if T == 280.7:
                Rs_list_mean280_7_low.append(mean(all_growing_season_1[:, 0])*0.003)
                Rj_list_mean280_7_low.append(mean(all_growing_season_1[:, 1])*0.003)
                Ra_list_mean280_7_low.append(mean(all_growing_season_1[:, 2])*0.003)
                Sj_list_mean280_7_low.append(mean(all_growing_season_1[:, 3])*0.003)
                Sa_list_mean280_7_low.append(mean(all_growing_season_1[:, 4])*0.003)
                Sb_list_mean280_7_low.append(mean(all_growing_season_1[:, 5])*0.003)
                Sga_list_mean280_7_low.append(mean(all_growing_season_1[:, 6])*0.003)
                Sgb_list_mean280_7_low.append(mean(all_growing_season_1[:, 7])*0.003)
                Cj_list_mean280_7_low.append(mean(all_growing_season_1[:, 8])*0.003)
                Ca_list_mean280_7_low.append(mean(all_growing_season_1[:, 9])*0.003)
                Cb_list_mean280_7_low.append(mean(all_growing_season_1[:, 10])*0.003)
                Cga_list_mean280_7_low.append(mean(all_growing_season_1[:, 11])*0.003)
                Cgb_list_mean280_7_low.append(mean(all_growing_season_1[:, 12])*0.003)

                adult_clupeid280_7_low = (mean(all_growing_season_1[:,4]) + mean(all_growing_season_1[:,5]) + mean(all_growing_season_1[:,6]) + mean(all_growing_season_1[:, 7])) * 0.003
                clupeid_adult_list_mean280_7_low.append(adult_clupeid280_7_low)

                adult_cod280_7_low = (mean(all_growing_season_1[:,9]) + mean(all_growing_season_1[:,10]) + mean(all_growing_season_1[:,11]) + mean(all_growing_season_1[:, 12])) * 0.003
                cod_adult_list_mean280_7_low.append(adult_cod280_7_low)

            if T == 283:
                Rs_list_mean283_low.append(mean(all_growing_season_1[:, 0])*0.003)
                Rj_list_mean283_low.append(mean(all_growing_season_1[:, 1])*0.003)
                Ra_list_mean283_low.append(mean(all_growing_season_1[:, 2])*0.003)
                Sj_list_mean283_low.append(mean(all_growing_season_1[:, 3])*0.003)
                Sa_list_mean283_low.append(mean(all_growing_season_1[:, 4])*0.003)
                Sb_list_mean283_low.append(mean(all_growing_season_1[:, 5])*0.003)
                Sga_list_mean283_low.append(mean(all_growing_season_1[:, 6])*0.003)
                Sgb_list_mean283_low.append(mean(all_growing_season_1[:, 7])*0.003)
                Cj_list_mean283_low.append(mean(all_growing_season_1[:, 8])*0.003)
                Ca_list_mean283_low.append(mean(all_growing_season_1[:, 9])*0.003)
                Cb_list_mean283_low.append(mean(all_growing_season_1[:, 10])*0.003)
                Cga_list_mean283_low.append(mean(all_growing_season_1[:, 11])*0.003)
                Cgb_list_mean283_low.append(mean(all_growing_season_1[:, 12])*0.003)

                adult_clupeid283_low = (mean(all_growing_season_1[:,4]) + mean(all_growing_season_1[:,5]) + mean(all_growing_season_1[:,6]) + mean(all_growing_season_1[:, 7])) * 0.003
                clupeid_adult_list_mean283_low.append(adult_clupeid283_low)

                adult_cod283_low = (mean(all_growing_season_1[:,9]) + mean(all_growing_season_1[:,10]) + mean(all_growing_season_1[:,11]) + mean(all_growing_season_1[:, 12])) * 0.003
                cod_adult_list_mean283_low.append(adult_cod283_low)

            if T == 285:
                Rs_list_mean285_low.append(mean(all_growing_season_1[:, 0])*0.003)
                Rj_list_mean285_low.append(mean(all_growing_season_1[:, 1])*0.003)
                Ra_list_mean285_low.append(mean(all_growing_season_1[:, 2])*0.003)
                Sj_list_mean285_low.append(mean(all_growing_season_1[:, 3])*0.003)
                Sa_list_mean285_low.append(mean(all_growing_season_1[:, 4])*0.003)
                Sb_list_mean285_low.append(mean(all_growing_season_1[:, 5])*0.003)
                Sga_list_mean285_low.append(mean(all_growing_season_1[:, 6])*0.003)
                Sgb_list_mean285_low.append(mean(all_growing_season_1[:, 7])*0.003)
                Cj_list_mean285_low.append(mean(all_growing_season_1[:, 8])*0.003)
                Ca_list_mean285_low.append(mean(all_growing_season_1[:, 9])*0.003)
                Cb_list_mean285_low.append(mean(all_growing_season_1[:, 10])*0.003)
                Cga_list_mean285_low.append(mean(all_growing_season_1[:, 11])*0.003)
                Cgb_list_mean285_low.append(mean(all_growing_season_1[:, 12])*0.003)

                adult_clupeid285_low = (mean(all_growing_season_1[:,4]) + mean(all_growing_season_1[:,5]) + mean(all_growing_season_1[:,6]) + mean(all_growing_season_1[:, 7])) * 0.003
                clupeid_adult_list_mean285_low.append(adult_clupeid285_low)

                adult_cod285_low = (mean(all_growing_season_1[:,9]) + mean(all_growing_season_1[:,10]) + mean(all_growing_season_1[:,11]) + mean(all_growing_season_1[:, 12])) * 0.003
                cod_adult_list_mean285_low.append(adult_cod285_low)

            if T == 279.5:
                Rs_list_mean279_5_low.append(mean(all_growing_season_1[:, 0])*0.003)
                Rj_list_mean279_5_low.append(mean(all_growing_season_1[:, 1])*0.003)
                Ra_list_mean279_5_low.append(mean(all_growing_season_1[:, 2])*0.003)
                Sj_list_mean279_5_low.append(mean(all_growing_season_1[:, 3])*0.003)
                Sa_list_mean279_5_low.append(mean(all_growing_season_1[:, 4])*0.003)
                Sb_list_mean279_5_low.append(mean(all_growing_season_1[:, 5])*0.003)
                Sga_list_mean279_5_low.append(mean(all_growing_season_1[:, 6])*0.003)
                Sgb_list_mean279_5_low.append(mean(all_growing_season_1[:, 7])*0.003)
                Cj_list_mean279_5_low.append(mean(all_growing_season_1[:, 8])*0.003)
                Ca_list_mean279_5_low.append(mean(all_growing_season_1[:, 9])*0.003)
                Cb_list_mean279_5_low.append(mean(all_growing_season_1[:, 10])*0.003)
                Cga_list_mean279_5_low.append(mean(all_growing_season_1[:, 11])*0.003)
                Cgb_list_mean279_5_low.append(mean(all_growing_season_1[:, 12])*0.003)

                adult_clupeid279_5_low = (mean(all_growing_season_1[:,4]) + mean(all_growing_season_1[:,5]) + mean(all_growing_season_1[:,6]) + mean(all_growing_season_1[:, 7])) * 0.003
                clupeid_adult_list_mean279_5_low.append(adult_clupeid279_5_low)

                adult_cod279_5_low = (mean(all_growing_season_1[:,9]) + mean(all_growing_season_1[:,10]) + mean(all_growing_season_1[:,11]) + mean(all_growing_season_1[:, 12])) * 0.003
                cod_adult_list_mean279_5_low.append(adult_cod279_5_low)

            x = x + 1

    # iterating over 4 different temperatures, with high clupeid fishing mortality (Fs)
    for temperature in temperatureKelvin:
        T = temperature
        x = 0
        
        # iterating over the different cod fishing mortality values (Fc)
        for codMortalityPerYear in codMortality:
            # dividing the Fc value by 250, to get the Fc value per day
            Fc = codMortalityPerYear / 250 

            Delta = 0.1 * math.exp((Edelta * (T-T0))/(k*T*T0))
            Rsmax = 98 * math.exp((ERmax * (T-T0))/(k*T*T0))
            Rjmax = 1 * math.exp((ERmax * (T-T0))/(k*T*T0))
            Ramax = 0.75 * math.exp((ERmax * (T-T0))/(k*T*T0))

            # if it is the first year, the initial values will be reset 
            if x == 0 :
                Rs2 = 32.7816744728
                Rj2 = 0.4713946441
                Ra2 =  0.7175924292
                Sj2 = 10.4824915096
                Sa2 =  69.8208217446
                Sb2 =  89.7371125270
                Sga2 = 0
                Sgb2 = 0
                Cj2 =  6.2772314150
                Ca2 = 9.4197483660
                Cb2 = 0.6284148535
                Cga2 = 0
                Cgb2 = 0

            # if it is not the first year, the initial values will be the biomass values of the last reproductive season
            if x > 0:
                Rs2, Rj2, Ra2, Sj2, Sa2, Sb2, Sga2, Sgb2, Cj2, Ca2, Cb2, Cga2, Cgb2 = P02

            # iterating over 100 years for each Fc value
            for n in range(1, 100):

                # Fs = fishing mortality on clupeid (low = 0.2 per year = 0.0008 per day, high = 0.5 per year = 0.002 per day)
                Fs = 0.002

                tn_minus= (n-1) * Y
                tn_plus = n * Y 
                ts = np.linspace(tn_minus,tn_plus,Y)

                growing_season_2 = solve_ivp(function_growing_season, [tn_minus, tn_plus],  [Rs2, Rj2, Ra2, Sj2, Sa2, Sb2, Sga2, Sgb2, Cj2, Ca2, Cb2, Cga2, Cgb2], first_step = 0.01, max_step = 1.0, t_eval = ts, rtol = 1.0E-8, atol = 1.0E-10)

                # getting the useful data out of growing_season_1 and transposing it, so the rows and columns are placed in the right way
                data_2 = growing_season_2.y 
                data_2 = data_2.transpose()

                # the data of the first 50 years won't be stored
                if n > t_trans:
                    if n == t_trans + 1:
                        all_growing_season_2 = data_2
                    else:
                        all_growing_season_2 = np.vstack((all_growing_season_2, data_2))
                        
                if n == t_measure: 
                    break
                        
                # storing the data in lists 
                if n == 99 and T == 279.5:
                    Rs_list_mean279_5_high.append(mean(all_growing_season_2[:, 0])*0.003)
                    Rj_list_mean279_5_high.append(mean(all_growing_season_2[:, 1])*0.003)
                    Ra_list_mean279_5_high.append(mean(all_growing_season_2[:, 2])*0.003)
                    Sj_list_mean279_5_high.append(mean(all_growing_season_2[:, 3])*0.003)
                    Sa_list_mean279_5_high.append(mean(all_growing_season_2[:, 4])*0.003)
                    Sb_list_mean279_5_high.append(mean(all_growing_season_2[:, 5])*0.003)
                    Sga_list_mean279_5_high.append(mean(all_growing_season_2[:, 6])*0.003)
                    Sgb_list_mean279_5_high.append(mean(all_growing_season_2[:, 7])*0.003)
                    Cj_list_mean279_5_high.append(mean(all_growing_season_2[:, 8])*0.003)
                    Ca_list_mean279_5_high.append(mean(all_growing_season_2[:, 9])*0.003)
                    Cb_list_mean279_5_high.append(mean(all_growing_season_2[:, 10])*0.003)
                    Cga_list_mean279_5_high.append(mean(all_growing_season_2[:, 11])*0.003)
                    Cgb_list_mean279_5_high.append(mean(all_growing_season_2[:, 12])*0.003)

                    adult_clupeid279_5_high = (mean(all_growing_season_2[:,4]) + mean(all_growing_season_2[:,5]) + mean(all_growing_season_2[:,6]) + mean(all_growing_season_2[:, 7])) * 0.003
                    clupeid_adult_list_mean279_5_high.append(adult_clupeid279_5_high)

                    adult_cod279_5_high = (mean(all_growing_season_2[:,9]) + mean(all_growing_season_2[:,10]) + mean(all_growing_season_2[:,11]) + mean(all_growing_season_2[:, 12])) * 0.003
                    cod_adult_list_mean279_5_high.append(adult_cod279_5_high)

                
                if n == 99 and T == 280.7:
                    Rs_list_mean280_7_high.append(mean(all_growing_season_2[:, 0])*0.003)
                    Rj_list_mean280_7_high.append(mean(all_growing_season_2[:, 1])*0.003)
                    Ra_list_mean280_7_high.append(mean(all_growing_season_2[:, 2])*0.003)
                    Sj_list_mean280_7_high.append(mean(all_growing_season_2[:, 3])*0.003)
                    Sa_list_mean280_7_high.append(mean(all_growing_season_2[:, 4])*0.003)
                    Sb_list_mean280_7_high.append(mean(all_growing_season_2[:, 5])*0.003)
                    Sga_list_mean280_7_high.append(mean(all_growing_season_2[:, 6])*0.003)
                    Sgb_list_mean280_7_high.append(mean(all_growing_season_2[:, 7])*0.003)
                    Cj_list_mean280_7_high.append(mean(all_growing_season_2[:, 8])*0.003)
                    Ca_list_mean280_7_high.append(mean(all_growing_season_2[:, 9])*0.003)
                    Cb_list_mean280_7_high.append(mean(all_growing_season_2[:, 10])*0.003)
                    Cga_list_mean280_7_high.append(mean(all_growing_season_2[:, 11])*0.003)
                    Cgb_list_mean280_7_high.append(mean(all_growing_season_2[:, 12])*0.003)

                    adult_clupeid280_7_high = (mean(all_growing_season_2[:,4]) + mean(all_growing_season_2[:,5]) + mean(all_growing_season_2[:,6]) + mean(all_growing_season_2[:, 7])) * 0.003
                    clupeid_adult_list_mean280_7_high.append(adult_clupeid280_7_high)

                    adult_cod280_7_high = (mean(all_growing_season_2[:,9]) + mean(all_growing_season_2[:,10]) + mean(all_growing_season_2[:,11]) + mean(all_growing_season_2[:, 12])) * 0.003
                    cod_adult_list_mean280_7_high.append(adult_cod280_7_high)

                if n == 99 and T == 283:
                    Rs_list_mean283_high.append(mean(all_growing_season_2[:, 0])*0.003)
                    Rj_list_mean283_high.append(mean(all_growing_season_2[:, 1])*0.003)
                    Ra_list_mean283_high.append(mean(all_growing_season_2[:, 2])*0.003)
                    Sj_list_mean283_high.append(mean(all_growing_season_2[:, 3])*0.003)
                    Sa_list_mean283_high.append(mean(all_growing_season_2[:, 4])*0.003)
                    Sb_list_mean283_high.append(mean(all_growing_season_2[:, 5])*0.003)
                    Sga_list_mean283_high.append(mean(all_growing_season_2[:, 6])*0.003)
                    Sgb_list_mean283_high.append(mean(all_growing_season_2[:, 7])*0.003)
                    Cj_list_mean283_high.append(mean(all_growing_season_2[:, 8])*0.003)
                    Ca_list_mean283_high.append(mean(all_growing_season_2[:, 9])*0.003)
                    Cb_list_mean283_high.append(mean(all_growing_season_2[:, 10])*0.003)
                    Cga_list_mean283_high.append(mean(all_growing_season_2[:, 11])*0.003)
                    Cgb_list_mean283_high.append(mean(all_growing_season_2[:, 12])*0.003)

                    adult_clupeid283_high = (mean(all_growing_season_2[:,4]) + mean(all_growing_season_2[:,5]) + mean(all_growing_season_2[:,6]) + mean(all_growing_season_2[:, 7])) * 0.003
                    clupeid_adult_list_mean283_high.append(adult_clupeid283_high)

                    adult_cod283_high = (mean(all_growing_season_2[:,9]) + mean(all_growing_season_2[:,10]) + mean(all_growing_season_2[:,11]) + mean(all_growing_season_2[:, 12])) * 0.003
                    cod_adult_list_mean283_high.append(adult_cod283_high)

                if n == 99 and T == 285:
                    Rs_list_mean285_high.append(mean(all_growing_season_2[:, 0])*0.003)
                    Rj_list_mean285_high.append(mean(all_growing_season_2[:, 1])*0.003)
                    Ra_list_mean285_high.append(mean(all_growing_season_2[:, 2])*0.003)
                    Sj_list_mean285_high.append(mean(all_growing_season_2[:, 3])*0.003)
                    Sa_list_mean285_high.append(mean(all_growing_season_2[:, 4])*0.003)
                    Sb_list_mean285_high.append(mean(all_growing_season_2[:, 5])*0.003)
                    Sga_list_mean285_high.append(mean(all_growing_season_2[:, 6])*0.003)
                    Sgb_list_mean285_high.append(mean(all_growing_season_2[:, 7])*0.003)
                    Cj_list_mean285_high.append(mean(all_growing_season_2[:, 8])*0.003)
                    Ca_list_mean285_high.append(mean(all_growing_season_2[:, 9])*0.003)
                    Cb_list_mean285_high.append(mean(all_growing_season_2[:, 10])*0.003)
                    Cga_list_mean285_high.append(mean(all_growing_season_2[:, 11])*0.003)
                    Cgb_list_mean285_high.append(mean(all_growing_season_2[:, 12])*0.003)

                    adult_clupeid285_high = (mean(all_growing_season_2[:,4]) + mean(all_growing_season_2[:,5]) + mean(all_growing_season_2[:,6]) + mean(all_growing_season_2[:, 7])) * 0.003
                    clupeid_adult_list_mean285_high.append(adult_clupeid285_high)

                    adult_cod285_high = (mean(all_growing_season_2[:,9]) + mean(all_growing_season_2[:,10]) + mean(all_growing_season_2[:,11]) + mean(all_growing_season_2[:, 12])) * 0.003
                    cod_adult_list_mean285_high.append(adult_cod285_high)

                P02 = function_reproduction_season(data_2[-1, :])

                # setting the new initial values for next year to be the reproduction values of this year 
                Rs2, Rj2, Ra2, Sj2, Sa2, Sb2, Sga2, Sgb2, Cj2, Ca2, Cb2, Cga2, Cgb2 = P02

            x = x + 1

# calculating mean annual yield for both low and high Fs
# Fc is bifurcation parameter 
if biomass_or_yield == 3:

    # calculating mean annual yield, iterating over 4 different temperatures, with low clupeid fishing mortality (Fs)
    for temperature in temperatureKelvin:
        T = temperature
        x = 0

        for codMortalityPerYear in codMortality:
            # dividing the Fc value by 250, to get the Fc value per day
            Fc = codMortalityPerYear / 250 
            
            Delta = 0.1 * math.exp((Edelta * (T-T0))/(k*T*T0))
            Rsmax = 98 * math.exp((ERmax * (T-T0))/(k*T*T0))
            Rjmax = 1 * math.exp((ERmax * (T-T0))/(k*T*T0))
            Ramax = 0.75 * math.exp((ERmax * (T-T0))/(k*T*T0))
            
            # if it is the first year, the initial values will be reset 
            if x == 0:
                Rs1 = 32.7816744728
                Rj1 = 0.4713946441
                Ra1 =  0.7175924292
                Sj1 = 10.4824915096
                Sa1 =  69.8208217446
                Sb1 =  89.7371125270
                Sga1 = 0
                Sgb1 = 0
                Cj1 =  6.2772314150
                Ca1 = 9.4197483660
                Cb1 = 0.6284148535
                Cga1 = 0
                Cgb1 = 0

            # if it is not the first year, the initial values will be the biomass values of the last reproductive season
            if x > 0:
                Rs1, Rj1, Ra1, Sj1, Sa1, Sb1, Sga1, Sgb1, Cj1, Ca1, Cb1, Cga1, Cgb1 = P01

            # lists to store mean annual yield in for clupeid and cod (in total)
            totalYieldClupeids = []
            totalYieldCod = []

            for n in range(1, 100):
                
                # Fs = fishing mortality on clupeid (low = 0.2 per year = 0.0008 per day, high = 0.5 per year = 0.002 per day)
                Fs= 0.0008

                tn_minus= (n - 1) * Y
                tn_plus = n * Y 
                ts = np.linspace(tn_minus, tn_plus, Y)

                # ODE function , for 250 days, starting with the initial values 
                growing_season_1 = solve_ivp(function_growing_season, [tn_minus, tn_plus],  [Rs1, Rj1, Ra1, Sj1, Sa1, Sb1, Sga1, Sgb1, Cj1, Ca1, Cb1, Cga1, Cgb1], first_step = 0.01 ,max_step = 1.0, t_eval = ts, rtol = 1.0E-8, atol = 1.0E-10)
                
                # getting the useful data out of growing_season_1 and transposing it, so the rows and columns are placed in the right way
                data_1 = growing_season_1.y 
                data_1 = data_1.transpose()

                juvenileClup_yield = data_1[:,3]
                juvenileClupeidretentionAnnual = 0
                juvenileClupeidretentionDay = 0

                # calculating yield for juvenile clupeids
                for elem in juvenileClup_yield:
                    juvenileClupeidretentionDay = elem * RhoS * Fs 
                    juvenileClupeidretentionAnnual = juvenileClupeidretentionAnnual + juvenileClupeidretentionDay

                smallAdultClup_yield = data_1[:,[4,6]]
                largeAdultClup_yield = data_1[:,[5,7]]

                smallAdultClup_yield = np.sum(smallAdultClup_yield, axis = 1)
                largeAdultClup_yield = np.sum(largeAdultClup_yield, axis = 1)

                smallAdultClupretentionAnnual = 0
                smallAdultClupretentionDay = 0

                # calculating yield for small adult clupeids
                for elem in smallAdultClup_yield:
                    smallAdultClupretentionDay = elem * Fs
                    smallAdultClupretentionAnnual = smallAdultClupretentionAnnual + smallAdultClupretentionDay

                largeAdultClupretentionAnnual = 0
                largeAdultClupretentionDay = 0
                
                # calculating yield for large adult clupeids
                for elem in largeAdultClup_yield:
                    largeAdultClupretentionDay = elem * Fs
                    largeAdultClupretentionAnnual = largeAdultClupretentionAnnual + largeAdultClupretentionDay

                yieldClupeids = juvenileClupeidretentionAnnual + smallAdultClupretentionAnnual + largeAdultClupretentionAnnual
                
                totalYieldClupeids.append(yieldClupeids) 

                smallAdultCod_yield = data_1[:,[9,11]]
                largeAdultCod_yield = data_1[:,[10,12]]

                smallAdultCod_yield = np.sum(smallAdultCod_yield, axis = 1)
                largeAdultCod_yield = np.sum(largeAdultCod_yield, axis = 1)

                smallAdultCodretentionAnnual = 0
                smallAdultCodretentionDay = 0
                
                # calculating yield for small adult cods
                for elem in smallAdultCod_yield:
                    smallAdultCodretentionDay = elem * Fc
                    smallAdultCodretentionAnnual = smallAdultCodretentionAnnual + smallAdultCodretentionDay

                largeAdultCodretentionAnnual = 0
                largeAdultCodretentionDay = 0
                
                # calculating yield for large adult cods
                for elem in largeAdultCod_yield:
                    largeAdultCodretentionDay = elem * Fc
                    largeAdultCodretentionAnnual = largeAdultCodretentionAnnual + largeAdultCodretentionDay

                yieldCod = smallAdultCodretentionAnnual + largeAdultCodretentionAnnual
                
                totalYieldCod.append(yieldCod) 

                if n > t_trans: 
                    if n == t_trans + 1:
                        all_growing_season_1 = data_1     
                    else:
                        all_growing_season_1 = np.vstack((all_growing_season_1, data_1))
                            
                if n == t_measure: 
                    break

                # using the reproduction function with the values of the last day of the growing season 
                P01 = function_reproduction_season(data_1[-1, :])

                # setting the new initial values for next year to be the reproduction values of this year 
                Rs1, Rj1, Ra1, Sj1, Sa1, Sb1, Sga1, Sgb1, Cj1, Ca1, Cb1, Cga1, Cgb1 = P01
            
            # storing the data in lists 
            if T == 279.5:
                meanTotalYieldClup_temp1_low.append(mean(totalYieldClupeids))
                meanTotalYieldCod_temp1_low.append(mean(totalYieldCod))

            if T == 280.7:
                meanTotalYieldClup_temp2_low.append(mean(totalYieldClupeids))
                meanTotalYieldCod_temp2_low.append(mean(totalYieldCod))

            if T == 283:
                meanTotalYieldClup_temp3_low.append(mean(totalYieldClupeids))
                meanTotalYieldCod_temp3_low.append(mean(totalYieldCod))

            if T == 285:
                meanTotalYieldClup_temp4_low.append(mean(totalYieldClupeids))
                meanTotalYieldCod_temp4_low.append(mean(totalYieldCod))

            x = x + 1


    # calculation mean annual yield, for 4 different temperatures and high Fs
    for temperature in temperatureKelvin:
        T = temperature
        x = 0
        
        for codMortalityPerYear in codMortality:
            
            # dividing the Fc value by 250, to get the Fc value per day
            Fc = codMortalityPerYear / 250 

            Delta = 0.1 * math.exp((Edelta * (T-T0))/(k*T*T0))
            Rsmax = 98 * math.exp((ERmax * (T-T0))/(k*T*T0))
            Rjmax = 1 * math.exp((ERmax * (T-T0))/(k*T*T0))
            Ramax = 0.75 * math.exp((ERmax * (T-T0))/(k*T*T0))

            # if it is the first year, the initial values will be reset 
            if x == 0 :
                Rs2 = 32.7816744728
                Rj2 = 0.4713946441
                Ra2 =  0.7175924292
                Sj2 = 10.4824915096
                Sa2 =  69.8208217446
                Sb2 =  89.7371125270
                Sga2 = 0
                Sgb2 = 0
                Cj2 =  6.2772314150
                Ca2 = 9.4197483660
                Cb2 = 0.6284148535
                Cga2 = 0
                Cgb2 = 0

            # if it is not the first year, the initial values will be the biomass values of the last reproductive season
            if x > 0:
                Rs2, Rj2, Ra2, Sj2, Sa2, Sb2, Sga2, Sgb2, Cj2, Ca2, Cb2, Cga2, Cgb2 = P02

            # lists to store mean annual yield in for clupeid and cod (in total)
            totalYieldClupeids2 = []
            totalYieldCod2 = []

            for n in range(1, 100):

                # Fs = fishing mortality on clupeid (low = 0.2 per year = 0.0008 per day, high = 0.5 per year = 0.002 per day)
                Fs = 0.002

                tn_minus= (n-1) * Y
                tn_plus = n * Y 
                ts = np.linspace(tn_minus,tn_plus,Y)

                growing_season_2 = solve_ivp(function_growing_season, [tn_minus, tn_plus],  [Rs2, Rj2, Ra2, Sj2, Sa2, Sb2, Sga2, Sgb2, Cj2, Ca2, Cb2, Cga2, Cgb2], first_step = 0.01, max_step = 1.0, t_eval = ts, rtol = 1.0E-8, atol = 1.0E-10)

                # getting the useful data out of growing_season_1 and transposing it, so the rows and columns are placed in the right way
                data_2 = growing_season_2.y 
                data_2 = data_2.transpose()

                juvenileClup_yield2 = data_2[:,3]
                juvenileClupeidretentionAnnual2 = 0
                juvenileClupeidretentionDay2 = 0

                # calculating yield for juvenile clupeids 
                for elem in juvenileClup_yield2:
                    juvenileClupeidretentionDay2 = elem * RhoS * Fs 
                    juvenileClupeidretentionAnnual2 = juvenileClupeidretentionAnnual2 + juvenileClupeidretentionDay2

                smallAdultClup_yield2 = data_2[:,[4,6]]
                largeAdultClup_yield2 = data_2[:,[5,7]]

                smallAdultClup_yield2 = np.sum(smallAdultClup_yield2, axis = 1)
                largeAdultClup_yield2 = np.sum(largeAdultClup_yield2, axis = 1)

                smallAdultClupretentionAnnual2 = 0
                smallAdultClupretentionDay2 = 0

                # calculating yield for small adult clupeids
                for elem in smallAdultClup_yield2:
                    smallAdultClupretentionDay2 = elem * Fs
                    smallAdultClupretentionAnnual2 = smallAdultClupretentionAnnual2 + smallAdultClupretentionDay2

                largeAdultClupretentionAnnual2 = 0
                largeAdultClupretentionDay2 = 0
                
                # calculating yield for large adult clupeids
                for elem in largeAdultClup_yield2:
                    largeAdultClupretentionDay2 = elem * Fs
                    largeAdultClupretentionAnnual2 = largeAdultClupretentionAnnual2 + largeAdultClupretentionDay2

                yieldClupeids2 = juvenileClupeidretentionAnnual2 + smallAdultClupretentionAnnual2 + largeAdultClupretentionAnnual2
                
                totalYieldClupeids2.append(yieldClupeids2) 

                smallAdultCod_yield2 = data_2[:,[9,11]]
                largeAdultCod_yield2 = data_2[:,[10,12]]

                smallAdultCod_yield2 = np.sum(smallAdultCod_yield2, axis = 1)
                largeAdultCod_yield2 = np.sum(largeAdultCod_yield2, axis = 1)

                smallAdultCodretentionAnnual2 = 0
                smallAdultCodretentionDay2 = 0

                # calculating yield for small adult cods
                for elem in smallAdultCod_yield2:
                    smallAdultCodretentionDay2 = elem * Fc
                    smallAdultCodretentionAnnual2 = smallAdultCodretentionAnnual2 + smallAdultCodretentionDay2

                largeAdultCodretentionAnnual2 = 0
                largeAdultCodretentionDay2 = 0
                
                # calculating yield for large adult cods
                for elem in largeAdultCod_yield2:
                    largeAdultCodretentionDay2 = elem * Fc
                    largeAdultCodretentionAnnual2 = largeAdultCodretentionAnnual2 + largeAdultCodretentionDay2

                yieldCod2 = smallAdultCodretentionAnnual2 + largeAdultCodretentionAnnual2
                
                totalYieldCod2.append(yieldCod2) 

                if n > t_trans: 
                    if n == t_trans + 1:
                        all_growing_season_2 = data_2
                                                
                    else:
                        all_growing_season_2 = np.vstack((all_growing_season_2, data_2))

                if n == t_measure: 
                    break
                        
                P02 = function_reproduction_season(data_2[-1, :])

                # setting the new initial values for next year to be the reproduction values of this year 
                Rs2, Rj2, Ra2, Sj2, Sa2, Sb2, Sga2, Sgb2, Cj2, Ca2, Cb2, Cga2, Cgb2 = P02

            # storing the data in lists 
            if T == 279.5:
                meanTotalYieldClup_temp1_high.append(mean(totalYieldClupeids2))
                meanTotalYieldCod_temp1_high.append(mean(totalYieldCod2))
        
            if T == 280.7:
                meanTotalYieldClup_temp2_high.append(mean(totalYieldClupeids2))
                meanTotalYieldCod_temp2_high.append(mean(totalYieldCod2))

            if T == 283:
                meanTotalYieldClup_temp3_high.append(mean(totalYieldClupeids2))
                meanTotalYieldCod_temp3_high.append(mean(totalYieldCod2))

            if T == 285:
                meanTotalYieldClup_temp4_high.append(mean(totalYieldClupeids2))
                meanTotalYieldCod_temp4_high.append(mean(totalYieldCod2))

            x = x + 1


#########################################################################################################################################################################################
###                                                              storing the data in CSV files                                                                                                         ###
#########################################################################################################################################################################################

# checking if biomass with temperature as bifurcation parameter is calculated , from low to high T
if biomass_or_yield == 0:
    temperature_data_low = {"Rs": Rs_list_low, 'Rj': Rj_list_low, 'Ra': Ra_list_low, 'Sj': Sj_list_low, 'Sa': Sa_list_low, 'Sb': Sb_list_low, 'Sga': Sga_list_low, 'Sgb': Sgb_list_low, 'Cj': Cj_list_low, 'Ca': Ca_list_low, 'Cb': Cb_list_low, 'Cga': Cga_list_low, 'Cgb': Cgb_list_low, 'Clupeid adult': clupeid_adult_list_low, 'Cod adult': cod_adult_list_low}
    temperature_data_high = {"Rs": Rs_list_high, 'Rj': Rj_list_high, 'Ra': Ra_list_high, 'Sj': Sj_list_high, 'Sa': Sa_list_high, 'Sb': Sb_list_high, 'Sga': Sga_list_high, 'Sgb': Sgb_list_high, 'Cj': Cj_list_high, 'Ca': Ca_list_high, 'Cb': Cb_list_high, 'Cga': Cga_list_high, 'Cgb': Cgb_list_high, 'Clupeid adult': clupeid_adult_list_high, 'Cod adult': cod_adult_list_high}

    dataframe_low = pd.DataFrame(temperature_data_low)
    dataframe_high = pd.DataFrame(temperature_data_high)

    dataframe_low.to_csv('results_low_clupeid_mortality4')
    dataframe_high.to_csv('results_high_clupeid_mortality4')

# checking if temperature is used as bifurcation parameter, with high to low T values
if biomass_or_yield == 1:
    temperature_data_low_reverse = {"Rs": Rs_list_low_reverse, 'Rj': Rj_list_low_reverse, 'Ra': Ra_list_low_reverse, 'Sj': Sj_list_low_reverse, 'Sa': Sa_list_low_reverse, 'Sb': Sb_list_low_reverse, 'Sga': Sga_list_low_reverse, 'Sgb': Sgb_list_low_reverse, 'Cj': Cj_list_low_reverse, 'Ca': Ca_list_low_reverse, 'Cb': Cb_list_low_reverse, 'Cga': Cga_list_low_reverse, 'Cgb': Cgb_list_low_reverse, 'Clupeid adult': clupeid_adult_list_low_reverse, 'Cod adult': cod_adult_list_low_reverse}
    temperature_data_high_reverse = {"Rs": Rs_list_high_reverse, 'Rj': Rj_list_high_reverse, 'Ra': Ra_list_high_reverse, 'Sj': Sj_list_high_reverse, 'Sa': Sa_list_high_reverse, 'Sb': Sb_list_high_reverse, 'Sga': Sga_list_high_reverse, 'Sgb': Sgb_list_high_reverse, 'Cj': Cj_list_high_reverse, 'Ca': Ca_list_high_reverse, 'Cb': Cb_list_high_reverse, 'Cga': Cga_list_high_reverse, 'Cgb': Cgb_list_high_reverse, 'Clupeid adult': clupeid_adult_list_high_reverse, 'Cod adult': cod_adult_list_high_reverse}

    dataframe_low_reverse = pd.DataFrame(temperature_data_low_reverse)
    dataframe_high_reverse = pd.DataFrame(temperature_data_high_reverse)

    dataframe_low_reverse.to_csv('results_low_clupeid_mortality_reverse')
    dataframe_high_reverse.to_csv('results_high_clupeid_mortality_reverse')
    
# checking if biomass is calculated with cod fishing mortality as bifurcation parameter
if biomass_or_yield == 2:
    
    # first temp, low and high clupeid fishing mortality 
    results_dictionary_low_1 = {"Rs": Rs_list_mean279_5_low, 'Rj': Rj_list_mean279_5_low, 'Ra': Ra_list_mean279_5_low, 'Sj': Sj_list_mean279_5_low, 'Sa': Sa_list_mean279_5_low, 'Sb': Sb_list_mean279_5_low, 'Sga': Sga_list_mean279_5_low, 'Sgb': Sgb_list_mean279_5_low, 'Cj': Cj_list_mean279_5_low, 'Ca': Ca_list_mean279_5_low, 'Cb': Cb_list_mean279_5_low, 'Cga': Cga_list_mean279_5_low, 'Cgb': Cgb_list_mean279_5_low, 'Clupeid adult': clupeid_adult_list_mean279_5_low, 'Cod adult': cod_adult_list_mean279_5_low}
    results_dictionary_high_1 = {"Rs": Rs_list_mean279_5_high, 'Rj': Rj_list_mean279_5_high, 'Ra': Ra_list_mean279_5_high, 'Sj': Sj_list_mean279_5_high, 'Sa': Sa_list_mean279_5_high, 'Sb': Sb_list_mean279_5_high, 'Sga': Sga_list_mean279_5_high, 'Sgb': Sgb_list_mean279_5_high, 'Cj': Cj_list_mean279_5_high, 'Ca': Ca_list_mean279_5_high, 'Cb': Cb_list_mean279_5_high, 'Cga': Cga_list_mean279_5_high, 'Cgb': Cgb_list_mean279_5_high, 'Clupeid adult': clupeid_adult_list_mean279_5_high, 'Cod adult': cod_adult_list_mean279_5_high}

    dataframe_low_1 = pd.DataFrame(results_dictionary_low_1)
    dataframe_high_1 = pd.DataFrame(results_dictionary_high_1)

    dataframe_low_1.to_csv('results_low_clupeid_mortality_temp1')
    dataframe_high_1.to_csv('results_high_clupeid_mortality_temp1')

    # second temp, low and high clupeid fishing mortality 
    results_dictionary_low_2 = {"Rs": Rs_list_mean280_7_low, 'Rj': Rj_list_mean280_7_low, 'Ra': Ra_list_mean280_7_low, 'Sj': Sj_list_mean280_7_low, 'Sa': Sa_list_mean280_7_low, 'Sb': Sb_list_mean280_7_low, 'Sga': Sga_list_mean280_7_low, 'Sgb': Sgb_list_mean280_7_low, 'Cj': Cj_list_mean280_7_low, 'Ca': Ca_list_mean280_7_low, 'Cb': Cb_list_mean280_7_low, 'Cga': Cga_list_mean280_7_low, 'Cgb': Cgb_list_mean280_7_low, 'Clupeid adult': clupeid_adult_list_mean280_7_low, 'Cod adult': cod_adult_list_mean280_7_low}
    results_dictionary_high_2 = {"Rs": Rs_list_mean280_7_high, 'Rj': Rj_list_mean280_7_high, 'Ra': Ra_list_mean280_7_high, 'Sj': Sj_list_mean280_7_high, 'Sa': Sa_list_mean280_7_high, 'Sb': Sb_list_mean280_7_high, 'Sga': Sga_list_mean280_7_high, 'Sgb': Sgb_list_mean280_7_high, 'Cj': Cj_list_mean280_7_high, 'Ca': Ca_list_mean280_7_high, 'Cb': Cb_list_mean280_7_high, 'Cga': Cga_list_mean280_7_high, 'Cgb': Cgb_list_mean280_7_high, 'Clupeid adult': clupeid_adult_list_mean280_7_high, 'Cod adult': cod_adult_list_mean280_7_high}

    dataframe_low_2 = pd.DataFrame(results_dictionary_low_2)
    dataframe_high_2 = pd.DataFrame(results_dictionary_high_2)

    dataframe_low_2.to_csv('results_low_clupeid_mortality_temp2')
    dataframe_high_2.to_csv('results_high_clupeid_mortality_temp2')

    # third temp, low and high clupeid fishing mortality 
    results_dictionary_low_3 = {"Rs": Rs_list_mean283_low, 'Rj': Rj_list_mean283_low, 'Ra': Ra_list_mean283_low, 'Sj': Sj_list_mean283_low, 'Sa': Sa_list_mean283_low, 'Sb': Sb_list_mean283_low, 'Sga': Sga_list_mean283_low, 'Sgb': Sgb_list_mean283_low, 'Cj': Cj_list_mean283_low, 'Ca': Ca_list_mean283_low, 'Cb': Cb_list_mean283_low, 'Cga': Cga_list_mean283_low, 'Cgb': Cgb_list_mean283_low, 'Clupeid adult': clupeid_adult_list_mean283_low, 'Cod adult': cod_adult_list_mean283_low}
    results_dictionary_high_3 = {"Rs": Rs_list_mean283_high, 'Rj': Rj_list_mean283_high, 'Ra': Ra_list_mean283_high, 'Sj': Sj_list_mean283_high, 'Sa': Sa_list_mean283_high, 'Sb': Sb_list_mean283_high, 'Sga': Sga_list_mean283_high, 'Sgb': Sgb_list_mean283_high, 'Cj': Cj_list_mean283_high, 'Ca': Ca_list_mean283_high, 'Cb': Cb_list_mean283_high, 'Cga': Cga_list_mean283_high, 'Cgb': Cgb_list_mean283_high, 'Clupeid adult': clupeid_adult_list_mean283_high, 'Cod adult': cod_adult_list_mean283_high}

    dataframe_low_3 = pd.DataFrame(results_dictionary_low_3)
    dataframe_high_3= pd.DataFrame(results_dictionary_high_3)

    dataframe_low_3.to_csv('results_low_clupeid_mortality_temp3')
    dataframe_high_3.to_csv('results_high_clupeid_mortality_temp3')

    # fourth temp, low and high clupeid fishing mortality 
    results_dictionary_low_4 = {"Rs": Rs_list_mean285_low, 'Rj': Rj_list_mean285_low, 'Ra': Ra_list_mean285_low, 'Sj': Sj_list_mean285_low, 'Sa': Sa_list_mean285_low, 'Sb': Sb_list_mean285_low, 'Sga': Sga_list_mean285_low, 'Sgb': Sgb_list_mean285_low, 'Cj': Cj_list_mean285_low, 'Ca': Ca_list_mean285_low, 'Cb': Cb_list_mean285_low, 'Cga': Cga_list_mean285_low, 'Cgb': Cgb_list_mean285_low, 'Clupeid adult': clupeid_adult_list_mean285_low, 'Cod adult': cod_adult_list_mean285_low}
    results_dictionary_high_4 = {"Rs": Rs_list_mean285_high, 'Rj': Rj_list_mean285_high, 'Ra': Ra_list_mean285_high, 'Sj': Sj_list_mean285_high, 'Sa': Sa_list_mean285_high, 'Sb': Sb_list_mean285_high, 'Sga': Sga_list_mean285_high, 'Sgb': Sgb_list_mean285_high, 'Cj': Cj_list_mean285_high, 'Ca': Ca_list_mean285_high, 'Cb': Cb_list_mean285_high, 'Cga': Cga_list_mean285_high, 'Cgb': Cgb_list_mean285_high, 'Clupeid adult': clupeid_adult_list_mean285_high, 'Cod adult': cod_adult_list_mean285_high}

    dataframe_low_4 = pd.DataFrame(results_dictionary_low_4)
    dataframe_high_4 = pd.DataFrame(results_dictionary_high_4)

    dataframe_low_4.to_csv('results_low_clupeid_mortality_temp4')
    dataframe_high_4.to_csv('results_high_clupeid_mortality_temp4')

# checking if it calculated yield
if biomass_or_yield == 3:
    results_yield_low_cod = {'temp1': meanTotalYieldCod_temp1_low, 'temp2': meanTotalYieldCod_temp2_low, 'temp3': meanTotalYieldCod_temp3_low, 'temp4': meanTotalYieldCod_temp4_low}
    results_yield_high_cod = {'temp1': meanTotalYieldCod_temp1_high, 'temp2': meanTotalYieldCod_temp2_high, 'temp3': meanTotalYieldCod_temp3_high, 'temp4': meanTotalYieldCod_temp4_high}
    results_yield_low_clup = {'temp1': meanTotalYieldClup_temp1_low, 'temp2': meanTotalYieldClup_temp2_low, 'temp3': meanTotalYieldClup_temp3_low, 'temp4': meanTotalYieldClup_temp4_low}
    results_yield_high_clup = {'temp1': meanTotalYieldClup_temp1_high, 'temp2': meanTotalYieldClup_temp2_high, 'temp3': meanTotalYieldClup_temp3_high, 'temp4': meanTotalYieldClup_temp4_high}

    dataframe_yield_low_cod = pd.DataFrame(results_yield_low_cod)
    dataframe_yield_high_cod = pd.DataFrame(results_yield_high_cod)
    dataframe_yield_low_clup = pd.DataFrame(results_yield_low_clup)
    dataframe_yield_high_clup = pd.DataFrame(results_yield_high_clup)

    dataframe_yield_low_cod.to_csv('results yield low cod')
    dataframe_yield_high_cod.to_csv('results yield high cod')
    dataframe_yield_low_clup.to_csv('results yield low clup')
    dataframe_yield_high_clup.to_csv('results yield high clup')