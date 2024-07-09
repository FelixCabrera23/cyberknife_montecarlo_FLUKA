# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 09:38:22 2024

@author: Felix Cabrera

GAMMA ANALYSIS OF MONTECARLO DATA AGAINS OCR DATA

THIS SCRIPT IS MENT TO BE USED ONLY WITH DATA FROM THE ADECUATE SOURCE AND MAY
NOT WORK IF USED WITH DIFFERENT FORMAT DATA.

The analysis is using data located in the file named CK2_DANE.xlsx, this file
must be in the same directory as the script. This file must not be modified
or the results should be not trusted.


"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d


df_org = pd.read_excel("CK2_DANE.xlsx")

#Defining global variables
depths = [15,100,300] #depths in mm
Radious_OCR = []
Dose_OCR = []
Radious_simul = []
Dose_simul = []
Gamma_map = []
Result_matrix = []

#Cleaning the data frame from not necesary data
df_org = df_org.drop(index=[0]) #this line drops the row where there is the radious info
#df_org = df_org.dropna() #this line drops al the NaN values

###### function definition ######
def Data_Specific(coll,depth):
    """
    This func takes collimator option and depth option and returns
    radious and dose date in form of lists

    >>> a,b = Data_especific(10,depths[1])
    >>> print(b)
    >>> [1.0,
         0.9993,
         0.9983,
         0.9974,...
    """
    index_x = coll
    index_y = str(coll)
    if depth==15: index_y=index_y+".1"
    elif depth==100: index_y=index_y+".2"
    elif depth==300: index_y=index_y+".3"
    else: print("error: data not found")

    # rad_t = df_org[index_x].tolist()
    # dose_t = df_org[index_y].tolist()
    # This has been changed to drop the NaN values in the lists
    
    rad_t = df_org[index_x].dropna().tolist()
    dose_t = df_org[index_y].dropna().tolist()
    
    return rad_t,dose_t

def read_data_file (file_path):
    """
    This func taces the file path of the data from the FLUKA simulation and 
    returns it into a matrix form
    
    ignores commented lines that beggin with '#'
    
    works with FLUKA 
    """ 
    data_temp = []
    with open(file_path, 'r') as file:
        for line in file:
            data_temp.append(line.strip())
    
    Dat_mat = []

    for line in data_temp:
        list_temp = line.split()
        if list_temp[0] != '#':        
            float_temp = [float(string) for string in list_temp]
            Dat_mat.append(float_temp)
        
    return(Dat_mat)

def Data_simulated (Data_matrix,norm=True, scale =10, norm_val = 1):
    """
    This function takes the data matrix from FLUKA 
    and returns list of the radious, dose, and error data
    
    the 'scale' changes the scale of the data in the x coordinate
    the 'norme' normalizes the measurment result to 'norm_val'

    """
  
    x1_tmp,x2_tmp,meas_tmp,error_tmp = [],[],[],[]

    for line in Data_matrix:
        x1_tmp.append(line[0]*scale)
        x2_tmp.append(line[1]*scale)
        meas_tmp.append(line[2])
        error_tmp.append(line[3])
        
    if norm: meas_tmp = normalize(meas_tmp, norm_Val=norm_val)
        
    return(x1_tmp,x2_tmp,meas_tmp,error_tmp)

def normalize (list_temp,norm_Val = 1):
    """
    Esta funcion normaliza una lista a 1
    """
    norm_list = []
    max_list = max(list_temp)/norm_Val
    for i in list_temp: 
        norm_list.append(i/max_list)
        
    return(norm_list)

def interpolate (rad,dose,rad_n):
    """
    This func takes lists of original data and extrapolates it to specific
    radious needed to analize the data comparing it with specific data from
    the simulation
    rad = should be the radious from the OCR
    dose = should be the dose of the ORC
    rad_n = should be the new radious of where we dont have specific data
    
    the result dose_n is the dose in the specific new radious rad_n
    
    """
    f = interp1d(rad, dose,kind='linear', fill_value='extrapolate')
    dose_n = f(rad_n)
    return(dose_n)

def gamma_index (dose_meas, rad_meas, dose_simul, rad_simul,dta = 2, dd = 3, dose_threshold = 0):
    """
    This func generates the Gamma map for the distributions given by 
    Note that this func asumes normalized doses
    """
    Gamma_map = []
    
    for i in range(len(dose_simul)):
        if dose_simul[i] >= dose_threshold:
            Gamma_min = np.inf
            for j in range(len(rad_meas)):
                distance = np.abs(rad_simul[i]-rad_meas[j])
                dose_diff = np.abs(dose_simul[i]-dose_meas[j])*100
                Gamma_value = np.sqrt((distance/dta)**2 + (dose_diff/dd)**2)
                if Gamma_value < Gamma_min:
                    Gamma_min = Gamma_value
            Gamma_map.append(Gamma_min)
            
    return(Gamma_map)

def Score_GPR (gamma_list):
    """
    This takes the gamma list and returns score
    """
    i = 0
    for number in gamma_list:
        if number <= 1: i+=1
    
    score = (i/len(gamma_list))*100
    
    return(score)

def G_Analize_data_file (file = "dane.dat", coll = 5, depth = 100):
    """
    This func runs the gamma analysis using as imput a data file from FLUKA
    and compares with OCR data. Is important to know that details of the data
    must be included as imput
    coll = collimator diameter in mm
    depth = depth of the measurment
    Data_simulated 
    WARNING! failure to include the correct details or data file can give not
    real results
    """
    if file == "dane.dat": print("WARNING! EXAMPLE DATA IS IN USE")
    
    global Radious_OCR, Dose_OCR, Radious_simul, Dose_simul, Gamma_map
    
    # Reading imput data
    Raw_data_matrix = read_data_file(file)
    R1, R2, dose, dose_err = Data_simulated(Raw_data_matrix)
    Radious_simul = R1
    Dose_simul = dose
    
    # Retriving original profile data
    Radious_OCR, Dose_OCR = Data_Specific(coll, depth)
    
    # Gamma analisis and scoring
    Gamma_map = gamma_index(Dose_OCR, Radious_OCR, Dose_simul, Radious_simul)
    Score = Score_GPR(Gamma_map)
    
    print("Gamma analysis completed, GPR = "+str(Score))
    
    global Result_matrix    
    Result_matrix = [Radious_OCR, Dose_OCR, Radious_simul, Dose_simul, Gamma_map]
    
    return(Result_matrix,Score)
    
def Ploting_gamma_results (data_matrix = None, norm = True):
    """
    This function prints the analisis of the global variables
    it will not work if no analyzis has been runed
    
    Another data can be ploted if is stored in proper matrix form

    """
    if data_matrix is None:
        data_matrix = Result_matrix
        
    R1,Dose1,R2,Dose2,Gamma_tmp = [list(data) for data in data_matrix]
    
    if norm:
        Gamma_plot = normalize(Gamma_tmp)
    else:
        Gamma_plot = Gamma_tmp
    
    plt.plot(R2,Dose2, label = 'MonteCarlo')        
    plt.plot(R1,Dose1, label = 'OCR')
    plt.plot(R2, Gamma_plot, label = 'Gamma')
    plt.xlabel(r"Radious [mm]")
    plt.ylabel(r"normalized dose")
    plt.legend()
    plt.show()
    
def LRD_Analize (data_matrix = None, norm = True):
    """
    This func plots the LRD from a dataset provided in matrix form
    the data must be analized first by the gamma analysis
    
    """
    if data_matrix is None:
        data_matrix = Result_matrix
        
    R1,Dose1,R2,Dose2,Gamma_tmp = [list(data) for data in data_matrix]
    
    mean_dose = interpolate(R1,Dose1,R2)
    
    LRD = []
    
    for i in range(len(Dose2)):
        if mean_dose[i] == 0:
            LRD.append(Dose2[i])
        else:
            LRD_i = np.abs(Dose2[i]-mean_dose[i])/mean_dose[i]
        LRD.append(LRD_i)
        
    if norm:
        LRD = normalize(LRD)

    plt.plot(R2,Dose2, label = 'MonteCarlo')    
    plt.plot(R2,mean_dose, label = 'OCR')
    plt.plot(R2,LRD,label = 'LRD')
    plt.xlabel(r"Radious [mm]")
    plt.ylabel(r"normalized dose")
    plt.legend()
    plt.show()
    
def Energy_spectra_plot (file_path = 'Sample_data/Energy_flux_before.lis', x_log = False, y_log = False ,norm = True):
    """
    This func plots the energy spectrum provided the data from a FLUKA '.lis' file
    with adecuate format

    ParametersData_simulated 
    ----------
    file_path : flie path to the flux data
        DESCRIPTION. The default is 'Sample_data/Energy_flux_before.lis'.
    x_log : boolean sets logarithmic scale on x axis
    y_log : boolean sets logarithmic scale on y axis
    norm : boolean sets the

    """
    data_flux = read_data_file(file_path)
    E1,E2,Flx,Flx_err = Data_simulated(data_flux, norm=norm,scale=1000)
    
    if x_log:
        plt.xscale('log')
        
    if y_log:
        plt.yscale('log')
        
    plt.plot(E1,Flx)
    plt.xlabel(r"MeV")
    plt.ylabel(r"Photon fluence")
    plt.show()
    
def PDD_plot (file_path = 'Sample_data/Coll_60mm_PDD.dat',norm=True):
    
    
    pdd_data = read_data_file(file_path)
    
    Z1, Z2, Dose, Dose_err = Data_simulated(pdd_data,norm=norm,scale=10)
    
    plt.plot(Z1,Dose)
    plt.xlabel(r'Distance [mm]')
    plt.ylabel(r'Normalized Dose')
    plt.show()
    
    
    
    
    
    
    
    
            
# ######### Test area #########    
# G_Analize_data_file()
# Ploting_gamma_results()
# LRD_Analize()
# Energy_spectra_plot()
