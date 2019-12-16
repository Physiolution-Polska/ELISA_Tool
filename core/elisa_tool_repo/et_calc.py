#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 2019
@author:  Edward Maievskij, Marek Bawiec, Grzegorz Banach
The ELISA tool plugin: module contains basic calculating functions.
"""

__author__ = "Marek Bawiec, Grzegorz Banach, Edward Maievskij"
__copyright__ = "Copyright 2019, Physiolution Polska"
__credits__ = [""]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Grzegorz Banach"
__email__ = "g.banach@physiolution.pl"
__status__ = "Production"

import numpy as np
import pandas as pd
import copy
from scipy.optimize import differential_evolution
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
from RegscorePy import *

def max_concentration(res_sam_mat, std_mat):
    """
    The function implemented to return the maximal known concentration from the whole data set for appropriate plotting.
    """
    res_sam=res_sam_mat[1]
    #print("res_sam:", res_sam)

    std=std_mat[std_mat.columns[0]]
    #print("std:", std)
    
    max_res_sam=max(res_sam)
    max_std= max(std)
    
    max_conc=0.
    if max_res_sam > max_std:
        max_conc=max_res_sam
    else:
        max_conc=max_std
    """returns the value of maximal titer concentration"""    
    return max_conc

def min_concentration(res_sam_mat, std_mat):
    """
    The function implemented to return the minimal known concentration from the whole data set for appropriate plotting.
    """
    res_sam=res_sam_mat[1]
    
    std=std_mat[std_mat.columns[0]]
    
    min_res_sam=min(res_sam)
    min_std= min(std)
    
    min_conc=0.
    #i am not sure about > or < sign
    if min_res_sam < min_std:
        min_conc=min_res_sam
    else:
        min_conc=min_std
    """
    returns the value of minimal titer concentration
    """    
    return(min_conc)


def theor_X(res_sam_mat, std_mat):
    """
    The function implemented for theoretical curve X values generating.
    """
    maxconc = max_concentration(res_sam_mat=res_sam_mat, std_mat=std_mat)
    minconc = min_concentration(res_sam_mat=res_sam_mat, std_mat=std_mat)
    #1000 no  points for theory curves between mincon and maxconc
    delta = (maxconc - minconc)/1000
    X_t = np.arange(minconc, maxconc, delta)
    #returns list of titer concentration valuesfor further calculations of theoretical curve
    return X_t
    
def data_st_to_print (imp_data, imp_data_map, data_standards):
    """
    The function used to generate an ordered DataFrame with initial measuremnt results.
    """
    g1=imp_data_map.stack()
    g2=set(g1)
    for g in g2:
        vars()[g]=[]
        
    imp_columns=list(imp_data.columns.values)
    imp_index=list(imp_data.index.values)
    
    for c_name in imp_columns:
        for r_name in imp_index:
            vars()[imp_data_map.loc[r_name][c_name]].append(imp_data.loc[r_name][c_name])
    sample=[]        
    X_st=[]
    Y_st=[]
    
    for g in g2:
        vars()[g]=sum(vars()[g]) / float(len(vars()[g]))
        if ('std' in g):
            sample.append(g)
            X_st.append(data_standards.loc[g]['values'])
            Y_st.append(vars()[g])
            continue
        else:
            continue
    
    result=pd.DataFrame()
    result["name"]=sample
    result["conc"]=X_st
    result["absorbance"]=Y_st
    
    result.sort_values(by=["name"], inplace=True)
    result.set_index(['name'], inplace=True)
    """returns DataFrame with extended table of data about calibration standards: (std_name,concentration, abs_ave) """
    return (result)
        
#ln function stuff
def ln_func_par_opt(D_st, start_range):
    """
    The function implemented to calculate an optimal ln function parameters. The logarithmic function is given below. Returns list of tuples with theparameters.
    Uses tuple of two 1-D lists of calibration standards data of the same size as an input data.
    """
    X_st1 = np.array(D_st[0])
    Y_st1 = np.array(D_st[1])
    
    def ln_func_opt(X_param):
        Y_th1=X_param[0]*np.log(X_st1)+X_param[1]
        MNK=sum((Y_th1-Y_st1)**2)
        return(MNK)
    
    A_l = start_range['A'][0]
    A_p = start_range['A'][1]
    B_l = start_range['B'][0]
    B_p = start_range['B'][1]
    bounds=[(A_l, A_p),(B_l, B_p)]

    result = differential_evolution(ln_func_opt, bounds, maxiter=750000, popsize= 45, strategy='best2bin')    
   
    par=[result.x[0], result.x[1]]
    RSS=ln_func_opt(X_param=par)/len(X_st1)
    if(result.x[0] == A_l or result.x[0] == A_p or result.x[1] == B_l or result.x[1] == B_p or RSS > 0.3):
        """ at least one of the parameters is equal to the range limit, error!"""
        error_bound = 1
    else:
        """ the parameters in proper range, ok"""
        error_bound = 0
      
    parameters=[('A_par1',result.x[0]), ('B_par1',result.x[1]), ('RSS', RSS)]

    R_squ=R_squared_LN(D_st, parameters)
    if(R_squ < 0.6):  error_bound = 1
    y_theor = LN_func(X_data=X_st1, param=parameters)
    akaike_crit = aic.aic(y=Y_st1, y_pred=y_theor, p=2)
    bayes_crit = bic.bic(y=Y_st1, y_pred=y_theor, p=2)
    r, prob = pearsonr(Y_st1, y_theor)

    parameters.append(('R_squ', R_squ))
    parameters.append(('AIC_crit', akaike_crit))
    parameters.append(('BIC_crit', bayes_crit))
    parameters.append(('R_corre', r))    
    parameters.append(('Error', error_bound))

    return(parameters)


def ln_func_concentration(Y_data, param):
    """
    The function implemented to calculate titer concentration values. Returns list and uses list of absorbance values and parameters obtained from the optimization procedure from ln_func_par_opt.
    """
    #print("Y from ln_func_conc: ",Y_data)
    A_par = param[0][1]
    B_par = param[1][1]

    X = np.exp((Y_data - B_par)/A_par)

    return(X)

def LN_func(X_data, param):
    """
    Returns the values of absorbance using titer concentration and model parameters
    Function input:
    X_data - list of titer concentrations
    param - list of tuples of the optimal parameters obtained for the logarithmic model
    """
    A_par = param[0][1]
    B_par = param[1][1]

    Y_func = A_par*np.log(X_data) + B_par
    """returns list of calculated absorbance values"""
    return(Y_func)

def remove_abs_smaller_0(X_data, Y_data):
    """
    Cleans the ln data from the values <0.
    """
    dat = pd.DataFrame()
    dat['X'] = X_data
    dat['Y'] = Y_data
    dat = dat[dat.Y > 0]
    """returns DataFrame"""
    return dat
    
def R_squared_LN(D_st, par):
    """Returns the R^2 model parameter from the calibration standards values and ln model parameters
    Function input:
    D_st - a tuple of calibration standard values list (D_st[0]- titer concentration, D_st[1] - absorbance)
    """
    X_dat=D_st[0]
    Ytheoretical=LN_func(X_data= X_dat, param=par)
    resp=r2_score(y_true= D_st[1], y_pred= Ytheoretical)
    """returns R^2 value for the fitted model and previously chosen data set"""
    return resp
    
#logit 5PL function stuff
def logit_5PL_par_opt(D_st, start_range):
    """
    The function implemented to calculate an optimal 5PL function parameters. The 5PL function is given below. Returns list of tuples with the parameters.
    Uses tuple or a list of two 1-D lists of calibration standards data of the same size as an input data.
    """
    X_st = np.array(D_st[0])
    Y_st = np.array(D_st[1])
    
    def logit_func_opt(X_param):
        Y_th = X_param[0] + ((X_param[1] - X_param[0])/(1 + (X_st/X_param[3])**X_param[2])**X_param[4])
        MNK = sum((Y_th-Y_st)**2)
        return MNK 
   
    D_l = start_range['D'][0]
    D_p = start_range['D'][1]
    A_l = start_range['A'][0]
    A_p = start_range['A'][1]
    B_l = start_range['B'][0]
    B_p = start_range['B'][1]
    C_l = start_range['C'][0]
    C_p = start_range['C'][1]
    E_l = start_range['E'][0]
    E_p = start_range['E'][1]
    bounds=[(D_l, D_p),(A_l, A_p),(B_l, B_p),(C_l, C_p),(E_l, E_p)]

       
    result = differential_evolution(logit_func_opt, bounds, maxiter=750000, popsize= 45, strategy='best2bin')

    par=[result.x[0], result.x[1], result.x[2], result.x[3], result.x[4]]  
    RSS=logit_func_opt(X_param=par)/len(X_st)
    if(result.x[1] == A_l or result.x[1] == A_p or result.x[2] == B_l or result.x[2] == B_p or
       result.x[3] == C_l or result.x[3] == C_p or result.x[0] == D_l or result.x[0] == D_p or
       result.x[4] == E_l or result.x[4] == E_p or RSS > 0.3):
        """ at least one of the parameters is equal to the range limit, error!"""
        error_bound = 1
    else:
        """ the parameters in proper range, ok"""
        error_bound = 0
    
    parameters=[('D_par1',result.x[0]), ('A_par1',result.x[1]), ('B_par1', result.x[2]), ('C_par1', result.x[3]), 
                ('E_par1', result.x[4]), ('RSS', RSS)]

    R_squ=R_squared_5PL(D_st, parameters)
    if(R_squ < 0.6):  error_bound = 1   
    y_theor = logit_5PL_func(X_data=X_st, param=parameters)
    akaike_crit = aic.aic(y=Y_st, y_pred=y_theor, p=5)
    bayes_crit = bic.bic(y=Y_st, y_pred=y_theor, p=5)
    r, prob = pearsonr(Y_st, y_theor)

    parameters.append(('R_squ', R_squ))
    parameters.append(('AIC_crit', akaike_crit))
    parameters.append(('BIC_crit', bayes_crit))
    parameters.append(('R_corre', r))
    parameters.append(('Error', error_bound))   
    """returns the list of tuples with 5PL model parameters"""
    return (parameters)


def logit_5PL_concentration(Y_data, param):
    """
    The function was implemented to calculate titer concentration from absorbance and 5PL model parameters.
    Function inputs:
    Y_data - values of absorbance
    param - list of tuples with optimal 5PL-model parameters
    """
    D_par = param[0][1]
    A_par = param[1][1]
    B_par = param[2][1]
    C_par = param[3][1]
    E_par = param[4][1]

    X = C_par*(((((A_par - D_par)/(Y_data - D_par))**(1/E_par))-1)**(1/B_par))
        
    """returns list of titer concentration values"""
    return(X)


def logit_5PL_func(X_data, param):
    """
    Returns the values of absorbance using titer concentration and 5PL model parameters
    Function input:
    X_data - list of titer concentrations
    param - list of tuples of the optimal parameters obtained for the 5PL model
    """
    
    D_par = param[0][1]
    A_par = param[1][1]
    B_par = param[2][1]
    C_par = param[3][1]
    E_par = param[4][1]

    Y_func = D_par + ((A_par - D_par)/(1 + (X_data/C_par)**B_par)**E_par)
    """returns calculated absorbance values"""
    return(Y_func)


def R_squared_5PL(D_st, par):
    """
    The function returns the R^2 model parameter from the calibration standards values and 5PL model parameters.
    Function input:
    D_st - a tuple of calibration standard values list (D_st[0]- titer concentration, D_st[1] - absorbance)
    """
    X_dat=D_st[0]
    Ytheoretical=logit_5PL_func(X_data= X_dat, param=par)
    resp=r2_score(y_true= D_st[1], y_pred= Ytheoretical)
    """returns R^2 value for the fitted model and previously chosen data set"""
    return(resp)
    

#logit 4PL function stuff
def logit_4PL_par_opt(D_st, start_range):
    """
    The function implemented to calculate an optimal 4PL function parameters. The 5PL function is given below. Returns a list of tuples with the 4PL model parameters.
    Function input:
    Uses tuple or a list of two 1-D lists of calibration standards data of the same size as an input data.
    """
    X_st = D_st[0]
    Y_st = D_st[1]
    
    def logit_func_opt(X_param):
        Y_th = X_param[0] + (X_param[1] - X_param[0])/(1 + (X_st/X_param[3])**X_param[2])
        MNK = sum((Y_th-Y_st)**2)
        return MNK

    D_l = start_range['D'][0]
    D_p = start_range['D'][1]
    A_l = start_range['A'][0]
    A_p = start_range['A'][1]
    B_l = start_range['B'][0]
    B_p = start_range['B'][1]
    C_l = start_range['C'][0]
    C_p = start_range['C'][1]

    bounds=[(D_l, D_p),(A_l, A_p),(B_l, B_p),(C_l, C_p)]
    result = differential_evolution(logit_func_opt, bounds, maxiter=750000, popsize= 45, strategy='best2bin')

    par=[result.x[0], result.x[1], result.x[2], result.x[3]]
    RSS=logit_func_opt(X_param=par)/len(X_st)  
    if(result.x[1] == A_l or result.x[1] == A_p or result.x[2] == B_l or result.x[2] == B_p or
       result.x[3] == C_l or result.x[3] == C_p or result.x[0] == D_l or result.x[0] == D_p or RSS > 0.3):
        """ at least one of the parameters is equal to the range limit, error!"""
        error_bound = 1
    else:
        """ the parameters in proper range, ok"""
        error_bound = 0
    
    parameters=[('D_par1',result.x[0]), ('A_par1',result.x[1]), ('B_par1', result.x[2]), 
                ('C_par1', result.x[3]), ('RSS', RSS)]

    R_squ=R_squared_4PL(D_st, parameters)
    if(R_squ < 0.6):  error_bound = 1
    y_theor = logit_4PL_func(X_data=X_st, param=parameters)
    akaike_crit = aic.aic(y=Y_st, y_pred=y_theor, p=4)
    bayes_crit = bic.bic(y=Y_st, y_pred=y_theor, p=4)
    r, prob = pearsonr(Y_st, y_theor)

    parameters.append(('R_squ', R_squ))
    parameters.append(('AIC_crit', akaike_crit))
    parameters.append(('BIC_crit', bayes_crit))
    parameters.append(('R_corre', r))
    parameters.append(('Error', error_bound))
    """returns list of tuples with optimal 4PL-model parameters and R_squ"""
    return (parameters)


def logit_4PL_concentration(Y_data, param):
    """
    The function was implemented to calculate titer concentration from absorbance and 4PL model parameters (reverse 4PL-function).
    Function inputs:
    Y_data - values of absorbance
    param - list of tuples with optimal 4PL-model parameters
    """
    
    D_par = param[0][1]
    A_par = param[1][1]
    B_par = param[2][1]
    C_par = param[3][1]

    X = C_par*((((A_par - D_par)/(Y_data - D_par)) - 1)**(1/B_par))        
        
    """returns list of titer concentration values"""
    return(X)


def logit_4PL_func(X_data, param):
    """
    The function implemented in order to return the values of absorbance using titer concentration and 4PL model parameters.
    Function input:
    X_data - list of titer concentrations
    param - list of tuples of the optimal parameters obtained for the 4PL model
    """
    D_par = param[0][1]
    A_par = param[1][1]
    B_par = param[2][1]
    C_par = param[3][1]

    Y_func = D_par + (A_par - D_par)/(1 + (X_data/C_par)**B_par)
    """returns calculated absorbance values"""
    return(Y_func)


def R_squared_4PL(D_st, par):
    """
    Returns the R^2 model parameter from the calibration standards values and ln model parameters
    Function input:
    D_st - a tuple of calibration standard values list (D_st[0]- titer concentration, D_st[1] - absorbance)
    """
    X_dat = D_st[0]
    Ytheoretical = logit_4PL_func(X_data= X_dat, param=par)
    resp = r2_score(y_true= D_st[1], y_pred= Ytheoretical)
    """returns R^2 value for the fitted model and previously chosen data set"""
    return(resp)

def data_std_format(imp_data, imp_data_map, data_standards):
    """
    The function implemented to extract, clean and order calibration standards for further ln-model fitting.
    """
    g1=imp_data_map.stack()
    g2=set(g1)
    for g in g2:
        vars()[g]=[]
        
    imp_columns=list(imp_data.columns.values)
    imp_index=list(imp_data.index.values)
    
    for c_name in imp_columns:
        for r_name in imp_index:
            vars()[imp_data_map.loc[r_name][c_name]].append(imp_data.loc[r_name][c_name])
            
    X_st = []
    Y_st = []
    
    for g in g2:
        vars()[g]=sum(vars()[g]) / float(len(vars()[g]))
        if ('std' in g):
            X_st.append(data_standards.loc[g]['values'])
            Y_st.append(vars()[g])
            continue
        else:
            continue
    """ 
    Remove zero value from concentration - unphysical and denger for logharitm
    """        
    X_st1=copy.deepcopy(X_st)
    Y_st1=copy.deepcopy(Y_st)
    N=0
    lenX=len(X_st1)
    N_of_del=0
    while N< lenX:
        if X_st1[N-N_of_del]==0:
            del(X_st1[N-N_of_del])
            del(Y_st1[N-N_of_del])
            N_of_del=N_of_del+1
            N=N+1
            continue
        else:
            N=N+1
            continue
    """returns the tuple of two arrays"""    
    return(X_st1, Y_st1)

    
def samples_concenration(imp_data, imp_data_map, params, ordered, d_st):
    
    def func_concentration(Y, param):
        """ Fill std concentration by input data, but what for NC an PC?"""
        return 0
    
    """Refactoring parsing_sam_and_concenration """
    df_local = pd.DataFrame(columns=['Y_names', 'Y_abs', 'samples_nb', 'Y_ave_abs', 'Y_abs_std', 'X_concentration',\
                                     'X_s_conc_std_bottom', 'X_s_conc_std_upper', 'Errors'])
    df_local.set_index('Y_names', inplace=True)  # setup Y_s_names as index of DF
    for i in range(0, len(imp_data_map.index)):
        for j in range(0, len(imp_data_map.columns)):
            is_string = imp_data_map.values[i,j]
            is_concen = imp_data.values[i,j]
            if(df_local.index.isin([is_string]).any()):
                """If YES then add to this position"""
                tmp_help = df_local.Y_abs[is_string]
                tmp_help.append(is_concen) #? dlaczego tak?
                df_local.ix[is_string, 'samples_nb'] +=  1
            else:
                """If NO then append to DataFrame new record"""
                tmp_df = pd.DataFrame([[[is_concen], 1]], columns=[ 'Y_abs', 'samples_nb'], index=[is_string])
                df_local = df_local.append(tmp_df, ignore_index=False)

    """
    First time for std, to prepare value for LN limit condition
    """
    for index, row in df_local.iterrows():
        row["Y_ave_abs"] = sum(row["Y_abs"])/row["samples_nb"]
        row["Y_abs_std"] = np.std(row["Y_abs"])
        if ('std' in index):  # put standard concentration (comming from parsing part) to the main DF
            row["X_concentration"] = d_st.conc[index]
            row["Errors"] = 100
    
    Max_ind_of_std = df_local['Y_ave_abs'].astype('float64').idxmax(skipna=True)  
    Max_val_of_std = df_local.loc[Max_ind_of_std]['Y_ave_abs']
    Max_val_of_dev = df_local.loc[Max_ind_of_std]['Y_abs_std']
    """
    Second time for sam only
    """       
    for index, row in df_local.iterrows():
        if ('sam' in index):           
            if ((ordered=="LN") or (ordered=="cfLN")):
                if(((row["Y_ave_abs"] + row["Y_abs_std"]) >= 0.0) and (row["Y_ave_abs"] - row["Y_abs_std"]) <= (Max_val_of_dev + Max_val_of_std)):
                    ave_concentr = ln_func_concentration(Y_data=row["Y_ave_abs"], param=params)
                    row["X_concentration"] = ave_concentr
                    row["X_s_conc_std_upper"] = ln_func_concentration(Y_data=(row["Y_ave_abs"] + row["Y_abs_std"]), param=params) - ave_concentr
                    row["X_s_conc_std_bottom"] = ave_concentr - ln_func_concentration(Y_data=(row["Y_ave_abs"] - row["Y_abs_std"]), param=params)
                    row["Errors"] = 0
                else:
                    row["X_concentration"] = 0.0
                    row["Errors"] = 1
            if ((ordered=="4PL") or (ordered=="cf4PL")): 
                if((row["Y_ave_abs"]  + row["Y_abs_std"] >= params[1][1]) and (row["Y_ave_abs"] - row["Y_abs_std"] <= params[0][1])):
                    """ 1) bottom limit for 4PL (Y_abs+Y_std >= A_par)
                        2) upper limit for 4PL (Y_abs-Y_std <= D_par)
                    """
                    ave_concentr = logit_4PL_concentration(Y_data=row["Y_ave_abs"], param=params)
                    row["X_concentration"] = ave_concentr
                    row["X_s_conc_std_upper"] = logit_4PL_concentration(Y_data=(row["Y_ave_abs"] + row["Y_abs_std"]), param=params) - ave_concentr
                    row["X_s_conc_std_bottom"] = ave_concentr - logit_4PL_concentration(Y_data=(row["Y_ave_abs"] - row["Y_abs_std"]), param=params)
                    row["Errors"] = 0
                else:
                    row["X_concentration"] = 0.0
                    row["Errors"] = 1
            if ((ordered=="5PL") or (ordered=="cf5PL")):            
                if((row["Y_ave_abs"] + row["Y_abs_std"] >= params[1][1]) and (row["Y_ave_abs"] - row["Y_abs_std"] <= params[0][1])): 
                    """ 1) bottom limit for 5PL (Y_abs+Y_std >= A_par)
                        2) upper limit for 5PL (Y_abs-Y_std <= D_par)
                    """
                    ave_concentr = logit_5PL_concentration(Y_data=row["Y_ave_abs"], param=params)
                    row["X_concentration"] = ave_concentr
                    row["X_s_conc_std_upper"] = logit_5PL_concentration(Y_data=(row["Y_ave_abs"] + row["Y_abs_std"]), param=params) - ave_concentr
                    row["X_s_conc_std_bottom"] = ave_concentr - logit_5PL_concentration(Y_data=(row["Y_ave_abs"] - row["Y_abs_std"]), param=params)
                    row["Errors"] = 0
                else:
                    row["X_concentration"] = 0.0
                    row["Errors"] = 1
               
        """                    
        if ('zero' in index): 
            df_local["X_concentration"] = func_concentration(Y=df_local['Y_ave_abs'], param=params)
        if ('nc' in index):  
            df_local["X_concentration"] = func_concentration(Y=df_local['Y_ave_abs'], param=params)
        if ('pc' in index):
            df_local["X_concentration"] = func_concentration(Y=df_local['Y_ave_abs'], param=params)
        """   
    """returns sample names, titer concentration, lower and upper boundaries of concentration standard deviations and simmetrical value of absorbance standard deviation."""

    print(df_local)

    """returns:
       Y_s_names           - sample names, 
       X_s_concentration   - titer concentration, 
       Y_s_good            - absorpcja probki zgodna z modelem aproksymujacym (np. Y_s_abs-Y_s_std => A_par and Y_s_abs+Y_s_std>D_par)
       Y_s_bad             - absorpcja probki wypadajaca poza model aproksymujacy (np. Y_s_abs-Y_s_std < A_par) string!
       X_s_conc_std_bottom - lower  boundaries of concentration standard deviations 
       X_s_conc_std_upper  - upper boundaries of concentration standard deviations  
       Y_s_abs_std         - simmetrical value of absorbance standard deviation
       X_std_abs           - concentrations of standards 
       Y_std_abs           - ave. absorption of standards 
       Y_std_std           - std of absorption of standards.
    """
    Y_s_names = []
    X_s_concentration = []
    Y_s_good = []
    Y_s_bad = []
    Y_s_abs_std = []
    X_s_conc_std_bottom = []
    X_s_conc_std_upper = []
    Y_std_std = []
    Y_std_abs = []
    X_std_abs = []
    for index, row in df_local.iterrows():
        if (('sam' in index) and (row['Errors']==0)):
            Y_s_names.append(index)
            X_s_concentration.append(row["X_concentration"])
            Y_s_good.append(row["Y_ave_abs"])
            Y_s_abs_std.append(row["Y_abs_std"])
            X_s_conc_std_upper.append(row["X_s_conc_std_upper"])
            X_s_conc_std_bottom.append(row["X_s_conc_std_bottom"])
        if (('sam' in index) and (row['Errors']==1)): # develop for more error levels
            Y_s_bad.append(index)
        if (('std' in index) and (row['Errors']==100)):
            Y_std_std.append(row["Y_abs_std"])
            Y_std_abs.append(row["Y_ave_abs"])
            X_std_abs.append(row["X_concentration"])
          
    #          0              1               2        3          4                   5                  6             7          8          9
    return (Y_s_names, X_s_concentration, Y_s_good, Y_s_bad, X_s_conc_std_bottom, X_s_conc_std_upper, Y_s_abs_std, X_std_abs, Y_std_abs, Y_std_std)

        
def samples_concenration_new(imp_data, imp_data_map, params, ordered, d_st, df_local):
    
    def func_concentration(Y, param):
        """ Fill std concentration by input data, but what for NC an PC?"""
        return 0
    
    """Refactoring parsing_sam_and_concenration """
    for index, row in df_local.iterrows():
        if ('sam' in index):           
            if ((ordered=="LN") or (ordered=="cfLN")):
                if(row["Y_ave_abs"] >= params[1][1]): #temporary solution
                    ave_concentr = ln_func_concentration(Y_data=row["Y_ave_abs"], param=params)
                    row["X_concentration"] = ave_concentr
                    row["X_s_conc_std_upper"] = ln_func_concentration(Y_data=(row["Y_ave_abs"] + row["Y_abs_std"]), param=params) - ave_concentr
                    row["X_s_conc_std_bottom"] = ave_concentr - ln_func_concentration(Y_data=(row["Y_ave_abs"] - row["Y_abs_std"]), param=params)
                    row["Errors"] = 0
                else:
                    row["X_concentration"] = 0.0
                    row["Errors"] = 1
            if ((ordered=="4PL") or (ordered=="cf4PL")): 
                if((row["Y_ave_abs"]  + row["Y_abs_std"] >= params[1][1]) and (row["Y_ave_abs"] - row["Y_abs_std"] <= params[0][1])):
                    """ 1) bottom limit for 4PL (Y_abs+Y_std >= A_par)
                        2) upper limit for 4PL (Y_abs-Y_std <= D_par)
                    """
                    ave_concentr = logit_4PL_concentration(Y_data=row["Y_ave_abs"], param=params)
                    row["X_concentration"] = ave_concentr
                    row["X_s_conc_std_upper"] = logit_4PL_concentration(Y_data=(row["Y_ave_abs"] + row["Y_abs_std"]), param=params) - ave_concentr
                    row["X_s_conc_std_bottom"] = ave_concentr - logit_4PL_concentration(Y_data=(row["Y_ave_abs"] - row["Y_abs_std"]), param=params)
                    row["Errors"] = 0
                else:
                    row["X_concentration"] = 0.0
                    row["Errors"] = 1
            if ((ordered=="5PL") or (ordered=="cf5PL")):            
                if((row["Y_ave_abs"] + row["Y_abs_std"] >= params[1][1]) and (row["Y_ave_abs"] - row["Y_abs_std"] <= params[0][1])): 
                    """ 1) bottom limit for 5PL (Y_abs+Y_std >= A_par)
                        2) upper limit for 5PL (Y_abs-Y_std <= D_par)
                    """
                    ave_concentr = logit_5PL_concentration(Y_data=row["Y_ave_abs"], param=params)
                    row["X_concentration"] = ave_concentr
                    row["X_s_conc_std_upper"] = logit_5PL_concentration(Y_data=(row["Y_ave_abs"] + row["Y_abs_std"]), param=params) - ave_concentr
                    row["X_s_conc_std_bottom"] = ave_concentr - logit_5PL_concentration(Y_data=(row["Y_ave_abs"] - row["Y_abs_std"]), param=params)
                    row["Errors"] = 0
                else:
                    row["X_concentration"] = 0.0
                    row["Errors"] = 1
               
        """                    
        if ('zero' in index): 
            df_local["X_concentration"] = func_concentration(Y=df_local['Y_ave_abs'], param=params)
        if ('nc' in index):  
            df_local["X_concentration"] = func_concentration(Y=df_local['Y_ave_abs'], param=params)
        if ('pc' in index):
            df_local["X_concentration"] = func_concentration(Y=df_local['Y_ave_abs'], param=params)
        """   
    """returns sample names, titer concentration, lower and upper boundaries of concentration standard deviations and simmetrical value of absorbance standard deviation."""

    print(df_local)

    """returns:
       Y_s_names           - sample names, 
       X_s_concentration   - titer concentration, 
       Y_s_good            - absorpcja probki zgodna z modelem aproksymujacym (np. Y_s_abs-Y_s_std => A_par and Y_s_abs+Y_s_std>D_par)
       Y_s_bad             - absorpcja probki wypadajaca poza model aproksymujacy (np. Y_s_abs-Y_s_std < A_par) string!
       X_s_conc_std_bottom - lower  boundaries of concentration standard deviations 
       X_s_conc_std_upper  - upper boundaries of concentration standard deviations  
       Y_s_abs_std         - simmetrical value of absorbance standard deviation
       X_std_abs           - concentrations of standards 
       Y_std_abs           - ave. absorption of standards 
       Y_std_std           - std of absorption of standards.
    """
    Y_s_names = []
    X_s_concentration = []
    Y_s_good = []
    Y_s_bad = []
    Y_s_abs_std = []
    X_s_conc_std_bottom = []
    X_s_conc_std_upper = []
    Y_std_std = []
    Y_std_abs = []
    X_std_abs = []
    for index, row in df_local.iterrows():
        if (('sam' in index) and (row['Errors']==0)):
            Y_s_names.append(index)
            X_s_concentration.append(row["X_concentration"])
            Y_s_good.append(row["Y_ave_abs"])
            Y_s_abs_std.append(row["Y_abs_std"])
            X_s_conc_std_upper.append(row["X_s_conc_std_upper"])
            X_s_conc_std_bottom.append(row["X_s_conc_std_bottom"])
        if (('sam' in index) and (row['Errors']==1)): # develop for more error levels
            Y_s_bad.append(index)
        if (('std' in index) and (row['Errors']==100)):
            Y_std_std.append(row["Y_abs_std"])
            Y_std_abs.append(row["Y_ave_abs"])
            X_std_abs.append(row["X_concentration"])
          
    #          0              1               2        3          4                   5                  6             7          8          9
    return (Y_s_names, X_s_concentration, Y_s_good, Y_s_bad, X_s_conc_std_bottom, X_s_conc_std_upper, Y_s_abs_std, X_std_abs, Y_std_abs, Y_std_std)
#    return df_local
        
 
def recognition_experiment(imp_data, imp_data_map, d_st):
    
    def func_concentration(Y, param):
        """ Fill std concentration by input data, but what for NC an PC?"""
        return 0
    
    """Refactoring parsing_sam_and_concenration """
    df_local = pd.DataFrame(columns=['Y_names', 'Y_abs', 'samples_nb', 'Y_ave_abs', 'Y_abs_std', 'X_concentration',\
                                     'X_s_conc_std_bottom', 'X_s_conc_std_upper', 'Errors'])
    df_local.set_index('Y_names', inplace=True)  # setup Y_s_names as index of DF
    for i in range(0, len(imp_data_map.index)):
        for j in range(0, len(imp_data_map.columns)):
            is_string = imp_data_map.values[i,j]
            is_concen = imp_data.values[i,j]
            if(df_local.index.isin([is_string]).any()):
                """If YES then add to this position"""
                tmp_help = df_local.Y_abs[is_string]
                tmp_help.append(is_concen) #? dlaczego tak?
                df_local.ix[is_string, 'samples_nb'] +=  1
            else:
                """If NO then append to DataFrame new record"""
                tmp_df = pd.DataFrame([[[is_concen], 1]], columns=[ 'Y_abs', 'samples_nb'], index=[is_string])
                df_local = df_local.append(tmp_df, ignore_index=False)
    flag = 'STD'
    for index, row in df_local.iterrows():
        row["Y_ave_abs"] = sum(row["Y_abs"])/row["samples_nb"]
        row["Y_abs_std"] = np.std(row["Y_abs"])
        #if ('sam' in index):                       
        if ('std' in index):  # put standard concentration (comming from parsing part) to the main DF
            row["X_concentration"] = d_st.conc[index]
            row["Errors"] = 100
        if (('nc' in index) and ('pc' in index)):  
            flag = 'SPC'
               
        """                    
        if ('zero' in index): 
            df_local["X_concentration"] = func_concentration(Y=df_local['Y_ave_abs'], param=params)
        if ('nc' in index):  
            df_local["X_concentration"] = func_concentration(Y=df_local['Y_ave_abs'], param=params)
        if ('pc' in index):
            df_local["X_concentration"] = func_concentration(Y=df_local['Y_ave_abs'], param=params)
        """   
    """returns sample names, titer concentration, lower and upper boundaries of concentration standard deviations and simmetrical value of absorbance standard deviation."""

    print(df_local)

    """returns:
       Y_s_names           - sample names, 
       X_s_concentration   - titer concentration, 
       Y_s_good            - absorpcja probki zgodna z modelem aproksymujacym (np. Y_s_abs-Y_s_std => A_par and Y_s_abs+Y_s_std>D_par)
       Y_s_bad             - absorpcja probki wypadajaca poza model aproksymujacy (np. Y_s_abs-Y_s_std < A_par) string!
       X_s_conc_std_bottom - lower  boundaries of concentration standard deviations 
       X_s_conc_std_upper  - upper boundaries of concentration standard deviations  
       Y_s_abs_std         - simmetrical value of absorbance standard deviation
       X_std_abs           - concentrations of standards 
       Y_std_abs           - ave. absorption of standards 
       Y_std_std           - std of absorption of standards.

    Y_s_names = []
    X_s_concentration = []
    Y_s_good = []
    Y_s_bad = []
    Y_s_abs_std = []
    X_s_conc_std_bottom = []
    X_s_conc_std_upper = []
    Y_std_std = []
    Y_std_abs = []
    X_std_abs = []
    for index, row in df_local.iterrows():
        if (('sam' in index) and (row['Errors']==0)):
            Y_s_names.append(index)
            X_s_concentration.append(row["X_concentration"])
            Y_s_good.append(row["Y_ave_abs"])
            Y_s_abs_std.append(row["Y_abs_std"])
            X_s_conc_std_upper.append(row["X_s_conc_std_upper"])
            X_s_conc_std_bottom.append(row["X_s_conc_std_bottom"])
        if (('sam' in index) and (row['Errors']==1)): # develop for more error levels
            Y_s_bad.append(index)
        if (('std' in index) and (row['Errors']==100)):
            Y_std_std.append(row["Y_abs_std"])
            Y_std_abs.append(row["Y_ave_abs"])
            X_std_abs.append(row["X_concentration"])
    """
          
    #          0              1               2        3          4                   5                  6             7          8          9
    #return (Y_s_names, X_s_concentration, Y_s_good, Y_s_bad, X_s_conc_std_bottom, X_s_conc_std_upper, Y_s_abs_std, X_std_abs, Y_std_abs, Y_std_std)
    return df_local, flag

       
def LN_func_curve_fit(D_st, start_range):
    """
    Optimization on the base of curvr_fit function
    """
    X_st = D_st[0]
    Y_st = D_st[1]
    
    def logit_func_LN(X_st, A, B):
        Y_th = A*np.log(X_st) + B
        return Y_th
  
    A_l = start_range['A'][0]
    A_p = start_range['A'][1]
    B_l = start_range['B'][0]
    B_p = start_range['B'][1]
    
    par, pcov = curve_fit(logit_func_LN, X_st, Y_st, bounds=([A_l, B_l],[A_p, B_p]))

    def avarage_LN(X_st, par):
        Y_new=logit_func_LN(X_st=X_st, A=par[0], B=par[1])
        return  sum((Y_new-Y_st)**2)
    
    RSS=avarage_LN(X_st, par)/len(X_st)
    if(par[0] == A_l or par[0] == A_p or par[1] == B_l or par[1] == B_p or RSS > 0.3):
        """ at least one of the parameters is equal to the range limit, error!"""
        error_bound = 1
    else:
        """ the parameters in proper range, ok"""
        error_bound = 0
    
    parameters=[('A_par1',par[0]), ('B_par1', par[1]), ('RSS', RSS)]
    
    R_squ=R_squared_LN(D_st, parameters)
    if(R_squ < 0.6):  error_bound = 1
    y_theor = LN_func(X_data=X_st, param=parameters)
    akaike_crit = aic.aic(y=Y_st, y_pred=y_theor, p=2)
    bayes_crit = bic.bic(y=Y_st, y_pred=y_theor, p=2)
    r, prob = pearsonr(Y_st, y_theor)

    parameters.append(('R_squ', R_squ))
    parameters.append(('AIC_crit', akaike_crit))
    parameters.append(('BIC_crit', bayes_crit))
    parameters.append(('R_corre', r))  
    parameters.append(('Error', error_bound))
    """returns list of tuples with optimal 4PL-model parameters and R_squ"""
    return (parameters)


def logit_4PL_curve_fit(D_st, start_range):
    """
    Optimization on the base of curvr_fit function
    """
    X_st = D_st[0]
    Y_st = D_st[1]
    
    def logit_func_4PL(X_st, A, B, C, D):
        Y_th = D + (A - D)/(1 + (X_st/C)**B)        
        return Y_th
  
    D_l = start_range['D'][0]
    D_p = start_range['D'][1]
    A_l = start_range['A'][0]
    A_p = start_range['A'][1]
    B_l = start_range['B'][0]
    B_p = start_range['B'][1]
    C_l = start_range['C'][0]
    C_p = start_range['C'][1]
    
    par, pcov = curve_fit(logit_func_4PL, X_st, Y_st, bounds=([D_l, A_l, B_l, C_l], [D_p, A_p, B_p, C_p]))

    def avarage_4PL(X_st, par):
        Y_new=logit_func_4PL(X_st=X_st, A=par[1], B=par[2], C=par[3], D=par[0])
        return  sum((Y_new-Y_st)**2)
    
    RSS=avarage_4PL(X_st, par)/len(X_st)
    if(par[1] == A_l or par[1] == A_p or par[2] == B_l or par[2] == B_p or
       par[3] == C_l or par[3] == C_p or par[0] == D_l or par[0] == D_p or RSS > 0.3):
        """ at least one of the parameters is equal to the range limit, error!"""
        error_bound = 1
    else:
        """ the parameters in proper range, ok"""
        error_bound = 0
    
    parameters=[('D_par1',par[0]), ('A_par1',par[1]), ('B_par1', par[2]), 
                ('C_par1', par[3]), ('RSS', RSS)]

    R_squ = R_squared_4PL(D_st, parameters)
    if(R_squ < 0.6):  error_bound = 1
    y_theor = logit_4PL_func(X_data=X_st, param=parameters)
    akaike_crit = aic.aic(y=Y_st, y_pred=y_theor, p=4)
    bayes_crit = bic.bic(y=Y_st, y_pred=y_theor, p=4)
    r, prob = pearsonr(Y_st, y_theor)

    parameters.append(('R_squ', R_squ))  
    parameters.append(('AIC_crit', akaike_crit))
    parameters.append(('BIC_crit', bayes_crit))
    parameters.append(('R_corre', r))
    parameters.append(('Error', error_bound))   
    """returns list of tuples with optimal 4PL-model parameters and R_squ"""
    return (parameters)

def logit_5PL_curve_fit(D_st, start_range):
    """
    Optimization on the base of curvr_fit function
    """
    X_st = D_st[0]
    Y_st = D_st[1]
    def logit_func_5PL(X_st, A, B, C, D, E):
        Y_th=D + (A - D)/((1 + (X_st/C)**B)**E)
        return Y_th
  
    D_l = start_range['D'][0]
    D_p = start_range['D'][1]
    A_l = start_range['A'][0]
    A_p = start_range['A'][1]
    B_l = start_range['B'][0]
    B_p = start_range['B'][1]
    C_l = start_range['C'][0]
    C_p = start_range['C'][1]
    E_l = start_range['E'][0]
    E_p = start_range['E'][1]

    limits=([D_l, A_l, B_l, C_l, E_l], [D_p, A_p, B_p, C_p, E_p])
    par, pcov = curve_fit(logit_func_5PL, X_st, Y_st, bounds=limits)        

    def avarage_5PL(X_st, par):
        Y_new = logit_func_5PL(X_st=X_st, A=par[1], B=par[2], C=par[3], D=par[0], E=par[4])
        return  sum((Y_new-Y_st)**2)

    RSS=avarage_5PL(X_st, par)/len(X_st)
    if(par[1] == A_l or par[1] == A_p or par[2] == B_l or par[2] == B_p or
       par[3] == C_l or par[3] == C_p or par[0] == D_l or par[0] == D_p or 
       par[4] == E_l or par[4] == E_p or RSS > 0.3):
        """ at least one of the parameters is equal to the range limit, error!"""
        error_bound = 1
    else:
        """ the parameters in proper range, ok"""
        error_bound = 0

    parameters = [('D_par1',par[0]), ('A_par1',par[1]), ('B_par1', par[2]), ('C_par1', par[3]), 
                  ('E_par1', par[4]), ('RSS', RSS)]
    
    R_squ = R_squared_5PL(D_st, parameters)    
    if(R_squ < 0.6):  error_bound = 1
    y_theor = logit_5PL_func(X_data=X_st, param=parameters)
    akaike_crit = aic.aic(y=Y_st, y_pred=y_theor, p=5)
    bayes_crit = bic.bic(y=Y_st, y_pred=y_theor, p=5)
    r, prob = pearsonr(Y_st, y_theor)

    parameters.append(('R_squ', R_squ))
    parameters.append(('AIC_crit', akaike_crit))
    parameters.append(('BIC_crit', bayes_crit))
    parameters.append(('R_corre', r))
    parameters.append(('Error', error_bound))   
    """returns list of tuples with optimal 5PL-model parameters and R_squ"""
    return (parameters)