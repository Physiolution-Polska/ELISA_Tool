#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 2019
@author:  Marek Bawiec, Grzegorz Banach
The ELISA tool plugin: module for validation functions of approximation models.
"""

__author__ = "Marek Bawiec, Grzegorz Banach"
__copyright__ = "Copyright 2019, Physiolution Polska"
__credits__ = [""]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Grzegorz Banach"
__email__ = "g.banach@physiolution.pl"
__status__ = "Testing"

import elisa_tool_repo.et_calc as ecf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from RegscorePy import *


def plot_init():
    '''
    Prepear all for plots
    #####################
    input: none
    output: handle for figure
    '''
    paperheight = 8.4
    paperwidth = 16.8
    margin = 1.0
    unitName = 'a.u.'
    chartName = "Absorption"
    fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
    plt.grid(True)
    plt.xlabel('C [a.u.]', fontsize='x-small')
    plt.ylabel('%s [%s]' % (chartName, unitName), fontsize='x-small')
    return (fig)


def validation_LN(noise_par, point_par):    
    param_val=[('A_par1', 0.55), ('B_par1', 1.22), ('MSE', 0.0), ('R_squ', 0.0)]
    config={"A":[0.0,20.0],	"B":[0.0,20.0]}

    x_val40 = np.arange(0.00001, 1.0, 1.0/point_par)
    y_val = ecf.LN_func(X_data=x_val40, param=param_val)
    np.random.seed(1729)
    y_noise = noise_par * np.random.normal(size=x_val40.size)
    y_val40 = y_val + y_noise
    
    dat_model_val = [x_val40, y_val40]
    par_model_val = ecf.ln_func_par_opt(D_st=dat_model_val, start_range=config)
    y_val40res = ecf.LN_func(X_data=x_val40, param=par_model_val)

    data ={'x_val': x_val40, 'y_val': y_val40, 'y_valres': y_val40res}
    val_export = pd.DataFrame(data)
    
    return val_export,par_model_val

def validation_fcLN(noise_par, point_par):    
    param_val=[('A_par1', 0.55), ('B_par1', 1.22), ('MSE', 0.0), ('R_squ', 0.0)]
    config={"A":[0.0,20.0],	"B":[0.0,20.0]}

    x_val40 = np.arange(0.00001, 1.0, 1.0/point_par)
    y_val = ecf.LN_func(X_data=x_val40, param=param_val)
    np.random.seed(1729)
    y_noise = noise_par * np.random.normal(size=x_val40.size)
    y_val40 = y_val + y_noise
    
    dat_model_val = [x_val40, y_val40]
    par_model_val = ecf.LN_func_curve_fit(D_st=dat_model_val, start_range=config)
    y_val40res = ecf.LN_func(X_data=x_val40, param=par_model_val)

    data ={'x_val': x_val40, 'y_val': y_val40, 'y_valres': y_val40res}
    val_export = pd.DataFrame(data)
    
    return val_export,par_model_val


def validation_4PL(noise_par, point_par):    
    param_val=[('D_par1', 1.0), ('A_par1', 0.001), ('B_par1', 4.5), ('C_par1', 0.5), ('MSE', 0.0), ('R_squ', 0.0)]
    config={"A":[0.0,10.0],	"B":[0.0,10.0], "C":[0.0,30.0], "D":[0.0,20.0]}

    x_val40 = np.arange(0.0, 1.0, 1.0/point_par)
    y_val = ecf.logit_4PL_func(X_data=x_val40, param=param_val)
    np.random.seed(1729)
    y_noise = noise_par * np.random.normal(size=x_val40.size)
    y_val40 = y_val + y_noise
    
    dat_model_val = [x_val40, y_val40]
    par_model_val = ecf.logit_4PL_par_opt(D_st=dat_model_val, start_range=config)
    y_val40res = ecf.logit_4PL_func(X_data=x_val40, param=par_model_val)

    data ={'x_val': x_val40, 'y_val': y_val40, 'y_valres': y_val40res}
    val_export = pd.DataFrame(data)
      
    return val_export,par_model_val

def validation_fc4PL(noise_par, point_par):    
    param_val=[('D_par1', 1.0), ('A_par1', 0.001), ('B_par1', 4.5), ('C_par1', 0.5), ('MSE', 0.0), ('R_squ', 0.0)]
    config={"A":[0.0,10.0],	"B":[0.0,10.0], "C":[0.0,30.0], "D":[0.0,20.0]}

    x_val40 = np.arange(0.0, 1.0, 1.0/point_par)
    y_val = ecf.logit_4PL_func(X_data=x_val40, param=param_val)
    np.random.seed(1729)
    y_noise = noise_par * np.random.normal(size=x_val40.size)
    y_val40 = y_val + y_noise
    
    dat_model_val = [x_val40, y_val40]
    par_model_val = ecf.logit_4PL_curve_fit(D_st=dat_model_val, start_range=config)
    y_val40res = ecf.logit_4PL_func(X_data=x_val40, param=par_model_val)

    data ={'x_val': x_val40, 'y_val': y_val40, 'y_valres': y_val40res}
    val_export = pd.DataFrame(data)
    
    return val_export,par_model_val

def validation_5PL(noise_par, point_par):    
    param_val=[('D_par1', 1.0), ('A_par1', 0.001), ('B_par1', 4.5), ('C_par1', 0.5), ('E_par1', 2.4), ('MSE', 0.0), ('R_squ', 0.0)]
    config={"A":[0.0,10.0],	"B":[0.0,30.0], "C":[0.0,30.0], "D":[0.0,20.0], "E":[0.0,10.0]}   
    
    x_val40 = np.arange(0.0, 1.0, 1.0/point_par)
    y_val = ecf.logit_5PL_func(X_data=x_val40, param=param_val)
    np.random.seed(1729)
    y_noise = noise_par * np.random.normal(size=x_val40.size)
    y_val40 = y_val + y_noise
    
    dat_model_val = [x_val40, y_val40]
    par_model_val = ecf.logit_5PL_par_opt(D_st=dat_model_val, start_range=config)
    y_val40res = ecf.logit_5PL_func(X_data=x_val40, param=par_model_val)

    data ={'x_val': x_val40, 'y_val': y_val40, 'y_valres': y_val40res}
    val_export = pd.DataFrame(data)

    return val_export,par_model_val

def validation_fc5PL(noise_par, point_par):    
    param_val=[('D_par1',1.0), ('A_par1',0.001), ('B_par1', 4.5), ('C_par1', 0.5), ('E_par1', 2.4), ('MSE', 0.0), ('R_squ', 0.0)]
    config={"A":[0.0,20.0],	"B":[0.0,50.0], "C":[0.0,50.0], "D":[0.0,50.0], "E":[0.0,20.0]}
     
    x_val40 = np.arange(0.0, 1.0, 1.0/point_par)
    y_val = ecf.logit_5PL_func(X_data=x_val40, param=param_val)
    np.random.seed(1729)
    y_noise = noise_par * np.random.normal(size=x_val40.size)
    y_val40 = y_val + y_noise
    
    dat_model_val = [x_val40, y_val40]
    par_model_val = ecf.logit_5PL_curve_fit(D_st=dat_model_val, start_range=config)
    y_val40res = ecf.logit_5PL_func(X_data=x_val40, param=par_model_val)
    print("Validation par fc5PL: ",par_model_val)
    data ={'x_val': x_val40, 'y_val': y_val40, 'y_valres': y_val40res}
    val_export = pd.DataFrame(data)
    
    return val_export,par_model_val

def validation_eng(option, noise, points, val_report):
    
    x_thoer = np.arange(0.0, 1.0, (1.0/(20.0*points)))
    
    if option=="fc4PL-val":
       val_export_csv,par_model = validation_fc4PL(noise_par=noise, point_par=points)
       y_theor = ecf.logit_4PL_func(X_data=x_thoer, param=par_model)
       no_param = 4
    if option=="fc5PL-val":
       val_export_csv,par_model = validation_fc5PL(noise_par=noise, point_par=points)
       y_theor = ecf.logit_5PL_func(X_data=x_thoer, param=par_model)
       no_param = 5
    if option=="fcLN-val":
       val_export_csv,par_model = validation_fcLN(noise_par=noise, point_par=points)
       y_theor = ecf.LN_func(X_data=x_thoer, param=par_model)
       no_param = 2
    if option=="4PL-val":
       val_export_csv,par_model = validation_4PL(noise_par=noise, point_par=points)
       y_theor = ecf.logit_4PL_func(X_data=x_thoer, param=par_model)
       no_param = 4
    if option=="5PL-val":
       val_export_csv,par_model = validation_5PL(noise_par=noise, point_par=points)
       y_theor = ecf.logit_5PL_func(X_data=x_thoer, param=par_model)
       no_param = 5
    if option=="LN-val":
       val_export_csv,par_model = validation_LN(noise_par=noise, point_par=points)
       y_theor = ecf.LN_func(X_data=x_thoer, param=par_model)
       no_param = 2

    y_test = val_export_csv['y_valres']
    Y_st = val_export_csv['y_val']
    # Czy na pewno tutaj?  A moze w funkcji nizej?
    akaike_crit = aic.aic(y=Y_st, y_pred=y_test, p=no_param)
    bayes_crit = bic.bic(y=Y_st, y_pred=y_test, p=no_param)
    r, prob = pearsonr(Y_st, y_test)

    file_name = val_report + "/validation_" + option + str(noise) + "_" + str(points)

    fig = plot_init()
    sub_text = "Validation for " + option + " (noise=" + str(noise) + " , no_point= " + str(points) + ")"
    plt.plot(val_export_csv['x_val'], val_export_csv['y_val'], 'o', color='tab:orange', alpha=0.6, label='Exp points')
    plt.plot(val_export_csv['x_val'], val_export_csv['y_valres'], 'x', color='tab:green', alpha=0.6, linewidth=2, label='Calc points')
    plt.plot(x_thoer, y_theor, '-', color='tab:blue', alpha=0.6, linewidth=2, label='Calc func')     
    plt.legend()
    fig.suptitle(sub_text)
    plt.savefig(file_name + '.png', dpi=150)
    plt.close()               

    val_export_csv.to_csv(file_name + '.csv', index=False, sep=",")

    log_file = open(file_name + '.txt', "w")
    log_file.write('Validation parameters: \n  ')
    log_file.write('Noise = %f, Number points = %d \n' %(noise, points))
    log_file.write('\n')
    log_file.write('Model parameters: \n  ')
    if ((option=="fc5PL-val") or (option=="5PL-val")):
        log_file.write('Absorbance=D+((A-D)/(1+(Conc/C)^B)^E) \n')
        log_file.write('A=%f, B=%f, C=%f, D=%f, E=%f \n' %(par_model[1][1], par_model[2][1], par_model[3][1], par_model[0][1], par_model[4][1]))
        RSS = par_model[5][1]
        RSQ = par_model[6][1]
    if ((option=="fc4PL-val") or (option=="4PL-val")):
        log_file.write('Absorbance=D+((A-D)/(1+(Conc/C)^B) \n')
        log_file.write('A=%f, B=%f, C=%f, D=%f \n' %(par_model[1][1], par_model[2][1], par_model[3][1], par_model[0][1]))
        RSS = par_model[4][1]
        RSQ = par_model[5][1]
    if ((option=="fcLN-val") or (option=="LN-val")):
        log_file.write('Absorbance=A*ln(Conc)+B \n ')
        log_file.write('A=%f, B=%f \n' %(par_model[0][1], par_model[1][1]))
        RSS = par_model[2][1]
        RSQ = par_model[3][1]
    
    log_file.write('\n')
    log_file.write('Model diagnostics: \n')
    log_file.write('  The Residual Sum of Squares RSS    = %f; \n' %(RSS))
    log_file.write('  Coefficient of Determination R^2   = %f; \n' %(RSQ))
    log_file.write('  Akaike Information Criterion AIC   = %f; \n' %(akaike_crit))
    log_file.write('  Bayesian Information Criterion BIC = %f; \n' %(bayes_crit))
    log_file.write('  Coefficient of Correlation r       = %f; \n' %(r))
    log_file.close()
    #data = {'X_theor': x_thoer, 'Y_theor': y_theor, 'X_std': val_export_csv['x_val'], 'Y_std': val_export_csv['y_val'],'X_sam': val_export_csv['x_val'], 'Y_sam': val_export_csv['y_valres']}
    #result_export = pd.DataFrame(data)

    return x_thoer, y_theor, val_export_csv['x_val'], val_export_csv['y_val'], val_export_csv['x_val'],  val_export_csv['y_valres']