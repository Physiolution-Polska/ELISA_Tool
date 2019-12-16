#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 2019
@author:  Edward Maievskij, Marek Bawiec, Grzegorz Banach
The ELISA tool plugin: module for graphics presentation of results.
"""

__author__ = "Marek Bawiec, Grzegorz Banach, Edward Maievskij"
__copyright__ = "Copyright 2019, Physiolution Polska"
__credits__ = [""]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Grzegorz Banach"
__email__ = "g.banach@physiolution.pl"
__status__ = "Production"

import matplotlib.pyplot as plt

def draw_and_save_plot(X_fun, Y_fun, X_std, Y_std, X_sam, Y_sam, title_string, order, X_sb_err, X_su_err, Y_s_err, Y_std_err):
    """
    The function was implemented to draw a simple plot in linear scale with a theoretical 4PL curve, calibration standards and measured samples.
    X_fun, Y_fun - theoretical curve data points
    X_std, Y_std - calibration standards data points
    X_sam, Y_sam - sample measurement points
    """
    Fpath= 'calc_results_common_' + order + '.png'

    #plt.figure(figsize=(7,4), dpi=150)
    
    if order=="lin":
        plt.plot(X_fun, Y_fun, label='Theor calibration')
        plt.errorbar(X_std, Y_std, yerr=Y_std_err, fmt='o', color='r', ecolor='black', label='Calibration standards')
        plt.errorbar(X_sam, Y_sam, xerr=[X_sb_err, X_su_err], yerr=Y_s_err, fmt='v', color='orange', ecolor='lightgray', label='Samples')
        plt.legend(loc='lower right', ncol=2)
    if order=="semilog":
        plt.semilogx(X_fun, Y_fun, label='Theor calibration')
        plt.semilogx(X_std, Y_std, 'o', color='r',label='Calibration standards')
        plt.semilogx(X_sam, Y_sam, 'v', label='Samples')
        plt.legend(loc='upper left', ncol=2)
        
    plt.grid(True) 
    plt.tick_params(axis='both', which='both', direction= 'in')
    plt.xlabel("Concentration")
    plt.ylabel("Absorbance (A.U.)")
    plt.title(title_string)
#    fig = plt.figure(figsize=(4,3))
    plt.savefig(Fpath, format='png', bbox_inches = 'tight', dpi=300)
#    fig.savefig(Fpath, dpi=150, format='png', bbox_inches = 'tight')
    plt.close()
    """Returns the filepath to obtained plot."""
    return(Fpath)
