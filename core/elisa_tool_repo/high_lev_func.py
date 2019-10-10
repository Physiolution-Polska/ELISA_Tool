#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on January 2019
@author:  Edward Maievskij, Marek Bawiec, Grzegorz Banach
The ELISA tool plugin: module contains high-level functions which combine all of the data flows in the program.
"""

__author__ = "Marek Bawiec, Grzegorz Banach, Edward Maievskij"
__copyright__ = "Copyright 2019, Physiolution Polska"
__credits__ = [""]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Grzegorz Banach"
__email__ = "g.banach@physiolution.pl"
__status__ = "Testing"

#import sys
import pandas as pd
import elisa_tool_repo.et_calc as ecf
import elisa_tool_repo.et_parse_func as epf
import elisa_tool_repo.et_plot_func as edpf
import elisa_tool_repo.pdf_creator as pdf

def parse_input_data(Fpath_tekan, Fpath_config):
    """
    The function used to parse, collect and process data for ln-model report generation.
    Function inputs Fpath_tekan and Fpath_config should contain filepaths to TEKAN raw data xlsx-file and standard xlsx config file respectively.
    """

    dat = pd.read_excel(Fpath_tekan)
    meas_res = epf.extract_meas_res(dat=dat)
    """
    means_res: read measurement table <> as DataFrame from tekan xls file, returns a DataFrame with initial measurement results
    """

    dat = pd.read_excel(Fpath_config)
    data_map = epf.extract_data_map(dat=dat)
    """
    data_map: read mapping table <_> as DataFrame from config xls file, returns a DataFrame with multiwell plate map
    """

    dat = pd.read_excel(Fpath_config)
    data_standards = epf.extract_data_standards(dat_config=dat)
    """
    data_standards: reduce number of std and read std concentration table <|> as DataFrame from config xls file:
    df include (std_name, concentration), returns calibration standards DataFrame from the excel_config file filled by the laboratory assistant
    """

    data_standards_for_rep = ecf.data_st_to_print(imp_data=meas_res, imp_data_map=data_map, data_standards=data_standards)
    """
    data_standards_for_rep: generate DataFrame(result) of data calibration standard
    result["name"]=sample
    result["conc"]=X_st
    result["absorbance"]=Y_st, abs_ave
    """

    pdf_leg = epf.extract_pdf_legend(Fpath_config)
    """
    pdf_leg: experiment information extracting from config xls file as dictionary for final pdf report;
    """

    specification = epf.get_full_specification(Fpath_config)
    """
    specification: extend sample specification (e.g. with the samples positions in the plate)
    """   

    dat_model = ecf.data_std_format(imp_data=meas_res, imp_data_map=data_map, data_standards=data_standards) 
    """
    dat_model: implemented to extract, clean and order calibration standards for further numerical-models fitting. Remove zero value from concentration.
    """   
    
    return(meas_res, data_map, data_standards, data_standards_for_rep, pdf_leg, specification, dat_model)
    #['Y_names', 'Y_abs', 'samples_nb', 'Y_ave_abs', 'Y_abs_std', 'X_concentration', 'X_s_conc_std_bottom', 'X_s_conc_std_upper', 'Errors']

def interpolation_eng(dat_model, order, param_range):
    """
    The function used to parse, collect and process data for ln-model report generation.
    Function inputs Fpath_tekan and Fpath_config should contain filepaths to TEKAN raw data xlsx-file and standard xlsx config file respectively.
    """
    print("param_range in ENG: ",param_range)
    if order=="LN":
        par_model=ecf.ln_func_par_opt(D_st=dat_model, start_range=param_range)
    if order=="4PL":
        par_model=ecf.logit_4PL_par_opt(D_st=dat_model, start_range=param_range)
    if order=="5PL":
        par_model=ecf.logit_5PL_par_opt(D_st=dat_model, start_range=param_range) 
    if order=="fcLN":
        par_model=ecf.LN_func_curve_fit(D_st=dat_model, start_range=param_range)
    if order=="fc4PL":
        par_model=ecf.logit_4PL_curve_fit(D_st=dat_model, start_range=param_range)
    if order=="fc5PL":
        par_model=ecf.logit_5PL_curve_fit(D_st=dat_model, start_range=param_range)
    return(par_model)

"""
output: (meas_res, data_map, data_standards, data_standards_for_rep, pdf_leg, specification, dat_model)
"""
def make_report(std_mat, X_std, Y_std, param, res_folder, r_name, meas_res, data_map, 
                data_standards, data_standards_for_rep, specification, pdf_leg, licence_notice, order_model, cal_time):
    #################################   Tutaj szukaj: data_standards_for_rep vs data_standards
    """
    The function used to generate pdf and csv reports for 5PL model.
    Function inputs Fpath_tekan and Fpath_config should contain filepaths to TEKAN raw data xlsx-file and standard xlsx config file respectively.
    Function input res_folder is used to provide information about destination folder for report generation opertions.
    Function input r_name is used to provide suffix for genrated reports (i.e. first part of the reports name)
    """
    # to jest do rozszycia! parsing_sam_and_concenration zawiera parsowanie sam i obliczenia koncentracji (!)  
    result_model=ecf.samples_concenration(imp_data=meas_res, imp_data_map=data_map, params=param, ordered=order_model, d_st=data_standards_for_rep)
    # dataframe['Y_names', 'Y_abs', 'samples_nb', 'Y_ave_abs', 'Y_abs_std', 'X_concentration', 'Errors']
    #print('z HL result_model: ', result_model)
    #print('Param approx for',order_model,': ', param)
    #sys.exit()
    """        
    result_model=ecf.parsing_sam_and_concenration(imp_data=meas_res, imp_data_map=data_map, params=param, ordered=order_model)
                0             1               2         3               4                5                 6            7        8          9
    output: (Y_s_names, X_s_concentration, Y_s_good, Y_s_bad, X_s_conc_std_bottom, X_s_conc_std_upper, Y_s_abs_std, X_std_abs, Y_std_abs, Y_std_std)
    """
    
    theor_model_X_axis=ecf.theor_X(res_sam_mat=result_model, std_mat=std_mat)
    if ((order_model=="LN") or (order_model=="fcLN")):
        theor_model_Y_axis=ecf.LN_func(X_data=theor_model_X_axis, param=param)
    if ((order_model=="4PL") or (order_model=="fc4PL")):
        theor_model_Y_axis=ecf.logit_4PL_func(X_data=theor_model_X_axis, param=param)
    if ((order_model=="5PL") or (order_model=="fc5PL")):
        theor_model_Y_axis=ecf.logit_5PL_func(X_data=theor_model_X_axis, param=param)

    theor_val = ecf.remove_abs_smaller_0(X_data=theor_model_X_axis, Y_data=theor_model_Y_axis)
    
    plot_name = edpf.draw_and_save_plot(X_fun=theor_val['X'], Y_fun=theor_val['Y'], 
                                       X_std=result_model[7], Y_std=result_model[8], X_sam=result_model[1], Y_sam=result_model[2], 
                                       title_string='Results for '+order_model+' logit curve', order="lin", 
                                       X_sb_err=result_model[4], X_su_err=result_model[5], Y_s_err=result_model[6], Y_std_err=result_model[9])

    plot_name1 = edpf.draw_and_save_plot(X_fun=theor_val['X'], Y_fun=theor_val['Y'], 
                                        X_std=result_model[7], Y_std=result_model[8], X_sam=result_model[1], Y_sam=result_model[2], 
                                        title_string='Results for '+order_model+' logit curve (semilog)', order="semilog",
                                        X_sb_err=result_model[4], X_su_err=result_model[5], Y_s_err=result_model[6], Y_std_err=result_model[9])
   
    pdf.pdf_report_generation(rep_name=r_name, p_folder=res_folder, plot_name=plot_name, D_samples=result_model, 
                              parameters=param, i_data=meas_res, i_data_map=data_map, d_st=data_standards_for_rep, legend=pdf_leg, 
                              plot_name1=plot_name1, spec=specification, notice=licence_notice, order=order_model, cal_time=cal_time)

    pdf.rep_csv_order(rep_name=r_name, p_folder=res_folder, D_sam=result_model, 
                      parameters=param, i_data=meas_res, i_data_map=data_map, d_st=data_standards_for_rep, order=order_model, cal_time=cal_time)
    
    """ 
    X_fun=theor_val['X'], Y_fun=theor_val['Y']     - Theor calibration
    X_std=run_step_1[6][0], Y_std=run_step_1[6][1] - Calibration standards
    X_std=result_model[7], Y_std=result_model[8]   - Calibration standards better ordered
    X_sam=result_model[1], Y_sam=result_model[2]   - Samples
    """
    #data = {'X_theor': theor_val['X'], 'Y_theor': theor_val['Y'], 'X_std': X_std, 'Y_std': Y_std, 'X_sam': result_model[1], 'Y_sam': result_model[2]}
    #result_export = pd.DataFrame(data)
    #return result_export

    return theor_val['X'], theor_val['Y'], X_std, Y_std, result_model[1], result_model[2]