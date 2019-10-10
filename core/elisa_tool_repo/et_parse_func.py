#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on January 2019
@author:  Edward Maievskij, Marek Bawiec, Grzegorz Banach
The ELISA tool plugin: module containing low- and mid- level functions for data parsing.
"""

__author__ = "Marek Bawiec, Grzegorz Banach, Edward Maievskij"
__copyright__ = "Copyright 2019, Physiolution Polska"
__credits__ = [""]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Grzegorz Banach"
__email__ = "g.banach@physiolution.pl"
__status__ = "Testing"

import pandas as pd
import re
from string import digits

def tekan_data_check(filepath):
    """
    Function implemented to check the validity of initial TEKAN output raw data_file.
    #Function input:
    filepath - filepath to the file from TEKAN output xlsx measurement file
    """
    pars_dat=pd.read_excel(filepath)
    filetype=filepath.endswith('.xlsx')
    pos_0=pars_dat[pars_dat =='<>'].stack().index.tolist()
    
    exc=((not filetype) or (not pos_0))
    """returns TRUE if bad data were given"""
    return (exc)

def config_data_check_0(filepath):
    """
    Function implemented to check the data formatting essential for calculations in the xlsx-config file.
    Function input:
    filepath - filepath to excel_config file filled by the laboratory assistant
    """
    pars_dat=pd.read_excel(filepath)
    filetype=filepath.endswith('.xlsx')
    pos_0=pars_dat[pars_dat =='<||>'].stack().index.tolist()
    pos_1=pars_dat[pars_dat =='<_>'].stack().index.tolist()
    pos_2=pars_dat[pars_dat =='<|>'].stack().index.tolist()
    exc=((not filetype) or (not pos_1)) or ((not pos_0) or (not pos_2))
    
    """Returns TRUE if bad data were given"""
    return (exc)

def config_data_check_1(filepath):
    """
    Function implemented  to check the data formatting essential for pdf-creation in the xlsx-config file.
    Function input:
    filepath - filepath to excel_config file filled by the laboratory assistant
    """
    pars_dat=pd.read_excel(filepath)
    pos_0=pars_dat[pars_dat =='SOP_name'].stack().index.tolist()
    pos_1=pars_dat[pars_dat =='SOP_version'].stack().index.tolist()
    pos_2=pars_dat[pars_dat =='RESEARCH_name'].stack().index.tolist()
    pos_3=pars_dat[pars_dat =='TEMPLATE_version'].stack().index.tolist()
    pos_4=pars_dat[pars_dat =='EFFECTIVE'].stack().index.tolist()
    pos_5=pars_dat[pars_dat =='DATE_of_exp'].stack().index.tolist()
    pos_6=pars_dat[pars_dat =='TEST_No'].stack().index.tolist()
    pos_7=pars_dat[pars_dat =='PROJECT_NAME'].stack().index.tolist()
    pos_8=pars_dat[pars_dat =='PERFORMED_BY'].stack().index.tolist()
# To do przepisania!!!
#    pos_9=pars_dat[pars_dat =='ABS_MAX'].stack().index.tolist()
#    pos_10=pars_dat[pars_dat =='ABS_MIN'].stack().index.tolist()
    exc=((not pos_3) or (not pos_4) or ((not pos_1)) or (((not pos_0) or (not pos_2)) or ((not pos_5) or (not pos_6))) or ((not pos_7) or (not pos_8)))
    
    """returns TRUE if bad data were given"""
    return (exc)

def extract_pdf_legend(filepath):
    """
    Function implemented for supporting information extracting.
    Function input:
    filepath - filepath to excel_config file filled by the laboratory assistant
    """
    pars_dat=pd.read_excel(filepath)
    
    pdf_pars=['SOP_name','SOP_version','RESEARCH_name','TEMPLATE_version','EFFECTIVE', 'DATE_of_exp','TEST_No','PROJECT_NAME','PERFORMED_BY']
    pdf_legend={}
    
    for nam in pdf_pars:
        pos=pars_dat[pars_dat ==nam].stack().index.tolist()
        start_pos=(pos[0][0],pars_dat.columns.get_loc(pos[0][1]))
        nam_val=str(pars_dat.iloc[start_pos[0],start_pos[1]+1])
        if nam_val =='nan':
            nam_val=' '
            pdf_legend[nam]=str(nam_val)
        else:
            pdf_legend[nam]=str(nam_val)
        
    """returns output in the dictionary format"""  
    return(pdf_legend)
    
def extract_meas_res(dat):
    """
    The function implemented for measurement results extracting from the data parsed from a TEKAN raw xlsx-datafile.
    Function input:
    dat - DataFrame with the data imported from TEKAN output measurement file
    """
    pos_0=dat[dat=='<>'].stack().index.tolist()
    start_pos=(pos_0[0][0],dat.columns.get_loc(pos_0[0][1]))
    imp_data=pd.DataFrame()
    N=0
    while ((N<13) and ((start_pos[1]+N) < (dat.shape[1]))) and (not(pd.isna(dat.iloc[start_pos[0],start_pos[1]+N]))) :
        if N==0:
            imp_data['raw_names']=dat.iloc[start_pos[0]+1:start_pos[0]+9,start_pos[1]].values.tolist()
            N=N+1
        else:
            imp_data[str(int(dat.iloc[start_pos[0],start_pos[1]+N]))]=dat.iloc[start_pos[0]+1:start_pos[0]+9,start_pos[1]+N].values.tolist()
            N=N+1
    imp_data.set_index('raw_names',inplace=True)
    """data output in the DataFrame format"""
    return(imp_data)


def extract_data_map(dat):
    """
    The function implemented for a multiwell plate map extracting from the data parsed from a xlsx-config file.
    """
    pos_1=dat[dat=='<_>'].stack().index.tolist()
    start_pos1=(pos_1[0][0],dat.columns.get_loc(pos_1[0][1]))
    imp_data_map=pd.DataFrame()
    """
    Building of map table <_> from full workshit table
    """
    N=0
    while ((N<13) and ((start_pos1[1]+N) < (dat.shape[1]))) and (not(pd.isna(dat.iloc[start_pos1[0],start_pos1[1]+N]))):
        if N==0:
            imp_data_map['raw_names']=dat.iloc[start_pos1[0]+1:start_pos1[0]+9,start_pos1[1]].values.tolist()
            N=N+1
        else:
            imp_data_map[str(int(dat.iloc[start_pos1[0],start_pos1[1]+N]))]= \
                                 dat.iloc[start_pos1[0]+1:start_pos1[0]+9, start_pos1[1]+N].values.tolist()
            N=N+1
            
    imp_data_map.set_index('raw_names',inplace=True)
        
    imp_columns_map=list(imp_data_map.columns.values)
    """ 
    Replace NaN from empty cell by valu=Empty
    """
    imp_data_map=imp_data_map.fillna(value='Empty')
    """
    Copy of the string in which all case-based characters have been lowercased
    """
    regularizer = lambda x: x.lower()
    """
    Remove spaces before word, and maybe after word?
    """
    ss_sub = lambda x: re.sub(r'[^\w]', '', x)
    """
    Copy of the string in which the occurrences of old have been replaced with new
    """
    sp_replacer = lambda x: x.replace(' ','')
        
    for name in imp_columns_map:
        imp_data_map[name]=imp_data_map[name].apply(ss_sub)
        imp_data_map[name]=imp_data_map[name].apply(regularizer)
        imp_data_map[name]=imp_data_map[name].apply(sp_replacer)
    """
    Returns a multiwell plate map in the DataFrame format
    """
    return(imp_data_map)

def data_map_check(imp_data_map):
    """The function implemented for multiwell plate map validation.
    Function input:
    imp_data_map - DataFrame with multiwell plate map
    """
    imp_columns_map=list(imp_data_map.columns.values)
    remove_digits = str.maketrans('', '', digits)
    translator= lambda x: x.translate(remove_digits)
    ret=False
    for name in imp_columns_map:
        g_check=imp_data_map[name].apply(translator)
        
        for ii in g_check:
            if (ii=='sam')or(ii=='std'):
                continue
            else:
                ret=True
                
    """returns TRUE if bad data were given (e.g. there are names which differ from the standard std and sam notations)"""
    return(ret)
            
        
def extract_data_standards(dat_config):
    """
    The function implemented to extract and process mean values of calibration standards.
    """
    epsilon = 0.0001
    """
    Workeround to avoid zero value in Std concentration, it should not be exist, but user ...
    """
    #counting of the postulated standards
    imp_data_map=extract_data_map(dat=dat_config)
    g1=imp_data_map.stack()
    g2=set(g1)
    N_of_std=0
    for val in g2:
        if ("std" in val):
            N_of_std=N_of_std+1
            continue
        else:
            continue
    print("Number of standards=" + str(N_of_std))
    
    #extracting and processing calibration standards values
    pos_2=dat_config[dat_config=='<|>'].stack().index.tolist()
    start_pos2=(pos_2[0][0],dat_config.columns.get_loc(pos_2[0][1]))

    data_standards=pd.DataFrame()
    regularizer = lambda x: x.lower()
    sp_replacer = lambda x: x.replace(' ','')
    
    N=1
    std_list=[]
    
    while N<=N_of_std:
        std_list.append(dat_config.iloc[start_pos2[0]+N,start_pos2[1]])
        N=N+1
#    print(std_list)   
    """
    Name of  the standards without free spaces and lowercase
    std_name - column for name of the standards
    value - column for value of concentration of standard
    """
    data_standards["std_names"]=std_list
    data_standards["values"]=dat_config.iloc[start_pos2[0]+1:start_pos2[0]+N,start_pos2[1]+1].values.tolist()
    data_standards["std_names"]=data_standards["std_names"].apply(regularizer)
    data_standards["std_names"]=data_standards["std_names"].apply(sp_replacer)
    data_standards.set_index('std_names',inplace=True)

    print("data_standards: ", data_standards)   
    
    data_standards=data_standards.replace({'values': {0: epsilon}}) 
    print("Test replace data_standards: ", data_standards)   
    """
    returns DataFrame with the shortened version of calibration standards table
    """
    return(data_standards)

def extract_specification(dat):
    """
    The function implemented to extract sample specification from the parsed excel_config file.
    Function input:
    dat - DataFrame parsed from excel_config file
    """
    pos_1=dat[dat=='<||>'].stack().index.tolist()
    start_pos1=(pos_1[0][0]+1,dat.columns.get_loc(pos_1[0][1]))
    imp_data_spec=pd.DataFrame()
    
    regularizer = lambda x: x.lower()
    ss_sub = lambda x: re.sub(r'[^\w]', '', x)
    sp_replacer = lambda x: x.replace(' ','')
    
    abbr=[]
    descr=[]
    
    l_extr_dat=len(dat)
    N=0
    while (((start_pos1[0]+N)<l_extr_dat) and (not(pd.isna(dat.iloc[start_pos1[0]+N,start_pos1[1]])))):
        abbr.append(dat.iloc[start_pos1[0]+N,start_pos1[1]])
        descr.append(dat.iloc[start_pos1[0]+N,start_pos1[1]+1])
        
        N=N+1
        
    imp_data_spec['Abbr']=abbr
    imp_data_spec['Abbr']=imp_data_spec['Abbr'].apply(ss_sub)
    imp_data_spec['Abbr']=imp_data_spec['Abbr'].apply(regularizer)
    imp_data_spec['Abbr']=imp_data_spec['Abbr'].apply(sp_replacer)
    imp_data_spec['Description']=descr
    imp_data_spec.sort_values(by=['Abbr'])
    """
    returns DataFrame with the extended samples specification
    """
    return (imp_data_spec)
        
def get_full_specification(Fpath_config):
    """The function implemented to extend sample specification (e.g. with the samples positions).
    Fpath_config - filepath to some excel_config file filled by the laboratory assistant
    """
    extr_data = pd.read_excel(Fpath_config)
    spec_short = extract_specification(dat=extr_data)  
    dat = pd.read_excel(Fpath_config)
    data_map = extract_data_map(dat=dat)
    
    spec_full=spec_short.fillna("-")
    
    N=0
    
    positions=[]
    while N<len(spec_short.iloc[:,0]):
        position=data_map[data_map==spec_short.iloc[N,0]].stack().index.tolist()
        
        lp=len(position)
        p_n=0
        position_wr=[]
        
        while p_n<lp:
            position_wr.append("-".join(position[p_n]))
            p_n=p_n+1
            
        position_wr= ", ".join(position_wr)
        
        positions.append(position_wr)
        N=N+1
    spec_full['Position']=positions
    """returns DataFrame with extended specification"""
    return(spec_full)
    

def check_function(Fpath_tekan, Fpath_config):
    """
    The function implemented for "parse and check" button.
    Function input:
    Fpath_tekan, Fpath_config - filepaths to raw data from TEKAN measuremnt system and excel_config file respectively
    """
    text_output="Appropriate data set"
    control_1=tekan_data_check(filepath=Fpath_tekan)
    if control_1 is True:
        text_output="Initial data was not marked. Inappropriate data style or file format." 
    else:
        control_2=config_data_check_0(filepath=Fpath_config)
        if control_2 is True:
            text_output=text_output="Data presented in xlsx-config file was not marked. Inappropriate xlsx-configurational data style or file format."
        else:
            control_3=config_data_check_1(filepath=Fpath_config)
            if control_3 is True:
                text_output=text_output="Configurational template for footers and headers is absent or has inappropriate form."
    """returns text output for the info field of GUI"""        
    return(text_output)
    