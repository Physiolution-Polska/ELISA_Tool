#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on January 2019
@author:  Edward Maievskij, Marek Bawiec, Grzegorz Banach
The ELISA tool plugin: module for pdf export functions.
"""

__author__ = "Marek Bawiec, Grzegorz Banach, Edward Maievskij"
__copyright__ = "Copyright 2019, Physiolution Polska"
__credits__ = [""]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Grzegorz Banach"
__email__ = "g.banach@physiolution.pl"
__status__ = "Testing"

from weasyprint import HTML
import datetime
import pandas as pd
import textwrap
from os import remove

def rej_samples_html_table(D_sam):
    """
    The function was implemented to generate html-formatted table with the rejected sample names.
    """
    g_names=D_sam[3]
    dftable=pd.DataFrame()
    dftable["bad_sample_names"]=g_names
    dftable.sort_values(by=["bad_sample_names"], inplace=True)
        
    html_table=dftable.to_html()
    return(html_table)


def rep_html_table(D_sam):
    """
    The function was implemented to create a HTML-table named "Calculation results" from the list of list containing the data about sample names,
    titer concentration, bottom edge ot titer concentration standard deviation, top edge of titer concentration standard deviation, sample absorbance
    value and simmetrical value of absorbance standard deviation.
    """
    g_names=D_sam[0]
    concentration_X=D_sam[1]
    conc_st_dev_bottom= D_sam[5]
    conc_st_dev_upper= D_sam[6]
    intensity_Y=D_sam[2]
    abs_st_dev=D_sam[4]
    dftable=pd.DataFrame()
    dftable["good_sample_names"]=g_names
    dftable["concentration_X_axis"]=concentration_X
    dftable["C_SD_bot"]=conc_st_dev_bottom
    dftable["C_SD_up"]=conc_st_dev_upper
    dftable["Absorbance_Y_axis"]=intensity_Y
    dftable["Abs_st_dev"]=abs_st_dev
    dftable.sort_values(by=["good_sample_names"], inplace=True)
    dftable.set_index(['good_sample_names'], inplace=True)
    
    ln_html_table=dftable.to_html()
    """returns html-formatted "calculation results" table"""
    return(ln_html_table)


def rep_csv_order(D_sam, rep_name, i_data, i_data_map, d_st, p_folder, parameters, order, cal_time):
    """
    The function was implemented to create and save csv report for 5PL model.
    """
    g_names=D_sam[0]
    concentration_X=D_sam[1]
    conc_st_dev_bottom= D_sam[5]
    conc_st_dev_upper= D_sam[6]
    intensity_Y=D_sam[2]
    abs_st_dev=D_sam[4]
    ts=datetime.datetime.now()
    ts_name=ts.strftime("%d_%b_%Y_%H_%M")
    dftable=pd.DataFrame()
    dftable["good_sample_names"]=g_names
    dftable["concentration_X_axis"]=concentration_X
    dftable["c_st_dev_bo"]=conc_st_dev_bottom
    dftable["c_st_dev_up"]=conc_st_dev_upper
    dftable["Absorbance_Y_axis"]=intensity_Y
    dftable["Abs_st_dev"]=abs_st_dev
    dftable.sort_values(by=["good_sample_names"], inplace=True)
    dftable.set_index(['good_sample_names'], inplace=True)
    
    s_csv_res=dftable.to_csv(path_or_buf=None)
    
    g_names=D_sam[3]
    dftable1=pd.DataFrame()
    dftable1["bad_sample_names"]=g_names
    dftable1.sort_values(by=["bad_sample_names"], inplace=True)
        
    s_csv_rej=dftable1.to_csv(path_or_buf=None)
    
    s_csv_i_data=i_data.to_csv(path_or_buf=None)
    
    s_csv_i_data_map=i_data_map.to_csv(path_or_buf=None)
    
    s_csv_d_st=d_st.to_csv(path_or_buf=None)
    
    doc_structure='Initial measurement results'+ '\r\n'+ s_csv_i_data + '\r\n'+ 'Multiwell plate map' +'\r\n'+s_csv_i_data_map + \
        '\r\n'+ 'Calibration standards'+'\r\n'+ s_csv_d_st + '\r\n'+ 'Calculation results'+'\r\n'+s_csv_res + '\r\n'+ 'BLQ samples'+ \
        '\r\n'+s_csv_rej
    
    if ((order=="LN") or (order=="fcLN")):
        s=doc_structure + '\r\n'+'Coefficient of Determination R^2,'+ str(round(parameters[3][1],8))+ '\r\n'\
                        +'Akaike Information Criterion AIC,'+ str(round(parameters[4][1],8))+ '\r\n'\
                        +'Bayesian Information Criterion BIC,'+ str(round(parameters[5][1],8))+ '\r\n'\
                        +'Coefficient of Correlation r,'+ str(round(parameters[6][1],8))+ '\r\n'\
                        +'The Residual Sum of Squares RSS,' + str(round(parameters[2][1],8))+'\r\n'\
                        +'\r\n' + 'Absorbance=A*ln(Conc)+B' + '\r\n'\
                        +'A,' + str(round(parameters[0][1],8))+'\r\n'\
                        +'B,'+ str(round(parameters[1][1],8)) + '\r\n'\
                        +'Time of calculation,' + str(round(cal_time,6))

    if ((order=="4PL") or (order=="fc4PL")):
        s=doc_structure + '\r\n'+'Coefficient of Determination R^2,'+ str(round(parameters[5][1],8))+ '\r\n'\
                        +'Akaike Information Criterion AIC,'+ str(round(parameters[6][1],8))+ '\r\n'\
                        +'Bayesian Information Criterion BIC,'+ str(round(parameters[7][1],8))+ '\r\n'\
                        +'Coefficient of Correlation r,'+ str(round(parameters[8][1],8))+ '\r\n'\
                        +'The Residual Sum of Squares RSS,' + str(round(parameters[4][1],8))+'\r\n'\
                        + 'Absorbance=D+((A-D)/(1+(Conc/C)^B)'+'\r\n'\
                        + 'D,' + str(round(parameters[0][1],8))+'\r\n' + 'A,'+str(round(parameters[1][1],8)) + '\r\n' + 'B,' + str(round(parameters[2][1],8))+'\r\n'+'C,'+str(round(parameters[3][1],8)) + '\r\n'\
                        +'Time of calculation,' + str(round(cal_time,6))

    if ((order=="5PL") or (order=="fc5PL")):
        s=doc_structure + '\r\n'+'Coefficient of Determination R^2,'+ str(round(parameters[6][1],8)) + '\r\n'\
                        +'Akaike Information Criterion AIC,'+ str(round(parameters[7][1],8))+ '\r\n'\
                        +'Bayesian Information Criterion BIC,'+ str(round(parameters[8][1],8))+ '\r\n'\
                        +'Coefficient of Correlation r,'+ str(round(parameters[9][1],8))+ '\r\n'\
                        +'The Residual Sum of Squares RSS,' + str(round(parameters[5][1],8))+'\r\n'\
                        +'\r\n' + 'Absorbance=D+((A-D)/(1+(Conc/C)^B)^E)'+'\r\n'\
                        +'D,' + str(round(parameters[0][1],8)) + '\r\n' + 'A,'+ str(round(parameters[1][1],8)) + '\r\n'\
                        +'B,' + str(round(parameters[2][1],8)) + '\r\n' + 'C,'+ str(round(parameters[3][1],8)) + '\r\n' + 'E,' + str(round(parameters[4][1],8)) + '\r\n'\
                        +'Time of calculation,' + str(round(cal_time,6))
    
    file= open(p_folder + '/' + rep_name + '_' + order + '_' + ts_name + "_data.csv","w")
    file.write(s)
    file.close() 


    
def pdf_report_generation(rep_name, i_data, i_data_map, d_st, p_folder, plot_name, D_samples, 
                             parameters, legend, plot_name1, spec, notice, order, cal_time):
    """
    The function was implemented in order to create logarithmic model pdf-report with the help of the functions from weasyprint module.
    """
    #variables definitions
    ts=datetime.datetime.now()
    ts1=ts.strftime("%d.%m.%Y %H:%M")
    ts_name=ts.strftime("%d_%m_%Y_%H_%M")
    expdate=legend['DATE_of_exp']
    projname=legend['PROJECT_NAME']
    sopname=legend["SOP_name"]
    resname=legend["RESEARCH_name"]
    test_no=legend["TEST_No"]
    sopver=legend["SOP_version"]
    username=legend['PERFORMED_BY']
    tempver=legend['TEMPLATE_version']
    effectivetemplate=legend['EFFECTIVE']
    if len(resname)>65:
        s61=textwrap.wrap(text=resname, width=65)
        resname= ' \A '.join(s61)
    
    if len(notice)>99:
        s71=textwrap.wrap(text=notice, width=99)
        notice= ' \A '.join(s71)
        
    if ((order=="LN") or (order=="fcLN")):
        core_name="_" + order + "_out_report_"
        report_name="Logarithmic"
        formula="Absorbance=A*ln(Conc)+B"

    if ((order=="4PL") or (order=="fc4PL")):
        core_name="_logit_" + order + "_out_report_"
        report_name="Logit " + order
        formula="Absorbance=D+((A-D)/(1+(Conc/C)^B)"

    if ((order=="5PL") or (order=="fc5PL")):
        core_name="_logit_" + order + "_out_report_"
        report_name="Logit " + order
        formula="Absorbance=D+((A-D)/(1+(Conc/C)^B)^E)"
       
    i_data=i_data.to_html()
    i_data_map=i_data_map.to_html()
    d_st=d_st.to_html()
    
    
    g_sample_table=rep_html_table(D_sam= D_samples)
    rej_sam= rej_samples_html_table(D_sam= D_samples)
    spec_html_table=spec.to_html()
    
    #html-body and css-formatting
    s=("""<html>
        <head>
            <meta name="pdfkit-page-size" content="A4"/>
            <meta name="pdfkit-orientation" content="Portrait"/>
        </head>
        <body>
        
          <div class="firstpage">
          <p class="firstframe"><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>Project:  {proj_name}<br><br></p>
          <p class="firstframe1">Date:  {exdate}<br><br>Performed by:  {uname}<br><br>Approved by:<br><br>Test No.:  {tno}<br><br>on:  {timestamp1}<br><br>on:<br><br></p>
          <div class="elem">
          
              <h3>Sample pipetting scheme</h3>
              <p>{spec_html}</p>
          </div>
          <div class="elem">
              <h3>Multiwell plate map</h3>
              <p>{data_map}</p>
              <h3>Initial measurement results</h3>
              <p>{data_in}</p>          
              <h3>Calibration standards</h3>
              <p>{cal_st}</p>
          </div>
        """).format(spec_html=spec_html_table, exdate=expdate, proj_name=projname,  uname=username, tno= test_no, timestamp1=ts1, 
                     timestamp=ts1, data_in=i_data, data_map=i_data_map, cal_st=d_st)
    
    if ((order=="LN") or (order=="fcLN")):     
        sb=("""
            <div class="elem">
                <h2>Model: {R_name} fitting</h2>
                <p style="font-size:12px"> 
                    {form}<br>
                </p>
                
                <h3>Model parameters</h3>
                <p style="font-size:12px">
                    A={A:.6f}, B={B:.6f} <br>
                </p>

                <h3>Model diagnostics</h3>
                <p style="font-size:12px"> The Residual Sum of Squares RSS    ={rss:.6f}<br>
                                           Coefficient of Determination R^2   ={rsq:.6f}<br>
                                           Akaike Information Criterion AIC   ={AIC:.6f}<br>
                                           Bayesian Information Criterion BIC ={BIC:.6f}<br>
                                           Coefficient of Correlation r       ={rcc:.6f}<br>
                                           Time of calculatin                 ={cal_time:.6f} [s]<br>
                </p>
            </div>
            """).format(R_name=report_name, form=formula, rss=parameters[2][1], rsq=parameters[3][1], AIC=parameters[4][1], BIC=parameters[5][1], rcc=parameters[6][1], \
                                      A=parameters[0][1], B=parameters[1][1], cal_time=cal_time)

    if ((order=="4PL") or (order=="fc4PL")): 
        sb=("""
            <div class="elem">
                <h2>Model: {R_name} fitting</h2>
                <p style="font-size:12px"> 
                    {form}<br>
                </p>
                
                <h3>Model parameters</h3>
                <p style="font-size:12px">
                      A={A:.6f}, B={B:.6f}, C={C:.6f}, D={D:.6f}</p>
                </p>

                <h3>Model diagnostics</h3>
          
                <p style="font-size:12px"> The Residual Sum of Squares RSS    ={rss:.6f}<br>
                                           Coefficient of Determination R^2   ={rsq:.6f}<br>
                                           Akaike Information Criterion AIC   ={AIC:.6f}<br>
                                           Bayesian Information Criterion BIC ={BIC:.6f}<br>
                                           Coefficient of Correlation r       ={rcc:.6f}<br>
                                           Time of calculatin                 ={cal_time:.6f} [s]<br>
                </p>
          </div>
            """).format(R_name=report_name, form=formula, rss=parameters[4][1], rsq=parameters[5][1], AIC=parameters[6][1], BIC=parameters[7][1], rcc=parameters[8][1], \
                                      D=parameters[0][1], A=parameters[1][1], B=parameters[2][1], C=parameters[3][1], cal_time=cal_time)

    if ((order=="5PL") or (order=="fc5PL")): 
        sb=("""
          <div class="elem">
                <h2>Model: {R_name} fitting</h2>
                <p style="font-size:12px"> 
                    {form}<br>
                </p>              
                <h3>Model parameters</h3>
                <p style="font-size:12px">
                    A={A:.6f}, B={B:.6f}, C={C:.6f}, D={D:.6f}, E={E:.6f}</p>
                </p>
                <h3>Model diagnostics</h3>
                <p style="font-size:12px"> The Residual Sum of Squares RSS    ={rss:.6f}<br>
                                           Coefficient of Determination R^2   ={rsq:.6f}<br>
                                           Akaike Information Criterion AIC   ={AIC:.6f}<br>
                                           Bayesian Information Criterion BIC ={BIC:.6f}<br>
                                           Coefficient of Correlation r       ={rcc:.6f}<br>                 
                                           Time of calculatin                 ={cal_time:.6f} [s]<br>
          </div>
            """).format(R_name=report_name, form=formula, rss=parameters[5][1], rsq=parameters[6][1], AIC=parameters[7][1], BIC=parameters[8][1], rcc=parameters[9][1], \
                                      D=parameters[0][1], A=parameters[1][1], B=parameters[2][1], C=parameters[3][1], E=parameters[4][1], cal_time=cal_time)
                 
    sc=("""
          <div class="elem">
              <h3>Plots</h3>
              <img src="{plot}" alt="resulting plot">
              <br>
              <img src="{plot1}" alt="resulting plot">
          </div>
          <div class="elem">
              <h2 class="calcres">Calculation results</h2>
              <p>{table1}</p>
          </div>
          <div class="elem">
              <h2>BLQ samples</h2>
              <p>{table2}</p>
          </div>
          
        </body>
        </html>
        """).format(plot=plot_name, pathfolder=p_folder, plot1=plot_name1, table1=g_sample_table, table2=rej_sam)
    
    
    s1=(""" <!doctype html>
        <html lang="pl">
        <head>
            <meta charset="utf-8">
            <meta name="pdfkit-page-size" content="A4"/>
            <meta name="pdfkit-orientation" content="Portrait"/>
            <style>
            table, th, td{
                text-align: center;
                border: 1px solid black;
                border-collapse: collapse;
                width: auto;
                font-size: 10px;
                padding: 2px;}
            table {
                text-align: center;
                width: 100%;
                white-space: nowrap;
                page-break-inside: avoid;
                width: auto;}
            .elem {
            page-break-inside: avoid;
            }
            @page {
                   @top-left{
                               display: inline-block;
                               width: 165px;
                               height: 75px;
                               margin-right: 5px;
                              
                               content: "";
                               background: url("template/Logo.png") no-repeat 0 0;
                               background-size: 70%;
                   
                                text-align: left;}
                   
                   
                   @top-center {
                               content:" """ )
#     content: url("../template/Logo.png");
    s2=("{sop_name} \A {res_name}").format(sop_name=sopname, res_name=resname)
    s3=("""";
    font-size: 11px;
    text-align: left;
    white-space: pre;
    max-width: 90mm;
    overflow-wrap:normal;
    word-wrap:break-word;}
    @top-right {
               content:" """)
    s4=("Ver: {sop_ver} \A \A Page").format(sop_ver=sopver)
    s5=(""" " counter(page) " of " counter(pages);
               text-align: right;
               font-size: 11px;
               white-space: pre;
               }
        @bottom-center {
                              content: " """)
    s6=("{noti}").format(noti=notice)
    s7=("""";
                              font-size: 10px;
                              text-align: left;
                              white-space: pre;
                              max-width: 90mm;
                              overflow-wrap:normal;
                              }
        @bottom-left {
                        content:" """)
    s8=("   Template: \A     Version:{temp_version} \A     Effective:{eff}").format(temp_version=tempver, eff=effectivetemplate)
    s9=(""" ";
                        font-size: 10px;
                        text-align: left;
                        white-space: pre;}                      
    }
            img {
                max-width: 100%;
                height: 400px;
                text-align: center;
                }
            p.firstframe1{
                        column-count: 2;
                        font-size: 12px;
                        text-align: left;}
            .firstframe{
                        
                        font-size: 12px;
                        text-align: left;
                        white-space: pre;}
            h2.calcres{
                        page-break-before: always;}
            h3 {
                font-size: 12px;
            }
            h2 {font-size: 12px;
                }
            h1 {font-size: 12px;
                }
            
            
            </style>""")
    s= s1+s2+s3+s4+s5+s6+s7+s8+s9+s+sb+sc
    #temporary html-file generation
    Html_file= open("AAA_tempfile1.html","w")
    Html_file.write(s)
    Html_file.close()

    #conversion to pdf with further removing of temporary files
    HTML("AAA_tempfile1.html").write_pdf(p_folder + '/' + rep_name + core_name + ts_name + '.pdf')
    remove("AAA_tempfile1.html")
    #remove(plot_name)
    #remove(plot_name1)
    return("pdf report generated")

