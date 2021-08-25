
import os
from datetime import datetime as date
import time
import warnings # suppress warnings
warnings.filterwarnings("ignore")

## 3rd-party modules
import pandas as pd
import numpy as np

from GAR import wb
from GAR.globals import read_parameters_global, read_partition_groups, show_message, add_logsheet
from .historical_gen import historical_gen
from matplotlib.backends.backend_pdf import PdfPages
#from .historical import historicalseq
#from .condqgreg import condquant

###############################################################################
#%% Functions for step 4: historical test
###############################################################################
def do_historical(debug=False):
    '''
    Entry point function called when button for cenario fits is called.

    This function cannot take in any arguments and can not have
    any return values due to limitation of the RunPython VBA code
    in xlwings.
    '''
    t0 = time.time()
    
    # Make sure a wb exists
    if wb is None:
        print('historical_main: wb is None')
        print('This may be due to not calling set_mock_caller_file')
        print('and setting the caller Workbook')
        import sys
        sys.exit(-1)
    else:
        print (wb)


    dict_input_historical, df_historical = prerun_historical(debug=debug)
       
    dict_output_historical = run_historical(dict_input_historical, df_historical, debug=debug)
    
    postrun_historical(dict_output_historical, debug=debug)

    # End measurement of time
    t1 = time.time()

    # Total time for this operation (formatted string)
    tdiff = "{:.1f}".format(t1 - t0)
    
    sheetname = dict_output_historical['sheet_historical']
    message = 'Finished with historical test in ' + tdiff + ' sec,\n'
    message += 'output is in sheets ' + ', '+sheetname
    show_message(message,msgtype='info')
    
def prerun_historical(debug=False):
    '''
    Prerun function for step 2, quantfit.
    
    This function cannot take in any arguments due to limitations
    of the RunPython VBA code in xlwings.

    Check that the necessary steps beforehand have been done.
    Read in/check the input parameters and return a
    dict for input parameters and a df for data.
    '''

    if debug:
        print('=' * 30)
        print('start of prerun_quantfit')
        print('=' * 30)

    # Check that the necessary steps beforehand have been done.
    
    # --------------------------
    # Read in parameters
    # --------------------------
    dict_input_historical = read_parameters_historical()
    
    # --------------------------
    # Check parameter values
    # --------------------------
    check_parameters_historical(dict_input_historical)

    # Read in global parameters
    # (this also checks if values have changed since being initially set)
    dict_global_params = read_parameters_global()



    # Add each key, val from dict_global_params to dict_input_quantfit
    for key, val in dict_global_params.items():
        # Check that the keys do not clash
        if key in dict_input_historical:
            message = 'dict_input_quantfit should not have key ' + key + ' that is common with dict_global_params'
            show_message(message, halt=True)
        dict_input_historical[key] = val

                        
    # Create df for data
    sheetname = dict_input_historical['sheet_cond_quant']
    
    print (dict_input_historical)
    df_historical = read_data_historical(sheetname)

    # return a dict for input parameters and a df
    return dict_input_historical,df_historical

def read_parameters_historical():
    '''
    Read in parameters for historical.

    The cell positions for inputs are hardcoded in.
    Parameter value and range checking is done in the check_parameters_quantfit function.
    '''

    # Create dict for parameters
    dict_parameters_historical = dict()

    # ---------------------------#
    # Read in necessary values
    # ---------------------------#
   
    # Read in start_date
    cellpos = 'B86'
    dict_parameters_historical['start_date'] = wb.sheets['Input_parameters'].range(cellpos).value
    # Read in end_date
    cellpos = 'B87'
    dict_parameters_historical['end_date'] = wb.sheets['Input_parameters'].range(cellpos).value    
    # Read in time period incremental
    cellpos = 'B88'
    dict_parameters_historical['time_inc'] = int(wb.sheets['Input_parameters'].range(cellpos).value)
    
    # Read in t-skew fit parameters
    dict_parameters_historical['fit_params'] = dict()
    
    cellpos = 'B59'
    dict_parameters_historical['fit_params']['fittype'] = wb.sheets['Input_parameters'].range(cellpos).value
    dict_parameters_historical['fit_params']['mode'] = dict()
    dict_parameters_historical['fit_params']['mode']['constraint'] = 'Free'
    dict_parameters_historical['fit_params']['mode']['value'] = None
    # Advanced t-skew fit parameters
    # Process each parameter in order
    pos = 74
    for param in ['dof',  'var_low', 'var_high', 'skew_low', 'skew_high']:
        dict_parameters_historical['fit_params'][param] = dict()
        cellpos = 'B' + str(pos)
        dict_parameters_historical['fit_params'][param]['constraint'] = wb.sheets['Input_parameters'].range(cellpos).value
        cellpos = 'D' + str(pos)
        dict_parameters_historical['fit_params'][param]['value'] = wb.sheets['Input_parameters'].range(cellpos).value
        pos += 1   

    # Read in output sheet
    cellpos = 'B91'
    dict_parameters_historical['sheet_historical'] = wb.sheets['Input_parameters'].range(cellpos).value


    startrow = 51
    for isheetname, sheetname in enumerate(['sheet_quantreg', 'sheet_cond_quant']):
        colnum = startrow + isheetname
        cellpos = 'B' + str(colnum)
        dict_parameters_historical[sheetname] = wb.sheets['Input_parameters'].range(cellpos).value
        # checking of values will be done in check_parameters_quantfit

    return dict_parameters_historical
def check_parameters_historical(dict_input_historical):
    '''
    Check the input parameters for historical.
    '''
    # Check that all keys exist
    keys=dict_input_historical.keys()    

    input_sheets = ['Readme', 'Input_parameters', 'Partition_groups', 'Data', 'Processing_Log']
    for key in keys:
        val = dict_input_historical[key]                 
        if key.find('sheet_') != -1:
            if val is None:
                if key == 'sheet_quantreg':
                    dict_input_historical[key] = 'Quant reg coefficients'
                elif key == 'sheet_cond_quant':
                    dict_input_historical[key] = 'Conditional quantiles'
                elif key == 'sheet_historical':
                    dict_input_historical[key] = 'Historical distribution'
                else:
                    message = 'No sheet called ' + key + ' should exist'
                    show_message(message, halt=True)
            else:
                # Check that the specified sheetname is not one of the inputs
                # (it is OK that it is the same name as an existing sheet if
                # that sheet is an output)
                if val in input_sheets:
                    message = key + ' specified as ' + val + ', cannot be the same as necessary input sheet'
                    show_message(message, halt=True)
                    
def read_data_historical(sheetname):
    '''
    Read in the input data for historical.
    Checks for the sheetname should have been done in check_parameters_quantfit.
    '''
    dall = wb.sheets[sheetname].range('A1').options(pd.DataFrame,index=False,expand='table').value
    return dall

def run_historical(dict_input_historical, df_historical, debug=False):
    

    '''
    Main run function for step 2, quantfit.

    Takes in as arguments a dict for input parameters
    and a df for data. Outputs a dict for output parameters.

    Does quantile fits and returns a dict of output parameters.
    ** This function should be independent of any Excel input/output
    and be executable as a regular Python function independent of Excel. **
    '''
    
    warnings.filterwarnings("ignore")

    if debug:
        print('=' * 30)
        print('start of run_quantfit')
        print('=' * 30)

    # ------------------------
    # Create output dict
    # ------------------------
    dict_output_historical = dict()
    for key in dict_input_historical:
        if key.find('sheet_') != -1:
            dict_output_historical[key] = dict_input_historical[key]
    # ------------------------
    # Copy the output sheet names
    # from dict_input_tsfit
    # ------------------------
    sdate=dict_input_historical['start_date']
    edate=dict_input_historical['end_date']
    dates=df_historical[(df_historical['tau']==0.5) & (df_historical['date']>=sdate) & (df_historical['date']<=edate) ]['date']
    dates=dates.iloc[list(range(0,len(dates),dict_input_historical['time_inc']))].values
    
    cond_quants=[]
    realvalues=[]
    olsmeans=[]
    for d in dates:
        clst = df_historical[(df_historical['tau']!='mean') & (df_historical['date']==d)]['tau'].values
        qval = df_historical[(df_historical['tau']!='mean') & (df_historical['date']==d)]['conditional_quantile_mean'].values
        olsmeans.append(df_historical[(df_historical['tau']=='mean') & (df_historical['date']==d)]['conditional_quantile_mean'].values[0])
        cond_quants.append(dict(zip(clst,qval)))
        realvalues.append(df_historical[(df_historical['tau']=='mean') & (df_historical['date']==d)]['realized_value'].values[0])
    fitparam=dict_input_historical['fit_params']
    figs,res,chartpacks = historical_gen(cond_quants,fitparam,dates,realvalues,olsmeans)
    df=pd.DataFrame(res)
    df.index=dates
    dict_output_historical['figs'] = figs
    dict_output_historical['charts'] = chartpacks
    dict_output_historical['data'] = df
    #dict_output_historical['data']= df_output      
        
    return dict_output_historical

    
def postrun_historical(dict_output_historical, debug=False):
    '''
    Postrun function for step 2, historical.

    Takes as input dict from main run function.

    This function cannot return any values due to limitations
    of the RunPython VBA code in xlwings.
    
    '''

    if debug:
        print('=' * 30)
        print('start of postrun_historical')
        print('=' * 30)

    # Create DataFrame for log
    log_frame = pd.DataFrame(columns=['Time','Action'])
    
    sheetname = dict_output_historical['sheet_historical']
    
        # Get existing sheetnames
    sheetnames = [sheet.name for sheet in wb.sheets]

    try:
        # Clear the sheet if it already exists
        if sheetname in sheetnames:
            wb.sheets[sheetname].clear()
            action = 'Cleared sheet ' + sheetname
            # Otherwise add it after the "Data" sheet
        else:
            wb.sheets.add(sheetname, after='Data')
                # Set output sheet colors to blue
            wb.sheets[sheetname].api.Tab.ColorIndex = 23
            action = 'Created sheet ' + sheetname
    except:
        print('Unable to acess '+sheetname)
         
    for p in wb.sheets[sheetname].shapes:
        try: 
            p.delete()
        except Exception as e: 
            print(e)

    tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
    
    sheet = wb.sheets[sheetname]    


    # Write out historical results
    figs = dict_output_historical['figs']
    res = dict_output_historical['data']
    charts= dict_output_historical['charts']
    # Set the path of the output file to be in the same dir as the
    # calling Excel file
    fullpath = os.path.abspath(os.path.dirname(wb.fullname) + '/figures')
    if not os.path.isdir(fullpath):
        os.makedirs(fullpath)
        
        
    outfilename = fullpath+'\\historical_'+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.png'
    try:
        figs['res'].savefig(outfilename)
    except:
         print('Fail to save historical figure.')
    
    try:
        pic=sheet.pictures.add(figs['res'], name='MyPlot', update=True, left=sheet.range('B30').left, top=sheet.range('B30').top, height=1700, width=480)
        pic.height=1700
        pic.width=480
        action = 'historical figure saved'
    except:
        action = 'Unable to add figure to sheet ' + sheetname
    
    
    outfilename = fullpath+'\\pittest_'+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.png'
    try:
        figs['pit'].savefig(outfilename)
    except:
         print('Fail to save pit figure.')
         
 
         
    try:    
        pic=sheet.pictures.add(figs['pit'], name='MyPlot2', update=True, left=sheet.range('B3').left, top=sheet.range('B3').top, height=360, width=350)
        pic.height=360
        pic.width=350
        action = 'historical figure saved'
    except:
        action = 'Unable to add figure to sheet ' + sheetname    
        
    outfilename = fullpath+'\\logscore_'+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.png'
    try:
        figs['ls'].savefig(outfilename)
    except:
         print('Fail to save logscore figure.')     
         
    try:    
        pic=sheet.pictures.add(figs['ls'], name='MyPlot3', update=True, left=sheet.range('J3').left, top=sheet.range('J3').top, height=360, width=480)
        pic.height=360
        pic.width=480
        action = 'historical figure saved'
    except:
        action = 'Unable to add figure to sheet ' + sheetname    
        
    try:
        wb.sheets[sheetname].range('U1').value = res
    except:
        action='Unable to output historical result.'    
        
    try:
        pdfile=fullpath+'\\historicalcharts_'+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.pdf'
        pp = PdfPages(pdfile)
        for e in charts:
            pp.savefig(e)              
        pp.close()
    except:
        print('Unable to save PDF chartpacks.')    
    wb.sheets[sheetname].autofit()
    # Add to log
    tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
    log = pd.Series({'Time': tn, 'Action': action})
    log_frame = log_frame.append(log, ignore_index=True)

    # Write out log_frame
    add_logsheet(wb, log_frame, colnum=3)
