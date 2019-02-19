
import os
from datetime import datetime as date
import time
import warnings # suppress warnings
warnings.filterwarnings("ignore")

## 3rd-party modules
import pandas as pd
import numpy as np
import copy

from GAR import wb
from GAR.globals import read_parameters_global, read_partition_groups, show_message, add_logsheet
from .relation import gen_relation
from .scenario_compare import scenario_compare
#from .condqgreg import condquant

###############################################################################
#%% Functions for step 4: scenario test
###############################################################################
def do_scenario(debug=False):
    '''
    Entry point function called when button for cenario fits is called.

    This function cannot take in any arguments and can not have
    any return values due to limitation of the RunPython VBA code
    in xlwings.
    '''
    t0 = time.time()
    
    # Make sure a wb exists
    if wb is None:
        print('scenario_main: wb is None')
        print('This may be due to not calling set_mock_caller_file')
        print('and setting the caller Workbook')
        import sys
        sys.exit(-1)
    else:
        print (wb)


    dict_input_scenario, df_collection_scenario = prerun_scenario(debug=debug)
    

    dict_output_scenario = run_scenario(dict_input_scenario, df_collection_scenario, debug=debug)
    postrun_scenario(dict_output_scenario, debug=debug)

    # End measurement of time
    t1 = time.time()

    # Total time for this operation (formatted string)
    tdiff = "{:.1f}".format(t1 - t0)
    
    sheetname = dict_output_scenario['sheet_scenario']
    message = 'Finished with scenario test in ' + tdiff + ' sec,\n'
    message += 'output is in sheets ' + ', '+sheetname
    show_message(message,msgtype='info')
    
def prerun_scenario(debug=False):
    '''
    Prerun function for step 2, quantfit.
    
    This function cannot take in any arguments due to limitations
    of the RunPython VBA code in xlwings.

    Check that the necessary steps beforehand have been done.
    Read in/check the input parameters and return a
    dict for input parameters and a df for data.
    '''

    # Check that the necessary steps beforehand have been done.
    
    # --------------------------
    # Read in parameters
    # --------------------------
    dict_input_scenario = read_parameters_scenario()
    
    # --------------------------
    # Check parameter values
    # --------------------------
    check_parameters_scenario(dict_input_scenario)

    # Read in global parameters
    # (this also checks if values have changed since being initially set)
    dict_global_params = read_parameters_global()
    #print(dict_input_scenario)
    #print(dict_global_params)


    # Add each key, val from dict_global_params to dict_input_quantfit
    for key, val in dict_global_params.items():
        # Check that the keys do not clash
        if key in dict_input_scenario:
            message = 'dict_input_quantfit should not have key ' + key + ' that is common with dict_global_params'
            show_message(message, halt=True)
        dict_input_scenario[key] = val

                        
    # Create df for data
    input_sheetnames = [dict_input_scenario['sheet_partition'],dict_input_scenario['sheet_quantreg'],dict_input_scenario['sheet_cond_quant']]
    sheet_partition=dict_input_scenario['sheet_partition']
    df_scenario_collections = read_data_scenario(input_sheetnames)
    df_scenario_collections['Data']=df_scenario_collections['Data'][(df_scenario_collections['Data'].index>= df_scenario_collections[sheet_partition]['date'].values[0]) & (df_scenario_collections['Data'].index<= df_scenario_collections[sheet_partition]['date'].values[-1])]

    # return a dict for input parameters and a df
    return dict_input_scenario,df_scenario_collections

def read_parameters_scenario():
    '''
    Read in parameters for quantfit.

    The cell positions for inputs are hardcoded in.
    Parameter value and range checking is done in the check_parameters_quantfit function.
    '''

    # Create dict for parameters
    dict_parameters_scenario = dict()

    # ---------------------------#
    # Read in necessary values
    # ---------------------------#

    # Read in quantlist
    cellpos = 'F31'
    # Get a list of all values starting at cellpos going down
    dict_parameters_scenario['quantlist'] = wb.sheets['Input_parameters'].range(cellpos).expand('down').value

    # Read in info on regressors.
    dict_parameters_scenario['regressors'] = dict()
    # Start with the first cell that contains the regressors.
    startrow = 31
    cellpos = 'A' + str(startrow)
    # Read down and get a list of all regressors
    regressors = wb.sheets['Input_parameters'].range(cellpos).expand('down').value
    if not isinstance(regressors, (list, tuple)):
        regressors = [regressors]
    for iregressor, regressor in enumerate(regressors):
       

        # Get the cells for transformation and optional parameter
        # from columns B and C for this row
        colnum = startrow+iregressor
        transform = wb.sheets['Input_parameters'].range('B' + str(colnum)).value
        option    = wb.sheets['Input_parameters'].range('D' + str(colnum)).value
        dict_parameters_scenario['regressors'][regressor+'_trans_'+str(iregressor)+'_'+transform] = dict()
        # Set as values of dict
        dict_parameters_scenario['regressors'][regressor+'_trans_'+str(iregressor)+'_'+transform]['transform'] = transform
        dict_parameters_scenario['regressors'][regressor+'_trans_'+str(iregressor)+'_'+transform]['option'] = option
    
    # Read in latest_date
    cellpos = 'B61'
    dict_parameters_scenario['latest_date'] = wb.sheets['Input_parameters'].range(cellpos).value
    

    # Read in t-skew fit parameters
    dict_parameters_scenario['fit_params'] = dict()
    cellpos = 'B59'
    dict_parameters_scenario['fit_params']['fittype'] = wb.sheets['Input_parameters'].range(cellpos).value
    
    dict_parameters_scenario['fit_params']['mode'] = dict()
    cellpos = 'B64'
    dict_parameters_scenario['fit_params']['mode']['constraint'] = wb.sheets['Input_parameters'].range(cellpos).value
    cellpos = 'D64'
    dict_parameters_scenario['fit_params']['mode']['value'] = wb.sheets['Input_parameters'].range(cellpos).value
    dict_parameters_scenario['fit_params']['qsmooth'] = dict()
    cellpos = 'B67'
    dict_parameters_scenario['fit_params']['qsmooth']['option'] = wb.sheets['Input_parameters'].range(cellpos).value
    cellpos = 'D67'
    dict_parameters_scenario['fit_params']['qsmooth']['period'] = wb.sheets['Input_parameters'].range(cellpos).value
    
    # Advanced t-skew fit parameters
    # Process each parameter in order
    pos = 74
    for param in ['dof',  'var_low', 'var_high', 'skew_low', 'skew_high']:
        dict_parameters_scenario['fit_params'][param] = dict()
        cellpos = 'B' + str(pos)
        dict_parameters_scenario['fit_params'][param]['constraint'] = wb.sheets['Input_parameters'].range(cellpos).value
        cellpos = 'D' + str(pos)
        dict_parameters_scenario['fit_params'][param]['value'] = wb.sheets['Input_parameters'].range(cellpos).value
        pos += 1
    
    # Read in shocked t-skew fit parameters
    dict_parameters_scenario['fit_params_shocked'] = copy.deepcopy(dict_parameters_scenario['fit_params'])
    cellpos = 'B110'
    dict_parameters_scenario['fit_params_shocked']['mode']['constraint'] = wb.sheets['Input_parameters'].range(cellpos).value
    cellpos = 'D110'
    dict_parameters_scenario['fit_params_shocked']['mode']['value'] = wb.sheets['Input_parameters'].range(cellpos).value
    
    # Read in scenario parameters
    
    dict_parameters_scenario['shockvars'] = dict()
    startrow=99
    cellpos = 'A'+str(startrow)    
    shockvars=wb.sheets['Input_parameters'].range(cellpos).expand('down').value
    if not isinstance(shockvars, (list, tuple)):
        shockvars = [shockvars]
        
    for ind,shockvar in enumerate(shockvars):
        col=startrow+ind
        shocktype = wb.sheets['Input_parameters'].range('B' + str(col)).value
        if shocktype=='By +/- STD':
            shockvalue = wb.sheets['Input_parameters'].range('F' + str(col)).value
        else:
            shockvalue = wb.sheets['Input_parameters'].range('D' + str(col)).value
    
        dict_parameters_scenario['shockvars'][shockvar] = dict()
        dict_parameters_scenario['shockvars'][shockvar]['shocktype']=shocktype
        dict_parameters_scenario['shockvars'][shockvar]['shockvalue']=shockvalue
    

    # Read in output sheet
    cellpos = 'B112'
    dict_parameters_scenario['sheet_scenario'] = wb.sheets['Input_parameters'].range(cellpos).value



    startrow = 51
    for isheetname, sheetname in enumerate(['sheet_quantreg', 'sheet_cond_quant']):
        colnum = startrow + isheetname
        cellpos = 'B' + str(colnum)
        dict_parameters_scenario[sheetname] = wb.sheets['Input_parameters'].range(cellpos).value
        # checking of values will be done in check_parameters_quantfit

    # Read in partitions sheet
    cellpos = 'B24'
    dict_parameters_scenario['sheet_partition'] = wb.sheets['Input_parameters'].range(cellpos).value
    dict_parameters_scenario['partition_groups']= read_partition_groups()
    return dict_parameters_scenario

def check_parameters_scenario(dict_input_scenario):
    '''
    Check the input parameters for scenario.
    '''
    # Check that all keys exist
    keys=dict_input_scenario.keys()    

    input_sheets = ['Readme', 'Input_parameters', 'Partition_groups', 'Data', 'Processing_Log']
    for key in keys:
        val = dict_input_scenario[key]

        if key == 'quantlist':
            # Check that all values are between 0 and 1,
            # and that the necessary values of 0.10, 0.25, 0.50, 0.75, 0.90 exist
            vals_np = np.array(val) # create np.array for value checking
            if not (np.all(0 < vals_np) and np.all(vals_np < 1)):
                message = 'All values of quantlist must be between 0 and 1'
                message+= 'Given values: ' + str(val)
                show_message(message, halt=True)
            # Check that necessary values are present
            necessary_vals = [0.10, 0.25, 0.50, 0.75, 0.90]
            for _val in necessary_vals:
                if _val not in val:
                    message = 'Value of ' + str(_val) + ' must be included in quantlist'
                    message += 'Given values: ' + str(val)
                    show_message(message, halt=True)

        if key == 'regressors':
            # val is a dict of dicts with keys [regressor]['transform/option']
            for regressor in val:
                transform = val[regressor]['transform']
                option    = val[regressor]['option']

                # Check that transform is a valid value from the pulldown menu
                if transform not in ['None', 'Lagged', 'MVA','Power','Diff','ChangeRate']:
                    message = 'transform for ' + regressor + ' was not a valid option, given ' + transform
                    show_message(message, halt=True)

                # If 'No transformation' or 'Log' was chosen, make sure no option was given
                if transform in ['None']:
                    if option is not None:
                        message = 'option for regressor = ' + regressor + ' with transform of ' + transform + ' must not have option set'
                        show_message(message, halt=True)
                    
                # If 'Lagged' or 'Moving Average' was chosen, make sure lag exists and is an int
                if transform in ['Lagged', 'MVA','Power','Diff','ChangeRate']:
                    if type(option) != float or abs(int(option) - option) > 1E-5:
                        message = 'option for regressor = ' + regressor + ' with transform of ' + transform + ' must have option of int, given ' + str(option)
                        show_message(message, halt=True)
                    # Since the value is less than 1E-5 away from an int,
                    # convert to int so that there are no problems later
                    dict_input_scenario['regressors'][regressor]['option'] = int(option)
                    
        if key == 'shockvars':
        
            for shockvar in val:
                shocktype = val[shockvar]['shocktype']
                shockvalue= val[shockvar]['shockvalue']

                # Check that shocktype is a valid value from the pulldown menu
                if shocktype not in ['None', 'By +/- STD','By +/- percentage']:
                    message = 'Shock type for ' + shockvar + ' was not a valid option, given ' + shocktype
                    show_message(message, halt=True)

                if shocktype in ['By +/- STD','By +/- percentage']:
                    if abs(shockvalue)>10:
                        message = 'Shock value for variable = ' + shockvar + ' with shocktype of ' + shocktype + ' must have option of int, given ' + str(shockvalue)
                        show_message(message, halt=True)

            
        elif key.find('sheet_') != -1:
            if val is None:
                if key == 'sheet_quantreg':
                    dict_input_scenario[key] = 'Quant reg coefficients'
                elif key == 'sheet_cond_quant':
                    dict_input_scenario[key] = 'Conditional quantiles'
                elif key == 'sheet_local_proj':
                    dict_input_scenario[key] = 'Local projections'
                elif key == 'sheet_partition':
                    dict_input_scenario[key] = 'Output_partitions'
                elif key == 'sheet_scenario':
                    dict_input_scenario[key] = 'Senario test'
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
                    
def read_data_scenario(sheetnames):
    '''
    Read in the input data for scenario.
    Checks for the sheetname should have been done in check_parameters_quantfit.
    '''
    data_collection={}
    for sheetname in sheetnames:
        print (sheetname)
        data_collection[sheetname] = wb.sheets[sheetname].range('A1').options(pd.DataFrame,index=False,expand='table').value
        if 'date' in data_collection[sheetname].columns:
            data_collection[sheetname].set_index(data_collection[sheetname]['date'],inplace=True)
       
    
    data_collection['Data'] = wb.sheets['Data'].range('A1').options(pd.DataFrame,index=False,expand='table').value     
    data_collection['Data'].set_index('date', inplace=True, drop=False)
    data_collection['Data'].fillna(method='ffill',inplace=True)
    return data_collection

def run_scenario(dict_input_scenario, df_collection_scenario, debug=False):
    

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
    dict_output_scenario = dict()
    
    # ------------------------
    # Copy the output sheet names
    # from dict_input_tsfit
    # ------------------------
    for key in dict_input_scenario:
        if key.find('sheet_') != -1:
            dict_output_scenario[key] = dict_input_scenario[key]
   
    df_shockedvar,df_shockedgrp = gen_relation(dict_input_scenario['shockvars'],dict_input_scenario['partition_groups'],df_collection_scenario['Data'],df_collection_scenario[dict_input_scenario['sheet_partition']])             
    df_output=pd.DataFrame(index=df_shockedvar.index)

    for shockvar in dict_input_scenario['shockvars'].keys():
        if shockvar not in dict_input_scenario['partition_groups']:
            df_output[shockvar+'_raw']=df_collection_scenario['Data'][shockvar]
            df_output[shockvar+'_shocked']=df_shockedvar[shockvar+'_shocked']
    df_output['Groups : ']=np.nan

        
    regressors=list(dict_input_scenario['regressors'].keys())
    


    for reg_long in  regressors:
        reg_short=reg_long.split('_trans_')[0]
        df_output[reg_short]=df_shockedgrp[reg_short]
        df_output[reg_short+'_shocked']=df_shockedgrp[reg_short+'_shocked']
        if dict_input_scenario['regressors'][reg_long]['transform']=='None':
            df_shockedgrp[reg_long]=df_shockedgrp[reg_short]
            df_shockedgrp[reg_long+'_shocked']=df_shockedgrp[reg_short+'_shocked']
        elif dict_input_scenario['regressors'][reg_long]['transform']=='Lagged':
            df_shockedgrp[reg_long]=df_shockedgrp[reg_short].shift(dict_input_scenario['regressors'][reg_long]['option'])
            df_shockedgrp[reg_long+'_shocked']=df_shockedgrp[reg_short+'_shocked'].shift(dict_input_scenario['regressors'][reg_long]['option'])            
        elif dict_input_scenario['regressors'][reg_long]['transform']=='MVA':
            df_shockedgrp[reg_long]=df_shockedgrp[reg_short].rolling(window=dict_input_scenario['regressors'][reg_long]['option']).mean()
            df_shockedgrp[reg_long+'_shocked']=df_shockedgrp[reg_short+'_shocked'].rolling(window=dict_input_scenario['regressors'][reg_long]['option']).mean()
        elif dict_input_scenario['regressors'][reg_long]['transform']=='Power': 
            df_shockedgrp[reg_long]=df_shockedgrp[reg_short]**(dict_input_scenario['regressors'][reg_long]['option'])
            df_shockedgrp[reg_long+'_shocked']=df_shockedgrp[reg_short+'_shocked']**(dict_input_scenario['regressors'][reg_long]['option'])            
        elif dict_input_scenario['regressors'][reg_long]['transform']=='Diff': 
            df_shockedgrp[reg_long]=df_shockedgrp[reg_short].diff(dict_input_scenario['regressors'][reg_long]['option'])
            df_shockedgrp[reg_long+'_shocked']=df_shockedgrp[reg_short+'_shocked'].diff(dict_input_scenario['regressors'][reg_long]['option']) 
        elif dict_input_scenario['regressors'][reg_long]['transform']=='ChangeRate': 
            df_shockedgrp[reg_long]=df_shockedgrp[reg_short].pct_change(dict_input_scenario['regressors'][reg_long]['option'])
            df_shockedgrp[reg_long+'_shocked']=df_shockedgrp[reg_short+'_shocked'].pct_change(dict_input_scenario['regressors'][reg_long]['option']) 
            
    #print(df_shockedgrp.iloc[-1])
    #print(df_shockedgrp.iloc[-2])
    horizon = dict_input_scenario['horizon']
    fitdate = dict_input_scenario['latest_date']
    fitparam=dict_input_scenario['fit_params']
    fitparam_shocked=dict_input_scenario['fit_params_shocked']
       
    # Fitdat
    cond_quant_raw={}
    cond_quant_shocked={}
    if dict_input_scenario['fit_params']['qsmooth']['option']=='None':
        try:
            df_partition_fit=df_shockedgrp[df_shockedgrp['date']==fitdate]
        except:
            df_partition_fit=df_shockedgrp.iloc[-1,:]
            print ('The latest date in the data will be used.')
                           
    elif dict_input_scenario['fit_params']['qsmooth']['option']=='Median':
        per=int(dict_input_scenario['fit_params']['qsmooth']['period'])
        df_partition_fit=df_shockedgrp[df_shockedgrp['date']<=fitdate].drop(['date','date_shocked'],axis=1)       
        df_partition_fit=df_partition_fit.tail(per)
        df_partition_fit.loc['median']=df_partition_fit.median()
        df_partition_fit=df_partition_fit.loc[['median'],:]
        
    elif dict_input_scenario['fit_params']['qsmooth']['option']=='Mean':        
        per=int(dict_input_scenario['fit_params']['qsmooth']['period'])
        df_partition_fit=df_shockedgrp[df_shockedgrp['date']<=fitdate].drop(['date','date_shocked'],axis=1)
        df_partition_fit=df_partition_fit.tail(per)
        #print(df_partition_fit)
        df_partition_fit.loc['mean']=df_partition_fit.mean()
        df_partition_fit=df_partition_fit.loc[['mean'],:]

    df_quantcoef=df_collection_scenario[dict_input_scenario['sheet_quantreg']]
    for quant in dict_input_scenario['quantlist']:
        cond_quant_raw[quant]=df_quantcoef[(df_quantcoef['variable']=='Intercept') & (df_quantcoef['quantile']==quant)]['coeff_noscale'].values[0]    
        cond_quant_shocked[quant]=df_quantcoef[(df_quantcoef['variable']=='Intercept') & (df_quantcoef['quantile']==quant)]['coeff_noscale'].values[0]
        for reg_long in  regressors:

            cond_quant_raw[quant]+=df_quantcoef[(df_quantcoef['variable']==reg_long) & (df_quantcoef['quantile']==quant)]['coeff_noscale'].values[0]*df_partition_fit[reg_long].values[0]
            cond_quant_shocked[quant]+=df_quantcoef[(df_quantcoef['variable']==reg_long) & (df_quantcoef['quantile']==quant)]['coeff_noscale'].values[0]*df_partition_fit[reg_long+'_shocked'].values[0]        
    ols_raw=df_quantcoef[(df_quantcoef['variable']=='Intercept') & (df_quantcoef['quantile']=='mean')]['coeff_noscale'].values[0]    
    ols_shocked=df_quantcoef[(df_quantcoef['variable']=='Intercept') & (df_quantcoef['quantile']=='mean')]['coeff_noscale'].values[0] 
    for reg_long in  regressors:
        ols_raw+=df_quantcoef[(df_quantcoef['variable']==reg_long) & (df_quantcoef['quantile']=='mean')]['coeff_noscale'].values[0]*df_partition_fit[reg_long].values[0]
        ols_shocked+=df_quantcoef[(df_quantcoef['variable']==reg_long) & (df_quantcoef['quantile']=='mean')]['coeff_noscale'].values[0]*df_partition_fit[reg_long+'_shocked'].values[0]
    print("***********")
    print(ols_raw,ols_shocked)
    fig,res,dfpdf= scenario_compare(cond_quant_raw,cond_quant_shocked,fitparam,fitparam_shocked,horizon,fitdate,ols_raw,ols_shocked)
    dict_output_scenario['fig'] = fig
    dict_output_scenario['res'] = res
    dict_output_scenario['dfpdf'] = dfpdf
    dict_output_scenario['data']= df_output      
        
    return dict_output_scenario

    
def postrun_scenario(dict_output_scenario, debug=False):
    '''
    Postrun function for step 2, scenario.

    Takes as input dict from main run function.

    This function cannot return any values due to limitations
    of the RunPython VBA code in xlwings.
    
    '''

    if debug:
        print('=' * 30)
        print('start of postrun_scenario')
        print('=' * 30)

    # Create DataFrame for log
    log_frame = pd.DataFrame(columns=['Time','Action'])
    
    sheetname = dict_output_scenario['sheet_scenario']
    
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
            wb.sheets[sheetname].api.Tab.Colorindex = 23
            action = 'Created sheet ' + sheetname
    except:
        print('Unable to acess '+sheetname)
         
    try:
        wb.sheets[sheetname].pictures[0].delete()
    except:
        pass

    tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
    
    sheet = wb.sheets[sheetname]    

    pd.options.display.float_format = '{:,.4f}'.format
    df_pdf= dict_output_scenario['dfpdf']
    df_data = dict_output_scenario['data']
    # Write out scenario resul
    res = dict_output_scenario['res']
    wb.sheets[sheetname].range('A1').value = res 
    wb.sheets[sheetname].range('N1').options(index=False).value =df_pdf
    wb.sheets[sheetname].range('T1').options(index=True).value =df_data
    wb.sheets[sheetname].range('T1').value='Variables : '
    fig = dict_output_scenario['fig']

    # Set the path of the output file to be in the same dir as the
    # calling Excel file
    fullpath = os.path.abspath(os.path.dirname(wb.fullname) + '/figures')
    if not os.path.isdir(fullpath):
        os.makedirs(fullpath)
    outfilename = fullpath+'\\scenario_'+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.png'
    try:
        fig.savefig(outfilename)
    except:
         print('Fail to save scenario figure.')
    
    try:
        sheet.pictures.add(fig, name='MyPlot', update=True, left=sheet.range('B12').left, top=sheet.range('B12').top, height=250, width=500)
        action = 'Scenario figure saved'
    except:
        action = 'Unable to add figure to sheet ' + sheetname
        
    wb.sheets[sheetname].autofit()
    # Add to log
    tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
    log = pd.Series({'Time': tn, 'Action': action})
    log_frame = log_frame.append(log, ignore_index=True)

    # Write out log_frame
    add_logsheet(wb, log_frame, colnum=3)
