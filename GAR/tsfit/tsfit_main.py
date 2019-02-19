
import os
from datetime import datetime as date
from matplotlib.backends.backend_pdf import PdfPages
from glob import glob
import time
import warnings # suppress warnings
warnings.filterwarnings("ignore")

## 3rd-party modules
import pandas as pd
import numpy as np


from GAR import wb
from GAR.globals import read_parameters_global, read_partition_groups, show_message, add_logsheet

from .gen_skewt import gen_skewt

# Separate the plotting into a separate plotRemove these plot
import matplotlib.pyplot as plt                       ## Plotting

###############################################################################
#%% Plotting
###############################################################################
## Style of the charts
plt.style.use('seaborn-white')
plt.close('all')  

###############################################################################
#%% Functions for step 3: tsfit
###############################################################################
def do_tsfit(debug=False):
    '''
    Entry point function called when button for t-skew fits is called.

    This function cannot take in any arguments and can not have
    any return values due to limitation of the RunPython VBA code
    in xlwings.
    '''

    # Start measurement of time
    t0 = time.time()
    
    # Make sure a wb exists
    if wb is None:
        print('partition_main: wb is None')
        print('This may be due to not calling set_mock_caller_file')
        print('and setting the caller Workbook')
        import sys
        sys.exit(-1)
    
    if debug:
        print('+' * 40)
        print('start of do_tsfit')
        print('+' * 40)

    # Call prerun
    if debug:
        print('---- calling prerun_tsfit')
    dict_input_tsfit, df_tsfit_collections = prerun_tsfit(debug=debug)
        
    # Call main run
    if debug:
        print('---- calling run_tsfit')
    dict_output_tsfit = run_tsfit(dict_input_tsfit, df_tsfit_collections, debug=debug)

    # Call postrun
    if debug:
        print('---- calling postrun_tsfit')
    postrun_tsfit(dict_output_tsfit, debug=debug)

    # End measurement of time
    t1 = time.time()

    # Total time for this operation (formatted string)
    tdiff = "{:.1f}".format(t1 - t0)
    
    sheet = dict_output_tsfit['sheet_tsfit']
    message = 'Finished with skewed T distribution in ' + tdiff + ' sec,\n'
    message += 'output is in sheets ' + sheet
    show_message(message,msgtype='info')


    
def prerun_tsfit(debug=False):
    '''
    Prerun function for step 3, tsfit.
    
    Check that the necessary steps beforehand have been done.
    Read in/check the input parameters and return a
    dict for input parameters and a df for data.
    '''

    if debug:
        print('=' * 30)
        print('start of prerun_tsfit')
        print('=' * 30)

    # Keys for input parameter dict
    keys = ['latest_date', 'fit_params', 'sheet_tsfit']

    # Read in parameters
    dict_input_tsfit = read_parameters_tsfit()
    # Check parameter values
    check_parameters_tsfit(dict_input_tsfit, keys)

    # Read in global parameters
    # (this also checks if values have changed since being initially set)
    dict_global_params = read_parameters_global()

    # Add each key, val from dict_global_params to dict_input_tsfit
    for key, val in dict_global_params.items():
        # Check that the keys do not clash
        if key in dict_input_tsfit:
            message = 'dict_input_tsfit should not have key ' + key + ' that is common with dict_global_params'
            show_message(message, halt=True)
        dict_input_tsfit[key] = val
 
    # Create df for data
    input_sheetnames = [dict_input_tsfit['sheet_partition'],dict_input_tsfit['sheet_quantreg'],dict_input_tsfit['sheet_cond_quant']]
    df_tsfit_collections = read_data_tsfit(input_sheetnames)

    # return a dict for input parameters and a df
    return dict_input_tsfit, df_tsfit_collections

def read_parameters_tsfit():
    '''
    Read in parameters for tsfit.

    The cell positions for inputs are hardcoded in.
    Parameter value and range checking is done in the check_parameters_tsfit function.
    '''

    # Create dict for parameters
    dict_parameters_tsfit = dict()

    # ---------------------------#
    # Read in necessary values
    # ---------------------------#
    cellpos = 'F31'
    # Get a list of all values starting at cellpos going down
    dict_parameters_tsfit['quantlist'] = wb.sheets['Input_parameters'].range(cellpos).expand('down').value

    # Read in info on regressors.
    dict_parameters_tsfit['regressors'] = dict()
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
        dict_parameters_tsfit['regressors'][regressor+'_trans_'+str(iregressor)+'_'+transform] = dict()
        # Set as values of dict
        dict_parameters_tsfit['regressors'][regressor+'_trans_'+str(iregressor)+'_'+transform]['transform'] = transform
        dict_parameters_tsfit['regressors'][regressor+'_trans_'+str(iregressor)+'_'+transform]['option'] = option
    
    # Read in latest_date
    cellpos = 'B61'
    dict_parameters_tsfit['latest_date'] = wb.sheets['Input_parameters'].range(cellpos).value

    # Read in t-skew fit parameters
    dict_parameters_tsfit['fit_params'] = dict()
    cellpos = 'B59'
    dict_parameters_tsfit['fit_params']['fittype'] = wb.sheets['Input_parameters'].range(cellpos).value
    
    dict_parameters_tsfit['fit_params']['mode'] = dict()
    cellpos = 'B64'
    dict_parameters_tsfit['fit_params']['mode']['constraint'] = wb.sheets['Input_parameters'].range(cellpos).value
    cellpos = 'D64'
    dict_parameters_tsfit['fit_params']['mode']['value'] = wb.sheets['Input_parameters'].range(cellpos).value
    dict_parameters_tsfit['fit_params']['qsmooth'] = dict()
    cellpos = 'B67'
    dict_parameters_tsfit['fit_params']['qsmooth']['option'] = wb.sheets['Input_parameters'].range(cellpos).value
    cellpos = 'D67'
    dict_parameters_tsfit['fit_params']['qsmooth']['period'] = wb.sheets['Input_parameters'].range(cellpos).value
    
    
    #Read in plot parameters 
    obj1 = wb.sheets['Input_parameters'].api   
    for ind,c in enumerate(obj1.CheckBoxes()):        
        e=c.value
        print(e)
        if ind==0:
            dict_parameters_tsfit['fit_params']['plot_mode']= (e>0)
        elif ind==1:
            dict_parameters_tsfit['fit_params']['plot_median']= (e>0)
        elif ind==2:
            dict_parameters_tsfit['fit_params']['plot_mean']= (e>0)
   
    
   
    # Read in output sheet
    cellpos = 'B69'
    dict_parameters_tsfit['sheet_tsfit'] = wb.sheets['Input_parameters'].range(cellpos).value

    # Advanced t-skew fit parameters
    # Process each parameter in order
    pos = 74
    for param in ['dof',  'var_low', 'var_high', 'skew_low', 'skew_high']:
        dict_parameters_tsfit['fit_params'][param] = dict()
        cellpos = 'B' + str(pos)
        dict_parameters_tsfit['fit_params'][param]['constraint'] = wb.sheets['Input_parameters'].range(cellpos).value
        cellpos = 'D' + str(pos)
        dict_parameters_tsfit['fit_params'][param]['value'] = wb.sheets['Input_parameters'].range(cellpos).value
        pos += 1

    # The sheetname for the input is read in from what is in cell B52.
    # This is read in for do_quantfit() but since there is no way to connect
    # the output, we have to assume that the user has not changed this cell
    # since running the partitions.

    
    startrow = 51
    for isheetname, sheetname in enumerate(['sheet_quantreg', 'sheet_cond_quant']):
        colnum = startrow + isheetname
        cellpos = 'B' + str(colnum)
        dict_parameters_tsfit[sheetname] = wb.sheets['Input_parameters'].range(cellpos).value
        # checking of values will be done in check_parameters_quantfit

    # Read in partitions sheet
    cellpos = 'B24'
    dict_parameters_tsfit['sheet_partition'] = wb.sheets['Input_parameters'].range(cellpos).value
    dict_parameters_tsfit['partition_groups']= read_partition_groups()
    
    return dict_parameters_tsfit

def check_parameters_tsfit(dict_input_tsfit, keys):
    '''
    Check the input parameters for tsfit.
    '''

    # Check that all keys exist
    for key in keys:
        if key not in dict_input_tsfit:
            message = 'key ' + key + ' not found in dict_input_tsfit'
            show_message(message)

    input_sheets = ['Readme', 'Input_parameters', 'Partition_groups', 'Data', 'Processing_Log']
    for key in dict_input_tsfit:
        val = dict_input_tsfit[key]
        
        if key == 'latest_date':
            if type(val) != date:
                message = 'edate = ' + str(val) + ' was not a datetime.datetime object'
                show_message(message)
            # the range of the date is checked in Excel and is not checked here

        elif key.find('sheet_') != -1:
            if val is None:
                if key == 'sheet_quantreg':
                    dict_input_tsfit[key] = 'Quant reg coefficients'
                elif key == 'sheet_cond_quant':
                    dict_input_tsfit[key] = 'Conditional quantiles'
                elif key == 'sheet_local_proj':
                    dict_input_tsfit[key] = 'Local projections'
                elif key == 'sheet_partition':
                    dict_input_tsfit[key] = 'Output_partitions'
                elif key == 'sheet_tsfit':
                    if dict_input_tsfit['fit_params']['fittype']=='Asymmetric T':
                        dict_input_tsfit[key] = 'Asymmetric T fit'
                    else:
                        dict_input_tsfit[key] ='T-skew fit'
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
        
        elif key == 'regressors':
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
                    dict_input_tsfit['regressors'][regressor]['option'] = int(option)                                       
        
        elif key == 'fit_params':
            # This is a dict with keys for the variable name of each constraint
            # 'mode', 'dof', 'var', 'skewness', 'var_low', 'var_high', 'skew_low', 'skew_high'
            # and for each variable name there are 2 keys 'constraint' and 'value'
            dict_params = val
            
            for varname in dict_params:
                
                if varname in ['dof',  'var_low', 'var_high', 'skew_low', 'skew_high']:
                    constraint = dict_params[varname]['constraint']
                    value      = dict_params[varname]['value']

                # 'mode' will always be specified
                    if varname == 'mode':
                    
                    # Check that value is from pulldown
                        if constraint not in ['Free', 'Fixed', 'Median','Mean']:
                            message = 'constraint for ' + varname + ' was ' + constraint + ', not in pulldown values'
                            show_message(message, halt=True)
                    else:
                    # For all other varnames, only options are 'Fixed' and 'Free'

                    # Check that value is from pulldown
                        if constraint not in ['Free', 'Fixed', 'Default']:
                            message = 'constraint for ' + varname + ' was ' + constraint + ', not in pulldown values'
                            show_message(message, halt=True)

                # Check that no value is given when Free
                    if constraint in ['Free', 'Median'] and value is not None:
                        message = 'constraint for ' + varname + ' was ' + constraint + ' so value cannot be given as ' + str(value)
                        show_message(message, halt=True)

                # Check that value is a float
                    if constraint in ['Fixed'] and type(value) != float:
                        message = 'If constraint for ' + varname + ' is ' + constraint + ', value must be float, given as ' + str(value)
                        show_message(message, halt=True)
                        
                elif varname == 'qsmooth':
                    if dict_params[varname]['option']!='None' and dict_params[varname]['period'] is None:
                        message = 'Please provide period number for quantile smooth.'
                        show_message(message, halt=True)
                elif varname == 'fittype':
                    if dict_params[varname]!='Asymmetric T' and dict_params[varname]!='T-skew':
                        message = 'Not valid skewed T distribution option.'
                        show_message(message, halt=True)
                
                    

def read_data_tsfit(sheetnames):
    '''
    Read in the input data for tsfit.
    Checks for the sheetname should have been done in check_parameters_quantfit.
    '''
    data_collection={}
    for sheetname in sheetnames:
        print (sheetname)
        data_collection[sheetname] = wb.sheets[sheetname].range('A1').options(pd.DataFrame,index=False,expand='table').value
        if 'date' in data_collection[sheetname].columns:
            data_collection[sheetname].set_index(data_collection[sheetname]['date'],inplace=True)

    return data_collection

def run_tsfit(dict_input_tsfit, df_tsfit_collections, debug=False):  
    '''
    Main run function for step 3, tsfit.

    Takes in as arguments a dict for input parameters
    and a df for data. Outputs a dict for output parameters.

    Does quantile fits and returns a dict of output parameters.
    ** This function should be independent of any Excel input/output
    and be executable as a regular Python function independent of Excel. **
    '''
    
    warnings.filterwarnings("ignore")
    
    if debug:
        print('=' * 30)
        print('start of run_tsfit')
        print('=' * 30)
 
    # ------------------------
    # Create output dict
    # ------------------------
    dict_output_tsfit = dict()

    # ------------------------
    # Copy the output sheet names
    # from dict_input_tsfit
    # ------------------------
    for key in dict_input_tsfit:
        if key.find('sheet_') != -1:
            dict_output_tsfit[key] = dict_input_tsfit[key]
    
    # ------------------------
    # Get parameters from
    # dict_input_tsfit
    # ------------------------
    horizon = dict_input_tsfit['horizon']
    fitdate = dict_input_tsfit['latest_date']
    # TODO: get freq of data directly from data
    freq    = 'Quarterly'
    
    regressors=list(dict_input_tsfit['regressors'].keys())
    dict_output_tsfit['regressors']=dict_input_tsfit['regressors']
    fitparam=dict_input_tsfit['fit_params']
    # TODO: Don't do this here, drop the column so it is not redundant
    df_quantfit =df_tsfit_collections[dict_input_tsfit['sheet_partition']].set_index(df_tsfit_collections[dict_input_tsfit['sheet_partition']]['date'], drop=False) 
    for reg_long in  regressors:
        reg_short=reg_long.split('_trans_')[0]
        if dict_input_tsfit['regressors'][reg_long]['transform']=='None':
            df_quantfit[reg_long]=df_quantfit[reg_short]
        elif dict_input_tsfit['regressors'][reg_long]['transform']=='Lagged':
            df_quantfit[reg_long]=df_quantfit[reg_short].shift(dict_input_tsfit['regressors'][reg_long]['option'])
        elif dict_input_tsfit['regressors'][reg_long]['transform']=='MVA':
            df_quantfit[reg_long]=df_quantfit[reg_short].rolling(window=dict_input_tsfit['regressors'][reg_long]['option']).mean()
        elif dict_input_tsfit['regressors'][reg_long]['transform']=='Power': 
            df_quantfit[reg_long]=df_quantfit[reg_short]**(dict_input_tsfit['regressors'][reg_long]['option'])
        elif dict_input_tsfit['regressors'][reg_long]['transform']=='Diff': 
            df_quantfit[reg_long]=df_quantfit[reg_short].diff(dict_input_tsfit['regressors'][reg_long]['option'])
        elif dict_input_tsfit['regressors'][reg_long]['transform']=='ChangeRate': 
            df_quantfit[reg_long]=df_quantfit[reg_short].pct_change(dict_input_tsfit['regressors'][reg_long]['option'])
    # Fitdat
    if fitparam['qsmooth']['option']=='None':
        try:
            df_partition_fit=df_quantfit[df_quantfit['date']==fitdate]
        except:
            df_partition_fit=df_quantfit.iloc[-1,:]
            print ('The latest date in the data will be used.')
                           
    elif fitparam['qsmooth']['option']=='Median':
        per=int(fitparam['qsmooth']['period'])
        df_partition_fit=df_quantfit[df_quantfit['date']<=fitdate].drop(['date'],axis=1)       
        df_partition_fit=df_partition_fit.tail(per)
        df_partition_fit.loc['median']=df_partition_fit.median()
        df_partition_fit=df_partition_fit.loc[['median'],:]
        
    elif fitparam['qsmooth']['option']=='Mean':        
        per=int(fitparam['qsmooth']['period'])
        df_partition_fit=df_quantfit[df_quantfit['date']<=fitdate].drop(['date'],axis=1)       
        df_partition_fit=df_partition_fit.tail(per)
        #print(df_partition_fit)
        df_partition_fit.loc['mean']=df_partition_fit.mean()
        df_partition_fit=df_partition_fit.loc[['mean'],:]
        
    cond_quant={}
    df_quantcoef=df_tsfit_collections[dict_input_tsfit['sheet_quantreg']]
    for quant in dict_input_tsfit['quantlist']:
        cond_quant[quant]=df_quantcoef[(df_quantcoef['variable']=='Intercept') & (df_quantcoef['quantile']==quant)]['coeff_noscale'].values[0]        
        for reg_long in  regressors:        
            cond_quant[quant]+=df_quantcoef[(df_quantcoef['variable']==reg_long) & (df_quantcoef['quantile']==quant)]['coeff_noscale'].values[0]*df_partition_fit[reg_long].values[0]
    
    
    df_tsfit=df_tsfit_collections[dict_input_tsfit['sheet_cond_quant']]

    df_tsfit=df_tsfit[df_tsfit['date']==fitdate]
    olsmean=df_tsfit[df_tsfit['tau']=='mean']['conditional_quantile_mean'].values[0]
    res,cqlist,fig,fig2,dfpdf= gen_skewt(fitdate,fitparam,cond_quant,horizon,freq,olsmean)
    if fitparam['fittype']=='Asymmetric T':
        dfpdf=dfpdf[['AsymT_PDF_x','AsymT_CDF','AsymT_PDF_y']].round(decimals=5)
    elif fitparam['fittype']=='T-skew':
        dfpdf=dfpdf[['Tskew_PDF_x','Tskew_CDF','Tskew_PDF_y']].round(decimals=5)
    dict_output_tsfit['result'] = res
    dict_output_tsfit['data']   = cqlist
    dict_output_tsfit['fig']    = fig
    dict_output_tsfit['fig2']    = fig2
    dict_output_tsfit['dfpdf']    = dfpdf
    
    
    return dict_output_tsfit

def postrun_tsfit(dict_output_tsfit, debug=False):
    '''
    Postrun function for step 3, tsfit.

    Takes as input dict from main run function.

    This function cannot return any values due to limitations
    of the RunPython VBA code in xlwings.
    
    '''

    if debug:
        print('=' * 30)
        print('start of postrun_tsfit')
        print('=' * 30)

    # Create DataFrame for log
    log_frame = pd.DataFrame(columns=['Time','Action'])

    # Create the output sheets
    sheetname = dict_output_tsfit['sheet_tsfit']
    sheetnames = [sheet.name for sheet in wb.sheets]
    pd.options.display.float_format = '{:,.4f}'.format
    
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
        action = 'Unable to access sheet ' + sheetname

        # Add to log
    tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
    log = pd.Series({'Time': tn, 'Action': action})
    log_frame = log_frame.append(log, ignore_index=True)
        
    # end of loop over output sheetvars

    # Write out tsfit results
    cond_quant = dict_output_tsfit['data']
    df_pdf= dict_output_tsfit['dfpdf']


    res = dict_output_tsfit['result']
    sheetname = dict_output_tsfit['sheet_tsfit']
    try:
        wb.sheets[sheetname].pictures[0].delete()
    except:
        pass
    try:
        wb.sheets[sheetname].pictures[1].delete()
    except:
        pass
    try:
        wb.sheets[sheetname].range('A1').value = res      
        wb.sheets[sheetname].range('M1').value = df_pdf
        wb.sheets[sheetname].range('M:M').clear()
        wb.sheets[sheetname].range('R1').value = cond_quant
        action='T-skew fit saved succesfully.'
    except:
        action='Unable to output t-skew fit result.'
    tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
    log = pd.Series({'Time': tn, 'Action': action})
    log_frame=log_frame.append(log,ignore_index=True)

    # Write out figure
    sheetname = dict_output_tsfit['sheet_tsfit']
    sheet = wb.sheets[sheetname]
    fig = dict_output_tsfit['fig']
    fig2= dict_output_tsfit['fig2']
    # Set the path of the output file to be in the same dir as the
    # calling Excel file
    fullpath = os.path.abspath(os.path.dirname(wb.fullname) + '/figures')
    if not os.path.isdir(fullpath):
        os.makedirs(fullpath)
    outfilename = fullpath+'\\tskfit_'+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.png'
    try:
        fig.savefig(outfilename,dpi='figure')
    except:
        print('Fail to save t-skew fit figure.')

    try:
        partitionfiglist=glob(fullpath+'\\partition*.png')
        pimgfile=max(partitionfiglist)
        qautfiglist=glob(fullpath+'\\quantfit*.png')
        qimgfile=max(qautfiglist)
        print(pimgfile,qimgfile)
        pimgr=plt.imread(pimgfile)
        qimgr=plt.imread(qimgfile)
        pfig,pax= plt.subplots(1,1, figsize=(57,38))
        pax.imshow(pimgr)
        plt.axis('off')
        qfig,qax= plt.subplots(1,1, figsize=(35,14))
        qax.imshow(qimgr)
        plt.axis('off')
        print(pimgfile,qimgfile)
        pdfile=fullpath+'\\resultplots_'+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.pdf'
        print(pdfile)
        pp = PdfPages(pdfile)
        pp.savefig(pfig)
        pp.savefig(qfig)
        pp.savefig(fig)
        pp.close()
    except:
        print('Unable to save PDF plots.')
    try:
        sheet.pictures.add(fig, name='MyPlot', update=True, left=sheet.range('B14').left, top=sheet.range('B14').top, height=250, width=500)
        action = 'Skewed T distribution fit figure saved'
    except:
        action = 'Unable to add figure to sheet ' + sheetname
        
    try:
        sheet.pictures.add(fig2, name='MyPlot_2', update=True, left=sheet.range('B33').left, top=sheet.range('B33').top, height=250, width=500)
        action = 'Skewed T distribution fit figure saved'
    except:
        action = 'Unable to add figure to sheet ' + sheetname
    sheet.autofit()
    
    # Add to log
    tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
    log = pd.Series({'Time': tn, 'Action': action})
    log_frame=log_frame.append(log,ignore_index=True)

    # Write out log_frame
    add_logsheet(wb, log_frame, colnum=5)
