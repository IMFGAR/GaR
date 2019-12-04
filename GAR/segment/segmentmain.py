
import os
from datetime import datetime as date
import time
import warnings # suppress warnings
warnings.filterwarnings("ignore")

## 3rd-party modules
import pandas as pd
import numpy as np


from GAR import wb
from GAR.globals import read_parameters_global, read_partition_groups, read_partition_groupsPLS , show_message, add_logsheet
from .partition_retro import partition_retro
from .condqgreg import condquant
from .gen_seg_skewt import gen_seg_skewt
from .tskewfit import tskew_fit
from .asymtfit import asymt_fit
from .plot_quantfit import coeff_plot
from .term_plot import termstruct_plot


###############################################################################
#%% Functions for step 4: segment test
###############################################################################
def do_segment(debug=False):
    '''
    Entry point function called when button for cenario fits is called.

    This function cannot take in any arguments and can not have
    any return values due to limitation of the RunPython VBA code
    in xlwings.
    '''
    t0 = time.time()
    
    # Make sure a wb exists
    if wb is None:
        print('segment_main: wb is None')
        print('This may be due to not calling set_mock_caller_file')
        print('and setting the caller Workbook')
        import sys
        sys.exit(-1)
    else:
        print (wb)


    dict_input_segment, df_collection_segment = prerun_segment(debug=debug)
    dict_output_segment = run_segment(dict_input_segment, df_collection_segment, debug=debug)
    postrun_segment(dict_output_segment, debug=debug)

    # End measurement of time
    t1 = time.time()

    # Total time for this operation (formatted string)
    tdiff = "{:.1f}".format(t1 - t0)
    
    sheetname = dict_input_segment['sheet_segment']
    message = 'Finished with multiple horizon projections in ' + tdiff + ' sec,\n'
    message += 'output is in sheets ' + ', '+sheetname
    show_message(message,msgtype='info')
    
def prerun_segment(debug=False):
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
    dict_input_segment = read_parameters_segment()
    
    # --------------------------
    # Check parameter values
    # --------------------------
    check_parameters_segment(dict_input_segment)

    # Read in global parameters
    # (this also checks if values have changed since being initially set)
    dict_global_params = read_parameters_global()
    #print(dict_input_segment)
    #print(dict_global_params)


    # Add each key, val from dict_global_params to dict_input_quantfit
    for key, val in dict_global_params.items():
        # Check that the keys do not clash
        if key in dict_input_segment:
            message = 'dict_input_quantfit should not have key ' + key + ' that is common with dict_global_params'
            show_message(message, halt=True)
        dict_input_segment[key] = val

                        
    # Create df for data
    df_segment_collections = read_data_segment()

    # return a dict for input parameters and a df
    return dict_input_segment,df_segment_collections

def read_parameters_segment():
    '''
    Read in parameters for quantfit.

    The cell positions for inputs are hardcoded in.
    Parameter value and range checking is done in the check_parameters_quantfit function.
    '''

    # Create dict for parameters
    dict_parameters_segment = dict()

    # ---------------------------#
    # Read in necessary values
    # ---------------------------#



    # Read in partition parament
    pos = 17

    # Process each parameter in order
    for param in ['freq', 'sdate', 'edate', 'method', 'pcutoff', 'method_growth', 'retropolate']:
        cellpos = 'B' + str(pos)
        dict_parameters_segment[param] = wb.sheets['Input_parameters'].range(cellpos).value
        pos += 1
        
    # Read in quantlist
    cellpos = 'F31'
    # Get a list of all values starting at cellpos going down
    dict_parameters_segment['quantlist'] = wb.sheets['Input_parameters'].range(cellpos).expand('down').value

    # Read in info on regressors.
    dict_parameters_segment['regressors'] = dict()
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
        dict_parameters_segment['regressors'][regressor+'_trans_'+str(iregressor)+'_'+transform] = dict()
        # Set as values of dict
        dict_parameters_segment['regressors'][regressor+'_trans_'+str(iregressor)+'_'+transform]['transform'] = transform
        dict_parameters_segment['regressors'][regressor+'_trans_'+str(iregressor)+'_'+transform]['option'] = option
    

    

    # Read in t-skew fit parameters
    dict_parameters_segment['fit_params'] = dict()
    cellpos = 'B59'
    dict_parameters_segment['fit_params']['fittype'] = wb.sheets['Input_parameters'].range(cellpos).value
    
    dict_parameters_segment['fit_params']['qsmooth'] = dict()
    cellpos = 'B67'
    dict_parameters_segment['fit_params']['qsmooth']['option'] = wb.sheets['Input_parameters'].range(cellpos).value
    cellpos = 'D67'
    dict_parameters_segment['fit_params']['qsmooth']['period'] = wb.sheets['Input_parameters'].range(cellpos).value
    
    # Advanced t-skew fit parameters
    # Process each parameter in order
    pos = 74
    for param in ['dof',  'var_low', 'var_high', 'skew_low', 'skew_high']:
        dict_parameters_segment['fit_params'][param] = dict()
        cellpos = 'B' + str(pos)
        dict_parameters_segment['fit_params'][param]['constraint'] = wb.sheets['Input_parameters'].range(cellpos).value
        cellpos = 'D' + str(pos)
        dict_parameters_segment['fit_params'][param]['value'] = wb.sheets['Input_parameters'].range(cellpos).value
        pos += 1
    
    # Read in segment parameters

    # Start with the first cell that contains the horizonlist.
    startrow = 120
    cellpos = 'A' + str(startrow)    
    # Read down and get a list of all horizons
    dict_parameters_segment['horizonlist'] = wb.sheets['Input_parameters'].range(cellpos).expand('down').value
    
    # Start with the first cell that contains the fitdatelist.
    cellpos = 'B' + str(startrow)    
    # Read down and get a list of all horizons
    dict_parameters_segment['fitdatelist'] = wb.sheets['Input_parameters'].range(cellpos).expand('down').value
    
    n_hzs=len(dict_parameters_segment['fitdatelist'])
    dict_parameters_segment['fitconstrainlist']=[]
    dict_parameters_segment['fitconstrainvalues']=[]
    # Start with the first cell that contains the mode.
    for i in range(n_hzs):
        cellpos = 'D' + str(startrow+i)    
    # Read down and get a list of all horizons
        dict_parameters_segment['fitconstrainlist'].append(wb.sheets['Input_parameters'].range(cellpos).value)
    
        # Start with the first cell of values of the mode constrains.
        cellpos = 'F' + str(startrow+i)    
        # Read down and get a list of all horizons
        dict_parameters_segment['fitconstrainvalues'].append(wb.sheets['Input_parameters'].range(cellpos).value)
    
    # Read in output sheet
    cellpos = 'B130'
    dict_parameters_segment['sheet_segment'] = wb.sheets['Input_parameters'].range(cellpos).value

    cellpos = 'B132'
    dict_parameters_segment['sheet_term'] = wb.sheets['Input_parameters'].range(cellpos).value

    dict_parameters_segment['partition_groups']= read_partition_groups()
    
    if dict_parameters_segment['method']=="PLS":
        dict_groups, dict_PLS = read_partition_groupsPLS()        
        dict_parameters_segment['PLS_target']=dict_PLS
        dict_parameters_segment['partition_groups']=dict_groups
    else:
        dict_groups = read_partition_groups()
        dict_parameters_segment['partition_groups']=dict_groups
        dict_parameters_segment['PLS_target']=None
    return dict_parameters_segment

def check_parameters_segment(dict_input_segment):
    '''
    Check the input parameters for segment.
    '''
    # Check that all keys exist
    keys=dict_input_segment.keys()    

    input_sheets = ['Readme', 'Input_parameters', 'Partition_groups', 'Data', 'Processing_Log']
    for key in keys:
        val = dict_input_segment[key]
        
        if key == 'freq' and val not in ['Monthly', 'Quarterly', 'Yearly']:
            message = 'freq = ' + val + ' was not a valid value'
            show_message(message)

        if key == 'sdate':
            if type(val) != date:
                message = 'sdate = ' + str(val) + ' was not a datetime.datetime object'
                show_message(message)
            # the range of the date is checked in Excel and is not checked here

        if key == 'edate':
            if type(val) != date:
                message = 'edate = ' + str(val) + ' was not a datetime.datetime object'
                show_message(message)
            # the range of the date is checked in Excel and is not checked here
            
        if key == 'method' and val not in ['LDA', 'PCA','PLS']:
            message = 'method = ' + val + ' was not a valid value'
            show_message(message)

        if key == 'pcutoff' and not (0 < val and val < 1):
            print('pcutoff = ' + str(val))
            message = 'pcutoff = ' + str(val) + ' was not a valid sheet name'
            show_message(message)

        if key == 'real_GDP':
            message = 'benchmark = ' + val + ' needs to be checked'
            show_message(message, output_messagebox=False)

        if key == 'method_growth':
            if val not in ['cpd', 'yoy','level']:
                message = 'method_growth = ' + val + ' must be one of cpd/yoy'
                show_message(message, output_messagebox=False)
            
        if key == 'retropolate' and val not in ['Yes', 'No']:
            message = 'retropolate = ' + val + ' was not a valid value'
            show_message(message)
            
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
                    dict_input_segment['regressors'][regressor]['option'] = int(option)
                    
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
                    dict_input_segment[key] = 'Quant reg coefficients'
                elif key == 'sheet_cond_quant':
                    dict_input_segment[key] = 'Conditional quantiles'
                elif key == 'sheet_local_proj':
                    dict_input_segment[key] = 'Local projections'
                elif key == 'sheet_partition':
                    dict_input_segment[key] = 'Output_partitions'
                elif key == 'sheet_segment':
                    dict_input_segment[key] = 'Multiple_projections'    
                elif key == 'sheet_term':
                    dict_input_segment[key] = 'Term_Structure'    
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
                    
def read_data_segment():
    '''
    Read in the input data for segment.
    Checks for the sheetname should have been done in check_parameters_quantfit.
    '''       
    data_collection={}
    data_collection['Data'] = wb.sheets['Data'].range('A1').options(pd.DataFrame,index=False,expand='table').value     
    data_collection['Data'].set_index('date', inplace=True, drop=False)
    data_collection['Data'].fillna(method='ffill',inplace=True)
    return data_collection

def run_segment(dict_input_segment, df_collection_segment, debug=False):
    

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
    dict_output_segment = dict()
    
    # ------------------------
    # Copy the output sheet names
    # from dict_input_tsfit
    # ------------------------
    for key in dict_input_segment:
        if key.find('sheet_') != -1:
            dict_output_segment[key] = dict_input_segment[key]
    
    sdate   = dict_input_segment['sdate']
    edate   = dict_input_segment['edate']

    tdep    = dict_input_segment['target']+'_hz_'
    df_partition = df_collection_segment['Data']
    df_partition = df_partition[sdate:edate] 
    method=dict_input_segment['method']
    benchcutoff = dict_input_segment['pcutoff']
    rgdp =  dict_input_segment['target'] # column name for real GDP
    method_growth = dict_input_segment['method_growth']
    dict_groups=dict_input_segment['partition_groups']
    PLStarget=dict_input_segment['PLS_target']
    regressors=list(dict_input_segment['regressors'].keys())
    fitparam=dict_input_segment['fit_params']
    # TODO: get freq of data directly from data
    freq    = 'Quarterly'
    skewtlist = []
    horizonlist = []
    medianlist = []
    loclist=[]
    inputdatelist=[]
    qlist=dict_input_segment['quantlist']
    qlist.sort()
    cqs=[qlist]

    hset=set()
    # ------------------------
    # Run the partition
    # ------------------------
    for indh,hz in enumerate (dict_input_segment['horizonlist']):
        horizon=int(hz)
        horizonlist.append(horizon)
        fitdate=dict_input_segment['fitdatelist'][indh]
        inputdate=fitdate
        inputdatelist.append(inputdate)
        df_quantfit, retroload, logretro, exitcode = partition_retro(dall=df_partition, groups_dict=dict_groups, tdep=tdep+str(horizon), rgdp=rgdp, method_growth=method_growth, horizon=horizon, method=method, sdate=sdate, edate=edate, benchcutoff=benchcutoff,PLStarget=PLStarget)            

        df_quantfit = df_quantfit.set_index(df_quantfit['date'], drop=False) 
        for reg_long in  regressors:
            reg_short=reg_long.split('_trans_')[0]
            if dict_input_segment['regressors'][reg_long]['transform']=='None':
                df_quantfit[reg_long]=df_quantfit[reg_short]
            elif dict_input_segment['regressors'][reg_long]['transform']=='Lagged':
                df_quantfit[reg_long]=df_quantfit[reg_short].shift(dict_input_segment['regressors'][reg_long]['option'])
            elif dict_input_segment['regressors'][reg_long]['transform']=='MVA':
                df_quantfit[reg_long]=df_quantfit[reg_short].rolling(window=dict_input_segment['regressors'][reg_long]['option']).mean()
            elif dict_input_segment['regressors'][reg_long]['transform']=='Power': 
                df_quantfit[reg_long]=df_quantfit[reg_short]**(dict_input_segment['regressors'][reg_long]['option'])
            elif dict_input_segment['regressors'][reg_long]['transform']=='Diff': 
                df_quantfit[reg_long]=df_quantfit[reg_short].diff(dict_input_segment['regressors'][reg_long]['option'])
            elif dict_input_segment['regressors'][reg_long]['transform']=='ChangeRate': 
                df_quantfit[reg_long]=df_quantfit[reg_short].pct_change(dict_input_segment['regressors'][reg_long]['option'])

        #df_quantcoef, dcond_quantiles_all, loco_all, exitcode = condquant(df_quantfit, tdep+str(horizon), regressors, horizon,dict_input_segment['quantlist'])
        df_quantcoef, exitcode = condquant(df_quantfit, tdep+str(horizon), regressors, horizon,dict_input_segment['quantlist'])
        
        
        if indh==0:
            df_coefs=pd.DataFrame(index=df_quantcoef.index)
            df_coefs['quantile']=df_quantcoef['quantile']
        df_coefs['coeff_scale_PROJ'+str(indh+1)]=df_quantcoef['coeff_scale']
        #df_coefs['R2_in_sample_PROJ'+str(indh+1)]=df_quantcoef['R2_in_sample']
        #df_coefs['errors_PROJ'+str(indh+1)]=(df_quantcoef['upper']-df_quantcoef['lower'])/2
        if len(hset)==0:
            df_term=pd.DataFrame(index=df_quantcoef.index)
            df_term['quantile']=df_quantcoef['quantile']
        if horizon not in hset:
            hset.add(horizon)
            df_term['coeff_scale_hz'+str(horizon)]=df_quantcoef['coeff_scale']
            df_term['upper_hz'+str(horizon)]=df_quantcoef['upper']
            df_term['lower_hz'+str(horizon)]=df_quantcoef['lower']
            df_term['error_hz'+str(horizon)]=(df_quantcoef['upper']-df_quantcoef['lower'])/2
            df_term['R2_in_sample_hz'+str(horizon)]=df_quantcoef['R2_in_sample']
            
            
        if fitparam['qsmooth']['option']=='None':
            try:
                df_partition_fit=df_quantfit[df_quantfit['date']==inputdate]
            except:
                df_partition_fit=df_quantfit.iloc[-1,:]
                print ('The latest date in the data will be used.')
                           
        elif fitparam['qsmooth']['option']=='Median':
            per=int(fitparam['qsmooth']['period'])
            df_partition_fit=df_quantfit[df_quantfit['date']<=inputdate].drop(['date'],axis=1)       
            df_partition_fit=df_partition_fit.tail(per)
            df_partition_fit.loc['median']=df_partition_fit.median()
            df_partition_fit=df_partition_fit.loc[['median'],:]
        
        elif fitparam['qsmooth']['option']=='Mean':        
            per=int(fitparam['qsmooth']['period'])
            df_partition_fit=df_quantfit[df_quantfit['date']<=inputdate].drop(['date'],axis=1)       
            df_partition_fit=df_partition_fit.tail(per)
        #print(df_partition_fit)
            df_partition_fit.loc['mean']=df_partition_fit.mean()
            df_partition_fit=df_partition_fit.loc[['mean'],:]
            
        elif fitparam['qsmooth']['option']=='Median within horizon':
            per=horizon
            df_partition_fit=df_quantfit[df_quantfit['date']<=inputdate].drop(['date'],axis=1)       
            df_partition_fit=df_partition_fit.tail(per)
            df_partition_fit.loc['median']=df_partition_fit.median()
            df_partition_fit=df_partition_fit.loc[['median'],:]    
        
        elif fitparam['qsmooth']['option']=='Mean within horizon':        
            per=horizon
            df_partition_fit=df_quantfit[df_quantfit['date']<=inputdate].drop(['date'],axis=1)       
            df_partition_fit=df_partition_fit.tail(per)
        #print(df_partition_fit)
            df_partition_fit.loc['mean']=df_partition_fit.mean()
            df_partition_fit=df_partition_fit.loc[['mean'],:]
        
        
        
        cond_quant={}
        cqlist=[]
        for quant in dict_input_segment['quantlist']:
            cond_quant[quant]=df_quantcoef[(df_quantcoef['variable']=='Intercept') & (df_quantcoef['quantile']==quant)]['coeff_noscale'].values[0]        
            for reg_long in  regressors:        
                cond_quant[quant]+=df_quantcoef[(df_quantcoef['variable']==reg_long) & (df_quantcoef['quantile']==quant)]['coeff_noscale'].values[0]*df_partition_fit[reg_long].values[0]
        for e in qlist:
            cqlist.append(cond_quant[e])
       
        cqs.append(cqlist)
        fitparam['mode']={}
        fitparam['mode']['constraint']=dict_input_segment['fitconstrainlist'][indh]
        fitparam['mode']['value']=dict_input_segment['fitconstrainvalues'][indh]
        
        
        olsmean=df_quantcoef[(df_quantcoef['variable']=='Intercept') & (df_quantcoef['quantile']=='mean')]['coeff_noscale'].values[0]
        for reg_long in  regressors:        
            olsmean+=df_quantcoef[(df_quantcoef['variable']==reg_long) & (df_quantcoef['quantile']=='mean')]['coeff_noscale'].values[0]*df_partition_fit[reg_long].values[0]

        if fitparam['fittype']=='T-skew':
            tsfit=tskew_fit(cond_quant,fitparam)
            skewtlist.append(tsfit)
            if fitparam['mode']['constraint']=='Fixed':
                loc=fitparam['mode']['value']
            else:
                loc=cond_quant[0.5]

            if fitparam['mode']['constraint']=='Free':
                loc=tsfit['loc']/tsfit['scale']
                
        elif fitparam['fittype']=='Asymmetric T':
            asymtfit=asymt_fit(cond_quant,fitparam,olsmean)
            skewtlist.append(asymtfit)
            loc=asymtfit['loc']
        medianlist.append(cond_quant[0.5])
        loclist.append(loc)
    
    df_cqs=pd.DataFrame(cqs)
    df_cqs=df_cqs.transpose()
    cqsheader=['Conditional quantiles :']
    hheader=['Projection '+str(i) for i in range(1,len(horizonlist)+1)]
    cqsheader.extend(hheader)
    df_cqs.columns=cqsheader

    res,fig,fig2,dfpdf= gen_seg_skewt(dict_input_segment['fitdatelist'],fitparam,skewtlist,horizonlist,medianlist,loclist,inputdatelist,freq)
    
    nhz=len(dict_input_segment['horizonlist'])
    fig3=coeff_plot(df_coefs,regressors, qlist,nhz)
       
    hlist=list(hset)
    hlist.sort()
    
    termfigs=termstruct_plot(df_term,regressors, qlist,hlist)
    
    dict_output_segment['fig'] = fig
    dict_output_segment['fig2'] = fig2
    dict_output_segment['fig3'] = fig3
    dict_output_segment['termfigs'] = termfigs
    dict_output_segment['cqs'] = df_cqs
    dict_output_segment['res'] = res
    dict_output_segment['dfpdf'] = dfpdf
    dict_output_segment['dfterm'] = df_term
    dict_output_segment['nhz']=len(hlist)
    dict_output_segment['nrg']=len(regressors)
    #print(df_coefs)
   
    return dict_output_segment

    
def postrun_segment(dict_output_segment, debug=False):
    '''
    Postrun function for step 2, segment.

    Takes as input dict from main run function.

    This function cannot return any values due to limitations
    of the RunPython VBA code in xlwings.
    
    '''

    if debug:
        print('=' * 30)
        print('start of postrun_segment')
        print('=' * 30)

    # Create DataFrame for log
    log_frame = pd.DataFrame(columns=['Time','Action'])
    
    sheetname = dict_output_segment['sheet_segment']
    
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
         
    for p in wb.sheets[sheetname].shapes:
        try: 
            p.delete()
        except Exception as e: 
            print(e)
    
    tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
    
    sheet = wb.sheets[sheetname]    

    pd.options.display.float_format = '{:,.4f}'.format
    
    df_pdf= dict_output_segment['dfpdf'].round(decimals=5)
    df_cqs= dict_output_segment['cqs'].round(decimals=5)
    df_term= dict_output_segment['dfterm'].round(decimals=5)
    nhz=dict_output_segment['nhz']
    nrg=dict_output_segment['nrg']
    # Write out segment resul
    res = dict_output_segment['res']
    wb.sheets[sheetname].range('A1').value = res 
    wb.sheets[sheetname].range('O1').options(index=False).value =df_cqs

    wb.sheets[sheetname].range('N'+str(len(df_cqs)+4)).options(index=False).value =df_pdf

    fig = dict_output_segment['fig']

    fullpath = os.path.abspath(os.path.dirname(wb.fullname) + '/figures')
    if not os.path.isdir(fullpath):
        os.makedirs(fullpath)

    outfilename = fullpath+'\\segment_pdf'+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.png'
    try:
        fig.savefig(outfilename)
    except:
         print('Fail to save segment figure.')

    try:
        sheet.pictures.add(fig, name='MyPlot', update=True, left=sheet.range('B12').left, top=sheet.range('B12').top, height=250, width=500)
        action = 'segment figure saved'
    except:
        action = 'Unable to add figure to sheet ' + sheetname
        
    fig2 = dict_output_segment['fig2']

    outfilename = fullpath+'\\segment_cdf'+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.png'
    try:
        fig2.savefig(outfilename)
    except:
         print('Fail to save segment figure.')
    
    try:
        sheet.pictures.add(fig2, name='MyPlot_2', update=True, left=sheet.range('B32').left, top=sheet.range('B32').top, height=250, width=500)
        action = 'segment figure saved'
    except:
        action = 'Unable to add figure to sheet ' + sheetname


    fig3 = dict_output_segment['fig3']

    outfilename = fullpath+'\\segment_qcoef'+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.png'
    try:
        fig3.savefig(outfilename)
    except:
         print('Fail to save segment figure.')
    
    try:
        sheet.pictures.add(fig3, name='MyPlot_3', update=True, left=sheet.range('B52').left, top=sheet.range('B52').top, height=230*nrg, width=500)
        action = 'segment figure saved'
    except:
        action = 'Unable to add figure to sheet ' + sheetname
        
        
        
    wb.sheets[sheetname].autofit()
    # Add to log
    tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
    log = pd.Series({'Time': tn, 'Action': action})
    log_frame = log_frame.append(log, ignore_index=True)
    
    sheetname = dict_output_segment['sheet_term']
    
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
         
    for p in wb.sheets[sheetname].shapes:
        try: 
            p.delete()
        except Exception as e: 
            print(e)
            
    wb.sheets[sheetname].range('T1').options(index=True).value =df_term
    
    
    
    sheet = wb.sheets[sheetname]    
    wb.sheets[sheetname].autofit()
    termfigs=dict_output_segment['termfigs']
    rs=(nhz-1)//4+1
    
    cs=min(4,nhz)
    for i,fig in enumerate(termfigs):
            

        outfilename = fullpath+'\\termstruct_'+str(i+1)+'_'+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.png'
        try:
            fig.savefig(outfilename)
        except:
            print('Fail to save term structure figure.')
    
        try:
            sheet.pictures.add(fig, name='Termplot_'+str(i), update=True, left=sheet.range('B3').left, top=sheet.range('B'+str(3+i*rs*21)).top, height=300*rs-40, width=200*cs)
            action = 'Term structure figure saved'
        except:
            action = 'Unable to add figure to sheet ' + sheetname
            
        tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
        log = pd.Series({'Time': tn, 'Action': action})
        log_frame = log_frame.append(log, ignore_index=True)
        # Write out log_frame
    add_logsheet(wb, log_frame, colnum=3)
    
    
    
