   
import os
from datetime import datetime as date
import time
import warnings # suppress warnings
warnings.filterwarnings("ignore")

## 3rd-party modules
import pandas as pd
import numpy as np

from GAR import wb
from GAR.globals import read_parameters_global, read_partition_groups, read_partition_groupsPLS, show_message, add_logsheet
from .plot_partition import partition_plot
from .partition_retro import partition_retro

###############################################################################
#%% Functions for step 1: partition
###############################################################################
def do_partition(debug=False):
    '''
    Entry point function called when button for partitions is called.

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
        print('start of do_partition')
        print('+' * 40)

    # Call prerun
    if debug:
        print('---- calling prerun_partition')
    dict_input_partition, dict_groups, df_partition = prerun_partition(debug=debug)

    # Call main run
    if debug:
        print('---- calling run_partition')
    dict_output_partition = run_partition(dict_input_partition, dict_groups, df_partition, debug=debug)

    # Call postrun
    if debug:
        print('---- calling postrun_partition')
    postrun_partition(dict_output_partition, debug=debug)

    # End measurement of time
    t1 = time.time()

    # Total time for this operation (formatted string)
    tdiff = "{:.1f}".format(t1 - t0)
    
    sheets = [dict_output_partition[key] for key in dict_output_partition if key.find('sheet') != -1 and key != 'sheet_input']
    message = 'Finished with partition in ' + tdiff + ' sec,\n'
    message += 'output is in sheets ' + ', '.join(sheets) 
    show_message(message,msgtype='info')
    
def prerun_partition(debug=False):
    '''
    Prerun function for step 1, partition.

    This function cannot take in any arguments due to limitations
    of the RunPython VBA code in xlwings.

    Read in/check the input parameters and return a
    dict for input parameters and a df for data.
    '''

    if debug:
        print('=' * 30)
        print('start of prerun_partition')
        print('=' * 30)

    # Keys for input parameter dict
    keys = ['freq', 'sdate', 'edate', 'method', 'pcutoff', 'method_growth', 'retropolate', 'sheet_partitions', 'sheet_loadings']

    # --------------------------
    # Read in parameters
    # --------------------------
    dict_input_partition = read_parameters_partition()

    # --------------------------
    # Check parameter values
    # --------------------------
    check_parameters_partition(dict_input_partition, keys)

    # --------------------------
    # Read in global parameters
    # --------------------------
    # (this also checks if values have changed since being initially set)
    dict_global_params = read_parameters_global()

    # Add each key, val from dict_global_params to dict_input_partition
    for key, val in dict_global_params.items():
        # Check that the keys do not clash
        if key in dict_input_partition:
            message = 'dict_input_partition should not have key ' + key + ' that is common with dict_global_params'
            show_message(message, halt=True)
        dict_input_partition[key] = val


    # --------------------------
    # Create df for data
    # --------------------------
    df_partition = read_data_partition()
    
    
    
    # --------------------------
    # Create a dict for groups
    # --------------------------

    if dict_input_partition['method']=="PLS":
        dict_groups, dict_PLS = read_partition_groupsPLS()        
        dict_input_partition['PLS_target']=dict_PLS
    else:
        dict_groups = read_partition_groups()
        dict_input_partition['PLS_target']=None
        

    # --------------------------
    # Check df for partition
    # --------------------------
    check_data_partition(df_partition, dict_input_partition)

    # --------------------------
    # Set start and end dates and
    # fill missing values for data df
    # --------------------------
    df_partition = format_data_partition(df_partition, dict_input_partition['sdate'], dict_input_partition['edate'])

    #---------------------------
    # Check PLS target coverage
    #---------------------------
    
    if dict_input_partition['method']=="PLS":      
        PLSvar=set()
        for g in dict_PLS.values():
            for e in g:
                PLSvar.add(e)
        PLSvar=list(PLSvar)
        if df_partition.loc[:,PLSvar].isnull().sum().sum()>0:
            print( df_partition.loc[:,PLSvar].head())
            message = 'PLS target should cover all dates'
            show_message(message, halt=True)

    # --------------------------
    # Check partition groups
    # --------------------------
    check_partition_groups(dict_groups, df_partition)
    
    # --------------------------
    # Return a dict for input parameters,
    # a dict for groups,
    # and a df for the data
    # --------------------------
    return dict_input_partition, dict_groups, df_partition

def read_parameters_partition():
    '''
    Read in parameters for partition.

    The cell positions for inputs are hardcoded in.
    Parameter value and range checking is done in the check_parameters_partition function.
    '''

    # Create dict for parameters
    dict_parameters_partition = dict()

    # ---------------------------#
    # Read in necessary values
    # ---------------------------#

    # This is the row number where the parameters start
    pos = 17

    # Process each parameter in order
    for param in ['freq', 'sdate', 'edate', 'method', 'pcutoff', 'method_growth', 'retropolate', 'sheet_partitions', 'sheet_loadings']:
        cellpos = 'B' + str(pos)
        dict_parameters_partition[param] = wb.sheets['Input_parameters'].range(cellpos).value
        pos += 1

    return dict_parameters_partition

def check_parameters_partition(dict_input_partition, keys):
    '''
    Check the input parameters for partition.
    '''

    # Check that all keys exist
    for key in keys:
        if key not in dict_input_partition:
            message = 'key ' + key + ' not found in dict_input_partition'
            show_message(message)

    # Get the sheets in the wb
    sheetnames = [sheet.name for sheet in wb.sheets]

    # These are the sheets that are used as input and never overwritten
    input_sheets = ['Readme', 'Input_parameters', 'Partition_groups', 'Data', 'Processing_Log']

    # Check for necessary sheets
    for sheetname in ['Partition_groups', 'Data']:
        if sheetname not in sheetnames:
            message = 'sheet with name ' + sheetname + ' must be in input Excel file'
            show_message(message)

    # ------------------------------------------- #
    # Go through each value and check the values
    # ------------------------------------------- #

    for key in keys:
        val = dict_input_partition[key]
        
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
            
        if key in ['sheet_partitions', 'sheet_loadings']:
            # If a value was specified, check that it is not one of the
            # input sheet names and use it as the output sheet name.
            # Otherwise we will use the default 'Output_partitions'
            if val is None:
                dict_input_partition[key] = key.replace('sheet', 'Output')
            else:
                # Check that the specified sheetname is not one of the inputs
                # (it is OK that it is the same name as an existing sheet if
                # that sheet is an output)
                if val in input_sheets:
                    message = key + ' specified as ' + val + ', cannot be the same as necessary input sheet'
                    show_message(message)

def read_data_partition():
    '''
    Read in the input data for partition.

    For partition, all data should be in the sheet called "Data".
    Another sheet, "Partition_groups" should 
    '''

    # Get the sheets in the wb
    sheetnames = [sheet.name for sheet in wb.sheets]
    # Make sure that Data sheet exists
    if 'Data' not in sheetnames:
        message = 'Sheet named Data does not exist'
        show_message(message, halt=True)

    # Read in the Data sheet as a df
    df_partition = wb.sheets['Data'].range('A1').options(pd.DataFrame,index=False,expand='table').value
    colset=set()
    dupset=set()
    for e in df_partition.columns:
        if e not in colset:
            colset.add(e)
        else:
            dupset.add(e)
    if len(dupset)>0:
        dlist=list(dupset)
        dstr=','.join(dlist)
        message = 'Duplicate variables '+dstr+' in datasheet, please check.'
        show_message(message, halt=True)
    # Set index to date
    df_partition.index=df_partition['date']
    df_partition.index.name=None

    # TODO: set index to PeriodIndex
        
    return df_partition

def check_data_partition(df_partition, dict_input_partition):
    '''
    Check that necessary columns in df_partition are available,
    and also check that all columns are numeric.
    '''
    
    # --------------------------
    # The parameter
    # dict_input_partition['target']
    # need to be columns in df_partition
    # --------------------------
    for key in ['target']:
        col = dict_input_partition[key]
        if col not in df_partition.columns:
            message = 'col ' + col + ' for key ' + key + ' not in columns of df_partition'
            show_message(message, halt=True)

    # --------------------------
    # If any of the columns in
    # df_partition are not numeric
    # raise an error since the user
    # may have included invalid text
    # --------------------------
    # Global check of col types
    coltypes = set(df_partition.dtypes.values)
    if np.dtype('O') in coltypes:
        # If we find an object type in the columns,
        # find which ones they are
        obj_cols = []
        for col in df_partition.columns:
            if df_partition[col].dtype == np.dtype('O'):
                obj_cols.append(col)
        message = 'The following columns were not numeric types\n'
        message += 'This may be due to the data containing characters\n'
        message += 'Please remove all characters and run again\n'
        message += ', '.join(obj_cols)
        show_message(message, halt=True)

def format_data_partition(df, startdate, enddate, interpolate_method='linear', ffill_limit=None, fill_warning=None, debug=False):
    '''
    Fill missing values and latest values for partition data.

    Method to fill holes in the data will be specified with interpolate_method
    (defaults to linear interpolation), and number of ffill will be limited by
    ffill_limit. If fill_warning is set to a number, a warning will be given
    when more ffills are done than that value.
    '''

    # Slice off dates before and after startdate, enddate
    original_length = len(df)
    
    # Make a copy of the original df
    df = df[startdate:enddate]

    if debug:
        if len(df) != original_length:
            print('Original Data had length of ' + str(original_length) + ' now is ' + str(len(df)))

    # Interpolate intermediate missing values
    df = interpolate_missing_values(wb, df)

    # ffill latest values where necessary
    df = ffill_values(wb, df)
    
    return df

def interpolate_missing_values(wb, df, debug=False):
    '''
    Interpolate missing in-between values.
    Input:
        df : DataFrame to be filled
        wb : workbook to be modified
    '''

    if debug:
        print('-' * 20 + ' start of interpolate_missing_values')
    
    # Make a copy so we can compare to the original
    _df = df.copy()
    if debug:
        print('Total of ' + str(len(df.columns)) + ' cols')
        print('_df before interpolate_missing_values:')
        print(_df)

    # Since df has dates before dict_input_partition['sdate'] sliced off,
    # we need an offset for the row number
    daterange = wb.sheets['Data'].range('A1').expand('down').value # this is a list
    # Find the index where _df.index[0] is
    try:
        offset = daterange.index(_df.index[0])
    except ValueError:
        message = 'Could not find ' + str(_df.index[0]) + ' in range of sheet Data starting at A1'
        show_message(message, halt=True)
        
    if debug:
        print('_df.index[0] = ' + str(_df.index[0]))
        print('offset = ' + str(offset))
        print('_df.index[:20]:')
        print(_df.index[:20])
        print('daterange[:20]:')
        print(daterange[:20])

    # range for all values in sheet Data
    range = wb.sheets['Data'].range('A1').expand()
    # For each column fill missing values but not latest missing values

    for icol, col in enumerate(_df.columns):
        s  = df[col] # original
        _s = _df[col] # copy that gets filled

        # Skip dtypes that are not floats.
        # This is because the original data includes columns for isocodes
        if df[col].dtype not in [float, np.float64, np.float32]:
            if debug:
                print('skipping col ' + col + ' due to dtype being ' + df[col].dtype.name)
            continue
        
        # Get last valid index value
        last_index = _s.last_valid_index()

        # Do interpolate up to last valid index so we interpolate all missing values
        _s[:last_index] = _s[:last_index].interpolate(method='linear')
        
        # We don't want the initial missing values, only values missing in between.
        # Use the first_valid_index to restrict the range
        first_index = s.first_valid_index()

        # Get all values that have nan as the difference
        # between the original and the interpolated copy
        filled = _s[first_index:last_index][(s[first_index:last_index] - _s[first_index:last_index]).isnull() == True]

        if debug:
            if len(filled) > 0:
                print('values that have been filled in interpolate:')
                print(filled)

        # Use the index and col name to fill the Excel sheet with the interpolated values
        for ind in filled.index:
            # Get the index location
            irow = _df.index.get_loc(ind)

            # Get the cell corresponding to the value we want to fill.
            # Note that depending on whether the date is a column or not, we need to add 1 to the column number,
            # and that the row has an offset determined at the beginning of this function
            cell = range[irow+offset, icol]
            if debug:
                print('Filling in row ' + str(cell.row) + ' column ' + str(cell.column))
                print('with value ' + str(_s.loc[ind]))

            # Fill the sheet with the interpolated values
            cell.value = _s.loc[ind]

            # Set the font color to blue
            # 3 is red
            # 4 is green
            # 5 is blue
            cell.api.Font.ColorIndex = 5
            # end of loop over index in filled
        # end of loop over index of filled
    # end of loop over columns in _df

    # Return the updated data
    return _df

def ffill_values(wb, df, method='ffill', limit=None, debug=False):
    '''
    Forward fill (ffill) missing latest values.
    Input:
        df    : DataFrame for data
        wb    : workbook to be modified
        limit : Limit on consecutive values to fill
    '''

    if debug:
        print('-' * 20 + ' start of ffill_values')

    # Make a copy so we can compare to the original
    _df = df.copy()
    if debug:
        print('_df before ffill:')
        print(_df)

    # Since df has dates before dict_input_partition['sdate'] sliced off,
    # we need an offset for the row number
    daterange = wb.sheets['Data'].range('A1').expand('down').value # this is a list
    # Find the index where _df.index[0] is
    try:
        offset = daterange.index(_df.index[0])
    except ValueError:
        message = 'Could not find ' + str(_df.index[0]) + ' in range of sheet Data starting at A1'
        show_message(message, halt=True)
        
    if debug:
        print('_df.index[0] = ' + str(_df.index[0]))
        print('offset = ' + str(offset))
        print('_df.index[:20]:')
        print(_df.index[:20])
        print('daterange[:20]:')
        print(daterange[:20])
        
    # range for all values in sheet Data
    range = wb.sheets['Data'].range('A1').expand()
    # For each column fill missing values but not latest missing values
    for icol, col in enumerate(_df):
        s  = df[col] # original
        _s = _df[col] # copy that gets filled

        # Skip dtypes that are not floats.
        # This is because the original data includes columns for isocodes
        if df[col].dtype not in [float, np.float64, np.float32]:
            if debug:
                print('skipping col ' + col + ' due to dtype being ' + df[col].dtype.name)
            continue
        
        # If doing ffill, we don't want the initial missing values, only the final missing values
        # Use the first_valid_index to restrict the range
        first_index = s.first_valid_index()
        if debug:
            print('first_index = ' + str(first_index))

        # If doing bfill, we don't want the final missing values, only the initial missing values
        # Use the last_valid_index to restrict the range
        last_index = s.last_valid_index()
        if debug:
            print('last_index = ' + str(last_index))
                
        # Do ffill
        if method == 'ffill':
            _s[first_index:] = _s[first_index:].fillna(method=method, limit=limit)
            
            # Get all values that have nan as the difference
            # between the original and the interpolated copy
            filled = _s[first_index:][(s[first_index:] - _s[first_index:]).isnull() == True]
        elif method == 'bfill':
            _s[:last_index] = _s[:last_index].fillna(method=method, limit=limit)
            
            # Get all values that have nan as the difference
            # between the original and the interpolated copy
            filled = _s[:last_index][(s[:last_index] - _s[:last_index]).isnull() == True]
        else:
            message = 'Function ffill_values cannot take in method ' + method
            show_message(message, halt=True)

        #if debug:
        if len(filled) > 0:
            print('values that have been filled in ffill:')
            print(filled)

        # Use the index and col name to fill the Excel sheet with the interpolated values
        # Get the index location
        for ind in filled.index:
            irow = _df.index.get_loc(ind)

            # Get the cell corresponding to the value we want to fill.
            # Note that depending on whether the date is a column or not, we need to add 1 to the column number,
            # and that the row has an offset determined at the beginning of this function
            cell = range[irow+offset, icol]
            # if debug:
            print('Filling in row ' + str(cell.row) + ' column ' + str(cell.column))
            print('with value ' + str(_s.loc[ind]))
            
            # Fill the sheet with the interpolated values
            cell.value = _s.loc[ind]
        
            # Set the font color to red
            # 3 is red
            # 4 is green
            # 5 is blue
            # 6 is yellow
            # 7 is magenta
            cell.api.Font.ColorIndex = 3
            # end of loop over index in filled
        # end of loop over index of filled
    # end of loop over columns in _df

    # Return the updated data
    return _df

def check_partition_groups(dict_groups, df_partition):
    '''
    Check that for the list of variables in each group the column
    exists in df_partition.
    '''

    # Loop over groups
    for group in dict_groups:
        # Loop over variables
        for varname in dict_groups[group]:
            if varname not in df_partition.columns:
                message = 'variable ' + varname + ' was specified for group ' + group + ' but does not exist in df_partition'
                show_message(message, halt=True)

def run_partition(dict_input_partition, dict_groups, df_partition, debug=False):
    '''
    Main run function for step 1, partition.

    Takes in as arguments a dict for input parameters
    and a df for data. Outputs a dict for output parameters.

    Does partitioning and returns a dict of output parameters.
    ** This function should be independent of any Excel input/output
    and be executable as a regular Python function independent of Excel. **
    '''

    if debug:
        print('=' * 30)
        print('start of run_partition')
        print('=' * 30)

        # Show input parameters
        print('dict_input_partition:')
        for key in dict_input_partition:
            print(key.ljust(30) + ':' + str(dict_input_partition[key]))
        print('dict_groups:')
        for key in dict_groups:
            print(key.ljust(30) + ':' + str(dict_groups[key]))
        print('df_partition:')
        print(df_partition)

    warnings.filterwarnings("ignore")

    # ------------------------
    # Create DataFrame for log
    # ------------------------
    log_frame=pd.DataFrame(columns=['Time','Action'])

    # ------------------------
    # Create output dict
    # ------------------------
    dict_output_partition = dict()

    # ------------------------
    # Copy the output sheet names
    # from dict_input_partition
    # ------------------------
    for key in dict_input_partition:
        if key.find('sheet_') != -1:
            dict_output_partition[key] = dict_input_partition[key]
            #print(key, dict_output_partition[key] , dict_input_partition[key])

    # ------------------------
    # Get parameters from
    # dict_input_partition
    # ------------------------
    sdate   = dict_input_partition['sdate']
    edate   = dict_input_partition['edate']
    horizon = dict_input_partition['horizon']
    tdep    = dict_input_partition['target']+'_hz_'+str(horizon)
    df_partition  = df_partition.set_index(df_partition['date'], drop=False)
    method=dict_input_partition['method']
    benchcutoff = dict_input_partition['pcutoff']
    rgdp =  dict_input_partition['target'] # column name for real GDP
    method_growth = dict_input_partition['method_growth']
    PLStarget=dict_input_partition['PLS_target']
    # ------------------------
    # Run the partition
    # ------------------------
    retroframe, retroload, logretro, exitcode = partition_retro(dall=df_partition, groups_dict=dict_groups, tdep=tdep, rgdp=rgdp, method_growth=method_growth, horizon=horizon, method=method, sdate=sdate, edate=edate, benchcutoff=benchcutoff,PLStarget=PLStarget)
    log_frame = log_frame.append(logretro,ignore_index=True)
    
    if exitcode==-1:
        message = 'In the given time period some groups are complete empty. No feasible partition can be made. Please adjust partition groups or start date'
        show_message(message, halt=True)

    # Add return values
    
    figs={}
    #print(list(dict_groups.keys()))
    
    figs=partition_plot(df_partition,retroframe,retroload,list(dict_groups.keys()),PLStarget,tdep,method)
    dict_output_partition['frame']   = retroframe
    dict_output_partition['loading'] = retroload
    dict_output_partition['log']     = logretro
    dict_output_partition['figs']     = figs
    dict_output_partition['groups']=list(dict_groups.keys())
    dict_output_partition['method'] = method
    return dict_output_partition

def postrun_partition(dict_output_partition, debug=False):
    '''
    Postrun function for step 1, partition.

    Takes as input dict from main run function.

    This function cannot return any values due to limitations
    of the RunPython VBA code in xlwings.
    
    '''

    if debug:
        print('=' * 30)
        print('start of run_partition')
        print('=' * 30)

    # Create DataFrame for log
    log_frame = pd.DataFrame(columns=['Time','Action'])

    # Create the output sheets
    
    for sheetvar in [key for key in dict_output_partition if key.find('sheet') != -1]:
        # Check that sheetvar exists as a key in dict_output_partition
        if sheetvar not in dict_output_partition:
            message = 'sheetvar ' + sheetvar + ' is not a key for dict_output_partition'
            show_message(message, halt=True)

        # Get the actual sheet name
        sheetname = dict_output_partition[sheetvar]

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
            action = 'Unable to access sheet ' + sheetname

        # Add to log
        tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
        log = pd.Series({'Time': tn, 'Action': action})
        log_frame = log_frame.append(log, ignore_index=True)
        
    # end of loop over output sheetvars

    # Write out partition and loadings
    try:
        sheetname = dict_output_partition['sheet_partitions']
        wb.sheets[sheetname].range('A1').options(index=False).value = dict_output_partition['frame']
        wb.sheets[sheetname].autofit()
        sheetname = dict_output_partition['sheet_loadings']
        wb.sheets[sheetname].range('A1').options(index=False).value = dict_output_partition['loading']
        wb.sheets[sheetname].autofit()
        action='Partitions and loadings saved succesfully.'
    except:
        action='Unable to output partitions and loadings.'
    
    sheetname = dict_output_partition['sheet_partitions']
    
    for p in wb.sheets[sheetname].shapes:
        try: 
            p.delete()
        except Exception as e: 
            print(e)
    sheet = wb.sheets[sheetname]
    if dict_output_partition['method']=='PLS':
        for i,fig in enumerate(dict_output_partition['figs']):
            fullpath = os.path.abspath(os.path.dirname(wb.fullname) + '/figures')
            if not os.path.isdir(fullpath):
                os.makedirs(fullpath)
            group=dict_output_partition['groups'][i]
            outfilename = fullpath+'\\partition_PLS_'+group+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.png'
            fig.savefig(outfilename)
            try:
                X=str(2+i*38)
                sheet.pictures.add(fig, name='MyPlot_P'+str(i+1), update=True, left=sheet.range('M'+X).left, top=sheet.range('M'+X).top, height=500, width=750)
                action = 'Partition figure saved'
            except:
                action = 'Unable to add figure to sheet ' + sheetname
    else:
        fig = dict_output_partition['figs'][0]
        # Set the path of the output file to be in the same dir as the
        # calling Excel file
        fullpath = os.path.abspath(os.path.dirname(wb.fullname) + '/figures')
        if not os.path.isdir(fullpath):
            os.makedirs(fullpath)
        outfilename = fullpath+'\\partition_'+date.now().strftime('%Y_%m-%d@%H_%M-%S')+'.png'
        fig.savefig(outfilename)
        cr=len(dict_output_partition['groups'])
        try:
            sheet.pictures.add(fig, name='MyPlot_P1', update=True, left=sheet.range('M2').left, top=sheet.range('M2').top, height=720, width=cr*255)
            action = 'Partition figure saved'
        except:
            action = 'Unable to add figure to sheet ' + sheetname
        tn = date.now().strftime('%Y-%m-%d %H:%M:%S')
        log = pd.Series({'Time': tn, 'Action': action})
        log_frame = log_frame.append(log,ignore_index=True)
        
        fig1 = dict_output_partition['figs'][1]
        if  dict_output_partition['method']=='PLS':
            ht=700
            wd=480
        else:
            ht=320
            wd=320
        try:
            sheet.pictures.add(fig1, name='MyPlot_P2', update=True, left=sheet.range('M54').left, top=sheet.range('M54').top, height=ht, width=wd)
            action = 'Partition figure saved'
        except:
            action = 'Unable to add figure to sheet ' + sheetname
    
    # Write out log_frame
    add_logsheet(wb, log_frame, colnum=1)
