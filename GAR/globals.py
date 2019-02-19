# -*- coding: utf-8 -*-
"""

File containing global constants and functions.
"""

import datetime

import win32api # used in show_message
import win32con
import sys
import pandas as pd

from GAR import wb

def show_message(message, output_console=False, output_messagebox=True, halt=False, msgtype='error'):
    '''
    Take in an input message to display and if necessary terminate execution.

    Useful function to use for debugging when output_console is True,
    for actual warnings to Excel users set output_messagebox = True
    so that a warning message pops up in Excel.

    If output_console is set to True, the message is printed to a console.
    If output_messagebox is True, the message is printed in an Excel dialog box.
    The two are not mutually exclusive.

    If halt is set to True, the program is halted.
    Otherwise the message is shown and no other action is taken.
    '''

    # Show popup in Excel
    if output_messagebox:
        if msgtype=='error':
            win32api.MessageBox(wb.app.hwnd, message,'Error',win32con.MB_ICONERROR)
        elif msgtype=='info':
            win32api.MessageBox(wb.app.hwnd, message,'Info',win32con.MB_ICONINFORMATION)
        elif msgtype=='warning':
            win32api.MessageBox(wb.app.hwnd, message,'Warning',win32con.MB_ICONWARNING)
        else:
            win32api.MessageBox(wb.app.hwnd, message,'Info',win32con.MB_ICONINFORMATION)
                        
    # Print message to console
    if output_console:
        print(message)

    if halt:
        sys.exit(-1)

       
dict_global_params = dict()
# Initialize the dict for global parameters
global_params = ['target', 'horizon']
for param in global_params:
    dict_global_params[param] = None
    
def read_parameters_global():
    '''
    Read in parameters that are common to all operations.
    These should come from the the of the Input_parameters sheet
    and should be included in the global variable dict_global_params.

    This function should be read in at the beginning of each operation.
    If the parameters have changed, warn the users.

    Also checks that global parameters have been given a value that is not None.

    This is read in by each operation so is included here.
    '''

    # Use the dict_global_params defined above this function
    global dict_global_params

    # There should not be too many global parameters
    for cellpos, varname in [('B11', 'target'), ('B12', 'horizon')]:
        val = wb.sheets['Input_parameters'].range(cellpos).value

        # For horizon make sure that it is an int.
        # Since Excel numbers are read in as floats, convert to int
        # and check for difference with original value
        if varname == 'horizon':
            if type(val) != float or abs(int(val) - val) > 1E-5:
                message = 'global parameter horizon ' + str(val) + ' must be int'
                show_message(message, halt=True)

            # Since the value is less than 1E-5 away from an int,
            # convert to int so that there are no problems later
            val = int(val)

        # If no value is specified, tell user to fill it
        if val is None:
            message = 'Cell ' + cellpos + ' must have valid value for variable ' + varname
            show_message(message, halt=True)
    
        # If dict_global_params has not been set, set the value
        if dict_global_params[varname] is None:
            dict_global_params[varname] = val
        # Otherwise raise a warning to the user that the value has changed
        # and all operations should be run in order.
        elif val == dict_global_params[varname]:
            pass
        else:
            message = 'global parameter target from ' + cellpos + ' already has value of ' + dict_global_params[varname] + '\n'
            message += 'make sure that all operations are re-run in order to avoid inconsistencies'
            show_message(message)
            # Set the parameter to the specified value
            dict_global_params[varname] = val
    # end of reading in parameters from Excel

    # Check that all global_params have been set to a non-None value
    for key in global_params:
        if key not in dict_global_params:
            message = 'key of ' + key + ' is in global_params but not in dict_global_params'
            show_message(message, halt=True)

    return dict_global_params
    
def read_partition_groups():
    '''
    Read in groups for partition.
    This should come from a sheet called Partition_groups.
    '''

    # Get the sheets in the wb
    sheetnames = [sheet.name for sheet in wb.sheets]
    sheetname = 'Partition_groups'
    # Make sure that Partition_groups sheet exists
    if sheetname not in sheetnames:
        message = 'Sheet named ' + sheetname + ' does not exist'
        show_message(message, halt=True)

    # Read in the contents
    df_groups = wb.sheets[sheetname].range('B1:B29').options(pd.DataFrame,expand='right').value
    # Create dict for output
    dict_groups = dict()

    # The columns of this df are the keys
    cols = list(df_groups.columns)

    for col in cols:
        # Select the column and remove all N/A values,
        # and return a list of values
        dict_groups[col] = list(df_groups[col].dropna(axis=0).values)

    # For each column add 
    return dict_groups

def add_logsheet(wb, log_frame, colnum=1, logsheetname = 'Processing_log'):
    '''
    Add the contents of a DataFrame log_frame to a sheet named logsheetname in wb.

    The parameter colnum specifies which column to start the output.
    This is 1-index based, so cell A1 is column 1, not 0.
    '''

    # Get sheets in wb
    sheetnames = [sheet.name for sheet in wb.sheets]
    
    try:
        # Clear the sheet for colnum to colnum+2 if it already exists
        if logsheetname in sheetnames:
            wb.sheets[logsheetname].range((1,colnum), (1,colnum+2)).expand('down').clear()
            action = 'Cleared contents of ' + logsheetname + ' starting at cols ' + str(colnum)
        # Otherwise add it
        else:
            wb.sheets.add(logsheetname, after="Input_parameters")
            wb.sheets[logsheetname].api.Tab.Colorindex = 38 # pinkish color
            action = 'Created sheet ' + logsheetname
    except:
        action = 'Unable to access sheet ' + logsheetname

    # Add to log
    tn=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    #log = pd.Series({'Time': tn, 'Action': action})
    wb.sheets[logsheetname].range((1,colnum)).options(index=False).value = log_frame
    wb.sheets[logsheetname].autofit('c')
