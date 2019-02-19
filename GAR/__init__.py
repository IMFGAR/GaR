
'''
Dynamic GaR (Growth at Risk) module for Python
==============================================

GAR is a Python module for running Growth at Risk calculations.
The user will use the provided Excel template to fill the data
and choose options, then run the code by clicking on Excel buttons.

For assistance contact
XXX YYY (XXX.imf.org)
'''

__version__ = "0.1.1"

import os
import sys
import imp

import xlwings as xw

# Initialize the wb to None
wb = None

# Set the wb from the caller if this library has been imported through Excel.
try:
    # If the GAR library is being imported by an Excel file through
    # xlwings, this will work and we can set up the wb
    wb = xw.Book.caller()
except:
    # If being imported from a regular Python run environment
    # this will be raised
    print('-' * 50)
    print('GAR WARNING:')
    print('    When importing the GAR library from a non-Excel')
    print('    environmentit is necessary to call')
    print("    GAR.set_mock_caller_file(xlfilename)")
    print('-' * 50)

from GAR.partition.partition_main import do_partition
from GAR.quantfit.quantfit_main import do_quantfit
from GAR.tsfit.tsfit_main import do_tsfit
from GAR.scenario.scenariomain import do_scenario
from GAR.historical.historicalmain import do_historical
from GAR.segment.segmentmain import do_segment
from GAR import historical
from GAR import scenario
from GAR import segment
def set_mock_caller_file(xlfilename='GaR.xlsm'):
    '''
    Function to set up mock caller for xlwings.
    This function needs to be called with the input Excel filename
    when trying to run the GaR package from the terminal.

    **When running from the input Excel this should not be called.**
    '''

    # Set the wb set in this function to be global
    # so it overwrites the wb throughout this file.
    global wb

    # Get the full path of the specified file
    xlfullpath = os.path.abspath(xlfilename)
    if not os.path.isfile(xlfullpath):
        print('File ' + xlfullpath + ' does not exist')
        sys.exit(-1)

    xw.Book(xlfullpath).set_mock_caller()
    wb = xw.Book.caller()
    print('Set file ' + xlfullpath + ' to be mock caller')

    # Once this is called we need to re-import the functions,
    # otherwise the variable wb imported is still the previous value
    imp.reload(globals)
    
    imp.reload(partition.partition_main)
    from GAR.partition.partition_main import do_partition
    
    imp.reload(quantfit.quantfit_main)
    from GAR.quantfit.quantfit_main import do_quantfit

    imp.reload(tsfit.tsfit_main)
    from GAR.tsfit.tsfit_main import do_tsfit  
        
    imp.reload(scenario.scenariomain)
    from GAR.scenario.scenariomain import do_scenario

    imp.reload(historical.historicalmain)
    from GAR.historical.historicalmain import do_historical
    
    imp.reload(segment.segmentmain)
    from GAR.segment.segmentmain import do_segment
