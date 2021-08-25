
'''
Simple utility script to test GAR from a non-Excel environment.

'''

import sys

import GAR

def main(argv):
    '''
    Main function that gets called when running run_GAR.py
    Usage:
    python run_GAR.py [input Excel file name]
    '''

    if len(argv) != 1:
        print('run_GAR.py usage:')
        print('python run_GAR.py [input Excel file name]')
        print('Must specify input Excel file name')
        sys.exit(-1)

    # Get the input file
    xlfile = argv[0]
    print(xlfile)
    # Set the mock caller
    GAR.set_mock_caller_file(xlfile)
    
    # Run the partition
    print('Running partition...')
    #GAR.do_partition()

    # Run the quantfit
    print('Running quantfit...')
    
    #GAR.do_quantfit()


    # Run the tsfit
    print('Running tsfit...')
    #GAR.do_tsfit()
    
    
    # Run the scenario
    print('Running scenario...')    
    #GAR.do_scenario()
    
    # Run the historical
    print('Running historical...')    
    #GAR.do_historical()
    
    # Run the segment
    print('Running segment...')    
    GAR.do_segment()
    
if __name__ == '__main__':
    #main(sys.argv[1:])
    main(['GaR.xlsm'])
