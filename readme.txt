* This version: 13 February 2019; Romain Lafarguette, rlafarguette@imf.org*

This folder contains the GaR excel tool developed at the International Monetary Fund. Specifically:

- gar.xlsx is the main Excel spreadsheet embedding the macros and the user-interface to run GaR estimation

- "run_GAR.py" is a Python file necessary to make the connection between the Excel spreadsheet and the GaR folder and should not be removed. 

- the folder /GAR contains the Python code, called by the Excel spreadsheet to run GaR. There is no need to know Python to use GaR, however Python should be installed on the computer

- the folder /figures contains the output in both .pdf and .png format of the different steps of the GaR analysis. The charts are also available in the excel spreadsheet

- the folder /Documentation contains the GaR project documentation, including:
	(i) HowToUseGaR.docx is an "how to file", explaining step by step how to use the Excel spreadhseet
	(ii) "GaR - Technical Appendix.docx" is a technical appendix providing the statistical background on GaR
	(iii) "Some examples of IMF GaR applications.docx" provides a list of references of IMF publications which have used the GaR excel tool on real country cases
