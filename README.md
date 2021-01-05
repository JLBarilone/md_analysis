# md_analysis
a suite of python and TCL scripts that can be used for trajectory analysis of molecular dynamic simulations


Contributing Authors: 
J. Prusa        email: prusa.jiri@gmail.com
J. Barilone

Authors of other codes that are implemented within some scripts are recognized within the respective script.


Python code was written to handle data from molecular dynamics simulations from polarizable force fields, but is capable of other functions as well.

#################### Relaxation process ####################
    # Example 1 - read one dipole moment evolution calculate it's autocorrelation function, plot it and save to text file
    # example_1()
    # Example 2 - read multiple correlation functions make mean and fit it by mono-exponential fit
    # example_2()
    # Example 3 - read multiple dipole moment evolution files, calculate correlation, calculate mean and calculate
    # the spectra as Fourier transform of correlation function
    # example_3()
    # Example 4 - Do all previous jobs (read multiple dipole moment evolution files, calculate correlation, calculate mean)
    # and make the fit on final correlation function
    # example_4()
    # Example 4.2 - Load mean autocorrelation function of all molecules in the box and fit it by biexponential function
    # example_4_2()
    # Example 4.3 - read multiple dipole moment evolution files, calculate correlation, calculate mean and fit
    # the correlation functioncalculate. Then calculate the spectra as a sum of peptide, water and crosscorelation part.
   # example_4_3()
    #################### Resonance process ####################
    # Example 5 - Calculate power spectral density from dipole moment evolution by single channel maximum entropy method (MEM)
    # Example_5()
    # Example 6 - Calculate power spectral density as Fourier transform of dipole moment time correlation function
    # Example_6()
