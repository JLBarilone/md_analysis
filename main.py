# This is a sample Python script.
"""
This code can optionally employ 'Spectrum: a Spectral Analysis Library in Python' written by Thomas Cokelaer.
Spectrum program can be found at : http://thomas-cokelaer.info/software/spectrum/html/contents.html
"""

from dipole_lib import importDipoleData, monoExp, numFLap, biExp
from MEM_burg import burg_AR
from fortran_dip import calc_corr, calc_two_side_corr, fl_transform
from spectrum import arma2psd
import numpy as np
from math import ceil
from numpy.fft import fft
from scipy.constants import k, epsilon_0, c
import matplotlib.pyplot as plt
import scipy.optimize
import pandas as pd

def example_1():
    # Define Constants
    conversion_factor = 3.33564 * 1e-30  # convert Debye to Cm
    T = 300  # thermodynamic temperature in Kelvin
    psPerFrame = 0.002  # ps between frames in dip file
    max_lag_ps = 3  # maximum evaluated point for correlation function in ps
    dipFileName = "../trajectories/protein_rot_trans_0.dat"
    #################
    # import dipole file
    dip_data, V = importDipoleData(dipFileName)
    dip_array = np.array(dip_data[['dip_x', 'dip_y', 'dip_z']])
    dip_array = np.transpose(dip_array)
    # calculate correlation function
    max_lag_frames = ceil(max_lag_ps / psPerFrame)
    corr_fun = calc_corr(dip_array, dip_array, max_lag_frames) # perform F-L transform
    # Convert to proper units
    constant = conversion_factor ** 2 / (V * 3 * k * T * epsilon_0)
    corr_fun = corr_fun * constant
    # plot autocorrelation
    x_ticks = np.linspace(0, max_lag_ps, len(corr_fun))
    plt.plot(x_ticks, corr_fun)
    plt.xlabel('Time (ps)')
    plt.ylabel('AC (A.U.)')
    plt.show()
    # Export new dataframes to csv files
    np.savetxt("autocorr.csv", np.vstack((x_ticks, corr_fun)).T, header='time(ps)  corr_fun', delimiter="\t", comments='')
    print('Autocorrelation function saved to autocorr.csv')


def example_2():
    num_of_traj = 10
    #################
    # import correlation function data calculate mean and plot it
    meanCorr = 0
    for i in range(num_of_traj):
        corrFileName = "autocorr_%i.csv" % i
        # import correlation function
        df = pd.read_csv(corrFileName, header=0, delimiter="\s+")
        corr_time = np.array(df['time(ps)'])
        corr_fun = np.array(df['corr_fun'])
        meanCorr += corr_fun
        plt.plot(corr_time, corr_fun, 'b', linewidth=1.0)
    meanCorr /= num_of_traj

    # plot autocorrelation
    plt.plot(corr_time, meanCorr, 'r', linewidth=2.0)
    plt.xlabel('Time (ps)')
    plt.ylabel('AC (A.U.)')
    plt.show()
    # Export new dataframes to csv files
    np.savetxt("mean_autocorr.csv", np.vstack((corr_time, corr_fun)).T, header='time(ps)  corr_fun', delimiter="\t")

def example_3():
    # Define Constants
    conversion_factor = 3.33564 * 1e-30  # convert Debye to Cm
    T = 300  # thermodynamic temperature in Kelvin
    psPerFrame = 0.5  # ps between frames in dip file
    max_lag_ps = 100  # maximum evaluated point for correlation function in ps
    num_of_traj = 10  # number of trajectories
    #################
    meanCorr = 0
    max_lag_frames = ceil(max_lag_ps / psPerFrame)
    corr_time = np.linspace(0, max_lag_ps, max_lag_frames )
    for i in range(num_of_traj):
        dipFileName = "../trajectories/relaxation/ala_10ns_p_%i.dat" % i
        # import dipole file
        dip_data, V = importDipoleData(dipFileName)
        dip_array = np.array(dip_data[['dip_x', 'dip_y', 'dip_z']])
        dip_array = np.transpose(dip_array)
        # calculate correlation function
        corr_fun = calc_corr(dip_array, dip_array, max_lag_frames)
        # Convert to proper units
        constant = conversion_factor ** 2 / (V * 3 * k * T * epsilon_0)
        corr_fun = corr_fun * constant
        meanCorr += corr_fun

    meanCorr /= num_of_traj
    diffCorr = -np.diff(meanCorr)
    ################ dfft of numpy ################
    denf = fft(diffCorr, 4096)
    freq = np.fft.fftfreq(len(denf), d=psPerFrame*1e-12)
    # strip the left (negative) side
    denf = denf[1:len(denf) // 2]
    freq = freq[1:len(freq) // 2]
    ##############################################
    ########### Our impementation fortan or python ########
    freq = np.logspace(8, 11.5, 1048)
    denf= numFLap(diffCorr, corr_time*1e-12, freq)
    #denf = fl_transform(diffCorr, corr_time*1e-12, freq)
    #######################################################
    plt.plot(freq, denf.real,'b')
    plt.plot(freq, denf.imag,'r')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Electric susceptibility')
    plt.legend(['real part', 'imaginary part'])
    plt.xscale('log')
    plt.show()

def example_4():
    # Define Constants
    conversion_factor = 3.33564 * 1e-30  # convert Debye to Cm
    T = 300  # thermodynamic temperature in Kelvin
    psPerFrame = 0.5  # ps between frames in dip file
    max_lag_ps = 500  # maximum evaluated point for correlation function in ps
    num_of_traj = 10  # number of trajectories
    #################
    meanCorr = 0
    max_lag_frames = ceil(max_lag_ps / psPerFrame)
    corr_time = np.linspace(0, max_lag_ps, max_lag_frames + 1)
    for i in range(num_of_traj):
        dipFileName = "../trajectories/relaxation/ala_10ns_p_%i.dat" % i
        # import dipole file
        dip_data, V = importDipoleData(dipFileName)
        dip_array = np.array(dip_data[['dip_x', 'dip_y', 'dip_z']])
        dip_array = np.transpose(dip_array)
        # calculate correlation function
        corr_fun = calc_corr(dip_array, dip_array, max_lag_frames)
        # Convert to proper units
        constant = conversion_factor ** 2 / (V * 3 * k * T * epsilon_0)
        corr_fun = corr_fun * constant
        meanCorr += corr_fun
        # output individual correlation functions to file
        out_fn = "autocorr_%i.csv" % i
        np.savetxt(out_fn, np.vstack((corr_time, corr_fun)).T, header='time(ps)  corr_fun', delimiter="\t", comments='')
    meanCorr /= num_of_traj
    # Export mean correlation function to csv files
    np.savetxt("mean_autocorr.csv", np.vstack((corr_time, corr_fun)).T, header='time(ps)  corr_fun', delimiter="\t", comments='')
    # Normalize correlation function
    A = meanCorr[0]
    norm_corr = meanCorr / A
    # Perform fit
    p0 = (10)  # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(monoExp, corr_time, norm_corr, p0)
    tau = params
    fig, axs = plt.subplots(2)
    axs[0].plot(corr_time, meanCorr, 'r', linewidth=2.0)
    axs[0].plot(corr_time, A*monoExp(corr_time, tau))
    axs[0].set_title('Dipole moment correlation function')
    axs[0].set_xlabel('Time (ps)')
    axs[0].legend(['corr.fun','monoexp fit'])
    print('Static susceptibility of sample equals to %.2f' % A)
    print('Relaxation time for Debye process is equal to %.2f ps' % tau)
    # Calculate the FT of monoexp decay function by analytical formula
    tau = tau*1e-12
    freq = np.logspace(8, 11.5, 1048)
    omega = 2*np.pi*freq
    susc = np.zeros(len(freq), dtype='complex_')
    for i in range(len(freq)):
        susc[i] = A - omega[i] * (A * tau**2 * omega[i] / (1 + tau**2 * omega[i]**2))
        susc[i] = susc[i] + 1j * (omega[i] * A * tau / (1 + tau**2 * omega[i]**2))
    # Plot the results
    axs[1].plot(freq, susc.real)
    axs[1].plot(freq, susc.imag)
    axs[1].set_xlabel('Angular frequency [Hz]')
    axs[1].set_ylabel('Susceptibility')
    axs[1].set_xscale('log')
    axs[1].legend(['real part', 'imaginary part'])
    plt.show()


def example_4_2():
    # Define Constants
    conversion_factor = 3.33564 * 1e-30  # convert Debye to Cm
    T = 300  # thermodynamic temperature in Kelvin
    psPerFrame = 0.5  # ps between frames in dip file
    max_lag_ps = 500  # maximum evaluated point for correlation function in ps
    num_of_traj = 10  # number of trajectories
    #################
    # import correlation function
    corrFileName = "mean_autocorr_all.csv"
    df = pd.read_csv(corrFileName, header=0, delimiter="\s+")
    corr_time = np.array(df['time(ps)'])
    corr_fun = np.array(df['corr_fun'])
    # Normalize correlation function
    A = corr_fun[0]
    norm_corr = corr_fun / A
    # Perform fit
    p0 = (0.8, 10, 10)  # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(biExp, corr_time, norm_corr, p0, bounds=(0, [1., 1000., 100.]))
    A1, tau1, tau2 = params
    fig, axs = plt.subplots(2)
    axs[0].plot(corr_time, corr_fun, 'r', linewidth=2.0)
    axs[0].plot(corr_time, A*biExp(corr_time, A1, tau1, tau2))
    axs[0].set_title('Dipole moment correlation function')
    axs[0].set_xlabel('Time (ps)')
    axs[0].legend(['corr.fun','biexp fit'])
    print('Total static susceptibility of sample equals to %.2f (%.2f and %.2f for first and second components respectively)' % (A, A1*A,(1-A1)*A))
    print('Relaxation time for first and second Debye process is equal to %.2f ps and %.2f ps respectively.' % (tau1, tau2))
    # Calculate the FT of biexp decay function by analytical formula
    tau1 = tau1 * 1e-12
    tau2 = tau2 * 1e-12
    A2 = (1 - A1) * A
    A1 = A1 * A
    freq = np.logspace(8, 11.5, 1048)
    omega = 2*np.pi*freq
    susc = np.zeros(len(freq), dtype='complex_')
    for i in range(len(freq)):
        susc[i] = A - omega[i] * (A1 * tau1**2 * omega[i] / (1 + tau1**2 * omega[i]**2) + A2 * tau2**2 * omega[i] / (1 + tau2**2 * omega[i]**2))
        susc[i] = susc[i] + 1j * omega[i] * (A1 * tau1 / (1 + tau1**2 * omega[i]**2) + A2 * tau2 / (1 + tau2**2 * omega[i]**2))

    # Plot the results
    axs[1].plot(freq, susc.real)
    axs[1].plot(freq, susc.imag)
    axs[1].set_xlabel('Angular frequency [Hz]')
    axs[1].set_ylabel('Susceptibility')
    axs[1].set_xscale('log')
    axs[1].legend(['real part', 'imaginary part'])
    plt.show()


def example_4_3():
    # Define Constants
    conversion_factor = 3.33564 * 1e-30  # convert Debye to Cm
    T = 300  # thermodynamic temperature in Kelvin
    psPerFrame = 0.5  # ps between frames in dip file
    max_lag_ps = 500  # maximum evaluated point for correlation function in ps
    num_of_traj = 10  # number of trajectories
    #################
    meanCorr = 0
    max_lag_frames = ceil(max_lag_ps / psPerFrame)
    corr_time = np.linspace(0, max_lag_ps, max_lag_frames + 1)

    f_name_list = ['ala_10ns_p_', 'ala_10ns_w_']
    nof_components = len(f_name_list)
    meanCorr = np.zeros((nof_components, nof_components, len(corr_time)))
    for i in range(num_of_traj):
        dip_part = []
        # import peptide and water part of dipole data
        for file_name in f_name_list:
            dipFileName = "../trajectories/relaxation/" + file_name + "%i.dat" % i
            # import dipole file
            dip_data, V = importDipoleData(dipFileName)
            dip_array = np.array(dip_data[['dip_x', 'dip_y', 'dip_z']])
            dip_array = np.transpose(dip_array)
            dip_part.append(dip_array)
        # calculate auto-correlation and cross-correlation functions
        for m in range(nof_components):
            for n in range(nof_components):
                corr_fun = calc_corr(dip_part[m], dip_part[n], max_lag_frames)
                # Convert to proper units
                constant = conversion_factor ** 2 / (V * 3 * k * T * epsilon_0)
                corr_fun = corr_fun * constant
                meanCorr[m, n] += corr_fun

    meanCorr /= num_of_traj
    freq = np.logspace(8, 11.5, 1048)
    omega = 2 * np.pi * freq
    total_susci = 0
    corr_fig, axs_corr = plt.subplots(nof_components**2)
    suscii_fig, axs_suscii = plt.subplots(2)
    axs_corr[0].set_title('Dipole moment correlation functions')
    axs_corr[nof_components**2 - 1].set_xlabel('Time [ps]')
    leg_text = []
    for m in range(nof_components):
        for n in range(nof_components):
            # Normalize correlation function
            A = meanCorr[m, n, 0]
            norm_corr = meanCorr[m, n] / A
            # Perform fit
            p0 = (10)  # start with values near those we expect
            params, cv = scipy.optimize.curve_fit(monoExp, corr_time, norm_corr, p0)
            tau = params
            print('Static susceptibility of sample equals to %.2f' % A)
            print('Relaxation time of Debye process is for process %i, %i equal to %.2f ps' % (m, n, tau))
            # plot individual correlation functions
            idx = (m * nof_components + n)
            axs_corr[idx].plot(corr_time, meanCorr[m, n], 'r', linewidth=2.0)
            axs_corr[idx].plot(corr_time, A*monoExp(corr_time, tau))
            axs_corr[idx].legend(['corr.fun component %i %i' % (m, n),'monoexp fit'])
            # Calculate the FT of monoexp decay function by analytical formula
            tau = tau*1e-12
            susc = np.zeros(len(freq), dtype='complex_')
            for i in range(len(freq)):
                susc[i] = A - omega[i] * (A * tau**2 * omega[i] / (1 + tau**2 * omega[i]**2))
                susc[i] = susc[i] + 1j * (omega[i] * A * tau / (1 + tau**2 * omega[i]**2))
            total_susci += susc
            # Plot resulting susceptibility
            color_code = np.random.rand(3,)
            axs_suscii[0].plot(freq, susc.real, c=color_code)
            axs_suscii[1].plot(freq, susc.imag, c=color_code)
            axs_suscii[1].set_xlabel('Angular frequency [Hz]')
            axs_suscii[0].set_ylabel('Real susceptibility')
            axs_suscii[1].set_ylabel('Imag. susceptibility')
            axs_suscii[0].set_xscale('log')
            axs_suscii[1].set_xscale('log')
            leg_text.append("component %i, %i" % (m, n))

    axs_suscii[0].plot(freq, total_susci.real)
    axs_suscii[1].plot(freq, total_susci.imag)
    leg_text.append("Total susceptibility")
    axs_suscii[1].legend(leg_text)
    plt.show()



def Example_5():
    timestep = 50 * 1e-15  # in sec
    nfft = 4096  # number of points to calculate FT for
    p = 8  # order of Burg
    num_of_traj = 10  # number of trajectories
    plot_style = 1  # 0 == Frequency, 1 == wavelength

    # Create a for loop to read in each csv file and append the results to a list
    data_list = []

    # import multiple CSV files
    for i in range(num_of_traj):
        dipFileName = "../trajectories/protein_rot_trans_%i.dat" % i
        # import dipole file
        dip_data, V = importDipoleData(dipFileName)
        dip_array = np.array(dip_data[['dip_x', 'dip_y', 'dip_z']])
        dip_array = np.transpose(dip_array)
        data_list.append(dip_array)

    dip_data = np.vstack(data_list)
    print("Available number of time series: %i with %i frames" % np.shape(dip_data))
    # uncomment next line to shrink the input data (Don't forget to change timestep if you do so)
    #dip_data = np.array(dip_data)[0:5, 0:5000:5]
    print("You calculate spectra from %i time series of %i length." % np.shape(dip_data))
    print("With timestep of %.4e sec between frames." % timestep)
    # Create an array of derivatives
    array_of_derivatives = np.diff(dip_data)
    # feed array into burg_AR class
    coeffs, sigma = burg_AR(p, array_of_derivatives)
    print(coeffs)

    ###############################################################
    # calculate power spectral density as FT of Burg coeffitients utilizing spectra package
    psd = arma2psd(A=coeffs, rho=sigma, T=1/timestep, NFFT=nfft, sides='default', norm=False)
    freq = np.fft.fftfreq(len(psd), d=timestep)
    # strip the left (negative) side
    psd = psd[1:len(psd) // 2]
    freq = freq[1:len(freq) // 2]
    ##################################################################
    if plot_style == 0:
        x_ticks = freq
    else:
        x_ticks = freq / c / 100

    fig, axs = plt.subplots(2)
    axs[0].plot(x_ticks, psd)
    axs[1].plot(coeffs)
    plt.show()

def Example_6():
    # Define Constants
    conversion_factor = 3.33564 * 1e-30  # convert Debye to Cm
    T = 300  # thermodynamic temperature in Kelvin
    psPerFrame = 0.05  # ps between frames in dip file
    max_lag_ps = 100  # maximum evaluated point for correlation function in ps
    num_of_traj = 10  # number of trajectories
    nfft = 1028
    plot_style = 1  # 0 == Frequency, 1 == wavelength
    #################
    meanCorr = 0
    max_lag_frames = ceil(max_lag_ps / psPerFrame)
    corr_time = np.linspace(0, max_lag_ps, max_lag_frames + 1)
    for i in range(num_of_traj):
        dipFileName = "../trajectories/protein_rot_trans_%i.dat" % i
        # import dipole file
        dip_data, V = importDipoleData(dipFileName)
        dip_array = np.array(dip_data[['dip_x', 'dip_y', 'dip_z']])
        dip_array = np.transpose(dip_array)
        # calculate correlation function
        dip_der = np.diff(dip_array)
        corr_fun = calc_two_side_corr(dip_der, dip_der, max_lag_frames)
        # Convert to proper units
        constant = conversion_factor ** 2 / (V * 3 * k * T * epsilon_0)
        corr_fun = corr_fun * constant
        meanCorr += corr_fun

    meanCorr /= num_of_traj
    denf = fft(meanCorr, nfft)
    freq = np.fft.fftfreq(len(denf), d=psPerFrame*1e-12)
    denf = denf[1:len(denf) // 2]
    freq = freq[1:len(freq) // 2]
    ##################################################################
    if plot_style == 0:
        x_ticks = freq
    else:
        x_ticks = freq / c / 100
    plt.plot(x_ticks, denf)
    plt.show()

if __name__ == '__main__':
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
    # Example 5 - Calculate power spectral density from dipole moment evolution by single channel MEM method
    # Example_5()
    # Example 6 - Calculate power spectral density as Fourier transform of dipole moment time correlation function
    # Example_6()
