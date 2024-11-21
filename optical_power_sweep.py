"""
Script to create a power_sweep

With this script, we are sweeping the power for the input RF powers
"""
import importlib.util

import pandas as pd

#default path for current release
spec_lin = importlib.util.spec_from_file_location('lumapi', "/opt/lumerical/v241/api/python/lumapi.py")
#Functions that perform the actual loading
lumapi = importlib.util.module_from_spec(spec_lin)
spec_lin.loader.exec_module(lumapi)

import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.io import savemat, loadmat
import pandas as pd
import os
from matplotlib.backends.backend_pdf import PdfPages
from init_constants import init_constants

# Open an INTERCONNECT session
# ic = lumapi.INTERCONNECT()
# ic.switchtodesign()

# First, execute the prerequisite script with open(
# "/home/maria/PycharmProjects/Solvers/Lumerical-python/Lum-python/interconnect/transmitter_MB/transmitter.py") as f:
# exec(f.read())


# (P_opt_dBm, SOA_gain, SOA_NF, SOA_NFactor, f_Signal, f_Spur, f_LO, c, lambda_central, f_Optical, RIN_dB,
#  BalancedDetection, V_pi, V_RF, V_LO, h, T0, kB, q, power_laser, sample_rate, R_out, R_in) = init_constants()

# Constant
P_opt_dBm = 10  # Input optical power [dBm]
SOA_gain = 15  # Gain dB of the SOA
SOA_NF = 6  # NF [dB] of the SOA
SOA_NFactor = 10 ** (SOA_NF / 10)

f_Signal = 30e9  # Signal frequency [GHz]
f_Spur = 30.1e9  # Spur frequency [GHz]
f_LO = 25  # LO frequency [GHz]

c = 3e8  # Speed of light (m/s)
lambda_central = 1.55e-6  # Central wavelength (meters)
f_Optical = c / lambda_central
RIN_dB = -160  # Laser RIN [dBc/Hz]
BalancedDetection = 0  # 0 -> single-ended detection, 1 -> balanced detection
V_pi = 4  # Half-wave voltage [V]
V_RF = 0.316  # RF voltage [V]
V_LO = 1.2  # LO voltage [V]

h = 6.626e-34  # Plank constant
T0 = 290  # System temperature [K]
kB = 1.38e-23  # Boltzman constant [J/K]
q = 1.6e-19  # Elementary electric charge

impedance = 50  # 50 Ohms
power_laser = 10  # 10 dBm
sample_rate = 300e9  # 300 GHz
R_out = impedance  # Output impedance [Ohm]
R_in = impedance  # Input impedance [Ohm]

H_pd = 1 / 2

pts = 27  # Number of points for the sweep in loop
pts_out = 9


def add_sweep(ic, P_opt_dBm, SOA_gain, SOA_NF, f_Signal, f_Spur, f_Optical,
              RIN_dB, Balanced_detection):
    # Create inner loop
    # Delete any existing sweep with the name "optical power sweep
    # ic.deletesweep("optical_power_sweep")

    # Initialize the inner loop sweep
    ic.addsweep()
    optPower = "optPower"
    ic.setsweep("sweep", "name", optPower)
    ic.setsweep(optPower, "type", "Values")
    ic.setsweep(optPower, "number of points", pts)

    # Define the values to be swept for P_opt
    # Generate the vector of values here
    power_values = [
        0.0001, 0.000125893, 0.000158489, 0.000199526, 0.000251189, 0.000316228,
        0.000398107, 0.000501187, 0.000630957, 0.000794328, 0.001, 0.00125893,
        0.00158489, 0.00199526, 0.00251189, 0.00316228, 0.00398107, 0.00501187,
        0.00630957, 0.00794328, 0.01, 0.0125893, 0.0158489, 0.0199526, 0.0251189,
        0.0316228, 0.0398107
    ]

    # Add parameter "P_opt" to the sweep with specified values
    para = {"Name": "P_opt",
            "Parameter": "::Root Element::CW_Laser::power",
            "Type": "Number"}

    for i, value in enumerate(power_values, start=1):
        para[f"Value_{i}"] = value

    ic.addsweepparameter(optPower, para)

    # Define results outputs for the inner sweep
    result_1 = {"Name": "Electrical_Spectrum",
                "Result": "::Root Element::ESA_output::signal"}

    result_2 = {"Name": "Optical_Spectrum",
                "Result": "::Root Element::OSA_PD::sum/signal"}

    result_3 = {"Name": "Optical_Power_PD",
                "Result": "::Root Element::PowerMeter_PD::sum/power"}

    result_4 = {"Name": "PowerMeter_SOAin",
                "Result": "::Root Element::PowerMeter_SOAin::sum/power"}

    results = [result_1, result_2, result_3, result_4]

    # Add all these results to the sweep
    for result in results:
        ic.addsweepresult(optPower, result)

    #############################################################################
    # Insert an outer loop sweep
    ic.insertsweep("optPower")
    opt_p_sweep = "optical_power_sweep"
    ic.setsweep("sweep", "name", opt_p_sweep)
    ic.setsweep(opt_p_sweep, "type", "Ranges")
    ic.setsweep(opt_p_sweep, "Number of points", pts_out)  #9 points

    # Add the parameters to be swept through
    para2 = {"Name": "RF_DC_bias",
             "Parameter": "::Root Element::RF_signal::amplitude",
             "Type": "Number",
             "Start": 0.25,
             "Stop": 2.25}
    ic.addsweepparameter(opt_p_sweep, para2)

    V_RF = "RF_voltage"
    para3 = {"Name": V_RF,
             "Parameter": "::Root Element::MZM_RF::bias voltage 1",
             "Type": "Number",
             "Start": 0.7 / np.pi,
             "Stop": 0.7 / np.pi * 9}
    ic.addsweepparameter(opt_p_sweep, para3)

    para4 = {"Name": "Vpi_RF_DC",
             "Parameter": "::Root Element::MZM_RF::pi dc voltage",
             "Type": "Number",
             "Start": 1,
             "Stop": 9}
    ic.addsweepparameter(opt_p_sweep, para4)

    Vpi_RF = "Vpi_RF"
    para5 = {"Name": Vpi_RF,
             "Parameter": "::Root Element::MZM_RF::pi rf voltage",
             "Type": "Number",
             "Start": 1,
             "Stop": 9}
    ic.addsweepparameter(opt_p_sweep, para5)

    para6 = {"Name": "LO_DC_bias",
             "Parameter": "::Root Element::MZM_LO::bias voltage 1",
             "Type": "Number",
             "Start": 0.25,
             "Stop": 2.25}
    ic.addsweepparameter(opt_p_sweep, para6)

    para7 = {"Name": "Vpi_LO_DC",
             "Parameter": "::Root Element::MZM_LO::pi dc voltage",
             "Type": "Number",
             "Start": 1,
             "Stop": 9}
    ic.addsweepparameter(opt_p_sweep, para7)

    para8 = {"Name": "Vpi_LO",
             "Parameter": "::Root Element::MZM_LO::pi rf voltage",
             "Type": "Number",
             "Start": 1,
             "Stop": 9}
    ic.addsweepparameter(opt_p_sweep, para8)

    V_LO = "LO_voltage"
    para9 = {"Name": V_LO,
             "Parameter": "::Root Element::LO_signal::amplitude",
             "Type": "Number",
             "Start": 0.942 / np.pi * 1,
             "Stop": 0.942 / np.pi * 9}
    ic.addsweepparameter(opt_p_sweep, para9)

    # Define outer loop parameters to sweep with ranges
    result_11 = {"Name": "Electrical_Spectrum",
                 "Result": "Electrical_Spectrum"}

    result_22 = {"Name": "Optical_Spectrum",
                 "Result": "Optical_Spectrum"}

    result_33 = {"Name": "Optical_Power_PD",
                 "Result": "Optical_Power_PD"}

    result_44 = {"Name": "PowerMeter_SOAin",
                 "Result": "PowerMeter_SOAin"}

    results2 = [result_11, result_22, result_33, result_44]

    # Add all these results to the sweep
    for result in results2:
        ic.addsweepresult(opt_p_sweep, result)

    para_values = [para2, para3, para4, para5, para6, para7, para8, para9]

    for result in results2:
        ic.addsweepresult(opt_p_sweep, result)

    # # ic.delete()
    # for value in para_values:
    #     ic.clear(value)
    # for result in results2:
    #     ic.clear(result)

    # Run the outer sweep
    ic.runsweep("optical_power_sweep")


def export_results_Popt_sweeps(ic, V_pi, V_RF, V_LO, P_opt_dBm, SOA_gain, SOA_NF, f_Signal, f_Spur, f_LO,
                               f_Optical, RIN_dB, Balanced_detection):
    # def export_results_Popt_sweeps(ic, P_opt_dBm=P_opt_dBm, SOA_gain=SOA_gain, SOA_NF=SOA_NF, f_Signal=f_Signal,
    #                                   f_Spur=f_Spur, f_LO=f_LO, f_Optical=f_Optical,
    #                                   RIN_dB=RIN_dB, Balanced_detection=BalancedDetection):
    """    Export the results from the sweep of the Optical power with Vpi, V_Rf and V_LO constant form the
    initialization of the simulation

    Args:
        ic:
        V_pi:
        V_RF:
        V_LO
        ic:
        P_opt_dBm:
        SOA_gain:
        SOA_NF:
        f_Signal:
        f_Spur:
        f_LO:
        f_Optical:
        RIN_dB:
        Balanced_detection:

    Returns: Results from the optical power sweep

    Note:   In this case, we are using constant values for Vpi, Vlo and Vrf, so that is why we export them like that
    In the case you want to get the data used for getting the sweeps, you can use the method described with the
    command ic.getsweepdata
    # Get the parameter values
    V_LO_ = "LO_voltage"
    Vpi_RF_ = "Vpi_RF"
    V_RF_ = "RF_voltage"

    V_LO = ic.getsweepdata(swpName_2, V_LO_)
    Vpi_RF = ic.getsweepdata(swpName_2, Vpi_RF_)
    V_RF = ic.getsweepdata(swpName_2, V_RF_)

    """
    # Get sweep results
    swpName_2 = "optical_power_sweep"

    # Optical power
    swp_2_data_Optical = ic.getsweepresult(swpName_2, "Optical_Spectrum")
    swp_2_data_Electrical = ic.getsweepresult(swpName_2, "Electrical_Spectrum")
    swp_2_data_PowerPD = ic.getsweepresult(swpName_2, "Optical_Power_PD")
    swp_2_data_PowerInSOA = ic.getsweepresult(swpName_2, "PowerMeter_SOAin")

    # In this case, we are using constant values for Vpi, Vlo and Vrf, so that is why we export them like that
    # In the case you want to get the data used for getting the sweeps, you can use the method described with the
    # command ic.getsweepdata
    # # Get the parameter values
    # V_LO_ = "LO_voltage"
    # Vpi_RF_ = "Vpi_RF"
    # V_RF_ = "RF_voltage"
    #
    # V_LO = ic.getsweepdata(swpName_2, V_LO_)
    # Vpi_RF = ic.getsweepdata(swpName_2, Vpi_RF_)
    # V_RF = ic.getsweepdata(swpName_2, V_RF_)

    # Create a dictionary to export the data
    # Check if it is better to save it asa numpy array or as a panda dataframe
    # And how will it manage the data if the columns are not of equal length

    data_export = {"Vpi_RF": V_pi,
                   "V_RF": V_RF,
                   "V_LO": V_LO,
                   "P_opt_dBm": P_opt_dBm,
                   "SOA_gain": SOA_gain,
                   "SOA_NF": SOA_NF,
                   "f_Signal": f_Signal,
                   "f_Spur": f_Spur,
                   "f_LO": f_LO,
                   "f_Optical": f_Optical,
                   "RIN_dB": RIN_dB,
                   "Balanced_detection": Balanced_detection,
                   "OS_Optical_power_sweep": swp_2_data_Optical,
                   "ES_Optical_power_sweep": swp_2_data_Electrical,
                   "PowerPD_Optical_power_sweep": swp_2_data_PowerPD,
                   "PowerInSOA_Optical_power_sweep": swp_2_data_PowerInSOA}

    # savemat("SimulationData_Optical_sweep.mat", data_export)
    # ic.write("SimulationData_Optical_RF.csv", data_export)

    return data_export


def trial_export(ic, data_export, file_name="SimulationData_Optical.csv"):
    ic.write(file_name, data_export)


def load_icp_file(ic, file_path):
    """
    To load and open an .icp file, to get the sweep data from there
    Args:
        ic:
        file_path:

    Returns:

    """
    # load the .icp file
    try:
        ic.load(file_path)
        print(f"File '{file_path}' successfully loaded")
    except Exception as e:
        print(f"Error loading file '{file_path}': {e} ")

    return ic


# file_path = "~/transmitter_MB/transmitter_optical_power_sweep.icp"
def load_data_from_mat(file_path):
    """
    Load data from a .mat file and return it as a dictionary.
    """
    data = loadmat(file_path)
    # Clean up the keys by removing MATLAB metadata keys
    cleaned_data = {key: value for key, value in data.items() if not key.startswith('__')}
    return cleaned_data

def export_results_ic(file_path):
    """
    To export the data from an .icp file

    Args:
        file_path:

    Returns:

    """
    # Start interconnect session
    ic = lumapi.INTERCONNECT()
    load_icp_file(ic, file_path)
    data_export = export_results_Popt_sweeps(ic, V_pi, V_RF, V_LO, P_opt_dBm, SOA_gain, SOA_NF, f_Signal, f_Spur, f_LO,
                                             f_Optical, RIN_dB, BalancedDetection)
    # trial_export(ic, data_export)
    return data_export


def data_to_csv(data_export):
    os_opt_sweep_data = np.squeeze(data_export["OS_Optical_power_sweep"])  # optical spectrum
    es_opt_sweep_data = np.squeeze(data_export["ES_Optical_voltage_sweep"])  # electrical spectrum
    P_pd_opt_sweep_data = np.squeeze(data_export["PowerPD_Optical_power_sweep"])  # data power PD

    os_freq_RF_sweep = data_export["OS_RF_voltage_sweep"]

    # np.squeeze()


def process_data(mat_file_path, f_Optical, R_in=50, detectionScheme=1, H_pd=1 / 2, freq_RF=30e9, freq_LO=25e9,
                 freq_Optical=f_Optical, SOA_NFactor=SOA_NFactor, R_out=50):
    """
    In here, if you want to get the data directly from Lumerical everytime, the 1st argument will be data_export
    If not use the file path of the .mat file
    Args:
        mat_file_path:
        f_Optical:
        R_in:
        detectionScheme:
        H_pd:
        freq_RF:
        freq_LO:
        freq_Optical:
        SOA_NFactor:
        R_out:

    Returns:

    """
    #In this code I am having an error because I need to process the .mat files
    #as they are not flatten
    #In this case to avoid this, probably, I have to process the .mat, if not
    #use another type of data to avoid this

    #Load data from the .mat file
    data_export = load_data_from_mat(mat_file_path)

    # global SinglePD_DC_current
    h = 6.626e-34  # Plank constant
    T0 = 290  # System temperature [K]
    kB = 1.38e-23  # Boltzman constant [J/K]
    q = 1.6e-19  # Elementary electric charge


    # Rename the data from the dictionary
    os_opt_sweep_data = data_export["OS_Optical_power_sweep"]  # optical spectrum
    es_opt_sweep_data = data_export["ES_Optical_power_sweep"]  # electrical spectrum
    P_opt_sweep_data = data_export["PowerPD_Optical_power_sweep"]  # data power PD

    print(np.shape(es_opt_sweep_data['P_opt']))
    print(es_opt_sweep_data['P_opt'])
    print(data_export)

    # Squeeze the data obtained, because it has multiple dimensions
    os_sweep_data = np.squeeze(os_opt_sweep_data['power (dBm)'])  # optical spectrum [dBm]
    es_sweep_data = np.squeeze(es_opt_sweep_data['power (dBm)'])  # electrical spectrum [dBm]

    os_freq_sweep = (c / os_opt_sweep_data["Frequency"]) * 1e6  # Wavelength -> Optical frequency [um]
    es_freq_sweep = (es_opt_sweep_data['frequency'] * 1e-9)  # Electrical frequency [GHz]

    # Extract P_opt from the data and ensure numerical compatibility
    P_opt = es_opt_sweep_data['P_opt']

    # Access and flatten the nested structure
    if isinstance(P_opt, np.ndarray) and P_opt.size == 1:
        # If the array contains one nested element, extract it
        P_opt = P_opt[0][0]  # Access the innermost array
    elif isinstance(P_opt, list) and len(P_opt) == 1:
        P_opt = P_opt[0][0]  # Adjust for list-based cases
    else:
        raise ValueError("Unexpected structure in P_opt")

    # Now P_opt should be a clean 1D array
    # print(f"Flattened P_opt: {P_opt}")

    # Apply the logarithmic conversion
    powerLaser = 10 * np.log10( P_opt * 1e3)  # Laser power [dBm] , because the laser power
    # is constant
    print(P_opt_sweep_data['power (dBm)'].shape)
    powerPD = np.squeeze(P_opt_sweep_data['power (dBm)'][:, 3]).T  # Power at PD [dBm] #Check the annotation here

    # opticalGain = powerPD[1] - powerLaser[1]  # Optical gain from the laser to the PD [dB]

    # Transpose is `.T` in NumPy
    data2Export = np.column_stack((powerLaser.T, powerPD.T))

    # Export the file on .csv format
    dataExportLocation = ("/home/maria/PycharmProjects/Solvers/Lumerical-python/Lum-python/interconnect/transmitter_MB"
                          "/data")
    filename = 'Optical Powers.csv'
    export_csv(data2Export, filename, dataExportLocation)

    # Calculate the output current and normalize output power to 50 [Ohms]
    # Watch out the power from the Spectrum analyzer component in Interconnect
    # Because it plots the power normalized to impedance of 1 Ohm

    powerOut_dB = es_sweep_data  # take the powers of the selected frequencies
    powerOut = 10 ** (powerOut_dB / 10) * 1e-3  # Conversion from power in [dBm] to [W]
    currentOut = np.sqrt(powerOut / R_in)  # Output current in [A]

    # print(detectionScheme)
    if detectionScheme == 1:
        temp = data_export["PowerInSOA_Optical_power_sweep"]["power (dBm)"]
        temp = np.squeeze(temp)
        temp = 10 ** (temp / 10) * 1e-3

        # Single PD current [A] matrix [frequency, optical power, Vpi]
        current_SinglePD = np.sqrt(temp / R_in)
        SinglePD_DC_current = np.squeeze(current_SinglePD[0, :])
        # print(f'Yes')
        temp = None  # To clear the variable

    power50_ohms = 10 * np.log10(currentOut ** 2 * R_in * np.abs(H_pd) ** 2 * 1e3)

    # Calculating the power output at 50 [Ohms]
    Vpi = np.linspace(1, 9, 9)  # Vpi[V] sweep

    # Select the frequencies with reasonable power output and calculate the different parameters for the signal
    freq_IF = freq_RF - freq_LO
    DC_index = 0  # Assuming DC is at index 0

    freqsTemp = np.abs(es_freq_sweep - freq_IF)  # So if I have a zero it will be at minimum distance
    IF_index = np.argmin(freqsTemp)
    # IF_index = np.where(freqsTemp == np.min(freqsTemp))

    freqsTemp = np.abs(es_freq_sweep - freq_LO)
    LO_index = np.argmin(freqsTemp)
    # LO_index = np.where(freqsTemp == np.min(freqsTemp))

    freqsTemp = np.abs(es_freq_sweep - freq_RF)
    RF_index = np.argmin(freqsTemp)
    # RF_index = np.where(freqsTemp == np.min(freqsTemp))

    display_freq = [DC_index, IF_index, LO_index, RF_index]  # Select the frequency index of the IF (or RF/LO) component
    freqs_OptSweep = [0, freq_IF, freq_LO, freq_RF]

    # RF Voltage normalization
    rf_volt_opt = es_opt_sweep_data['RF_voltage']
    P0_RF_dBm = 10 * np.log10(rf_volt_opt / np.sqrt(2) ** 2 / R_in * 1e3)  # Average power extraction (A/sqrt(2))^2/ 50

    # Gain calculation [dB ] calculation at the selected frequencies
    # gain_OptSweep = es_sweep_data[display_freq, :, :] - P0_RF_dBm[:, np.newaxis, np.newaxis]
    gain_OptSweep = es_sweep_data[display_freq, :, :] - P0_RF_dBm

    # Print circuit parameters
    VpiIndex_plot = 3  #Python starts by 0
    Vpi_plot = Vpi[VpiIndex_plot]

    optPower_Index = 4
    # print(len(powerLaser))
    # powerLaser = powerL
    print(f'power laser: ', powerLaser)
    powerLaser_plot = powerLaser[0, optPower_Index]
    powerPD_OptSweep_plot = powerPD[optPower_Index]

    # Normalize output electrical spectrum
    power50_ohms_data_norm = (power50_ohms[:, optPower_Index, VpiIndex_plot] -
                              power50_ohms[ IF_index, optPower_Index, VpiIndex_plot])

    # Export normalized spectrum
    data2Export = power50_ohms_data_norm
    pd.DataFrame(data2Export).to_csv("Electrical Normalized Spectrum from Optical Sweep.csv", index=False)

    # N_out_dBm, NF_dB = noise_power_calculations(Vpi, powerLaser, gain_OptSweep, currentOut, freq_Optical,
    #                                             SOA_NFactor, H_pd, detectionScheme, SinglePD_DC_current, T0,
    #                                             R_out)
    #
    # pdf_filename = "Optical_power_sweep_graphs.pdf"
    # plot_and_save_to_pdf(pdf_filename, dataExportLocation, es_freq_sweep, power50_ohms, optPower_Index,
    #                      VpiIndex_plot, P0_RF_dBm, Vpi_plot, powerLaser_plot, os_freq_sweep, os_sweep_data,
    #                      powerPD, gain_OptSweep, freqs_OptSweep, currentOut, display_freq,
    #                      N_out_dBm, NF_dB, Vpi)

def process_data_icp(data_export, f_Optical, R_in=50, detectionScheme=1, H_pd=1 / 2, freq_RF=30e9, freq_LO=25e9,
                 freq_Optical=f_Optical, SOA_NFactor=SOA_NFactor, R_out=50):
    """
    In here, if you want to get the data directly from Lumerical everytime, the 1st argument will be data_export
    If not use the file path of the .mat file
    Args:
        mat_file_path:
        f_Optical:
        R_in:
        detectionScheme:
        H_pd:
        freq_RF:
        freq_LO:
        freq_Optical:
        SOA_NFactor:
        R_out:

    Returns:

    """
    # global SinglePD_DC_current
    h = 6.626e-34  # Plank constant
    T0 = 290  # System temperature [K]
    kB = 1.38e-23  # Boltzman constant [J/K]
    q = 1.6e-19  # Elementary electric charge

    # Rename the data from the dictionary
    os_opt_sweep_data = data_export["OS_Optical_power_sweep"]  # optical spectrum
    es_opt_sweep_data = data_export["ES_Optical_power_sweep"]  # electrical spectrum
    P_opt_sweep_data = data_export["PowerPD_Optical_power_sweep"]  # data power PD

    # Squeeze the data obtained, because it has multiple dimensions
    os_sweep_data = np.squeeze(os_opt_sweep_data['power (dBm)'])  # optical spectrum [dBm]
    es_sweep_data = np.squeeze(es_opt_sweep_data['power (dBm)'])  # electrical spectrum [dBm]

    os_freq_sweep = (c / os_opt_sweep_data["Frequency"]) * 1e6  # Wavelength -> Optical frequency [um]
    es_freq_sweep = (es_opt_sweep_data['frequency'] * 1e-9)  # Electrical frequency [GHz]

    powerLaser = 10 * np.log10(es_opt_sweep_data['P_opt'] * 1e3)  # Laser power [dBm] , because the laser power
    # is constant
    print(P_opt_sweep_data['power (dBm)'].shape)
    powerPD = np.squeeze(P_opt_sweep_data['power (dBm)'][:, 3]).T  # Power at PD [dBm] #Check the annotation here

    # opticalGain = powerPD[1] - powerLaser[1]  # Optical gain from the laser to the PD [dB]

    # Transpose is `.T` in NumPy
    data2Export = np.column_stack((powerLaser.T, powerPD.T))

    # Export the file on .csv format
    dataExportLocation = ("/home/maria/PycharmProjects/Solvers/Lumerical-python/Lum-python/interconnect/transmitter_MB"
                          "/data")
    filename = 'Optical Powers.csv'
    export_csv(data2Export, filename, dataExportLocation)

    # Calculate the output current and normalize output power to 50 [Ohms]
    # Watch out the power from the Spectrum analyzer component in Interconnect
    # Because it plots the power normalized to impedance of 1 Ohm

    powerOut_dB = es_sweep_data  # take the powers of the selected frequencies
    powerOut = 10 ** (powerOut_dB / 10) * 1e-3  # Conversion from power in [dBm] to [W]
    currentOut = np.sqrt(powerOut / R_in)  # Output current in [A]

    print(detectionScheme)
    if detectionScheme == 1:
        temp = data_export["PowerInSOA_Optical_power_sweep"]["power (dBm)"]
        temp = np.squeeze(temp)
        temp = 10 ** (temp / 10) * 1e-3

        # Single PD current [A] matrix [frequency, optical power, Vpi]
        current_SinglePD = np.sqrt(temp / R_in)
        SinglePD_DC_current = np.squeeze(current_SinglePD[0, :])
        print(f'Yes')
        temp = None  # To clear the variable

    power50_ohms = 10 * np.log10(currentOut ** 2 * R_in * np.abs(H_pd) ** 2 * 1e3)

    # Calculating the power output at 50 [Ohms]
    Vpi = np.linspace(1, 9, 9)  # Vpi[V] sweep

    # Select the frequencies with reasonable power output and calculate the different parameters for the signal
    freq_IF = freq_RF - freq_LO
    DC_index = 0  # Assuming DC is at index 0

    freqsTemp = np.abs(es_freq_sweep - freq_IF)  # So if I have a zero it will be at minimum distance
    IF_index = np.argmin(freqsTemp)
    # IF_index = np.where(freqsTemp == np.min(freqsTemp))

    freqsTemp = np.abs(es_freq_sweep - freq_LO)
    LO_index = np.argmin(freqsTemp)
    # LO_index = np.where(freqsTemp == np.min(freqsTemp))

    freqsTemp = np.abs(es_freq_sweep - freq_RF)
    RF_index = np.argmin(freqsTemp)
    # RF_index = np.where(freqsTemp == np.min(freqsTemp))

    display_freq = [DC_index, IF_index, LO_index, RF_index]  # Select the frequency index of the IF (or RF/LO) component
    freqs_OptSweep = [0, freq_IF, freq_LO, freq_RF]

    # RF Voltage normalization
    rf_volt_opt = es_opt_sweep_data['RF_voltage']
    P0_RF_dBm = 10 * np.log10(rf_volt_opt / np.sqrt(2) ** 2 / R_in * 1e3)  # Average power extraction (A/sqrt(2))^2/ 50

    # Gain calculation [dB ] calculation at the selected frequencies
    # gain_OptSweep = es_sweep_data[display_freq, :, :] - P0_RF_dBm[:, np.newaxis, np.newaxis]
    gain_OptSweep = es_sweep_data[display_freq, :, :] - P0_RF_dBm

    # Print circuit parameters
    VpiIndex_plot = 3  #Python starts by 0
    Vpi_plot = Vpi[VpiIndex_plot]

    optPower_Index = 4
    print(len(powerLaser))
    # powerLaser = powerL
    print(f'power laser: ', powerLaser)
    powerLaser_plot = powerLaser[0, optPower_Index]
    powerPD_OptSweep_plot = powerPD[optPower_Index]

    # Normalize output electrical spectrum
    power50_ohms_data_norm = (power50_ohms[:, optPower_Index, VpiIndex_plot] -
                              power50_ohms[ IF_index, optPower_Index, VpiIndex_plot])

    # Export normalized spectrum
    data2Export = power50_ohms_data_norm
    pd.DataFrame(data2Export).to_csv("Electrical Normalized Spectrum from Optical Sweep.csv", index=False)

    # N_out_dBm, NF_dB = noise_power_calculations(Vpi, powerLaser, gain_OptSweep, currentOut, freq_Optical,
    #                                             SOA_NFactor, H_pd, detectionScheme, SinglePD_DC_current, T0,
    #                                             R_out)
    #
    # pdf_filename = "Optical_power_sweep_graphs.pdf"
    # plot_and_save_to_pdf(pdf_filename, dataExportLocation, es_freq_sweep, power50_ohms, optPower_Index,
    #                      VpiIndex_plot, P0_RF_dBm, Vpi_plot, powerLaser_plot, os_freq_sweep, os_sweep_data,
    #                      powerPD, gain_OptSweep, freqs_OptSweep, currentOut, display_freq,
    #                      N_out_dBm, NF_dB, Vpi)


def noise_power_calculations(Vpi, powerLaser, gain_OptSweep, currentOut, freq_Optical,
                             SOA_NFactor, H_pd, detectionScheme, SinglePD_DC_current, T0=300, R_out=50):
    """
        Function to calculate noise power and noise figure based on the given parameters.
        Translated from MATLAB to Python. """

    # Plot ( os_freq_sweep, data2Export)
    # Noise power calculations ----------------------------------------------------------------------------------------
    numVpi = len(Vpi)
    numOpticalPowers = len(powerLaser)
    print(numOpticalPowers)

    I_DC = np.zeros(numOpticalPowers)
    N_out = np.zeros((numVpi, numOpticalPowers))
    NF = np.zeros((numVpi, numOpticalPowers))
    # I_DC = []

    # Ad-hoc solution for power at SOA
    powerIn_SOA = 10 ** ((powerLaser - 8) / 10) * 1e-3  # Convert dBm to Watts
    RIN_SOA = 2 * h * freq_Optical * SOA_NFactor / powerIn_SOA  # Using Planck constant for h

    if detectionScheme == 0:
        # Single-ended detection
        for ii in range(numVpi):
            # Signal gain [linear]
            g_signal_dB = np.squeeze(gain_OptSweep[1, :, ii])  # Gain in dB
            g_signal = 10 ** (g_signal_dB / 10)  # Convert to linear

            # LO gain [linear]
            g_LO_dB = gain_OptSweep[2, :, ii]  # Gain in dB
            g_LO = 10 ** (g_LO_dB / 10)  # Convert to linear

            # Output noise power spectral density [W/Hz]
            I_DC = currentOut[1, :, ii]
            print(f"I_DC shape: {I_DC.shape}")
            N_out[ii, :] = (   # This is giving me problems, bc of the shape of the vector
                    kB * T0 * g_signal +
                    kB * T0 * g_LO +
                    kB * T0 +
                    2 * q * I_DC * R_out * np.abs(H_pd) ** 2 +
                    RIN_SOA * I_DC ** 2 * R_out * np.abs(H_pd) ** 2
            )

            # Noise figure calculation
            NF[ii, :] = N_out[ii, :] / (g_signal * kB * T0)  # Link Noise Figure (NF) #prob error ii-1
    else:
        # Differential detection
        for ii in range(numVpi):
            # RF-to-RF/RF-to-IF and LO gain in the link to calculate the amplification of inout
            # Signal gain [linear]
            # Thermal noise
            g_signal_dB = np.squeeze(gain_OptSweep[1, :, ii])  # Gain in dB #prob error ii-1
            g_signal = 10 ** (g_signal_dB / 10)  # Convert to linear

            # LO gain [linear]
            g_LO_dB = gain_OptSweep[2, :, :]  # Gain in dB
            g_LO = np.squeeze(10 ** (g_LO_dB / 10))  # Convert to linear

            # Output noise power spectral density [W/Hz] for each incident optical power
            I_DC = SinglePD_DC_current
            print(len(I_DC))

            # Adjustments for differential detection
            N_out[ii, :] = (
                    kB * T0 * g_signal +
                    kB * T0 * g_LO[:, ii] +
                    2 * kB * T0 +
                    2 * q * (2 * I_DC[:, ii]) * R_out * abs(H_pd) ** 2
            )

            # Noise figure calculation
            NF[ii, :] = N_out[ii, :] / (g_signal * kB * T0)  # Link Noise Figure (NF)

        # Convert to dBm/Hz
    N_out_dBm = 10 * np.log10(N_out * 1e3)  # [dBm/Hz]
    NF_dB = 10 * np.log10(NF)  # [dB]

    return N_out_dBm, NF_dB


def plot_and_export(dataExportLocation, ES_freqs_OptSweep, power50Ohms, optPower_Index, VpiIndex_plot,
                    P0_RF_dBm, Vpi_plot, powerLaser_plot, OS_freqs_OptSweep, OS_OptSweep, powerPD_OptSweep,
                    gain_OptSweep, freqs_OptSweep, current_OptSweep, display_freqs_OptSweep,
                    N_out_dBm, NF_dB, Vpi_OptSweep):
    # Plot: Electrical Spectrum at link output
    x_axis = ES_freqs_OptSweep
    y_axis = power50Ohms[:, optPower_Index, VpiIndex_plot]

    plt.figure("Output electrical power spectrum")
    plt.plot(x_axis, y_axis)
    plt.xlabel("f (GHz)")
    plt.ylabel("P_out (dBm)")
    plt.ylim([-100, 10])
    plt.xlim([0, 50])
    plt.grid(True)
    plt.title("Output electrical power spectrum")
    plt.text(8, -15, f"P_RF, in. = {P0_RF_dBm:.1f} dBm", fontsize=8, fontweight='bold', color='red')
    plt.text(8, -10, f"V_pi = {Vpi_plot:.1f} V", fontsize=8, fontweight='bold', color='red')
    plt.text(8, -5, f"P_opt., in. = {powerLaser_plot:.1f} dBm", fontsize=8, fontweight='bold', color='red')

    filename = 'Electrical Spectrum from Optical Sweep.csv'
    np.savetxt(os.path.join(dataExportLocation, filename), np.column_stack((x_axis, y_axis)), delimiter=",")

    # Plot: Optical Spectrum at PD
    x_axis = OS_freqs_OptSweep
    y_axis = OS_OptSweep[:, optPower_Index, VpiIndex_plot]

    plt.figure("Optical power spectrum at PD")
    plt.plot(x_axis, y_axis)
    plt.xlabel("Optical offset frequency, f_off. (GHz)")
    plt.ylabel("P_opt. (dBm)")
    plt.ylim([-100, 20])
    plt.xlim([-150, 150])
    plt.grid(True)
    plt.title("Optical power spectrum at PD")
    plt.text(193.107, -5, f"V_pi = {Vpi_plot:.1f} V", fontsize=8, fontweight='bold', color='red')
    plt.text(193.107, -10, f"P_opt., in. = {powerLaser_plot:.1f} dBm", fontsize=8, fontweight='bold', color='red')

    filename = 'Optical Spectrum from Optical Sweep.csv'
    np.savetxt(os.path.join(dataExportLocation, filename), np.column_stack((x_axis, y_axis)), delimiter=",")

    # Plot: Output electrical gain @ specific frequency vs. Optical Power
    x_axis = powerPD_OptSweep
    y_axis = gain_OptSweep[1:, :, VpiIndex_plot]

    plt.figure("Electrical gain vs. Optical power at PD")
    plt.plot(x_axis, y_axis.T)
    plt.xlabel("P_laser + G_opt. (dBm)")
    plt.ylabel("G (dB)")
    plt.grid(True)
    plt.xlim([-17, 22])
    plt.ylim([-85, -15])
    plt.legend([f"{freq:.1f}" for freq in freqs_OptSweep[1:]], loc="northwest", title="Frequency (GHz)")
    plt.title("Electrical gain vs. Optical power at the PD")
    plt.text(10, -50, f"V_pi = {Vpi_plot:.1f} V", fontsize=8, fontweight='bold', color='red')

    # Export gain data
    np.savetxt(os.path.join(dataExportLocation, 'IF Gain from Optical Sweep.csv'),
               gain_OptSweep[1, :, :], delimiter=",")
    np.savetxt(os.path.join(dataExportLocation, 'LO Gain from Optical Sweep.csv'),
               gain_OptSweep[2, :, :], delimiter=",")
    np.savetxt(os.path.join(dataExportLocation, 'RF Gain from Optical Sweep.csv'),
               gain_OptSweep[3, :, :], delimiter=",")

    # Plot: Electrical current vs. Optical Power
    x_axis = powerPD_OptSweep
    y_axis = current_OptSweep[display_freqs_OptSweep, :, VpiIndex_plot] * 1e3  # Convert to mA

    plt.figure("Electrical current vs. Optical power at PD")
    plt.semilogy(x_axis, y_axis.T)
    plt.xlabel("Laser power, P_opt. (dBm)")
    plt.ylabel("I (mA)")
    plt.grid(True)
    plt.ylim([1e-3, 1e2])
    plt.legend(["DC"] + [f"{freq:.1f}" for freq in freqs_OptSweep[1:]], loc="northwest", title="Frequency (GHz)")
    plt.title("Electrical current vs. Optical power at the PD")
    plt.text(10, -50, f"V_pi = {Vpi_plot:.1f} V", fontsize=8, fontweight='bold', color='red')

    filename = 'Electrical Current from Optical Sweep.csv'
    np.savetxt(os.path.join(dataExportLocation, filename), y_axis.T, delimiter=",")

    # Plot: Output noise PSD vs. Optical Power
    x_axis = powerPD_OptSweep
    y_axis = N_out_dBm[VpiIndex_plot, :]

    plt.figure("Output noise power spectral density")
    plt.plot(x_axis, y_axis)
    plt.xlabel("P_opt. (dBm)")
    plt.ylabel("N_out (dBm/Hz)")
    plt.grid(True)
    plt.xlim([-17, 22])
    plt.title("Output noise power spectral density")

    filename = 'Output Noise PSD from Optical Sweep.csv'
    np.savetxt(os.path.join(dataExportLocation, filename), N_out_dBm.T, delimiter=",")

    # Plot: NF vs. Optical Power
    x_axis = powerPD_OptSweep
    y_axis = NF_dB[[1, 3, 5, 7], :]
    label_strings = [f"{Vpi:.1f}" for Vpi in Vpi_OptSweep[[1, 3, 5, 7]]]

    plt.figure("Link noise figure")
    plt.plot(x_axis, y_axis.T)
    plt.xlabel("P_opt. + G (dBm)")
    plt.ylabel("NF (dB)")
    plt.grid(True)
    plt.xlim([-17, 22])
    plt.ylim([35, 95])
    plt.legend(label_strings, loc="northeast", title="V_pi (V)")
    plt.title("Link noise figure")

    filename = 'NF from Optical Sweep.csv'
    np.savetxt(os.path.join(dataExportLocation, filename), NF_dB.T, delimiter=",")


def plot_and_save_to_pdf(pdf_filename, dataExportLocation, ES_freqs_OptSweep, power50Ohms, optPower_Index,
                         VpiIndex_plot, P0_RF_dBm, Vpi_plot, powerLaser_plot, OS_freqs_OptSweep, OS_OptSweep,
                         powerPD_OptSweep, gain_OptSweep, freqs_OptSweep, current_OptSweep, display_freqs_OptSweep,
                         N_out_dBm, NF_dB, Vpi_OptSweep):
    # Create a PDF file for saving all plots
    with PdfPages(pdf_filename) as pdf:
        # Plot: Electrical Spectrum at link output
        x_axis = ES_freqs_OptSweep
        y_axis = power50Ohms[:, optPower_Index, VpiIndex_plot]

        plt.figure("Output electrical power spectrum")
        plt.plot(x_axis, y_axis)
        plt.xlabel("f (GHz)")
        plt.ylabel("P_out (dBm)")
        plt.ylim([-100, 10])
        plt.xlim([0, 50])
        plt.grid(True)
        plt.title("Output electrical power spectrum")
        plt.text(8, -15, f"P_RF, in. = {P0_RF_dBm:.1f} dBm", fontsize=8, fontweight='bold', color='red')
        plt.text(8, -10, f"V_pi = {Vpi_plot:.1f} V", fontsize=8, fontweight='bold', color='red')
        plt.text(8, -5, f"P_opt., in. = {powerLaser_plot:.1f} dBm", fontsize=8, fontweight='bold', color='red')
        pdf.savefig()
        plt.close()

        # Plot: Optical Spectrum at PD
        x_axis = OS_freqs_OptSweep
        y_axis = OS_OptSweep[:, optPower_Index, VpiIndex_plot]

        plt.figure("Optical power spectrum at PD")
        plt.plot(x_axis, y_axis)
        plt.xlabel("Optical offset frequency, f_off. (GHz)")
        plt.ylabel("P_opt. (dBm)")
        plt.ylim([-100, 20])
        plt.xlim([-150, 150])
        plt.grid(True)
        plt.title("Optical power spectrum at PD")
        plt.text(193.107, -5, f"V_pi = {Vpi_plot:.1f} V", fontsize=8, fontweight='bold', color='red')
        plt.text(193.107, -10, f"P_opt., in. = {powerLaser_plot:.1f} dBm", fontsize=8, fontweight='bold', color='red')
        pdf.savefig()
        plt.close()

        # Plot: Output electrical gain vs. Optical Power
        x_axis = powerPD_OptSweep
        y_axis = gain_OptSweep[1:, :, VpiIndex_plot]

        plt.figure("Electrical gain vs. Optical power at PD")
        plt.plot(x_axis, y_axis.T)
        plt.xlabel("P_laser + G_opt. (dBm)")
        plt.ylabel("G (dB)")
        plt.grid(True)
        plt.xlim([-17, 22])
        plt.ylim([-85, -15])
        plt.legend([f"{freq:.1f}" for freq in freqs_OptSweep[1:]], loc="northwest", title="Frequency (GHz)")
        plt.title("Electrical gain vs. Optical power at the PD")
        plt.text(10, -50, f"V_pi = {Vpi_plot:.1f} V", fontsize=8, fontweight='bold', color='red')
        pdf.savefig()
        plt.close()

        # Plot: Electrical current vs. Optical Power
        x_axis = powerPD_OptSweep
        y_axis = current_OptSweep[display_freqs_OptSweep, :, VpiIndex_plot] * 1e3  # Convert to mA

        plt.figure("Electrical current vs. Optical power at PD")
        plt.semilogy(x_axis, y_axis.T)
        plt.xlabel("Laser power, P_opt. (dBm)")
        plt.ylabel("I (mA)")
        plt.grid(True)
        plt.ylim([1e-3, 1e2])
        plt.legend(["DC"] + [f"{freq:.1f}" for freq in freqs_OptSweep[1:]], loc="northwest", title="Frequency (GHz)")
        plt.title("Electrical current vs. Optical power at the PD")
        plt.text(10, -50, f"V_pi = {Vpi_plot:.1f} V", fontsize=8, fontweight='bold', color='red')
        pdf.savefig()
        plt.close()

        # Plot: Output noise PSD vs. Optical Power
        x_axis = powerPD_OptSweep
        y_axis = N_out_dBm[VpiIndex_plot, :]

        plt.figure("Output noise power spectral density")
        plt.plot(x_axis, y_axis)
        plt.xlabel("P_opt. (dBm)")
        plt.ylabel("N_out (dBm/Hz)")
        plt.grid(True)
        plt.xlim([-17, 22])
        plt.title("Output noise power spectral density")
        pdf.savefig()
        plt.close()

        # Plot: NF vs. Optical Power
        x_axis = powerPD_OptSweep
        y_axis = NF_dB[[1, 3, 5, 7], :]
        label_strings = [f"{Vpi:.1f}" for Vpi in Vpi_OptSweep[[1, 3, 5, 7]]]

        plt.figure("Link noise figure")
        plt.plot(x_axis, y_axis.T)
        plt.xlabel("P_opt. + G (dBm)")
        plt.ylabel("NF (dB)")
        plt.grid(True)
        plt.xlim([-17, 22])
        plt.ylim([35, 95])
        plt.legend(label_strings, loc="northeast", title="V_pi (V)")
        plt.title("Link noise figure")
        pdf.savefig()
        plt.close()


def export_csv(data2Export, filename, dataExportLocation):
    # Specify filename and file path

    fullPathFile = os.path.join(dataExportLocation, filename)

    # Save to CSV
    np.savetxt(fullPathFile, data2Export, delimiter=",", fmt="%.6f")

    print(f"Data successfully exported to {fullPathFile}")
