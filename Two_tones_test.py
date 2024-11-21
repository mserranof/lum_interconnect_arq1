"""
Script to create a Two tones test

This script should be compiled in order 3.
"""
import importlib.util

#default path for current release
spec_lin = importlib.util.spec_from_file_location('lumapi', "/opt/lumerical/v241/api/python/lumapi.py")
#Functions that perform the actual loading
lumapi = importlib.util.module_from_spec(spec_lin)
spec_lin.loader.exec_module(lumapi)

import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.io import savemat, loadmat
import optical_power_sweep as ops
import os

# Open an INTERCONNECT session
# ic = lumapi.INTERCONNECT()
# ic.switchtodesign()

# Constant
P_opt_dBm = 10  # Input optical power [dBm]
SOA_gain = 15  # Gain dB of the SOA
SOA_NF = 6  # NF [dB] of the SOA

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
R_in = 50


def add_sweep(ic):
    """
    Script created for the 2 tones test
    Args:
        ic:

    Returns:

    """
    pts = 46  # Number of points for the sweep in loop
    pts_out = 9

    # Initialize the inner loop sweep
    ic.addsweep()
    rf_volt = "RFVolt"
    ic.setsweep("sweep", "name", rf_volt)
    ic.setsweep(rf_volt, "type", "Values")
    ic.setsweep(rf_volt, "number of points", pts)

    # Add parameter RF_voltage signal to the sweep with specified values
    para = {"Name": "RF_voltage_signal",
            "Parameter": "::Root Element::RF_signal::amplitude",
            "Type": "Number"}

    # Define the values to be swept for RF_voltage_signal
    # Generate the vector of values here

    volt_values = [
        0.01, 0.0125893, 0.0158489, 0.0199526, 0.0251189, 0.0316228, 0.0398107, 0.0501187, 0.0630957,
        0.0794328, 0.1, 0.112202, 0.125893, 0.141254, 0.158489, 0.177828, 0.199526, 0.223872, 0.251189,
        0.281838, 0.316228, 0.354813, 0.398107, 0.446684, 0.501187, 0.562341, 0.630957, 0.707946, 0.794328,
        0.891251, 1., 1.12202, 1.25893, 1.41254, 1.58489, 1.77828, 1.99526, 2.23872, 2.51189, 2.81838, 3.16228,
        3.54813, 3.98107, 4.46684, 5.01187, 5.62341
    ]

    for i, value in enumerate(volt_values, start=1):
        para[f"Value_{i}"] = value

    ic.addsweepparameter(rf_volt, para)

    para2 = {"name": "RF_voltage_spur",
             "Parameter": "::Root Element::RF_spur::amplitude",
             "Type": "Number"}

    amplitude_values = [
        0.01, 0.0125893, 0.0158489, 0.0199526, 0.0251189, 0.0316228, 0.0398107, 0.0501187, 0.0630957, 0.0794328, 0.1,
        0.112202, 0.125893, 0.141254, 0.158489, 0.177828, 0.199526, 0.223872, 0.251189, 0.281838, 0.316228, 0.354813,
        0.398107, 0.446684, 0.501187, 0.562341, 0.630957, 0.707946, 0.794328, 0.891251, 1., 1.12202, 1.25893, 1.41254,
        1.58489, 1.77828, 1.99526, 2.23872, 2.51189, 2.81838, 3.16228, 3.54813, 3.98107, 4.46684, 5.01187, 5.62341
    ]

    for i, value in enumerate(amplitude_values, start=1):
        para2[f"Value_{i}"] = value
    ic.addsweepparameter(rf_volt, para2)

    # Define the result outputs
    result_1 = {"Name": "Electrical Spectrum",
                "Result": "::Root Element::ESA_output::signal"}

    result_2 = {"Name": "Optical_Spectrum",
                "Result": "::Root Element::OSA_PD::sum/signal"}

    result_3 = {"Name": "Optical_Power_PD",
                "Result": "::Root Element::PowerMeter_PD::sum/power"}

    # Add result to sweep
    results = [result_1, result_2, result_3]

    for result in results:
        ic.addsweepresult(rf_volt, result)

    ##################################################
    # Insert the outer loop (Vpi sweep)

    #Adding the inner loop
    ic.insertsweep(rf_volt)
    vpi_sweep = "Vpi_sweep"

    ic.setsweep("sweep", "name", vpi_sweep)
    ic.setsweep(vpi_sweep, "type", "Ranges")
    ic.setsweep(vpi_sweep, "Number of points", pts_out)

    # Add the parameters to be swept through
    para11 = {"Name": "RF_DC_bias",
              "Parameter": "::Root Element::MZM_RF::bias voltage 1",
              "Type": "Number",
              "Start": 0.25,
              "Stop": 2.25}
    ic.addsweepparameter(vpi_sweep, para11)

    para22 = {"Name": "Vpi_RF_DC",
              "Parameter": "::Root Element::MZM_RF::pi dc voltage",
              "Type": "Number",
              "Start": 1,
              "Stop": 9}
    ic.addsweepparameter(vpi_sweep, para22)

    para33 = {"Name": "Vpi_RF",
              "Parameter": "::Root Element::MZM_RF::pi rf voltage",
              "Type": "Number",
              "Start": 1,
              "Stop": 9}
    ic.addsweepparameter(vpi_sweep, para33)

    para44 = {"Name": "LO_DC_bias",
              "Parameter": "::Root Element::MZM_LO::bias voltage 1",
              "Type": "Number",
              "Start": 0.25,
              "Stop": 2.25}
    ic.addsweepparameter(vpi_sweep, para44)

    para55 = {"Name": "Vpi_LO_DC",
              "Parameter": "::Root Element::MZM_LO::pi dc voltage",
              "Type": "Number",
              "Start": 1,
              "Stop": 9}
    ic.addsweepparameter(vpi_sweep, para55)

    para66 = {"Name": "Vpi_LO",
              "Parameter": "::Root Element::MZM_LO::pi rf voltage",
              "Type": "Number",
              "Start": 1,
              "Stop": 9}
    ic.addsweepparameter(vpi_sweep, para66)

    para77 = {"Name": "LO_voltage",
              "Parameter": "::Root Element::LO_signal::amplitude",
              "Type": "Number",
              "Start": 0.942 / np.pi * 1,
              "Stop": 0.942 / np.pi * 1}
    ic.addsweepparameter(vpi_sweep, para77)

    #Define the results outputs
    result_11 = {"Name": "Electrical_Spectrum",
                 "Result": "Electrical_Spectrum"}

    result_22 = {"Name": "Optical_Spectrum",
                 "Result": "Optical Spectrum"}

    result_33 = {"Name": "Optical_Power_PD",
                 "Result": "Optical_Power_PD"}

    results2 = [result_11, result_22, result_33]

    #add result to sweep
    for result in results2:
        ic.addsweepresult("Vpi_sweep", result_2)

    # Insert an outer loop (Optical power sweep)
    # Insert the inner sweep loop
    ic.insertsweep(vpi_sweep)

    pts_out3 = 17
    tones2_sweep = "two_tones_sweep"

    ic.setsweep("sweep", "name", tones2_sweep)
    ic.setsweep(tones2_sweep, "type", "Values")
    ic.setsweep(tones2_sweep, "Number of points", pts_out3)

    para111 = {"Name": "P_out",
               "Parameter": "::Root Element::CW_Laser::power",
               "Type": "Number"}

    power_values = [0.001, 0.00125893, 0.00158489, 0.00199526, 0.00251189, 0.00316228, 0.00398107,
                    0.00501187, 0.00630957, 0.00794328, 0.01, 0.0125893, 0.0158489, 0.0199526, 0.0251189,
                    0.0316228, 0.0398107

                    ]

    for i, value in enumerate(power_values, start=1):
        para111[f"Value_{i}"] = value

    ic.addsweepparameter(tones2_sweep, para111)

    # Define the result outputs
    result_111 = {"Name": "Electrical_Spectrum",
                  "Result": "Electrical_Spectrum"}

    result_222 = {"Name": "Optical_Spectrum",
                  "Result": "Optical_Spectrum"}

    result_333 = {"Name": "Optical_Power_PD",
                  "Result": "Optical_Power_PD"}

    result_444 = {"Name": "PowerMeter_SOAin",
                  "Result": "PowerMeter_SOAin"}

    results3 = [result_111, result_222, result_333, result_444]

    #Add result to sweep
    for result in results3:
        ic.addsweepresult(tones2_sweep, result)

    para_values = [para, para11, para22, para33, para44, para55, para66, para77, para111, para2]
    #
    # for value in para_values:
    #     ic.clear(value)
    # for result in results3:
    #     ic.clear(result)

    #Run the sweep
    ic.runsweep(tones2_sweep)
    export_results_2tones(ic)


def export_results_2tones(ic):
    """

    Args:
        ic:

    Returns: The data from the 2 tones sweep

    """
    # Get sweep
    tones2_sweep = "two_tones_sweep"

    data_Optical = ic.getsweepresult(tones2_sweep, "Optical_Spectrum")
    data_Electrical = ic.getsweepresult(tones2_sweep, "Electrical_Spectrum")
    data_PowerPD = ic.getsweepresult(tones2_sweep, "Optical_Power_PD")

    # Create a struct to export the data

    data_export_2tones = {"OS_two_tones_sweep": data_Optical,
                          "ES_two_tones_sweep": data_Electrical,
                          "PowerPD_two_tones_sweep": data_PowerPD}

    savemat("SimulationData_2tones_sweep.mat", data_export_2tones)

    return data_export_2tones


def export_results_icp(file_path):
    """
    To export the data from an .icp file

    Args:
        file_path:

    Returns:

    """
    # Start interconnect session
    ic = lumapi.INTERCONNECT()
    ops.load_icp_file(ic, file_path)
    data_export = export_results_2tones(ic)


def process_data_icp(data_export, f_Optical, R_in=50):
    # Rename the data from the dictionary
    os_2tones_data = data_export["OS_two_tones_sweep"]  # optical spectrum
    es_2tones_data = data_export["ES_two_tones_sweep"]  # electrical spectrum
    P_pd_2tones_data = data_export["PowerPD_two_tones_sweep"]  # data power PD

    # Squeeze the data obtained, because it has multiple dimensions
    os_sweep_data = np.squeeze(data_export["OS_two_tones_sweep"])  # optical spectrum
    es_sweep_data = np.squeeze(data_export["ES_two_tones_sweep"])  # electrical spectrum
    P_pd_sweep_data = np.squeeze(data_export["PowerPD_two_tones_sweep"])  # data power PD

    os_freq_sweep = (c / os_2tones_data["Frequency"]) * 1e6  # Wavelength -> Optical frequency [um]
    es_freq_sweep = (es_2tones_data['frequency'] * 1e-9)  # Electrical frequency [GHz]

    voltageRF = es_2tones_data['RFVolt']  # RF voltage [V]
    powerRF = 10 * np.log10(voltageRF ** 2 / (2 * R_in) * 1e3)  # RF power [dBm]
    replicate_dims = (2, 1, 7, 4)
    powerRF_replicated = np.tile(powerRF, replicate_dims)

    powerLaser = 10 * np.log10(es_2tones_data['P_opt_dBm'])  # Laser power [dBm] , because the laser power is constant
    powerPD = np.squeeze(os_2tones_data['power (dBm)'][0, 0, 0, 0, :]).T  # Power at PD [dBm] #Check the annotation here

    Vpi = np.linspace(2, 8, 7)

    # Transpose is `.T` in NumPy
    data2Export = np.column_stack((powerLaser.T, powerPD.T))

    # Export the file on .csv format
    dataExportLocation = ("/home/maria/PycharmProjects/Solvers/Lumerical-python/Lum-python/interconnect/transmitter_MB"
                          "/data")
    filename = 'Optical Powers_TwoTonesSweep.csv'
    ops.export_csv(data2Export, filename, dataExportLocation)



# After this line 27 of the code