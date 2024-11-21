"""
Script to determine the RF voltage sweep for calculations
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
import optical_power_sweep as ops
from scipy.io import savemat, loadmat

# Open an INTERCONNECT session
# ic = lumapi.INTERCONNECT()
# ic.switchtodesign()

# with open("/home/maria/PycharmProjects/Solvers/Lumerical-python/Lum-python/interconnect/transmitter_MB/transmitter.py") as f:
#  exec (f.read())

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


def run_sweep(ic):
    # Create inner loop
    ic.addsweep()

    # Constants
    pts = 46  #inner loop
    pts_out = 9  #out loop

    # General setting of the sweep
    rf_volt = "RFVolt"
    ic.setsweep("sweep", "name", rf_volt)
    ic.setsweep(rf_volt, "type", "Values")
    ic.setsweep(rf_volt, "number of points", pts)

    # Add the parameters to be swept through
    para = {"Name": rf_volt,
            "Parameter": "::Root Element::RF_signal::amplitude",
            "Type": "Number"}

    #Parameters
    rf_voltages = [0.01, 0.0125893, 0.0158489, 0.0199526, 0.0251189,
                   0.0316228, 0.0398107, 0.0501187, 0.0630957, 0.0794328, 0.1,
                   0.112202, 0.125893, 0.141254, 0.158489, 0.177828, 0.199526, 0.223872,
                   0.251189, 0.281838, 0.316228, 0.354813, 0.398107, 0.446684, 0.501187,
                   0.562341, 0.630957, 0.707946, 0.794328, 0.891251, 1, 1.12202,
                   1.25893, 1.41254, 1.58489, 1.77828, 1.99526, 2.23872, 2.51189, 2.81838,
                   3.16228]

    for i, value in enumerate(rf_voltages, start=1):
        para[f"Value_{i}"] = value

    ic.addsweepparameter(rf_volt, para)

    # Define the results outputs

    result_1 = {"Name": "Electrical_Spectrum",
                "Result": "::Root Element::ESA_output::signal"}

    result_2 = {"Name": "Optical_Spectrum",
                "Result": "::Root Element::OSA_PD::sum/signal"}

    result_3 = {"Name": "Optical_Power_PD",
                "Result": "::Root Element::PowerMeter_PD::sum/power"}

    # Add result to sweep
    results = [result_1, result_2, result_3]

    for i in results:
        ic.addsweepresult(rf_volt, i)

    ########################################################################
    # Insert the outer loop
    ic.insertsweep(rf_volt)

    # Continue with the outer loop
    rf_volt_sweep = "RF_voltage_sweep"
    ic.setsweep("sweep", "name", rf_volt_sweep)
    ic.setsweep(rf_volt_sweep, "type", "Ranges")
    ic.setsweep(rf_volt_sweep, "number of points", pts_out)

    # Add the parameters to be swept through
    para2 = {"Name": "RF_DC_bias",
             "Parameter": "::Root Element::MZM_RF::bias voltage 1",
             "Type": "Number",
             "Start": 0.25,
             "Stop": 2.25}
    ic.addsweepparameter(rf_volt_sweep, para2)

    para3 = {"Name": "Vpi_RF_DC",
             "Parameter": "::Root Element::MZM_RF::pi dc voltage",
             "Type": "Number",
             "Start": 1,
             "Stop": 9}
    ic.addsweepparameter(rf_volt_sweep, para3)

    para4 = {"Name": "Vpi_RF",
             "Parameter": "::Root Element::MZM_RF::pi rf voltage",
             "Type": "Number",
             "Start": 1,
             "Stop": 9}
    ic.addsweepparameter(rf_volt_sweep, para4)

    para5 = {"Name": "LO_DC_bias",
             "Parameter": "::Root Element::MZM_LO::bias voltage 1",
             "Type": "Number",
             "Start": 0.25,
             "Stop": 2.25}
    ic.addsweepparameter(rf_volt_sweep, para5)

    para6 = {"Name": "Vpi_LO_DC",
             "Parameter": "::Root Element::MZM_LO::pi dc voltage",
             "Type": "Number",
             "Start": 1,
             "Stop": 9}
    ic.addsweepparameter(rf_volt_sweep, para6)

    para7 = {"Name": "Vpi_LO",
             "Parameter": "::Root Element::MZM_LO::pi rf voltage",
             "Type": "Number",
             "Start": 1,
             "Stop": 9}
    ic.addsweepparameter(rf_volt_sweep, para7)

    para8 = {"Name": "LO_voltage",
             "Parameter": "::Root Element::LO_signal::amplitude",
             "Type": "Number",
             "Start": 0.942 / np.pi * 1,
             "Stop": 0.942 / np.pi * 9}
    ic.addsweepparameter(rf_volt_sweep, para8)

    # Define the result outputs
    result_11 = {"Name": "Electrical_Spectrum",
                 "Result": "Electrical_Spectrum"}

    result_22 = {"Name": "Optical_Spectrum",
                 "Result": "Optical_Spectrum"}

    result_33 = {"Name": "Optical_Power_PD",
                 "Result": "Optical_Power_PD"}

    results2 = [result_11, result_22, result_33]

    # Add result to sweep
    for result in results2:
        ic.addsweepresult(rf_volt_sweep, result)

    para_values = [para, para2, para3, para4, para5, para6, para7, para8]

    # for value in para_values:
    #     ic.clear(value)
    # for result in results2:
    #     ic.clear(result)

    # Run the sweep
    ic.runsweep(rf_volt_sweep)


def export_results_ic_RF_sweep(ic, V_pi, V_RF, V_LO, P_opt_dBm=P_opt_dBm, SOA_gain=SOA_gain, SOA_NF=SOA_NF,
                            f_Signal=f_Signal,
                            f_Spur=f_Spur, f_LO=f_LO, f_Optical=f_Optical,
                            RIN_dB=RIN_dB, Balanced_detection=BalancedDetection):
    """    Export the results from the sweep of the RF voltage sweep

    Args:
        ic:
        V_pi:
        V_RF:
        V_LO
        P_opt_dBm:
        SOA_gain:
        SOA_NF:
        f_Signal:
        f_Spur:
        f_LO:
        f_Optical:
        RIN_dB:
        Balanced_detection:

    Returns: Results from the optical power sweep and the RF voltage sweep

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
    swpName_1 = "RF_voltage_sweep"

    # RF voltage
    swp_1_data_Optical = ic.getsweepresult(swpName_1, "Optical_Spectrum")
    swp_1_data_Electrical = ic.getsweepresult(swpName_1, "Electrical_Spectrum")
    swp_1_data_PowerPD = ic.getsweepresult(swpName_1, "Optical_Power_PD")

    # # Get the parameter values
    # V_LO_ = "LO_voltage"
    # Vpi_RF_ = "Vpi_RF"
    # # V_RF_ = "RF_voltage"
    #
    # V_LO = ic.getsweepdata(swpName_1, V_LO_)
    # Vpi_RF = ic.getsweepdata(swpName_1, Vpi_RF_)
    # # V_RF = ic.getsweepdata(swpName_1, V_RF_)

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
                   "OS_RF_voltage_sweep": swp_1_data_Optical,
                   "ES_RF_voltage_sweep": swp_1_data_Electrical,
                   "PowerPD_RF_Optical_power_sweep": swp_1_data_PowerPD}

    savemat("SimulationData_RF_sweep.mat", data_export)
    # ic.write("SimulationData_Optical_RF.csv", data_export)

    return data_export


def export_results_ic(file_path):
    """
    To export the data from an .icp file

    Args:
        file_path:

    Returns:

    """
    # Start interconnect session
    ic = lumapi.INTERCONNECT()
    ops.load_icp_file(ic, file_path)
    data_export = export_results_ic_RF_sweep(ic, V_pi, V_RF, V_LO, P_opt_dBm, SOA_gain, SOA_NF,
                            f_Signal,
                            f_Spur, f_LO, f_Optical,
                            RIN_dB, BalancedDetection)
    # ops.trial_export(ic, data_export, "SimulationData_RF.csv")
