"""
Script to create a sweep for the voltage of the local oscillator

With this script, we are sweeping the power for the input RF powers
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
import pandas as pd


# Open an INTERCONNECT session
# ic = lumapi.INTERCONNECT()
# ic.switchtodesign()

def V_LO_sweep(ic):
    """
    Function to add a sweep for the local oscillator (LO) sweep
    Args:
        ic:

    Returns:

    """
    start_v = 0.1
    final_v = 4.1
    step_v = 0.1
    pts = 40

    sweepName = "V_LO_sweep"
    # ic.deletesweep(sweepName)

    ic.addsweep()

    ic.setsweep("sweep", "name", sweepName)
    ic.setsweep(sweepName, "type", "Values")
    ic.setsweep(sweepName, "number of points", pts)

    para = {}
    para["Name"] = "V_LO"
    para["Parameter"] = "::Root Element::LO_signal::amplitude"
    para["Type"] = "Number"

    values = np.arange(start_v, final_v, step_v).tolist()

    for i, value in enumerate(values, start=1):
        para[f"Value_{i}"] = value

    ic.addsweepparameter(sweepName, para)

    #Define the result outputs
    result_1 = {}
    result_1["Name"] = "Electrical_Spectrum"
    result_1["Result"] = "::Root Element::ESA_output::signal"

    ic.addsweepresult(sweepName, result_1)

    # Run the sweep
    ic.runsweep(sweepName)


    # ic.clear(para, result_1)




