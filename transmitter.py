"""
Architecture 1

In here we describe a transmitter system composed by a laser, MZM, SOA, another MZM,
a photodiode, and a low pass filter.

The purpose is to modulate 2 RF frequencies, into an optical signal provided by the
Continuous Laser, modulated at a central wavelength of 1.55um (balanced-single drived - single push pull ).

Balanced single-drive MZMs indeed operate in a push-pull configuration,
where the two arms of the interferometer are driven in opposite directions.
This configuration enhances modulation efficiency, speed, and linearity, making it
suitable for high-speed optical communication applications.
The push-pull design also helps in achieving high dynamic range and suppressing unwanted
distortions, ensuring high-quality signal transmission.

A gain has been included to amplify the signal.
Another MZM, is used with a Modulation signal from the Local Oscillator (LO), because
we want to modulate the signal to an IF signal
The PD will detect the signal
And the LPF, will filer only he signal with the wanted frequency



In this case the simulations will show a transmission that should go up to 0.50, as the component is divided in 2 outputs
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
import V_LO_sweep
import RF_voltage_sweep  #RF voltage sweep script
import optical_power_sweep
import Two_tones_test

# Open an INTERCONNECT session
ic = lumapi.INTERCONNECT()
ic.switchtodesign()
ic.deleteall()  # Clear previous elements
# ic.closeall()

# Constants
c = 3e8  # Speed of light (m/s)
lambda_central = 1.55e-6  # Central wavelength (meters)
lambda_range = 100e-9  # Wavelength range (meters)
lambda_points = 100  # Number of points in the wavelength sweep
NeffTE_waveguide = 3.2  # Effective index
time_window = 1e-7  # s Time for the simulation

impedance = 50  # 50 Ohms
power_laser = 10  # 10 dBm
sample_rate = 300e9  # 300 GHz

V_pi = 4  # Half-wave voltage [V]
V_RF = 0.316  # RF voltage [V]
V_LO = 1.2  # LO voltage [V]
P_opt_dBm = 10  # Input optical power [dBm]
SOA_gain = 15  # Gain dB of the SOA
SOA_NF = 6  # NF [dB] of the SOA
SOA_NFactor = 10 ** (SOA_NF / 10)

f_Signal = 30e9  # Signal frequency [GHz]
f_Spur = 30.1e9  # Spur frequency [GHz]
f_LO = 25  # LO frequency [GHz]
f_Optical = c / lambda_central
RIN_dB = -160  # Laser RIN [dBc/Hz]
BalancedDetection = 0  # 0 -> single-ended detection, 1 -> balanced detection

R_out = impedance  # Output impedance [Ohm]
R_in = impedance  # Input impedance [Ohm]

h = 6.626e-34  # Plank constant
T0 = 290  # System temperature [K]
kB = 1.38e-23  # Boltzman constant [J/K]
q = 1.6e-19  # Elementary electric charge


def transmitter(ic):
    # (f_Optical = f_Optical, time_window=time_window,
    # f_Signal = f_Signal, f_Spur = f_Spur, impedance=impedance, P_opt_dBm = P_opt_dBm,
    # sample_rate=sample_rate):
    """
    Function that sets the architecture of the transmitter. This transmitter does not include
    any specific technology
    Returns: Various graphs

    """
    # General root element
    root = "::Root Element"
    # ic.setnamed(root, 'name', root)
    # Set up simulation to match time window and number of samples
    # The time windows should match because the INTERCONNECT simulation
    # will stop once the simulation time exceeds the time window.
    # The number of samples defines the INTERCONNECT time step, dt, by dt = time_window/(Nsamples+1).
    # The time steps do NOT have to match, although in this example they do. Indeed,
    # the time step of an external simulator can be variable
    ic.setnamed(root, "time window", time_window)
    ic.setnamed(root, "sample rate", sample_rate)
    # ic.setnamed(root, "output signal mode", "sample") 3Block is giving errors
    ic.setnamed(root, "monitor data", "save to memory")
    ic.setnamed(root, "temperature", T0)  #K

    # Add the laser
    ic.addelement("CW Laser")
    laser = "CW_Laser"
    ic.set("name", laser)
    ic.set("x position", 0)
    ic.set("y position", 0)
    ic.set("frequency", f_Optical)
    ic.set("power", P_opt_dBm)  #dBm for it to give me 10 dBm
    ic.set("enable RIN", 1)
    ic.set("RIN", RIN_dB)
    # ic.set("linewidth, 5 ") #MHz
    # ic.set("phase", 0) #rad

    # Add the RF signal
    ic.addelement("Sine Wave")
    rf1 = "RF_signal"
    ic.set("name", rf1)
    ic.set("x position", -100)
    ic.set("y position", -200)
    ic.set("amplitude", 1)
    ic.set("bias", 0)
    ic.set("frequency", f_Signal)
    ic.set("phase", 0)  #rad

    ic.addelement("Sine Wave")
    rf2 = "RF_spur"
    ic.set("name", rf2)
    ic.set("x position", -100)
    ic.set("y position", -400)
    ic.set("amplitude", 1)
    ic.set("bias", 0)
    ic.set("frequency", f_Spur)
    ic.set("phase", 0)  #rad

    # Add the sum of the signals
    ic.addelement("Electrical Adder")
    sum1 = "SUM_1"
    ic.set("name", sum1)
    ic.set("x position", 200)
    ic.set("y position", -300)
    ic.set("run diagnostic", 0)

    # Add the Modulator for the 2 RF signals
    ic.addelement("Mach-Zehnder Modulator")
    mzm1 = "MZM_RF"
    ic.set("name", mzm1)
    ic.set("x position", 200)
    ic.set("y position", 0)
    ic.set("modulator type", "balanced single drive")
    # ic.set("modulator type", "dual drive")
    # ic.set("modulator type", "unbalanced single drive")
    ic.set("dc bias source", "internal")  #internal because I do not need extra components
    ic.set("bias voltage 1", 1)  #Why is it 1? Because is cuadrature and the it divides by 2, and then
    #Lumerical counts it as you would have a bias voltage in 2 arms, and that is why
    #you have to divide it by 2 again
    ic.set("pi dc voltage", V_pi)
    ic.set("pi rf voltage", V_RF)
    # ic.set("extinction ratio", 100) #dB
    ic.set("insertion loss", 5)  #dB
    ic.set("phase shift", 0)  #rad

    # Optical Amplifier
    ic.addelement("Optical Amplifier")
    gain1 = "AMP_1"
    ic.set("name", gain1)
    ic.set("x position", 500)
    ic.set("y position", 0)
    ic.set("gain", SOA_gain)
    ic.set("noise figure", SOA_NF)
    ic.set("enable noise", False)

    # Add the Modulator for the Local Oscillator
    ic.addelement("Mach-Zehnder Modulator")
    mzm_lo = "MZM_LO"
    ic.set("name", mzm_lo)
    ic.set("x position", 700)
    ic.set("y position", 0)
    ic.set("modulator type", "balanced single drive")
    # ic.set("modulator type", "dual drive")
    # ic.set("modulator type", "unbalanced single drive")
    ic.set("dc bias source", "internal")  #internal because I do not need extra components
    ic.set("bias voltage 1", 1)  #Why is it 1? Because is cuadrature and the it divides by 2, and then
    #Lumerical counts it as you would have a bias voltage in 2 arms, and that is why
    #you have to divide it by 2 again
    ic.set("pi dc voltage", V_pi)
    ic.set("pi rf voltage", V_LO)
    # ic.set("extinction ratio", 100) #dB
    ic.set("insertion loss", 5)  #dB
    ic.set("phase shift", 0)  #rad

    # Add the LO signal
    ic.addelement("Sine Wave")
    lo1 = "LO_signal"
    ic.set("name", lo1)
    ic.set("x position", 500)
    ic.set("y position", -300)
    ic.set("amplitude", 1.2)
    ic.set("bias", 0)
    ic.set("frequency", f_LO)
    ic.set("phase", 0)  #rad

    # Add the photodiode
    ic.addelement("PIN Photodetector")
    pd1 = "PIN_1"
    ic.set("name", pd1)
    ic.set("x position", 1000)
    ic.set("y position", 0)
    ic.set("frequency at max power", 1)  #0 meaning false
    # ic.set("frequency", c / lambda_central)
    ic.set("input parameter", "constant")
    ic.set("responsivity", 0.85)  # 0.85 A/W
    ic.set("dark current", 2.5e-8)  #A
    ic.set("enable power saturation", True)
    ic.set("saturation power", 15)  #dBm
    ic.set("enable thermal noise", False)
    ic.set("enable shot noise", False)
    ic.set("convert noise bins", False)
    ic.set("automatic seed", True)

    # Low-pass filter
    ic.addelement("LP Butterworth Filter")
    lpf1 = "LPF_1"
    ic.set("name", lpf1)
    ic.set("x position", 1150)
    ic.set("y position", 0)
    ic.set("order", 3)
    # ic.set("cutoff frequency", 30e9)

    # Add the OSA
    ic.addelement("Optical Spectrum Analyzer")
    osa1 = "OSA_1"
    ic.set("name", osa1)
    ic.set("x position", 400)
    ic.set("y position", -150)
    ic.set("limit frequency range", 0)
    # ic.set("sensitivity", -100) #dBm
    ic.set("limit time range", 0)
    ic.set("resolution", "Gaussian function")
    ic.set("bandwidth", 1e6)

    # Add the Electrical Spectrum Analyzer
    ic.addelement("Spectrum Analyzer")
    rf_sa1 = "RF_SA1"
    ic.set("name", rf_sa1)
    ic.set("x position", 1050)
    ic.set("y position", -200)
    ic.set("limit frequency range", False)
    ic.set("remove dc", 0)
    ic.set("limit time range", 0)
    ic.set("resolution", "rectangular function")
    ic.set("bandwidth", 1e6)
    # ic.set("sensitivity", -130)

    # Power meter
    ic.addelement("Power Meter")
    pw1 = "PWM_1"
    ic.set("name", pw1)
    ic.set("x position", 50)
    ic.set("y position", -130)
    ic.set("input kind", "voltage")
    ic.set("impedance", impedance)
    ic.set("power unit", "dBm")

    # Power meter
    ic.addelement("Power Meter")
    pw_lo = "PWM_LO"
    ic.set("name", pw_lo)
    ic.set("x position", 700)
    ic.set("y position", -330)
    ic.set("input kind", "voltage")
    ic.set("impedance", impedance)
    ic.set("power unit", "dBm")

    # Add the OSA for MZM LO
    ic.addelement("Optical Spectrum Analyzer")
    osa_lo = "OSA_PD"
    ic.set("name", osa_lo)
    ic.set("x position", 900)
    ic.set("y position", -200)
    ic.set("limit frequency range", 0)
    # ic.set("sensitivity", -130) #dBm
    ic.set("limit time range", 0)
    ic.set("resolution", "Gaussian function")
    ic.set("bandwidth", 1e6)

    # Optical power meter for the MZM-RF
    ic.addelement("Optical Power Meter")
    opw_mzm1 = "PowerMeter_SOAin"
    ic.set("name", opw_mzm1)
    ic.set("x position", 200)
    ic.set("y position", 200)

    # Optical power meter for the PD
    ic.addelement("Optical Power Meter")
    opw_pd = "PowerMeter_PD"
    ic.set("name", opw_pd)
    ic.set("x position", 1050)
    ic.set("y position", 200)

    # Add the Electrical Spectrum Analyzer
    ic.addelement("Spectrum Analyzer")
    rf_sa2 = "ESA_output"
    ic.set("name", rf_sa2)
    ic.set("x position", 1200)
    ic.set("y position", -200)
    ic.set("limit frequency range", False)
    ic.set("remove dc", 0)
    ic.set("limit time range", 0)
    ic.set("resolution", "rectangular function")
    ic.set("bandwidth", 1e6)
    # ic.set("sensitivity", -130)

    # Connect the elements
    # MZM1 Connections
    ic.connect(laser, "output", mzm1, "input")
    ic.connect(sum1, "output", mzm1, "modulation 1")
    ic.connect(mzm1, "output", gain1, "input")
    ic.connect(mzm1, "output", opw_mzm1, "input")
    ic.connect(mzm1, "output", osa1, "input")

    # RF connections with adder
    ic.connect(rf1, "output", sum1, "input 2")
    ic.connect(rf2, "output", sum1, "input 1")
    ic.connect(rf1, "output", pw1, "input")

    # MZM_LO Connections
    ic.connect(gain1, "output", mzm_lo, "input")
    ic.connect(mzm_lo, "output", pd1, "input")
    ic.connect(mzm_lo, 'output', osa_lo, "input")
    ic.connect(mzm_lo, "modulation 1", lo1, "output")
    ic.connect(lo1, "output", pw_lo, "input")

    # PD to LPF
    ic.connect(pd1, "output", lpf1, "input")
    ic.connect(lpf1, "output", rf_sa2, "input")
    ic.connect(pd1, "output", rf_sa1, "input")
    ic.connect(mzm_lo, "output", opw_pd, "input")

    # Run the simulation
    ic.run()

    # Retrieve the simulation results
    # Optical Spectrum Analyzer
    # OSA 1
    # P = ic.getresult(osa1, "mode 1/signal")
    # # print(P)
    # p1 = P.get("X power (dBm)")
    # ## Extract the wavelength
    # freq = P['Frequency']
    #
    # # OSA LO
    # P_lo = ic.getresult(osa_lo, "mode 1/signal")
    # # print(P)
    # p1_lo = P_lo.get("X power (dBm)")
    # ## Extract the wavelength
    # freq_lo = P_lo['Frequency']
    #
    # ##Electrical Spectrum Analyzer
    # # RF_SA1
    # P_rfsa = ic.getresult(rf_sa1,"signal")
    # p1_rf = P_rfsa.get("power (dBm)")
    # # P_rfsa = ic.getresult(rf_sa1,"spectrum")
    # freq_rf = P_rfsa['frequency']
    #
    # # RF_SA_LPF
    # P_rfsa_lpf = ic.getresult(rf_sa2,"signal")
    # p1_rf_lpf = P_rfsa_lpf.get("power (dBm)")
    # # P_rfsa = ic.getresult(rf_sa1,"spectrum")
    # freq_rf_lpf= P_rfsa_lpf['frequency']
    #
    # # Power meter
    # P_pm1 = ic.getresult(pw1, "total power")
    # print(f'Power meter pw1 \n', P_pm1, "dBm")
    # # p1_pm = P_pm1.get("")
    # # Power meter
    # P_pm2 = ic.getresult(pw_lo, "total power")
    # print(f'Power meter pw_lo \n', P_pm2, "dBm")
    #
    # # Optical power meter
    # P_opm1 = ic.getresult(opw_mzm1, "sum/power")
    # print(f'Optical Power meter OPM_MZM1  \n', P_opm1, "dBm")
    #
    # P_opm2 = ic.getresult(opw_mzm1, "sum/power")
    # print(f'Optical Power meter OPM_PD  \n', P_opm2, "dBm")
    #
    #
    #
    # # Plot Frequency vs Power (dBm)
    # # OSA 1
    # plt.figure()
    # plt.plot(c / freq * 1e9, p1, label="Input 1 Power")
    # plt.xlabel("Wavelength (nm)")
    # plt.ylabel("Power (dBm)")
    # plt.title("Power vs Wavelength (OSA 1)")
    # plt.legend()
    # plt.grid(True)
    # plt.show()
    #
    # # OSA LO
    # plt.figure()
    # plt.plot(c / freq * 1e9, p1_lo, label="Input 1 Power")
    # plt.xlabel("Wavelength (nm)")
    # plt.ylabel("Power (dBm)")
    # plt.title("Power vs Wavelength (OSA_LO)")
    # plt.legend()
    # plt.grid(True)
    # plt.show()
    #
    #
    # #
    # # Plot Frequency vs Power (dBm)
    # # RF_SA1
    # plt.figure()
    # plt.plot(freq_rf * 1e-9, p1_rf, label="Input rf 1 Power")
    # plt.xlabel("Frequency (GHz)")
    # plt.ylabel("Power (dBm)")
    # plt.title("Power vs Frequency RF (RF_SA1)")
    # plt.legend()
    # plt.grid(True)
    # plt.show()
    #
    # # RF_SA_LPF
    # plt.figure()
    # plt.plot(freq_rf_lpf * 1e-9, p1_rf_lpf, label="Input rf 1 Power")
    # plt.xlabel("Frequency (GHz)")
    # plt.ylabel("Power (dBm)")
    # plt.title("Power vs Frequency RF (RF_SA_LPF)")
    # plt.legend()
    # plt.grid(True)
    # plt.show()


transmitter(ic)
ic.switchtodesign()
# V_LO_sweep.V_LO_sweep(ic)
optical_power_sweep.add_sweep(ic, P_opt_dBm, SOA_gain, SOA_NF, f_Signal, f_Spur, f_Optical,
                              RIN_dB, BalancedDetection)
# RF_voltage_sweep.run_sweep(ic)
# Two_tones_test.add_sweep(ic)


pdb.set_trace()  # Constants
