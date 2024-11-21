"""
Code that has the values for initialization and constant values that will be needed for the code
"""

def init_constants():
    """

    Returns:
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
        H_pd = 1/2


    """

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

    H_pd = 1/2

    return (P_opt_dBm, SOA_gain, SOA_NF, SOA_NFactor, f_Signal, f_Spur, f_LO, c, lambda_central, f_Optical, RIN_dB,
            BalancedDetection, V_pi, V_RF, V_LO, h, T0, kB, q, power_laser, sample_rate, R_out, R_in, H_pd)