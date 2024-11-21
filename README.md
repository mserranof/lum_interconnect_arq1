# lum_interconnect_arq1
For using this code, first add the lumapi path according to your operative system. In here, I worked with Linux, so this will arise an error if you have Windows (also for the filepaths)

Transmitter architecture simulation w/lumapi &amp; data process

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

For this different codes have been defined
- transmitter.py
- init_constants.py
- V_LO_sweep.py
- optical_power_sweep.py
- RF_voltage_sweep.py
- Two_tones_test.py
- export_results_test.py

## Definition of sweeps
In here various sweeps in interconnect have been defined

### Local Oscillator sweep
In here you can find the functions:
-  V_LO_sweep.py: to define the sweep for the local oscillator voltage

### Optical Power sweep:
In here you can find the functions:
-  add_sweep: To define the sweep
-  export_results_Popt_sweeps: To export the data & save the .mat file
-  export_results_ic: To export the data from a .icp file
-  process_data: to process the data from a .mat file (unfinished)
-  process_data_icp: to process the data from a .icp file (debugging)
-  noise_power_calculations: to calculate the noise (debugging)
-  plot_and_export: to plot graphs and save the sata as .csv
-  plot_and_save_to_pd: to plot graphs directly on a .pdf file

### RF voltage sweep
To vary the RF voltage from a signal, you cna find the functions
-  add_sweep: To define the sweep
-  export_results_ic_RF_sweep: To export the data from a .icp file
