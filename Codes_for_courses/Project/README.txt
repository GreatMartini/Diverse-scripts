Python code for the identification of a gravitational wave signal. The data consist in three detections of a supposed gravitational wave signal, one from LIGO Hanford, one from LIGO Livingstone and one from VIRGO. The objective of this code was to determine if the detected signal was from a real gravitational wave based on test masses and the time between detections.

"strain_plot(strain)" Is a function that plot the strain signal.

"psd_calculate(name_strain)" calculates the power spectral density of the signal and "psd_plot(name_strain)" plots them.

"snr_calculate(strain, mass, name_strain)" calculates the signal-noise ratio and "snr_plot(strain, mass, name_strain)" plots it.

"test_mass" calculates the maximal signal over noise ratio, the GPS at the maximal SNr and the corresponding masses of the sources.