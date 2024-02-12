
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pylab
import numpy as np
import pycbc.filter
import pycbc.waveform
from pycbc.psd import welch, interpolate
from pycbc.detector import Detector
from pycbc import frame
from pycbc.catalog import Merger
from pycbc.types.timeseries import FrequencySeries, TimeSeries
from pycbc.filter import highpass_fir, matched_filter


# Load the strain data
strain_h1 = frame.read_frame('data/group4/H1_1359826218_800.gwf',
                             'H1:LOSC-STRAIN')
strain_l1 = frame.read_frame('data/group4/L1_1359826218_800.gwf', 
                             'L1:LOSC-STRAIN')
strain_v1 = frame.read_frame('data/group4/V1_1359826218_800.gwf', 
                             'V1:LOSC-STRAIN')

# Variables
flow = 20.0
delta_f = 1.0 / 800.0
flen = int(4096.0 / delta_f) + 1


def strain_plot(strain):
    # Plot the time series of h(t)
    pylab.plot(strain.sample_times, strain)
    pylab.title('h(t)')
    pylab.xlabel('Time (s)')
    pylab.ylabel('h(t)')
    pylab.show()
    pylab.clf()
    
    
def psd_calculate(name_strain):
    #psd = interpolate(welch(strain_h1), 1.0 / strain_h1.duration)
    psd = []
    if "L" in name_strain :
        psd = pycbc.psd.read.from_txt('psd_data/aLIGO_O4_high_psd.txt', 
                                    length=flen, delta_f=delta_f,
                                    low_freq_cutoff=flow, is_asd_file=False)
    elif "V" in name_strain :
        psd = pycbc.psd.read.from_txt('psd_data/AdV_psd.txt', 
                                    length=flen, delta_f=delta_f,
                                    low_freq_cutoff=flow, is_asd_file=False)
    
    return psd


def psd_plot(name_strain):
    
    psd = psd_calculate(name_strain)
    
    pylab.loglog(psd.sample_frequencies, psd)
    pylab.title('PSD')
    pylab.xlabel('Hz')
    pylab.ylabel('1/Hz')
    pylab.xlim(10, 5e3)
    pylab.show()
    pylab.clf()
    

def snr_calculate(strain, mass, name_strain):
    
    psd = psd_calculate(name_strain)

    # Generate a template to filter with
    hp, hc = pycbc.waveform.get_fd_waveform(approximant="IMRPhenomD", 
                                            mass1=mass, mass2=mass,
                                            f_lower=flow, 
                                            delta_f=delta_f)
    stilde = strain.to_frequencyseries()
    hp.resize(len(stilde))
    
    snr = pycbc.filter.matched_filter(hp, stilde, psd=psd, low_frequency_cutoff=flow)

    snr = abs(snr)
    gps = snr.sample_times
    
    # Find the maximum of the SNR
    arg_max = np.argmax(snr)

    return gps, snr, arg_max


def test_mass(strain, mass_list, name_strain):

    snrmaxtot = 0;
    gpsmaxtot = 0;
    massmaxtot = 0;
    for mass in mass_list:
        snrc = snr_calculate(strain, mass, name_strain)
        snrmax = snrc[1][snrc[2]]
        gpsmax = snrc[0][snrc[2]]
        
        if snrmax > snrmaxtot :
            snrmaxtot = snrmax
            gpsmaxtot = gpsmax
            massmaxtot = mass
        else: continue
        
    return gpsmaxtot, snrmaxtot, massmaxtot


def snr_plot(strain, mass_list, name_strain):
    mass = test_mass(strain, mass_list, name_strain)[2]
    snrc = snr_calculate(strain, mass, name_strain)
    
    pylab.plot(snrc[0], snrc[1])
    pylab.title('SNR')
    pylab.xlabel('GPS Time (s)')
    pylab.ylabel('signal-to-noise')
    pylab.show()
    pylab.clf()

if __name__ == '__main__':
    times = np.zeros(3)#En ordre HLV
    delta_times = np.zeros(3)#
    c = 299792458.
    c_calc = np.zeros(3)
    deltaHL = 3001.e3
    deltaHV = 8181.e3
    deltaLV = deltaHV-deltaHL
    masses = [1.4, 10, 30, 50]
    strains = [strain_h1, strain_l1, strain_v1]
    names = ["LIGO Hanford", "LIGO Livingston", "Virgo"]
    
    for i in range(len(strains)) :
        s = strains[i]
        n = names[i]
        
        print(f"\nDétecteur considéré {n}")
        
        result = test_mass(s, masses, n)
        
        print(f"\nSNR maximal trouvé = {result[1]}")
        print(f"GPS au SNR maximal = {result[0]}")
        print(f"Masse correspondante = {result[2]}\n")   
        times[i]=result[0]        
        answer1 = str(input("Voulez-vous voir les time series de h(t) (y/n) ?"))
        if answer1 == "y":
            strain_plot(s)
        elif answer1 == "n":
            pass
        else: 
            print("ERROR : Not a valid entry !")
            break
        
        answer2 = str(input("Voulez-vous voir la psd (y/n) ?"))
        if answer2 == "y":
            psd_plot(n)
        elif answer2 == "n":
            pass
        else: 
            print("ERROR : Not a valid entry !")
            break
        
        answer3 = str(input("Voulez-vous voir le SNR max (y/n) ?"))
        if answer3 == "y":
            snr_plot(s,masses,n)
        elif answer3 == "n":
            pass
        else: 
            print("ERROR : Not a valid entry !")
            break
    #Ordre HLV
    delta_times[0] = abs(times[0]-times[1])
    delta_times[1] = abs(times[1]-times[2])
    delta_times[2] = abs(times[2]-times[0])
    c_calc[0]=deltaHL/delta_times[0]
    c_calc[1]=deltaLV/delta_times[1]
    c_calc[2]=deltaHV/delta_times[2]
    print("Rapport vitesses HL:" c_calc[0]/c)
    print("Rapport vitesses LV:" c_calc[1]/c)
    print("Rapport vitesses HL:" c_calc[2]/c)
