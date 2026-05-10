#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created in 2025

@author: huilian Ma

Function: This script implements the generation of corresponding templates 
          based on P-waves at each station to ultimately determine polarity results

    
"""


import pandas as pd
import fnmatch, os, csv, re
from obspy import read
from obspy.geodetics.base import kilometer2degrees
import numpy as np
from obspy.taup import TauPyModel, taup_create
from obspy.core import UTCDateTime,Stream
import matplotlib.pyplot as plt
import math
from scipy.signal import hilbert, find_peaks
import obspy
from obspy.signal.trigger import classic_sta_lta, trigger_onset
import gc
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
plt.rcParams["font.family"] = "Times New Roman"
    
     

#%%---- Key parameters given by user
 
# Cut the time window range of noise and signal segments
windowlength = 20
# The signal-to-noise ratio(SNR) threshold parameter is used to screen data that meets 
# the SNR threshold and further determine polarity
SNRthreshold = 3
# Directory for storing downloaded data files
src_dir = './2025-03-28-mw7.7-Myanmar/'
# Event location information
ev_lat, ev_lon = 22, 95.92



#%%---- Determine polarity
def pol(st, P_pick):
    
    # fig = plt.figure(constrained_layout=True, figsize=(8, 4))
    # ax0 = fig.add_subplot(111)
    
    # waveform data
    Z = st  

    starttime = Z.stats.starttime
    delta = Z.stats.delta
    S_R = Z.stats.sampling_rate
    
    network = Z.stats.network
    station = Z.stats.station
    
    # vertical component
    stZ = Z   
    stZno = stZ.copy() # Copy data to extract noise segment and protect the original data
    stZsi = stZ.copy() # Copy data to extract signal segment and protect the original data
    stZ2 = stZ.copy() # Copy data for drawing, protect original data
    stZ4 = stZ.copy() # Copy data to pick polarity and protect original data
    
    # P-wave arrival
    ENZonsetP = P_pick
            
    # Calculate SNR
    Znosiebeg  = ENZonsetP-windowlength
    Znosieend  = ENZonsetP
    Zsignalbeg = ENZonsetP
    Zsignalend = ENZonsetP+windowlength

    Znoisewave  = stZno.trim( starttime+Znosiebeg, starttime+Znosieend )
    Zsignalwave = stZsi.trim( starttime+Zsignalbeg, starttime+Zsignalend )
    yZ = Zsignalwave.data / max(np.fabs(stZ2.data)) #归一化用于后续计算周期
    
    noiseZ  = np.max( np.fabs( Znoisewave.data ) )
    signalZ = np.max( np.fabs( Zsignalwave.data ) )
    snrZ = round(signalZ/noiseZ, 2)
    
    '''
    # Drawing is used for inspection             
    dataZ = stZ2.data / max(np.fabs(stZ2.data)) # Z分量  
    xZ = np.arange(0, len(dataZ), 1)*delta
                               
    ax0.plot(xZ, dataZ, lw=0.9, label='Z', c='black') 
    ax0.axvline( ENZonsetP, linestyle='--', c='red', lw=1) #phasenet    
    ax0.axvspan( Znosiebeg, Znosieend, alpha=0.1, color='red')
    ax0.axvspan( Zsignalbeg, Zsignalend, alpha=0.1, color='blue')    
    ax0.legend(fontsize=16, loc='upper right')
    
    ax0.text(ENZonsetP-1500*delta+0.05, 0.8, str(station), va="center", fontsize=20)
    ax0.text(ENZonsetP-1500*delta+0.05, 0.6, 'SNR: '+str(snrZ), va="center", fontsize=20)
    
    ax0.set_xlim(ENZonsetP-1500*delta, ENZonsetP+3000*delta)
    ax0.set_ylim(-1, 1)
    
    ax0.set_xlabel("Time(s)", fontsize=20)
    ax0.set_ylabel("Amp.", fontsize=20)
    
    plt.tick_params(axis='x', labelsize=20)  
    plt.tick_params(axis='y', labelsize=20)  
    plt.tight_layout()    
    figurepath = './figure_ENZ/'
    if not os.path.exists(figurepath):
        os.makedirs(figurepath )
    
    fig.savefig(figurepath+str(network)+'_'+str(station)+'.png', dpi=300)
    plt.close(fig)  
    
    ax0.clear()   
    '''
    
    # Pick up the polarity of station data that meets the SNR ratio threshold
    if snrZ >= SNRthreshold :        

        from obspy import read_inventory
        # Read the instrument response file
        inv = read_inventory(src_dir+network+'.'+station+'.xml')  
        
        netName = inv.select(network=network)
        staName = netName.select(station=station)[0][0]
        
        staLat = staName.latitude
        staLon = staName.longitude
        
        # Calculate azimuth and epicenter distance               
        from obspy.geodetics import gps2dist_azimuth, locations2degrees
        dist_m, azimuthEachStation, baz = gps2dist_azimuth(ev_lat, ev_lon, staLat, staLon)
        recDisInKm = dist_m / 1000
        recDisInDeg = locations2degrees(staLat, staLon, ev_lat, ev_lon)
        
        model = TauPyModel(model="ak135")
        arrivals = model.get_ray_paths(10, recDisInDeg, phase_list=["P"])
        
        # Calculate the takeoff_angle
        if arrivals:
            arrival = arrivals[0]
            takeoff_angle = arrival.takeoff_angle
            offAngP = takeoff_angle
        else:
            offAngP = 0
    
        # Use a sine wave with one period as a template, and its polarity is positive
        tempolZ = 1
        
        # Use the maximum absolute amplitude of the noise segment as the threshold for peak filtering
        height_waveZ = Znoisewave.data / max(np.fabs(stZ2.data)) 
        heightZ = np.max( np.fabs( height_waveZ ) )
    
        # Find the first and second peaks from the intercepted signal segment and calculate the period
        try:
            posPeakIdxZ, _ = find_peaks( yZ, heightZ, distance=int(delta/delta) )
            negPeakIdxZ, _ = find_peaks( yZ*-1, heightZ, distance=int(delta/delta) )
            periodZ = 2.0*abs( (posPeakIdxZ[0]-negPeakIdxZ[0])*delta )
            
            # Generate sine wave function  
            frequency = 1 / periodZ  
            t = np.linspace(0, periodZ, int(S_R * periodZ))  
            temWaveZ = 1 * np.sin(2 * np.pi * frequency * t)  
            total_nptsZ = len(temWaveZ)
            
            # Extracting waveform segments based on peak values for subsequent cross-correlation calculations
            if posPeakIdxZ[0] < negPeakIdxZ[0]:
                aZ = [posPeakIdxZ[0], negPeakIdxZ[0]]
            else:
                aZ = [negPeakIdxZ[0], posPeakIdxZ[0]]
            
            BZ = int( (Zsignalbeg+aZ[0]*delta-periodZ/4)/delta )
            EZ = BZ+total_nptsZ
            
            traWaveZ = stZ4.data[BZ: EZ]  
            
            '''
            # Drawing is used for inspection                                      
            xSinZ = (np.arange(0, len(temWaveZ), 1)+(ENZonsetP/delta+aZ[0]-(periodZ/4)/delta))*delta
        
            ax0.plot(xZ, dataZ, lw=0.9, label='Z', c='black') 
            ax0.plot(xSinZ, temWaveZ, lw=0.9, label='Template', c='blue') 
            ax0.axvline( Zsignalbeg+aZ[0]*delta, linestyle='--', c='blue', lw=1) 
            ax0.axvline( Zsignalbeg+aZ[1]*delta, linestyle='--', c='blue', lw=1)     
            ax0.axvline( ENZonsetP, linestyle='--', c='red', lw=1) #phasenet            
            ax0.legend(fontsize=16, loc='upper right')
            
            ax0.set_xlim(ENZonsetP-1500*delta, ENZonsetP+3000*delta)
            ax0.set_ylim(-1, 1)
            
            ax0.text(ENZonsetP-1500*delta+0.02, 0.8, str(station), va="center", fontsize=20)
            ax0.text(ENZonsetP-1500*delta+0.02, 0.6, 'SNR: '+str(snrZ), va="center", fontsize=20)
            '''    
            
            # Calculate the cross-correlation value between the template and P-wave
            bZ = temWaveZ - np.mean(temWaveZ)
            aZ = traWaveZ.data - np.mean(traWaveZ.data)
            stdevZ = (np.sum(aZ**2)) ** 0.5 * (np.sum(bZ**2)) ** 0.5
            
            if stdevZ != 0:
                corrZ = np.sum(aZ*bZ)/stdevZ
            else:
                corrZ = 0
            
            if corrZ >= 0:
                polZ = tempolZ
                pZ = 'U'
            else:
                polZ = tempolZ*(-1) 
                pZ = 'D'
                
            '''
            ax0.text(ENZonsetP-0.05, 0.2, str(pZ), va="center", fontsize=20)
        
            ax0.set_xlabel("Time(s)", fontsize=20)
            ax0.set_ylabel("Amp.", fontsize=20)
            plt.tick_params(axis='x', labelsize=20)  
            plt.tick_params(axis='y', labelsize=20)  
            plt.tight_layout()
            
            corrfigurepath = './figure_polENZ/'
            if not os.path.exists(corrfigurepath):
                os.makedirs(corrfigurepath )

            fig.savefig(corrfigurepath+str(network)+'_'+str(station)+'.png', dpi=300)
            plt.close(fig)

            ax0.clear()   
            '''   
            
        except:
            polZ = 'nan'
            azimuthEachStation=0
            recDisInKm = 0
            offAngP = 0
                                
    else:
        polZ = 'nan'
        azimuthEachStation=0
        recDisInKm = 0
        offAngP = 0
    
    return polZ, round(azimuthEachStation,2), round(recDisInKm,2), round(offAngP,2)
 
       
    
#%%---- Parallel work
def process_single_mseed(iFile, dataPath, pick_dict):

    try:
        # Read waveform
        st_ENZ = obspy.read(os.path.join(dataPath, iFile))
        
        Net = st_ENZ[0].stats.network
        Sta = st_ENZ[0].stats.station
        key = f"{Net}.{Sta}"

        # If there is no corresponding P-wave arrival, skip it
        if key not in pick_dict:
            return None
            
        # vertical component
        st_Z = st_ENZ[2]
        
        # P-wave arrival  
        P_pick = pick_dict[key]
        
        # Call polarity calculation function 'pol'
        pol_Z, az, DisInKm, AngP = pol(st_Z, P_pick)
        
        # Return the result for the main process to write
        return (str(Net), str(Sta), az, DisInKm, AngP, pol_Z)

    except Exception as e:
        print(f"{iFile} : {e}")
        return None



#%%---- Main program
if __name__ == '__main__':
    start_time = time.perf_counter()
    
    # Create a CSV file for storing information
    resultcsv = './'
    outFile = str(resultcsv) + 'StationinformationSTALTA.csv'
    
    with open(outFile, mode='w', newline='') as fp:
        writer = csv.writer(fp, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        writer.writerow( [ 'NetName', 'StName', 'Az(deg)', 'EpDis(km)', 'p_offAng(deg)', 'p_POL(Z)' ] )

    # Read CSV file of first arrival time information   
    filePath = './templates.csv' 
    data = pd.read_csv(filePath)
    net = data['network']
    sta = data['station']
    onsetP = data['onsetP']

    # Read the mseed data file containing three components
    dataPath = './dataENZ/'
    if os.path.exists(dataPath):
        nFiles = fnmatch.filter(sorted(os.listdir(dataPath)), '*.mseed') 
    else:
        print(f"error: catalog {dataPath} is not exist")
        nFiles = []

    # Build a dictionary to store station and arrival information
    pick_dict = {}
    for i in range(len(data)):
        key = f"{net[i]}.{sta[i]}"
        pick_dict[key] = onsetP[i]

    # Parallel computing
    with open(outFile, mode='a', newline='') as fp:
        writer = csv.writer(fp, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        
        # Automatic scheduling, but also the option to set maxw_workers
        with ProcessPoolExecutor() as executor:

            # pick_dict is passed as a parameter to avoid global variable lookup issues
            futures = [
                executor.submit(process_single_mseed, f, dataPath, pick_dict) 
                for f in nFiles ]
            
            # Polarity results
            count = 0
            for future in as_completed(futures):
                result = future.result()
                if result is not None:
                    #---- 获取返回的 (Sta, pol_Z) 并写入
                    row_to_write = result
                    writer.writerow(row_to_write)
                    count += 1
                    if count % 10 == 0:
                        print(f"finished: {count} ")

    # Running program time
    end_time = time.perf_counter()
    duration = end_time - start_time
    
    print("\n" + "="*30)
    print(f"Total time taken: {duration:.2f} s")
    if duration > 60:
        mins = duration // 60
        secs = duration % 60
        print(f"Format Time: {int(mins)} min {secs:.2f} sec")
    print("="*30) 
    
    
                        
                

                
                
                


    







        
    







