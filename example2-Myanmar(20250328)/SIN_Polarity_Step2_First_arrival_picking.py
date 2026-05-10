#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created in 2025

@author: huilianma

Function: This script can be used for data processing 
           and utilizes STA/LTA(Allen, 1978) to obtain P-wave arrival time

 
"""



import fnmatch, os, csv
import obspy
import shutil
from obspy import Stream
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
from obspy import read_inventory
from obspy.signal.trigger import classic_sta_lta
import numpy as np



#%%---- Key parameters given by user
 
# Directory for storing downloaded data files
src_dir = './2025-03-28-mw7.7-Myanmar/'
# Event location information
ev_lat, ev_lon, ev_depth_km = 22, 95.92, 10  
# filter parameter
fqmin = 0.01
fqmax = 0.2
# Time window range, used to capture waveform segments and use STALTA to pick first-arrival times
winbefore = 30
winafter = 150
# Long term window and short term window range
longlen = 15
shortlen = 3



#%%---- Transfer the waveforms to the corresponding folder
#-- folder path
dst_dir1 = './data1'
dst_dir2 = './data2'
dst_dirZ = './dataZ'
os.makedirs(dst_dir1, exist_ok=True)
os.makedirs(dst_dir2, exist_ok=True)
os.makedirs(dst_dirZ, exist_ok=True)

for fname in os.listdir(src_dir):
    if fname.endswith(".mseed") and ("BH1" in fname or "BHE" in fname):
        src_path = os.path.join(src_dir, fname)
        dst_path1 = os.path.join(dst_dir1, fname)
        shutil.move(src_path, dst_path1)
    
    if fname.endswith(".mseed") and ("BH2" in fname or "BHN" in fname):
        src_path = os.path.join(src_dir, fname)
        dst_path2 = os.path.join(dst_dir2, fname)
        shutil.move(src_path, dst_path2)

    if fname.endswith(".mseed") and ("BHZ" in fname):
        src_path = os.path.join(src_dir, fname)
        dst_pathZ = os.path.join(dst_dirZ, fname)
        shutil.move(src_path, dst_pathZ)
    
print("finished!")



#%%---- Calculate the theoretical time and capture the time window
resultcsv = './dataENZ'
os.makedirs(resultcsv, exist_ok=True)

#-- Create a CSV file to save component information
outFile = str(resultcsv)+'/mseedNameList.csv'
with open( outFile, mode='w', newline='' ) as fp:
    writer = csv.writer( fp, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    writer.writerow( [ 'fname', 'E', 'N', 'Z' ] )
    
#-- Create a CSV file to save arrival time information
outFile2 = './templates.csv'
with open( outFile2, mode='w', newline='' ) as fp2:
    writer = csv.writer( fp2, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    writer.writerow( [ 'network', 'station', 'onsetP' ] )

#-- Read data files
dataPathE = './data1/'
dataPathN = './data2/'
dataPathZ = './dataZ/'
nFilesE = fnmatch.filter( sorted(os.listdir(dataPathE)), '*.mseed') 
nFilesN = fnmatch.filter( sorted(os.listdir(dataPathN)), '*.mseed') 
nFilesZ = fnmatch.filter( sorted(os.listdir(dataPathZ)), '*.mseed') 

#-- Merge three component waveform data
for i in range(len(nFilesE)):

    stE = obspy.read(dataPathE + nFilesE[i])
    stN = obspy.read(dataPathN + nFilesN[i])
    stZ = obspy.read(dataPathZ + nFilesZ[i])
    ENZ = Stream(traces=[stE[0], stN[0], stZ[0]])

    # Remove linear trend, mean, and bandpass filtering
    ENZ.detrend( type='demean')
    ENZ.detrend( type='simple')
    ENZ = ENZ.filter('bandpass', freqmin=fqmin, freqmax=fqmax, corners=4, zerophase=False)  
    
    ENZ_P = ENZ.copy()
    ENZ_ALL = ENZ.copy()
    
    delta = ENZ[0].stats.delta
    starttime = ENZ[0].stats.starttime
    orgTimeUTC = ENZ[0].stats.starttime
    
    network = ENZ[0].stats.network
    station = ENZ[0].stats.station
    channelE = ENZ[0].stats.channel
    channelN = ENZ[1].stats.channel
    channelZ = ENZ[2].stats.channel
    
    # Instrument response file    
    inv = read_inventory(src_dir+network+'.'+station+'.xml')      
    netName = inv.select(network=network)
    staName = netName.select(station=station)[0][0]
    
    # Station location
    st_lat = staName.latitude
    st_lon = staName.longitude    
                           
    # Calculate spherical distance
    gcarc_deg = locations2degrees(ev_lat, ev_lon, st_lat, st_lon)
    
    # Calculating P-wave theoretical travel time based on ak135
    model = TauPyModel(model="ak135")
    arrivals = model.get_travel_times(source_depth_in_km=ev_depth_km,
                                      distance_in_degree=gcarc_deg,
                                      phase_list=["P"])   
    
    if arrivals:
        p_arrival = arrivals[0].time     
        p_abs_time = orgTimeUTC + p_arrival  
        print("gcarc (deg):", gcarc_deg)
        print("P travel time (s):", p_arrival)
        print("Predicted P arrival (UTC):", p_abs_time.isoformat())
    else:
        print("No P arrival found for this geometry/model.")
    
    # Capture the time window for subsequent STA/LTA use
    ENZ1 = ENZ_P.trim(p_abs_time-winbefore, p_abs_time+winafter)
    
    # vertical component
    Zsac = ENZ1[2]
    
    nstaZ = int(shortlen * Zsac.stats.sampling_rate)   
    nltaZ = int(longlen * Zsac.stats.sampling_rate)   
    
    # STA/LTA
    cftZ = classic_sta_lta(Zsac.data, nstaZ, nltaZ)
    imaxZ = np.argmax(cftZ)
    onsetP = imaxZ / Zsac.stats.sampling_rate - 5.5
    
    # Write as sac file
    ENZ_ALL = ENZ_ALL.trim(p_abs_time-winbefore, p_abs_time+winafter)
    
    filenameE = str(network)+'.'+str(station)+'.'+str(channelE)+'.sac'
    ENZ_ALL[0].write(resultcsv+'/'+filenameE, format="SAC")
    
    filenameN = str(network)+'.'+str(station)+'.'+str(channelN)+'.sac'
    ENZ_ALL[1].write(resultcsv+'/'+filenameN, format="SAC")
    
    filenameZ = str(network)+'.'+str(station)+'.'+str(channelZ)+'.sac'
    ENZ_ALL[2].write(resultcsv+'/'+filenameZ, format="SAC")
    
    # Write as mseed file
    filename = str(network)+'.'+str(station)+'.mseed'
    ENZ_ALL.write(resultcsv+'/'+filename, format="MSEED")
    
    # Write component information into CSV file
    with open( outFile, mode='a', newline='' ) as fp:
        writer = csv.writer( fp, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        writer.writerow( [ str(filename), channelE, channelN, channelZ ] ) 
    
    # Save the P-wave arrival times of each station to CSV file
    with open( outFile2, mode='a', newline='' ) as fp2:
        writer = csv.writer( fp2, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        writer.writerow( [ str(network), str(station), onsetP] )


     
                



    
   