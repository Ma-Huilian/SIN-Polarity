#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created in 2024

@author: Jianlong Yuan and Huilian Ma

Function: This script can be used for data processing 
           and calling PhaseNet(Zhu and Beroza, 2019) to pick P-wave arrival time

 
"""

from obspy.core import read
import fnmatch, os, csv
import obspy
import pandas as pd
import shutil
from obspy import Stream
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
from obspy import read_inventory



#%%---- Key parameters given by user
 
# Directory for downloaded data
src_dir = './2014-10-10-mb4.5-Oklahoma/'
# Event location information
ev_lat, ev_lon, ev_depth_km = 35.9677, -96.7344, 15.9 
# filter parameter
fqmin = 0.5
fqmax = 15
# Time window range, used to capture waveform segments and input them into PhaseNet to pick first-arrival times
winbefore = 30
winafter = 60
# Threshold parameter, used to obtain more reliable first arrival results 
threshold = 0.5



#%%---- Transfer the waveforms to the corresponding folder
#-- folder path
dst_dir1 = './data1'
dst_dir2 = './data2'
dst_dirZ = './dataZ'
os.makedirs(dst_dir1, exist_ok=True)
os.makedirs(dst_dir2, exist_ok=True)
os.makedirs(dst_dirZ, exist_ok=True)

for fname in os.listdir(src_dir):
    if fname.endswith(".mseed") and ("H1" in fname or "HE" in fname):
        src_path = os.path.join(src_dir, fname)
        dst_path1 = os.path.join(dst_dir1, fname)
        shutil.move(src_path, dst_path1)
        
    if fname.endswith(".mseed") and ("H2" in fname or "HN" in fname):
        src_path = os.path.join(src_dir, fname)
        dst_path2 = os.path.join(dst_dir2, fname)
        shutil.move(src_path, dst_path2)
    
    if fname.endswith(".mseed") and ("HZ" in fname):
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

#-- Read data files
dataPathE = dst_dir1 + '/'
dataPathN = dst_dir2 + '/'
dataPathZ = dst_dirZ + '/'
nFilesE = fnmatch.filter( sorted(os.listdir(dataPathE)), '*.mseed') 
nFilesN = fnmatch.filter( sorted(os.listdir(dataPathN)), '*.mseed') 
nFilesZ = fnmatch.filter( sorted(os.listdir(dataPathZ)), '*.mseed') 

#-- Merge three component waveform data
for i in range(len(nFilesE)):

    stE = obspy.read(dataPathE + nFilesE[i])
    stN = obspy.read(dataPathN + nFilesN[i])
    stZ = obspy.read(dataPathZ + nFilesZ[i])
    ENZ = Stream(traces=[stE[0],stN[0],stZ[0]])

    # Remove linear trend, mean, and bandpass filtering
    ENZ.detrend( type='demean')
    ENZ.detrend( type='simple')
    ENZ = ENZ.filter('bandpass', freqmin=fqmin, freqmax=fqmax, corners=4, zerophase=False)  

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
    model = TauPyModel(model= 'ak135' )
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
    
    # Capture the time window for subsequent PhaseNet use
    ENZ = ENZ.trim(p_abs_time-winbefore, p_abs_time+winafter)

    # Write as sac file
    filenameE = str(network)+'.'+str(station)+'.'+str(channelE)+'.sac'
    ENZ[0].write( resultcsv+'/'+filenameE, format="SAC" )
    
    filenameN = str(network)+'.'+str(station)+'.'+str(channelN)+'.sac'
    ENZ[1].write( resultcsv+'/'+filenameN, format="SAC" )
    
    filenameZ = str(network)+'.'+str(station)+'.'+str(channelZ)+'.sac'
    ENZ[2].write( resultcsv+'/'+filenameZ, format="SAC" )
    
    # Write as mseed file
    filename = str(network)+'.'+str(station)+'.mseed'
    ENZ.write( resultcsv+'/'+filename, format="MSEED" )

    # Write component information to CSV file
    with open( outFile, mode='a', newline='' ) as fp:
        writer = csv.writer( fp, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        writer.writerow( [ str(filename), channelE, channelN, channelZ ] ) 



#%%---- PhaseNet        
print('Now processing (PhaseNet picking):')
os.system( f"python phasenet/run.py --mode=pred --model_dir=model/190703-214543\
                --data_dir=dataENZ --data_list=dataENZ/mseedNameList.csv\
                --output_dir=output --batch_size=1 --input_mseed\
                --plot_figure --save_result --input_length=14000" )
     
     
           
#%%---- Select actual P-wave arrival time
#-- Create a CSV file to save arrival time information
outFile2 = './templates.csv'
with open( outFile2, mode='w', newline='' ) as fp2:
    writer = csv.writer( fp2, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    writer.writerow( [ 'network', 'station', 'onsetP' ] )

#-- Read P and S picked up by Phasenet
phaseNetPickedPath = './output/picks.csv'
picks = pd.read_csv(phaseNetPickedPath)
numPicks = len(picks)
print('numPicks =', numPicks)

fname = picks['fname']
itp   = picks['itp']
its   = picks['its']
tp_prob = picks['tp_prob']
ts_prob = picks['ts_prob']

for i in range(numPicks):
  if i >= 0:
    print('Now processing (picking):', i+1, '/', numPicks)
    ifname = str(fname[i])
    iitp   = str(itp[i])
    iits   = str(its[i])
    itp_prob  = str(tp_prob[i])
    its_prob  = str(ts_prob[i])
    # Filter the first arrival of the picking window (i.e. only data with a suffix of * mseed_0)
    if 'mseed_0' in ifname:
        # Get the current event ID
        refNe  = ifname.split(".")[0]
        refSt  = ifname.split(".")[1]
        print('refNe, refSt =', refNe, refSt)
        mseedFileName = resultcsv+"/{0}.{1}.mseed".format( str(refNe), refSt )
        infileENZ = open(mseedFileName)
        ENZ = read(infileENZ.name, debug_headers=True)
        
        orgTimeUTC = ENZ[0].stats.starttime        
        station = ENZ[0].stats.station
        delta= ENZ[0].stats.delta
        
        # Retrieve the initial index value within square brackets
        leftIndexP  = iitp.find("[")
        rightIndexP = iitp.find("]")
        leftIndexS  = iits.find("[")
        rightIndexS = iits.find("]")
        leftIndex_tp_prob  = itp_prob.find("[")
        rightIndex_tp_prob = itp_prob.find("]")
        leftIndex_ts_prob  = its_prob.find("[")
        rightIndex_ts_prob = its_prob.find("]")
        print('leftIndexP, rightIndexP =', leftIndexP, rightIndexP)
        print('leftIndexS, rightIndexS =', leftIndexS, rightIndexS)
    
        # If no P-wave is found (i.e. [] has no content), 
        # it is considered that the signal-to-noise ratio is not high 
        # or the initial arrival is not obvious, 
        # and the waveform data of the current event is not selected
        if leftIndexP+1 == rightIndexP:
            continue
        
        if leftIndexS+1 == rightIndexS:
            iitp  = iitp.split(" ")[0]
            iitp_prob  = itp_prob.split(" ")[0]
            iitp = int(iitp[ leftIndexP + 1:rightIndexP ])
            iitp_prob = float(iitp_prob[ leftIndex_tp_prob + 1:rightIndex_tp_prob ])
            print('iitp =', iitp)
            print('itp_prob =', itp_prob)
            
            # Due to the fact that the initial arrival of picking is generally later, 
            # the results are placed one sampling point in advance
            onsetP = (iitp-1)*delta           
            print('onsetP  =', onsetP)
            
            if iitp_prob < threshold:
                onsetP = 'NA'               
                print('P picking error!\n')
                
            # continue
        
        else:
            # If there are multiple first arrival indexes, 
            # divide each first arrival index by a space symbol, 
            # and then select the first arrival as the P-wave first arrival index
            iitp  = iitp.split(" ")[0]
            iits  = iits.split(" ")[0]
            iitp_prob  = itp_prob.split(" ")[0]
            iits_prob  = its_prob.split(" ")[0]
            iitp = int(iitp[ leftIndexP + 1:rightIndexP ])
            iits = int(iits[ leftIndexS + 1:rightIndexS ]) 
            iitp_prob = float(iitp_prob[ leftIndex_tp_prob + 1:rightIndex_tp_prob ])
            iits_prob = float(iits_prob[ leftIndex_ts_prob + 1:rightIndex_ts_prob ]) 
            print('iitp =', iitp)
            print('itp_prob =', itp_prob)
            
            onsetP = (iitp-1)*delta
            print('onsetP  =', onsetP)
            
            # If the initial arrival of P-wave is not less than the initial arrival of S-wave, 
            # it is considered that there is an error in picking up and 
            # the waveform data of the current event is not selected. 
            #Alternatively, if the probability of picking up P-wave initial arrival 
            # is less than a certain threshold, the station is not selected.
            if iitp >= iits or iitp_prob < threshold:
                onsetP = 'NA'               
                print('P picking error!\n')
                
                # continue
                                    
        # Save the P-wave arrival times of each station to CSV file
        with open( outFile2, mode='a', newline='' ) as fp2:
            writer = csv.writer( fp2, delimiter=',', quoting=csv.QUOTE_MINIMAL)
            writer.writerow( [ str(refNe), str(station), onsetP] )

    
   
