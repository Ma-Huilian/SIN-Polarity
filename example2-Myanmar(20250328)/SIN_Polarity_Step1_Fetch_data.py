#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 14:53:55 2019

@author: jianlongyuan

Function: This script can be used to obtain Miniseed waveform data and
          instrument response XML files from IRIS DMC directly for DSA

reference；
    https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.mass_downloader.html
"""
import os
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.clients.fdsn.mass_downloader import CircularDomain, Restrictions, MassDownloader
    
#-------- Key parameters given by user
#-- Directory defined for storing data 
eventPath = './2025-03-28-mw7.7-Myanmar/'
#-- Event origin time (refer to the report of IRIS, ISC, USGS, etc )
origin_time = UTCDateTime(2025, 3, 28, 6, 20, 54) # from IRIS
#-- epicenter
lat = 22    # from IRIS
lon = 95.92  # from IRIS
#-- epicenter distance range for requesting stations
minDisInDeg, maxDisInDeg = 0, 90#-- unit: Degree
#-- service of the client
clientName = "IRIS"


#-------- In general, the following code does not need to be changed
if not os.path.exists(str(eventPath)):
    os.mkdir(str(eventPath))
#-- fetch the event parameters and save as a xml format file
client      = Client(clientName)
catalog = client.get_events( starttime = origin_time - 60,
                             endtime   = origin_time + 60,
                             latitude  = lat,
                             longitude = lon,
                             minradius = minDisInDeg,
                             maxradius = maxDisInDeg )

#-- print catalog to screen and save it to file in QUAKEML format 
print( catalog )
catalog.write( eventPath+'eventParameters.xml', format="QUAKEML" )
    
#-------- download data
# Circular domain around the epicenter. This will download all data between
# 'minDisInDeg' and 'maxDisInDeg' degrees distance from the epicenter. This module also offers
# rectangular and global domains. More complex domains can be defined by
# inheriting from the Domain class.
domain = CircularDomain(latitude =lat,
                        longitude=lon,
                        minradius=minDisInDeg,
                        maxradius=maxDisInDeg )

restrictions = Restrictions(
    # Get data from 30 minutes before the event to 60 minutes after the
    # event. This defines the temporal bounds of the waveform data.
    starttime=origin_time + 0*60, 
    endtime  =origin_time + 30*60,
    # You might not want to deal with gaps in the data. If this setting is
    # True, any trace with a gap/overlap will be discarded.
    reject_channels_with_gaps = 0,
    # And you might only want waveforms that have data for at least 95 % of
    # the requested time span. Any trace that is shorter than 95 % of the
    # desired total duration will be discarded.
    minimum_length = 0.95,
    # No two stations should be closer than 0 km to each other. This is
    # useful to for example filter out stations that are part of different
    # networks but at the same physical station. Settings this option to
    # zero or None will disable that filtering.
    minimum_interstation_distance_in_m = 0,
    network = '*',
    station = '*',
    channel = 'BH*',
    # channel = 'BH*,HH*,EH*',
    # Location codes are arbitrary and there is no rule as to which
    # location is best. Same logic as for the previous setting.
    location_priorities = ["*"] )

# No specified providers will result in all known ones being queried.
mdl = MassDownloader()
# The data will be downloaded to the 'eventPath'
# folders with automatically chosen file names.
# Control how many threads are used to download data in parallel
# per data center - 3 is a value in agreement with some data centers.
mdl.download(domain, restrictions, threads_per_client=3,
             mseed_storage=eventPath,
             stationxml_storage=eventPath)


