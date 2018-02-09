#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 15:21:51 2018

@author: milesbodmer
"""

from obspy.core import read

#   Read in all of the stations for the first day and plot them
morechannels = read('/Volumes/research/users/mbodmer/TwoTowers_data/Two_Towers/SAC_Data/*/20160614*Z*')
morechannels.decimate(10)
#   pad the stations with zeros so the show up correctly
for i in range(0,morechannels.count()):
    morechannels[i].stats.station = morechannels[i].stats.station.zfill(2)
    
morechannels.plot(size=(1000, 1500))

#   Grab the last day of data and plot a line of stations recording a single 
#   shot at L01

import obspy

#   set the time to look at
dt = obspy.UTCDateTime('2016-06-21T19:00:00.0Z')

linechannels = read('/Volumes/research/users/mbodmer/TwoTowers_data/Two_Towers/SAC_Data/*/20160621*Z*')
linechannels.decimate(10)

for i in range(0,linechannels.count()):
    linechannels[i].stats.station = linechannels[i].stats.station.zfill(2)
    
linechannels.plot(starttime=dt+58, endtime=dt+65, size=(800,1500),equal_scale=False)


#   try to set up the automated triggering

from obspy.signal.trigger import trigger_onset
from obspy.signal.trigger import recursive_sta_lta
from obspy.signal.trigger import plot_trigger
import numpy as np
import matplotlib.pyplot as plt

#   load traces and decimate 
event_trace = read('/Volumes/research/users/mbodmer/TwoTowers_data/Two_Towers/SAC_Data/*/20160621*Z*', starttime=dt+59, endtime=dt+62)
event_trace.decimate(2)

#   plot a spectrogram of the source shot
event_trace[0].spectrogram()

#   pad the station identifier with 0 so they can be sorted
for i in range(0,event_trace.count()):
    event_trace[i].stats.station = event_trace[i].stats.station.zfill(2)
event_trace.sort()    
    
#   plot all traces
event_trace.plot(size=(800,1500), equal_scale=False)


#   create time vector
npts = event_trace[0].stats.npts
df = event_trace[0].stats.sampling_rate
t = np.arange(npts, dtype=np.float32) / df



#   loop over each trace and run sta lta analysis to auto pic the arrival
for i in range(0,linechannels.count()):
    df = event_trace[i].stats.sampling_rate

    cft = recursive_sta_lta(event_trace[i].data, 5, 200)
    #plot_trigger(event_trace[0], cft, 12, 5)

    trig = trigger_onset(cft,np.max(cft)*2/4, np.max(cft)*1/3)
    
    #pick = event_trace[0].stats.starttime + event_trace[0].stats.delta*trig[0,0]
    fig = plt.figure()
    fig.set_size_inches(15, 2)
    plt.xlim([t[trig[0,0]]-1,t[trig[0,0]]+1])
    plt.plot(t,event_trace[i].data,'k')
    plt.axvline(t[trig[0,0]])
    pick[0,i] =  trig[0,0]


#   look at all the shots for station 1 and pick them
shotpick_trace = read('/Volumes/research/users/mbodmer/TwoTowers_data/Two_Towers/SAC_Data/1/20160621*Z*', starttime=dt+35, endtime=dt+68)
shotpick_trace.plot()

df = shotpick_trace[0].stats.sampling_rate

#   sta lta picker
cft = recursive_sta_lta(shotpick_trace[0].data, 5, 200)
plot_trigger(shotpick_trace[0], cft, 25, 5)
trig_picks = trigger_onset(cft,np.max(cft)*3/4, np.max(cft)*1/3)

npts = shotpick_trace[0].stats.npts
t = np.arange(npts, dtype=np.float32) / df

print(t[trig_picks[:,0]])


