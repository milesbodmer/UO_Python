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
morechannels.plot(size=(1000, 1500))

#   Grab the last day of data and plot a line of stations recording a single 
#   shot at L01

import obspy

#   set the time to look at
dt = obspy.UTCDateTime('2016-06-21T19:00:00.0Z')

linechannels = read('/Volumes/research/users/mbodmer/TwoTowers_data/Two_Towers/SAC_Data/1*/20160621*Z*')
linechannels.decimate(10)
linechannels[1:11].plot(starttime=dt+58, endtime=dt+65, size=(800,400))