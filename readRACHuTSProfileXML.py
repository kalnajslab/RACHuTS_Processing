#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri August 2nd 13:25:44 2019

@author: kalnajs

Simple routine to read in Profiler binary data packets from the RACHuTS
profiler and convert to a time stamped CSV file.   The binary format is defined
in RACHUTS_PU_V2.2.ino


"""
import struct
import csv
import os
from datetime import datetime
from numpy import log
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits import mplot3d
#from netCDF4 import Dataset
import glob
from matplotlib.transforms import Transform
from matplotlib.ticker import (
    AutoLocator, AutoMinorLocator)

matplotlib.use('Qt5Agg')

def get_distance(lat_1, lng_1, lat_2, lng_2): 
    d_lat = lat_2 - lat_1
    d_lng = lng_2 - lng_1 

    temp = (  
         np.sin(d_lat / 2) ** 2 
       + np.cos(lat_1) 
       * np.cos(lat_2) 
       * np.sin(d_lng / 2) ** 2
    )

    return 6373.0 * (2 * np.arctan2(np.sqrt(temp), np.sqrt(1 - temp)))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def TSENCalT(T_counts):

    a = 4120.33771
    b = -114.69175
    k = 62941.9437
    
    a0 = -0.001508480319
    a1 = 0.001586963877
    a2 = -0.0002605665956
    a3 = 0.00002627958344
    a4 = -0.000001287017349
    a5 = 2.512381256e-08
    
    R = k*(T_counts-a)/(b-T_counts)
    
    T=1.0/(a0+a1*log(R)+a2*(log(R)**2)+a3*(log(R)**3)+a4*(log(R)**4)+a5*(log(R)**5))
    
    return T

def TSENCalP(Count1, Count2):
    c1 = 41040
    c2 = 37196
    c3 = 24195
    c4 = 23747
    c5 = 32904
    c6 = 28471
    
    dT = Count2 - c5*2**8
    Off = c2*2**16 + (c4*dT/2**7)
    Sense = c1*2**15 + (c3*dT/2**8)
    
    T = (2000 + dT*c6/2**23)/100
    P = (((Count1*Sense/2**21)-Off)/2**15)/100
    
    return T,P

def TMAppend(path):
    
    files = [f for f in glob.glob(path + "*.ready_tm")]
    
    files.sort()
    
    for f in files:
        print(f)
    
    with open(files[0], "rb") as binary_file:
           fileOneData = binary_file.read()  # Read the whole file at once
           start = fileOneData.find(b'START') + 5  # Find the 'START' string that mark the start of the binary section
           end = fileOneData.find(b'END')
           if end == -1 :
               print('END not found in first file, exiting')
               exit()
           else:
               data = fileOneData[:end-2]
               print('Appending to; ' + files[0])
               
    for indx in range(len(files)-1):
        with open(files[indx+1], "rb") as binary_file:
           fileData = binary_file.read()  # Read the whole file at once
           start = fileData.find(b'START') + 5  # Find the 'START' string that mark the start of the binary section
           end = fileData.find(b'END')
           if start == -1 or end == -1:
               print('START and END not found, exiting')
               exit()
           dataToAppend = fileData[start:(end-2)]
           print('Appending ' + str(len(dataToAppend)/30.0) + ' to ouput')
           data = data + dataToAppend
    
    data = data + b'END'
    OutFile = os.path.splitext(files[0])[0] + '.TMA'
    with open(OutFile, "wb") as output:
        output.write(data)
    
    return OutFile

def plotProfile(filename, **kwargs):
    
    PUdata = np.genfromtxt(filename, skip_header = 3, delimiter = ',')
    time = PUdata[:,0]
   
    folder = os.path.basename(os.path.dirname(filename))
    path = os.path.dirname(filename)
    
    eplasedTime = PUdata[:,1]
    alt = PUdata[:,2]
    lat = PUdata[:,3]
    lon = PUdata[:,4]
    FLASH_Flour = PUdata[:,5]
    FLASH_bkg = PUdata[:,6]
    TSEN_T = PUdata[:,7]
    TSEN_P = PUdata[:,8]
    TSEN_3 = PUdata[:,9]
    OPC_300 = PUdata[:,10]
    OPC_500 = PUdata[:,11]
    OPC_700 = PUdata[:,12]
    OPC_1000 = PUdata[:,13]
    OPC_2000 = PUdata[:,14]
    OPC_3000 = PUdata[:,15]
    OPC_5000 = PUdata[:,16]
    OPC_10000 = PUdata[:,17]
    OPC_fps = PUdata[:,18]
    FLASH_T = PUdata[:,19]
    FLASH_pmtV = PUdata[:,20]
    FLASH_lampI = PUdata[:,21]
    FLASH_lampV = PUdata[:,22]
    FLASH_lampT = PUdata[:,23]
    FLASH_V = PUdata[:,24]
    FLASH_uC_T = PUdata[:,25]
    V_Battery = PUdata[:,26]
    V_Pump = PUdata[:,27]
    I_Tsen = PUdata[:,28]
    I_Flash = PUdata[:,29]
    I_Opc = PUdata[:,30]
    T_Pump = PUdata[:,31]
    T_Battery = PUdata[:,32]
    T_Chassis = PUdata[:,33]
    Battery_Heater_State = PUdata[:,34]
    Chassis_Heater_State = PUdata[:,35]
    
    GPSFile = '/Users/kalnajs/Documents/Strateole/Python/Flight/RACHuTS/ST2_C0_03_TTL3_gps.nc'
    Balloon_GPS = Dataset(GPSFile, mode='r')
    Ball_lon = Balloon_GPS.variables['lon'][:] - 360.0
    Ball_lat = Balloon_GPS.variables['lat'][:]
    Ball_alt = Balloon_GPS.variables['alt'][:]
    Ball_time = Balloon_GPS.variables['time'][:] + 1569888000 # convert from seconde since 10/1/2019 to POSIX
    
    start = find_nearest(Ball_time, time[0])
    end = find_nearest(Ball_time, time[-1])
    
    Ball_time= Ball_time[start:end]
    Ball_lon = np.interp(time, Ball_time, Ball_lon[start:end])
    Ball_lat = np.interp(time, Ball_time, Ball_lat[start:end])
    Ball_alt = np.interp(time, Ball_time, Ball_alt[start:end])

    
    print(Ball_lon)

    
    min_alt = np.argmin(alt)
    date_time = datetime.fromtimestamp(float(time[0]))
    d = date_time.strftime("%m/%d/%Y, %H:%M:%S")
    title = 'RACHuTS Profile ' + d+ ' ' + str(lat[0]) + '$^\circ$S, '+ str(lon[0]) + '$^\circ$E'
        
    plt.close('all')
    if 'vsTime' in kwargs:
        
        fps = OPC_fps
        plt.plot(time,smooth(OPC_300/fps,10),'r-')
        plt.plot(time[OPC_500 < 9990],OPC_500[OPC_500 < 9990]/fps[OPC_500 < 9990],'b-', label = '0.5um')
        plt.plot(time[OPC_700 < 9990],OPC_700[OPC_700 < 9990]/fps[OPC_700 < 9990],'g-', label = '0.7um')
        plt.plot(time[OPC_1000 < 9990],OPC_1000[OPC_1000 < 9990]/fps[OPC_1000 < 9990],'c-', label = '1.0um')
        plt.plot(time[OPC_2000 < 9990],OPC_2000[OPC_2000 < 9990]/fps[OPC_2000 < 9990],'m-', label = '2.0um')
        plt.plot(time[OPC_3000 < 9990],OPC_3000[OPC_3000 < 9990]/fps[OPC_3000 < 9990],'y-', label = '3.0um')
        plt.plot(time[OPC_5000 < 9990],OPC_5000[OPC_5000 < 9990]/fps[OPC_5000 < 9990],'k-', label = '5.0um')
        plt.ylabel('Concentration #/cc')
        plt.xlabel('Time [S]')
        plt.yscale('log')
        plt.legend(loc = 'upper right', fontsize='xx-small')
                   
        fig1, ax1 = plt.subplots()
        ax3 = plt.subplot2grid((2,2), (0, 0))
        ax3.plot(time[T_Pump < 9990],T_Pump[T_Pump < 9990], 'r-', label = 'Pump')
        ax3.plot(time[T_Battery < 9990],T_Battery[T_Battery < 9990], 'b-', label = 'Battery')
        ax3.plot(time[T_Chassis < 9990],T_Chassis[T_Chassis < 9990], 'g-', label = 'Chassis')
        ax3.legend(loc='upper left')
        ax3.set_ylabel('Temperature [C]')
        ax3.set_xlabel('Elapsed Time [s]')
        
        ax4 = plt.subplot2grid((2,2), (0, 1))
        ax4.plot(time[I_Opc < 9990],I_Opc[I_Opc < 9990], 'r-', label = 'ROPC')
        ax4.plot(time[I_Flash < 9990],I_Flash[I_Flash < 9990], 'b-', label = 'FLASHB')
        ax4.plot(time[I_Tsen < 9990],I_Tsen[I_Tsen < 9990], 'g-', label = 'TSEN')
        ax4.legend(loc='upper left')
        ax4.set_ylabel('Current [A]')
        ax4.set_xlabel('Elapsed Time [s]')
        
        ax5 = plt.subplot2grid((2,2), (1, 0))
        ax5.plot(time[V_Battery < 9990],V_Battery[V_Battery < 9990], 'g-', label = 'Battery')
        ax5.plot(time[V_Pump < 9990],V_Pump[V_Pump < 9990], 'm-', label = 'Pump')
        ax5.plot(time[FLASH_V < 9990],FLASH_V[FLASH_V < 9990], 'c-', label = 'FLASHB')
        ax5.legend(loc='upper left')
        ax5.set_ylabel('Voltage [V]')
        ax5.set_xlabel('Elapsed Time [s]')
        
        ax6 = plt.subplot2grid((2,2), (1, 1))
        ax6.plot(time, TSEN_T, 'r-')
        ax7 = ax6.twinx() 
        ax7.plot(time, TSEN_P, 'b-')
        ax6.set_ylabel('TSEN T and P')
        ax6.set_xlabel('Pressure [hPa]')
    
    if 'GPS' in kwargs: 
        fig2 = plt.figure( figsize = (9,9))
        fig2.suptitle (title,size = 'large', weight = 'bold')
        ax5 = plt.axes(projection="3d")
        ax5.plot3D(lat, lon, alt/1000.0, 'blue', label = 'bottom dwell + ascent')
        ax5.plot3D(lat[:min_alt], lon[:min_alt], alt[:min_alt]/1000.0, 'red',label = 'flight level + descent' )
        ax5.plot3D(Ball_lat,Ball_lon,Ball_alt/1000.0, 'green', label = 'Balloon')
        ax5.set_xlabel('Latitude [deg]')
        ax5.set_ylabel('Longitude [deg]')
        ax5.set_zlabel('Altitude [km]')
        ax5.set_title(title,size = 'large', weight = 'bold')
        ax5.legend(loc = 'lower right')
    
    if 'vsAlt' in kwargs:
    
        #fix vertical scale limits
        alt_min = 17000
        alt_max = 19200 
                  
        fig3, (ax1,ax3,ax2) = plt.subplots(1,3,figsize = (12,9), sharey = True)
        fig3.suptitle (title,size = 'large', weight = 'bold')
        ax1.plot(TSEN_T[:min_alt], alt[:min_alt], 'r.',label = 'descent')
        ax1.plot(TSEN_T[min_alt:], alt[min_alt:], 'b.',label = 'ascent')
        ax1.set_xlabel('Air Temperature [K]',size='large')
        ax1.set_ylabel('Altitude [m]',size='large')
        ax1.grid('on')
        ax1.set_ylim(alt_min, alt_max)
        ax1.legend(loc = 'lower right')
        ax1.set_title('Air Temperature',size='large')
        locs = ax1.get_yticks()
        
        ax2.plot(smooth(FLASH_Flour[:min_alt]*0.012,10), alt[:min_alt], 'r.', label = 'descent')
        ax2.plot(smooth(FLASH_Flour[min_alt:]*0.012,10), alt[min_alt:], 'b.',label = 'ascent')
        ax2.set_xlim([1,8])
        ax2.set_xlabel('Water Vapor [ppmv]',size='large')
        ax2.grid('on')
        ax2.set_ylim(alt_min, alt_max)
        ax2.legend(loc = 'lower right')
        ax2.set_title('Water Vapor',size='large')
     
        fps = OPC_fps
        #fig4, ax4 = plt.subplots(figsize = (9,9))
        ax3.plot(smooth(OPC_300[:min_alt]/fps[:min_alt],10)*1000, alt[:min_alt],'r.', label = 'd>0.3um desc')
        ax3.plot(smooth(OPC_300[min_alt:]/fps[min_alt:],10)*1000, alt[min_alt:],'b.', label = 'd>0.3um asc')
#        ax3.plot(smooth(OPC_500[OPC_500 < 9990]/fps[OPC_500 < 9990],10),alt[OPC_500 < 9990],'g.', label = '0.5um')
#        ax3..plot(OPC_700[OPC_700 < 9990]/fps[OPC_700 < 9990],alt[OPC_700 < 9990],'g-', label = '0.7um')
#        ax3.plot(smooth(OPC_1000[OPC_1000 < 9990]/fps[OPC_1000 < 9990],10),alt[OPC_1000 < 9990],'c.', label = '1.0um')
#        ax3.plot(smooth(OPC_2000[OPC_2000 < 9990]/fps[OPC_2000 < 9990],10),alt[OPC_2000 < 9990],'m.', label = '2.0um')
        ax3.plot(smooth(OPC_3000[:min_alt]/fps[:min_alt],1)*1000, alt[:min_alt],'r*', markersize = 10, label = 'd>3.0um desc')
        ax3.plot(smooth(OPC_3000[min_alt:]/fps[min_alt:],1)*1000, alt[min_alt:],'b*', markersize = 10, label = 'd>3.0um asc')
#        ax3.plot(smooth(OPC_5000[OPC_5000 < 9990]/fps[OPC_5000 < 9990],10),alt[OPC_5000 < 9990],'g.', label = '5.0um')
        ax3.set_xlabel('Aerosol/Cloud Concentration [#/L]',size='large')
        ax3.set_xlim([10,10000])
        ax3.grid('on')
        ax3.set_ylim(alt_min, alt_max)
        ax3.set_xscale('log')
        ax3.legend(loc = 'lower left')
        ax3.set_title('Aerosol Concentration',size='large')

        #Uncomment for vertical velocity plot        
#        ax4.plot(smooth(np.diff(alt[:min_alt]),30), alt[1:min_alt], 'r.', label = 'descent')
#        ax4.plot(smooth(np.diff(alt[min_alt:]),30), alt[min_alt:-1], 'b.', label = 'ascent')
#        ax4.set_xlim([-2,2])
#        ax4.set_xlabel('Vertical Velocity [m/s]')
#        #ax2.set_ylabel('Altitude [m]')
#        ax4.grid('on')
#        ax4.set_yticklabels([])
#        #ax2.set_title(os.path.basename(filename))
#        ax2.legend(loc = 'lower right')
#        ax4.set_title('Vertical Velocity')
    
        plt.savefig(path + '/' + folder + '.png')
    if 'vsTheta' in kwargs:
       
        theta_min = 380
        theta_max = 420    
        theta = TSEN_T * (1000.0/TSEN_P)**0.286
    
        fig4, (ax1,ax3,ax2) = plt.subplots(1,3,figsize = (12,9), sharey = True)
        fig4.suptitle (title,size = 'large', weight = 'bold')
        ax1.plot(TSEN_T[:min_alt], theta[:min_alt], 'r.',label = 'descent')
        ax1.plot(TSEN_T[min_alt:], theta[min_alt:], 'b.',label = 'ascent')
        ax1.set_xlabel('Air Temperature [K]',size='large')
        ax1.set_ylabel('Theta [K]',size='large')
        ax1.grid('on')
        ax1.set_ylim(theta_min, theta_max)
        ax1.legend(loc = 'lower right')
        ax1.set_title('Air Temperature',size='large')
        locs = ax1.get_yticks()
        
        ax2.plot(smooth(FLASH_Flour[:min_alt]*0.012,10), theta[:min_alt], 'r.', label = 'descent')
        ax2.plot(smooth(FLASH_Flour[min_alt:]*0.012,10), theta[min_alt:], 'b.',label = 'ascent')
        ax2.set_xlim([1,8])
        ax2.set_xlabel('Water Vapor [ppmv]',size='large')
        ax2.grid('on')
        ax2.set_ylim(theta_min, theta_max)
        ax2.legend(loc = 'lower right')
        ax2.set_title('Water Vapor',size='large')
        ax2twin = ax2.twinx()
        ax2twin.set_yticklabels([17.0,18.0,18.2,18.30,18.5,19.0])
        ax2twin.set_ylabel('Approx. Altitude [km]',size='large')
        
        fps = OPC_fps
        ax3.plot(smooth(OPC_300[:min_alt]/fps[:min_alt],10)*1000, theta[:min_alt],'r.', label = 'd>0.3um desc')
        ax3.plot(smooth(OPC_300[min_alt:]/fps[min_alt:],10)*1000, theta[min_alt:],'b.', label = 'd>0.3um asc')
#        ax3.plot(smooth(OPC_500[OPC_500 < 9990]/fps[OPC_500 < 9990],10),theta[OPC_500 < 9990],'g.', label = '0.5um')
#        ax3.plot(OPC_700[OPC_700 < 9990]/fps[OPC_700 < 9990], theta[OPC_700 < 9990],'g-', label = '0.7um')
#        ax3.plot(smooth(OPC_1000[OPC_1000 < 9990]/fps[OPC_1000 < 9990],10), theta[OPC_1000 < 9990],'c.', label = '1.0um')
#        ax3.plot(smooth(OPC_2000[OPC_2000 < 9990]/fps[OPC_2000 < 9990],10), theta[OPC_2000 < 9990],'m.', label = '2.0um')
        ax3.plot(smooth(OPC_3000[:min_alt]/fps[:min_alt],1)*1000, theta[:min_alt],'r*', markersize = 10, label = 'd>3.0um desc')
        ax3.plot(smooth(OPC_3000[min_alt:]/fps[min_alt:],1)*1000, theta[min_alt:],'b*', markersize = 10, label = 'd>3.0um asc')
#        ax3.plot(smooth(OPC_5000[OPC_5000 < 9990]/fps[OPC_5000 < 9990],10),theta[OPC_5000 < 9990],'g.', label = '5.0um')
        ax3.set_xlabel('Aerosol/Cloud Concentration [#/L]',size='large')
        ax3.set_xlim([10,10000])
        ax3.grid('on')
        ax3.set_ylim(theta_min, theta_max)
        ax3.set_xscale('log')
        ax3.legend(loc = 'lower left')
        ax3.set_title('Aerosol Concentration',size='large')
        
    if 'HK' in kwargs:
        alt_min = 17000
        alt_max = 19200 
        fig5, (ax1,ax2,ax3) = plt.subplots(1,3,figsize = (12,9), sharey = True)
        fig5.suptitle (title,size = 'large', weight = 'bold')
        
        ax1.plot(V_Battery, alt, 'r.',label = 'Battery')
        ax1.plot(V_Pump, alt, 'b.',label = 'Pump')
        ax1.set_xlabel('Battery Voltage [V]',size='large')
        ax1.set_ylabel('Altitude [m]',size='large')
        ax1.grid('on')
        ax1.set_ylim(alt_min, alt_max)
        ax1.set_xlim([11.5,12.5])
        ax1.legend(loc = 'lower right')
        ax1.set_title('Battery Voltage',size='large')
        
        ax2.plot(T_Battery, alt, 'r.',label = 'Battery')
        ax2.plot(T_Chassis, alt, 'b.',label = 'Chassis')
        ax2.plot(T_Pump, alt, 'g.',label = 'OPC Pump')
        ax2.plot(FLASH_T, alt, 'c.',label = 'FLASH')
        ax2.set_xlabel('Temperature [C]',size='large')
        ax2.set_ylabel('Altitude [m]',size='large')
        ax2.grid('on')
        ax2.set_ylim(alt_min, alt_max)
        ax2.set_xlim([-15,10])
        ax2.legend(loc = 'lower right')
        ax2.set_title('Internal Temperatures',size='large')
        
        ax3.plot(FLASH_pmtV/100.0, alt, 'r.',label = 'PMT_V/100')
        ax3.plot(FLASH_lampI, alt, 'b.',label = 'Lamp I')
        ax3.plot(FLASH_lampV, alt, 'g.',label = 'Lamp V')
        ax3.plot(FLASH_lampT, alt, 'c.',label = 'LampT')
        ax3.set_xlabel('Temperature [C]',size='large')
        ax3.set_ylabel('Altitude [m]',size='large')
        ax3.grid('on')
        ax3.set_ylim(alt_min, alt_max)
        ax3.set_xlim(-30,20)
        ax3.legend(loc = 'lower right')
        ax3.set_title('FLASH HK',size='large')
        
 
        plt.savefig(path + '/' + folder + '_HK.png')
           

def parseFirstProfileDatatoCSV(InputFile):    
    
    OutFile = os.path.splitext(InputFile)[0] + '.csv'
    with open(InputFile, "rb") as binary_file:
        data = binary_file.read()  # Read the whole file at once
    
    start = data.find(b'START') + 5  # Find the 'START' string that mark the start of the binary section
    end = data.find(b'END')
    
    binData = data[start:end]
    
    ROPC_flow = 3.18 * (30.5 + 273.15)/273.15 * 1013.0/1006 #Voulme liters per minute measured in lab using a mass flow 
    #meter and adjusted from SLPM to volume using the ambient P and T.
    
    with open(OutFile, mode='w') as out_file:
        file_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
        #read the header
        startTime = "%.1f"%(struct.unpack_from('<I',binData,0)[0])
        date_time = datetime.fromtimestamp(float(startTime))
        d = date_time.strftime("%m/%d/%Y, %H:%M:%S")
        GPSStartLat = "%.6f"%(struct.unpack_from('<f',binData,4)[0])
        GPSStartLon = "%.6f"%(struct.unpack_from('<f',binData,8)[0])
        ZephyrAlt = "%.1f"%(struct.unpack_from('<H',binData,12)[0])
        ZephyrLat = "%.6f"%(struct.unpack_from('<f',binData,14)[0])
        ZephyrLon = "%.6f"%(struct.unpack_from('<f',binData,18)[0])
        V_3v3 = "%.3f"%(struct.unpack_from('<H',binData,22)[0] / 1000.0) #3v3 supply in volts
        V_Input = "%.3f"%(struct.unpack_from('<H',binData,24)[0] / 1000.0) #V from PIB supply in volts
        I_Charge = "%.3f"%(struct.unpack_from('<H',binData,26)[0] / 1000.0) #3v3 supply in volts
        
        #write the header data to csv
        header = ['Profile Start Time', 'Time Stamp', 'Initial Lat', 'Initial Lon', 'Zephyr Alt', 'Zephyr Lat', 'Zephyr Lon', '3.3V supply', 
                  'Input V', 'I Charge', 'ROPC Lab Flow [VLPM]']
        file_writer.writerow(header)
        header = [d,startTime,GPSStartLat, GPSStartLon, ZephyrAlt, ZephyrLat, ZephyrLon, V_3v3, V_Input, I_Charge, ROPC_flow]
        file_writer.writerow(header)
        
        header = (['Time POSIX','Time Elapsed [S]','Altitude [m]','Latitude','Longitude','FLASH_flour',
                                 'FLASH_bkg','TSEN Temp [K]','TSEN P [hPa]','TSEN P Sensor T [C]','OPC_300 [cummulative counts]','OPC_500 [cummulative counts]','OPC_700 [cummulative counts]',
                                 'OPC_1000[cummulative counts]', 'OPC_2000[cummulative counts]', 'OPC_3000[cummulative counts]', 'OPC_5000 [cummulative counts]', 'OPC_10000 [cummulative counts]', 'OPC cm3 per sample',
                                 'FLASH_T [C]','FLASH_pmtV [V]','FLASH_lampI [mA]','FLASH_lampV [V]','FLASH_lampT [C]',
                                 'FLASH_V','FLASH_uC_T','V_Battery','V_Pump','I_Tsen','I_Flash','I_Opc',
                                 'T_Pump [C]','Battery_T [C]','Chassis_T [C]', 'Battery_Heater_Status [1=ON]','Heater2_Status [1=ON]'])
        file_writer.writerow(header)
        
        
        Time = startTime
        ElapsedTime = 9999.9
        TSEN_1 = -75.0
        FLASH_T = 9999.9
        FLASH_pmtV = 9999.9
        FLASH_lampI = 9999.9
        FLASH_lampV = 9999.9
        FLASH_lampT = 9999.9
        FLASHuC_T = 9999.9
        FLASH_V = 9999.9
        OPC_500 = 9999
        OPC_700 = 9999
        OPC_1000 = 9999
        OPC_2000 = 9999
        OPC_3000 = 9999
        OPC_5000 = 9999
        OPC_10000 = 9999
        OPC_fps = 9999
        V_Battery = 9999.9
        V_Pump = 9999.9
        I_Tsen = 9999.9
        I_Flash = 9999.9
        I_Opc = 9999.9
        T_Pump = 10.0 #set this to something reasonable as it is used in ROPC volume flow calcultaion
        Heater1_T = 9999.9
        Heater2_T = 9999.9
        
        
        for y in range(int(len(binData)/30)-1):
           indx = (y+1)*30
           
           PriorSampleTime = float(Time)
           ElapsedTime = struct.unpack_from('<H',binData,indx)[0] 
           Time = ElapsedTime + float(startTime)  
           
           #Lat/long are differences from start position *10000 - 0.0001 degree resolution
           Latitude = "%.6f"%(struct.unpack_from('<h',binData,indx + 2)[0]/50000.0 + float(GPSStartLat))
           Longitude = "%.6f"%(struct.unpack_from('<h',binData,indx + 4)[0]/50000.0 + float(GPSStartLon))
           #GPS Altitude in meters
           Altitude = "%.1f"%(struct.unpack_from('<H',binData,indx + 6)[0]) 
           
           #FLASH-B Science Data
           FLASH_flour = struct.unpack_from('<H',binData,indx + 8)[0]
           FLASH_bkg = struct.unpack_from('<H',binData,indx + 10)[0]
           
           #Calculate the OPC flow per sampel based on pump T, sample time and lab flow - need to update to use TSEN ambient 
           OPC_fps = ROPC_flow * 1000.0 / 60.0 * (273.15-75.0)/(273.15+float(T_Pump)*(Time - PriorSampleTime))
           OPC_300 = "%.4f"%(struct.unpack_from('<H',binData,indx + 12)[0])
           
           #TSEN Science data as Unsigned 16 bit int
           TSEN_1 = struct.unpack_from('<H',binData,indx + 14)[0]
           TSEN_2 = struct.unpack_from('<H',binData,indx + 16)[0]
           TSEN_3 = struct.unpack_from('<H',binData,indx + 18)[0]
           TSEN2_MSB = struct.unpack_from('B',binData,indx + 21)[0]
           TSEN3_MSB = struct.unpack_from('B',binData,indx + 20)[0]
           TSEN_T_P,TSEN_P = TSENCalP(TSEN_2+TSEN2_MSB*65535, TSEN_3 + TSEN3_MSB*65535)
           TSEN_T = TSENCalT(TSEN_1)
           
           TSEN_2 = TSEN_2 + TSEN2_MSB*6535
           TSEN_3 = TSEN_3 + TSEN3_MSB*6535

           #Enumeration that tells which round-robbin house keeping variables we have
           HK_indx = struct.unpack_from('<H',binData,indx + 22)[0] // 255
           Heater1_Status = struct.unpack_from('<H',binData,indx + 22)[0] % 2
           Heater2_Status = struct.unpack_from('<H',binData,indx + 22)[0] & 0x0004
           
           #Depending on the HK enumeration the HK values are different.  
           #Each value is sent every 8 samples, including the larger OPC channels.
           if HK_indx == 0:
               FLASH_T = "%.2f"%(struct.unpack_from('<H',binData,indx + 24)[0]/100.0 - 273.15)
               OPC_500 = "%.4f"%(struct.unpack_from('<H',binData,indx + 26)[0]/8)
               V_Battery = "%.3f"%(struct.unpack_from('<H',binData,indx + 28)[0]/1000.0)
           elif HK_indx == 1:
               FLASH_pmtV = "%.3f"%(struct.unpack_from('<H',binData,indx + 24)[0]/10.0)
               OPC_700 = "%.4f"%(struct.unpack_from('<H',binData,indx + 26)[0]/8)
               V_Pump = "%.3f"%(struct.unpack_from('<H',binData,indx + 28)[0]/1000.0)
           elif HK_indx == 2:
               FLASH_lampI = "%.2f"%(struct.unpack_from('<H',binData,indx + 24)[0]/10.0)
               OPC_1000 = "%.4f"%(struct.unpack_from('<H',binData,indx + 26)[0]/8)
               I_Tsen = "%.2f"%(struct.unpack_from('<H',binData,indx + 28)[0]/1000.0)
           elif HK_indx == 3:
               FLASH_lampV = "%.2f"%(struct.unpack_from('<H',binData,indx + 24)[0]/10.0)
               OPC_2000 = "%.4f"%(struct.unpack_from('<H',binData,indx + 26)[0]/8)
               I_Flash = "%.2f"%(struct.unpack_from('<H',binData,indx + 28)[0]/1000.0)
           elif HK_indx == 4:
               FLASH_lampT = "%.2f"%(struct.unpack_from('<H',binData,indx + 24)[0]/100.0 - 273.15)
               OPC_3000 = "%.4f"%(struct.unpack_from('<H',binData,indx + 26)[0]/8)
               I_Opc = "%.2f"%(struct.unpack_from('<H',binData,indx + 28)[0]/1000.0)
           elif HK_indx == 5:
               FLASH_V = "%.2f"%(struct.unpack_from('<H',binData,indx + 24)[0]/1000.0)
               OPC_5000 = "%.4f"%(struct.unpack_from('<H',binData,indx + 26)[0]/8)
               T_Pump = "%.2f"%(struct.unpack_from('<H',binData,indx + 28)[0]/100.0 - 273.15)
           elif HK_indx == 6:
               FLASHuC_T = "%.2f"%(struct.unpack_from('<H',binData,indx + 24)[0]/100.0 - 273.15)
               OPC_10000 = "%.4f"%(struct.unpack_from('<H',binData,indx + 26)[0]/8)
               Heater1_T = "%.2f"%(struct.unpack_from('<H',binData,indx + 28)[0]/100.0 - 273.15)
           elif HK_indx == 7:
               Heater2_T = "%.2f"%(struct.unpack_from('<H',binData,indx + 28)[0]/100.0 - 273.15)
           else:
               print("HK Enumeration out of range")
               
           file_writer.writerow([Time,ElapsedTime,Altitude,Latitude,Longitude,FLASH_flour,
                                 FLASH_bkg,TSEN_T,TSEN_P,TSEN_T_P,OPC_300,OPC_500,OPC_700,
                                 OPC_1000, OPC_2000, OPC_3000, OPC_5000, OPC_10000, OPC_fps,
                                 FLASH_T,FLASH_pmtV,FLASH_lampI, FLASH_lampV, FLASH_lampT,
                                 FLASH_V,FLASHuC_T,V_Battery,V_Pump,I_Tsen, I_Flash, I_Opc,
                                 T_Pump, Heater1_T, Heater2_T,Heater1_Status, Heater2_Status])
    return OutFile




def main():
    
    Path = '/Users/kalnajs/Documents/Strateole/Python/2021/RACHuTS/Test_profile/'
    InputFile = TMAppend(Path)
    CSVFile = parseFirstProfileDatatoCSV(InputFile)
    #plotProfile(CSVFile,vsAlt=True, HK=True, GPS=True)
    
if __name__ == "__main__": 
  # calling main function 
  main() 