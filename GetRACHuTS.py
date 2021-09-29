#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
LPC download and processing script

This is based on Albert's automatic download script, customized for LPC.
The script uses sftp to mirror the LPC directories on the CCMz.  If there is a
a new LPC TM file, it processes that file using the readLPCXML module which generates
a csv file of that TM.   It then adds that data to two global CSV files, one 
containing every measurement from each TM file, the other containing a single 
average measurement from each TM file.  It also reads the state messages from 
XML portion of the TM file and appends those to the LPC_State_Log text file.  

This interface uses Python 3.x, pysftp, gzip, csv and readLPCXML modules.

"""

import datetime as dt
import csv
import os
import glob
import pysftp
import gzip


#https://pysftp.readthedocs.io/en/latest/pysftp.html

# Some constant values

default_local_target_dir="RACHuTS/TMs/" # directory where to store mirrored data on your local machine
Output_dir = "RACHuTS/"
csv_dir = "csv/" # dir where to put processesed csv files 
log_file_suffix = "_RACHuTS_Log.txt" #file to save log of XML messages
ccmz_url="sshstr2.ipsl.polytechnique.fr" # CCMz URL from where to download data
ccmz_user="lkalnajs" # Your login on the CCMz
ccmz_pass="Xm<]5D7j" # Your password on the CCMz
# ID of flights with RACHuTS in 2021 campaign
my_flights=['ST2_C1_04_TTL3','ST2_C1_09_TTL2','ST2_C1_19_TTL3'] #All the flights with RACHuTS

# ID of my instrument
my_instruments=['RACHUTS'] # Adapt according to your needs
flight_or_test='Flight'
tm_or_tc='TM'
raw_or_processed='Processed'


def mirror_ccmz_folder(instrument, ccmz_folder, local_target_dir=default_local_target_dir, show_individual_file=True):
   """
   Mirror one CCMz folder.
   Files are stored locally in local_target_dir/ccmz_path/to/ccmz_folder/
   Files already downloaded are not downloaded again.
   local_target_dir prescribes where CCMz files will be downloaded locally.
   show_individual_file controls whether the name of each downloaded file is displayed or not.
   """

   print('---------------------------------')
   print('Trying to mirror CCMz folder: \033[1m'+ccmz_folder+'\033[0m')
   
   downloaded_files = []

   # Create (if needed) the appropriate local directory
   local_folder=os.path.join(local_target_dir,ccmz_folder)
   if not os.path.exists(local_folder):
      os.makedirs(local_folder)

   # Connect to CCMz
   try:
       with pysftp.Connection(host=ccmz_url, username=ccmz_user, password=ccmz_pass) as sftp:
          print("\033[1mConnection to CCMz succesfully established\033[0m...")

          # Switch to the remote directory
          try:
              sftp.cwd(ccmz_folder)
          except IOError:
              print('\033[1m\033[91mNo such directory on CCMz: '+ccmz_folder+'\033[0m')
              return

          # Get file list in current directory, i.e. those that have been already downloaded from CCMz
          local_files=glob.glob(os.path.join(local_folder,'*')) # filenames with relative path
          local_filenames=[os.path.basename(f) for f in local_files] # filenames without

          # Get file list from the CCMz directory with file attributes
          ccmz_file_list = sftp.listdir_attr()

          # check wether CCMz files need to be downloaded
          n_downloads=0

          for ccmz_file in ccmz_file_list:
              # Get rid of directories in CCMZ folder (if any)
              if ccmz_file.longname[0] == '-':
                 ccmz_filename=ccmz_file.filename
                 # Check whether the file has already been downloaded (so as to not download it again)
                 if not ccmz_filename in local_filenames:
                    # file has to be downloaded
                    if show_individual_file == True:
                       print('Downloading \033[92m'+ccmz_filename+'\033[0m...') # display file name
                       downloaded_files.append(os.path.join(local_folder,ccmz_filename))
                    sftp.get(ccmz_filename,os.path.join(local_folder,ccmz_filename),preserve_mtime=True)
                    #print(ccmz_filename)
                    n_downloads=n_downloads+1
            # Create a ZipFile Object and load sample.zip in it
                    print (os.path.join(local_folder,ccmz_filename))

            
          # and print some statistics
          if n_downloads == 0:
              print('\nYour local repository \033[92m'+local_folder+'\033[0m looks\033[1m up do date\033[0m')
          else:
              print('\n\033[1m'+str(n_downloads)+ '\033[0m file(s) downloaded in \033[92m'+local_folder+'\033[0m')
              print('List of downloaded Files')
              return downloaded_files
              

   except:
          print('\033[1m\033[91mConnection to CCMz failed\033[0m: check your login/password')
          return

def loop_over_flights_and_instruments():
    """
    Get all data from CCMz for the input list of flights/instruments
    """
    for flight in my_flights:
        for instrument in my_instruments:
            ccmz_folder=os.path.join(flight,instrument,flight_or_test,tm_or_tc,raw_or_processed)
            #mirror_ccmz_folder(ccmz_folder)
            new_files = mirror_ccmz_folder(instrument,ccmz_folder, show_individual_file=True)
            if new_files != None:
                if os.path.exists(Output_dir+csv_dir + flight + '/') == False:
                    print('Creating Directory: ' + Output_dir + csv_dir + flight + '/')
                    os.makedirs(Output_dir + csv_dir + flight + '/')
                for file in new_files:
                    if file.endswith('.gz'):
                        path = os.path.dirname(file)
                        filename = os.path.basename(file)
                    
                    if 'RACHUTS' == instrument:
                        InputFile = path + '/' + filename
                        #print('Input Filename: ' + InputFile)
                        try:
                            readHeader(InputFile, Output_dir + flight + log_file_suffix)
                        except:
                            print("Unable to read header from: " + os.path.basename(file))
                        try:
                            #OutputFile = LPC_csv_dir + os.path.splitext(os.path.basename(file))[0]
                            #OutputFile = os.path.splitext(OutputFile)[0] + '.csv'
                            #print('Processing to: ' + OutputFile)
                            #csvFile = parseLCPdatatoCSV(InputFile,OutputFile)
                            #plotLPC(csvFile)
                        except:
                            print('Unable to Process Data From: ' + os.path.basename(file) )
                
                
                #master_csv(LPC_csv_dir + "*.csv",mean_file_name,master_file_name)
                        

def readHeader(InputFile,logFile):    
    
    
    if os.path.exists(logFile) == False:
        with open(logFile, 'w') as csvfile:
            csvwriter = csv.writer(csvfile)
            Header = ['MsgID','date','time','message_type','BinLength','XMLMsg','profile_num','profile_segment','reel_pos','pu_current','pu_battery_v','pu_battery_t','pu_time','initial_lat', 'initial_lon','initial_alt','filename']
            csvwriter.writerow(Header)
    
    with open(InputFile, "rb") as binary_file:
        data = binary_file.read()
    #/Users/kalnajs/Documents/Strateole/Python/2021/LPC_Test/ST2_C0_03_TTL3/LPC/Flight/TM/Processed/ST2_C0_03_TTL3.2447_20200228_02_37_07.LPC.dat.gz
    filename = os.path.basename(InputFile)
    #date = filename.split('_')[4]
    #time = filename.split('_')[5] + ':' + filename.split('_')[6] + ':' + filename.split('_')[7][0:2]
    #for the simulator:
    date = '20'+filename[8:14]
    time = filename[14:20]
    
    start = data.find(b'<StateMess1>')
    end = data.find(b'</StateMess1>')
    XMLMsg = data[start+12:end].decode()
    start = data.find(b'<Msg>')
    end = data.find(b'</Msg>')
    MsgID = int(data[start+5:end].decode())
    start = data.find(b'<Length>')
    end = data.find(b'</Length>')
    BinLength = int(data[start+8:end].decode())
    
    
    if XMLMsg.startswith('PU TM:'):    #This is a Profiler TM downloaded through docking connector
        #PU TM: 17.4, 6150, -2.8745, 170.9126, 18367.5
        message_type = 'PU TM'
        substr = XMLMsg.split(': ')     #split off the data from the message type
        substr = substr[1].split(', ')  #break the data up by commas, should return list of 5 strings
        profile_num = int(substr[0].split('.')[0]) #get the profile number (eg 17)
        profile_segment = int(substr[0].split('.')[1]) #get the profile segment (eg 4)
        pu_time = int(substr[1])            #get the time stamp
        initial_lat = float(substr[2])      #get the starting latitude for the profile
        initial_lon = float(substr[3])      #get the starting longitdue for the profile
        initial_alt = float(substr[4])      #get the altitude
        reel_pos = ""                       #no reel position - profiler is docked
        pu_current = ""                     #no profiler charge current in this message
        pu_battery_v = ""                   #no profiler battery voltage in this message
        pu_battery_t = ""                   #no profiler battery temperature in this message
        csvLine = [MsgID,date,time,message_type,BinLength,XMLMsg,profile_num,profile_segment,reel_pos,pu_current,pu_battery_v,pu_battery_t,pu_time,initial_lat, initial_lon,initial_alt,filename]
    
    elif XMLMsg.startswith('PU LoRa TM:'):  #This is a Profiler TM transmitted via LoRa
        #PU LoRa TM: 17.4, -4589.5, -2.8745, 170.9126, 18367.5
        message_type = 'PU LORA TM'
        substr = XMLMsg.split(': ')     #split off the data from the message type
        substr = substr[1].split(', ')  #break the data up by commas, should return list of 5 strings
        profile_num = int(substr[0].split('.')[0]) #get the profile number (eg 17)
        profile_segment = int(substr[0].split('.')[1]) #get the profile segment (eg 4)
        reel_pos = float(substr[1])         #get the reel position
        initial_lat = float(substr[2])      #get the starting latitude for the profile
        initial_lon = float(substr[3])      #get the starting longitdue for the profile
        initial_alt = float(substr[4])      #get the altitude
        pu_time = ""                        #no pu time stamp in this message
        pu_current = ""                     #no profiler charge current in this message
        pu_battery_v = ""                   #no profiler battery voltage in this message
        pu_battery_t = ""                   #no profiler battery temperature in this message
        csvLine = [MsgID,date,time,message_type,BinLength,XMLMsg,profile_num,profile_segment,reel_pos,pu_current,pu_battery_v,pu_battery_t,pu_time,initial_lat, initial_lon,initial_alt,filename]
    
    elif XMLMsg.startswith('PU ST:'):    #This is a PU LoRa status message
        #PU ST: 6150,18367,675.4,5.67,-9.87,-12.56,12.13,15.56,0.976,0,0,-4589.5,867.3"
        message_type = 'PU LORA STATUS'
        substr = XMLMsg.split(': ')         #split off the data from the message type
        substr = substr[1].split(',')      #break the data up by commas, should return list of 13 strings
        profile_num = ""                    #No profile number
        profile_segment = ""                #No profile segment
        reel_pos = float(substr[11])        #reel position (eg -4589.5)
        initial_lat = ""                    #no lat
        initial_lon = ""                    #no lon
        initial_alt = float(substr[1])      #get the altitude (eg 18367)
        pu_time = float(substr[0])          #pu time stamp (eg 6150)
        pu_current = float(substr[8])       #profiler charge current from PU (eg 0.976)
        pu_battery_v = float(substr[6])     #profiler battery voltage (eg 12.13)
        pu_battery_t = float(substr[3])     #profiler battery temperature (eg 5.67)
        csvLine = [MsgID,date,time,message_type,BinLength,XMLMsg,profile_num,profile_segment,reel_pos,pu_current,pu_battery_v,pu_battery_t,pu_time,initial_lat, initial_lon,initial_alt,filename]
    
    elif XMLMsg.startswith('PU TSEN:'): #this is a docked TSEN TM
        #PU TSEN: 6150, 12.13, 0.976, 5.67, -9.87, 0, -2.8745, 170.9126, 18367.5
        message_type = 'PU TSEN'
        substr = XMLMsg.split(': ')         #split off the data from the message type
        substr = substr[1].split(', ')      #break the data up by commas, should return list of 9 strings
        profile_num = ""                    #No profile number
        profile_segment = ""                #No profile segment
        reel_pos = ""                       #No reel position - profiler is docked
        initial_lat = float(substr[6])      #Latitude from Zephyr
        initial_lon = float(substr[7])      #Longitude from Zephyr
        initial_alt = float(substr[8])      #get the altitude (eg 18367)
        pu_time = float(substr[0])          #pu time stamp (eg 6150)
        pu_current = float(substr[2])       #profiler charge current from PU (eg 0.976)
        pu_battery_v = float(substr[1])     #profiler battery voltage (eg 12.13)
        pu_battery_t = float(substr[3])     #profiler battery temperature (eg 5.67)
        csvLine = [MsgID,date,time,message_type,BinLength,XMLMsg,profile_num,profile_segment,reel_pos,pu_current,pu_battery_v,pu_battery_t,pu_time,initial_lat, initial_lon,initial_alt,filename]
    
    else:
        message_type = 'OTHER'
        
        profile_num = ""                    #No profile number
        profile_segment = ""                #No profile segment
        reel_pos = ""                       #No reel position - profiler is docked
        initial_lat = ""      #Latitude from Zephyr
        initial_lon = ""      #Longitude from Zephyr
        initial_alt = ""      #get the altitude (eg 18367)
        pu_time = ""         #pu time stamp (eg 6150)
        pu_current = ""       #profiler charge current from PU (eg 0.976)
        pu_battery_v = ""     #profiler battery voltage (eg 12.13)
        pu_battery_t = ""     #profiler battery temperature (eg 5.67)
        csvLine = [MsgID,date,time,message_type,BinLength,XMLMsg,profile_num,profile_segment,reel_pos,pu_current,pu_battery_v,pu_battery_t,pu_time,initial_lat, initial_lon,initial_alt,filename]
    
    print(csvLine)
    
    with open(logFile, 'a') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(csvLine)
    
    
#    if os.path.basename(InputFile).startswith('ST2'):
#        #print(os.path.basename(InputFile))
#        #print(MsgID + ' ' + XMLMsg)
#        
#        with open(logFile, "a") as log:
#            log.write(os.path.basename(InputFile) + ': ' + MsgID + ' ' + XMLMsg + '\n')                        

if __name__ == '__main__':
    #loop_over_flights_and_instruments()
    files = glob.glob('/Users/kalnajs/Documents/Strateole/SIMULATEUR_SCI_02_08/LV/TM_FILES/*.ready_tm')
    files.sort()
    for file in files:
        readHeader(file,'RACHuTS_log.csv')
    