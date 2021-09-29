# RACHuTS_Processing
Python Code to download and process RACHuTS data

The current version of the code is untested alpha version.  It has the same dependencies and usage as the LPC download code.   When data becomes available on the CCMz, GetRACHuTS will download all the RACHuTS TM data to the specified local directory, extract the ascii headers from the TMs and generate a CSV log file with the status information, then if the file is from a profile, it will copy that file into a subdir for that profile.

If there is profile data that has been downloaded from the CCMz and copied into a subdir, it can then be processed using the routines in readRACHuTSProfileXML.py.   The TMAppend function will append all the TM files in a given directory into a single TM file.   This appended TM file can then processed into a csv file using the parseProfileDatatoCSV(directory) function.  Finally some rudimentry vizualizations are available by passing the profile CSV file to plotProfile(file).  For a docked profile (ie a profile at a fixed altitude) use the form: plotProfile(CSVFile,vsTime=True, HK=True, GPS=False).
