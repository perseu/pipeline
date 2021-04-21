# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 18:35:20 2021

This script will be used to attempt to gather photometric data from GALEX images.

@author: Joao Aguas
"""


import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
#import csv
import sys
import os
import csv
from sys import exit
from os import path
from sys import exit
from astroquery.mast import Observations
from astropy import units as u
from astropy.coordinates import SkyCoord


##############################################################################
#                   Config variables, if needed...                           #
##############################################################################

targetlist = [['Object','Band','Observation','Info','filepath','ra','dec']]
targetfail = []
TO_counts = []
radius= "0.02 deg"
descript = 'Background subtracted intensity map (J2000)'
skybkg = 'Sky background image (J2000)'
intmap = 'Intensity map (J2000)'
targetfilename = 'outputfile.csv'
failedfilename = 'failedtarget.csv'

##############################################################################
#                            The Functions                                   #
##############################################################################

def conv_dms_deg(degstr):
    # This function takes the declination in format DEG:MM:SS as a string and 
    # returns a degrees as a float. It's a very complex thing to do. :-\
    # This code will not be elegant, it'll be linear and functional. :-|
    # As dec goes between [-90,90], if it goes out of this interval returns -9999
    
    nsign = False
    ang=np.float32(degstr.split(':'))
    
    val = ang[0]
    
    if val < 0:
        nsign = True
        
    val = np.abs(val)
    val = val + (ang[1]/60)
    val = val + (ang[2]/3600)
            
    if nsign == True:
        val = -val
        
    if np.abs(val) > 90:
        val = -9999

    return val
    
##############################################################################

def conv_hms_deg(degstr):
    
    ## SGG: This function is NOT doing the right thing! Be careful!


    # This function takes the declination in format HH:MM:SS as a string and 
    # returns a degrees as a float. 
    # As RA goes between [0,24], if it goes out of this interval returns -9999
    
    ang=np.float32(degstr.split(':'))
    
    val = ang[0]
    
    if val < 0:
        val=-9999
    else:        
        val = val + (ang[1]/60)
        val = val + (ang[2]/3600)
        val = 360/24*val            
    if np.abs(val) > 24:
        val = -9999

    return val

##############################################################################

def check_exist_location(ra, dec,radius):
    

    # This function takes the columns of RA and DEC to check if there are any observations 
    # of that region, and how many are there in MAST. 
    # RA and DEC must be in decimal format, NOT in DEG:MM:SS.
    # print("\nChecking...\n")
    
    results=[]
    
    location = str(ra)+' '+str(dec)
        
    results = Observations.query_region(location,radius = radius)
    final = results[(results['obs_collection']=='GALEX') & 
                    (results['dataproduct_type']=='image')]
    
    data_products_by_obs = Observations.get_product_list(final)
     	
    ObjIDs = data_products_by_obs[((data_products_by_obs['description']==descript) | 
                                   (data_products_by_obs['description']==skybkg) |
                                   (data_products_by_obs['description']==intmap)) & 
                                  ((data_products_by_obs['project']=='MIS') | 
                                   (data_products_by_obs['project']=='AIS'))]
    
    counts = len(ObjIDs)
    if counts == 0:
        counts=-9999
        ObjIDs=-9999
        
    return counts, ObjIDs

##############################################################################

def strip_filename(filelocation):
    
    # This function will parse the file name to return information about the photo run,
    # and band.
    
    temp = filelocation.split('\\')
    expinfo=(temp[-1].split('-'))
    band = expinfo[-2]
    runid = (expinfo[-3]).split('_')[-1]
    einfo = (expinfo[-1]).split('.')[0]
    return einfo,runid,band


##############################################################################
##########                 The main!!!                            ############
##############################################################################

# Getting script inline parameters
args = sys.argv

# Debug arguments... Comment the next line when the program is running.
args = ['batch.py','t=snlistsample.txt']

# Parsing and interpreting the command line.
for ii in range(len(args)):
    argtemp = args[ii].split('=')
    if argtemp[0] == 't':
        filelist = argtemp[1]
        
if len(args) > 1:
    if os.path.isfile(filelist):
        if os.path.getsize(filelist) > 0:
            targets_df = pd.read_csv(filelist)
        else:
            print('\n\nTarget list file '+filelist+' it\'s empty.')
            exit()
            
    else:
        print('\n\nTarget list file '+filelist+' not found.')
        exit()

# Politely inquiring the MAST server for the targets. One has to be always polite.
from astroquery.mast import Observations
from astropy.coordinates import SkyCoord

for ii in range(targets_df.shape[0]):
    
    #print('\nObject '+str(targets_df['snname'][ii])+' ')
    if (str(targets_df['decsn'][ii]).split() != ['nan']) & (str(targets_df['rasn'][ii]).split() != ['nan']):
        
        RAstr = str(targets_df['rasn'][ii]).split()[0]
        #RAstr = str(conv_hms_deg(RAstr))
        DECstr = str(targets_df['decsn'][ii]).split()[0]
        #DECstr = str(conv_dms_deg(DECstr))
        pos = SkyCoord(RAstr,DECstr,frame='fk5',unit=(u.hourangle,u.deg))
        #location = str(RAstr)+' '+str(DECstr)
        count, ObjIDs = check_exist_location(pos.ra.value, pos.dec.value, radius)
        if count > 0:
            TO_counts.append([str(targets_df['snname'][ii]).split(), count])
            print('\n'+str(targets_df['snname'][ii])+' OK!!! Number of entries for area: '+str(count))
            print('Downloading Files...')
            data_products = Observations.download_products(ObjIDs,description=descript)
                
            print('\n\nFor the object ' + (targets_df['snname'][ii]).split()[0])
            print(data_products)  
            for kk in range(len(data_products)):
                einfo, runid, band = strip_filename(data_products['Local Path'][kk])
                insinfo = ((targets_df['snname'][ii]).split()[0]),band,runid,einfo,str(data_products['Local Path'][kk]),str(pos.ra.value),str(pos.dec.value)
                targetlist.append(insinfo)
                insinfo = []
        else:
            print(str(targets_df['snname'][ii])+' No dice!!! I could not find anything. Adding target to failed list!')
            targetfail.append((targets_df['snname'][ii]).split()[0])
            
    else:
        print(str(targets_df['snname'][ii])+' Missing Location!!!')
        targetfail.append((targets_df['snname'][ii]).split()[0])

file = open(targetfilename, 'w+', newline ='')
with file:
    write = csv.writer(file)
    write.writerows(targetlist)
    
file.close()

if len(targetfail) > 0:
    file1 = open(failedfilename, 'w+', newline ='')
    with file1:
        write1 = csv.writer(file1)
        write1.writerow(targetfail)
        
    file.close()
    print('Sadly there were some objects that had missing data. Could not help there.\nThe failed targets file is: '+failedfilename)
else:
    print('CONGRATULATIONS!!! We had no problem locating your targets. No target failure file written. ')
    
print('All that is good must end. My mission ends here.\nThe output file is called: '+ targetfilename)