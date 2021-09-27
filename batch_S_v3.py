# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 01:31:31 2021

This script will be used to attempt to gather photometric data from the surveys contained in MAST.

@author: Joao Aguas
"""

import numpy as np
import pandas as pd
import sys
import os
import csv
from os import path
from sys import exit
from time import sleep
from astroquery.mast import Observations
from astropy import units as u
from astropy.coordinates import SkyCoord



#############################################################################
# Global Variables                                                          #
#############################################################################

filelist = []
outputfile = []
count = 0
survey = ['GALEX', 'HST']
dict_bands = {'GALEX':['FUV','NUV'],'HST_WFC3':['F1','F2','F3']}
descript_params = ['Background subtracted intensity map (J2000)', 
                    'Sky background image (J2000)', 
                    'Intensity map (J2000)']

project_params = ['MIS', 'AIS']

#############################################################################
# Control variables                                                         #
#############################################################################

chkavail = True

#############################################################################
# Constants                                                                 #
#############################################################################

radius= "0.02 deg"

#############################################################################
# Functions                                                                 #
#############################################################################

def check_exist_location(ra, dec, radius, survey, descript, project):
    

    # This function takes the columns of RA and DEC to check if there are any observations 
    # of that region, and how many are there in MAST. 
    # RA and DEC must be in decimal format, NOT in DEG:MM:SS.
    # print("\nChecking...\n")
    
    results=[]
    
    location = str(ra)+' '+str(dec)
        
    results = Observations.query_region(location,radius = radius)
    if(survey == 'GALEX'):
        final = results[(results['obs_collection']==survey) & 
                        (results['dataproduct_type']=='image')]
    
        data_products_by_obs = Observations.get_product_list(final)
     	
        ObjIDs = data_products_by_obs[((data_products_by_obs['description']==descript_params[0]) | 
                                   (data_products_by_obs['description']==descript_params[1]) |
                                   (data_products_by_obs['description']==descript_params[2])) & 
                                  ((data_products_by_obs['project']==project_params[0]))] 
    
    if(survey == 'HST'):
        print('\nWIP!!!\n')
    
    counts = len(ObjIDs)
    if counts == 0:
        counts=-9999
        ObjIDs=-9999
        
    return counts, ObjIDs

#############################################################################

def check_avail(ra, dec, radius):
    
    location = str(ra)+' '+str(dec)
    avail_surveys=np.unique(np.array(Observations.query_region(location)['obs_collection']))
    return avail_surveys

#############################################################################
# The Main                                                                  #
#############################################################################

# Getting information from the command line.
args = sys.argv

# FOR DEBUGING ONLY!!! COMMENT THIS OR ERASE WHEN IN PRODUCTION.
args = ['batch.py','t=list_snlist.txt', 'o=test_output.txt']

# Parsing and interpreting the command line.
for ii in range(len(args)):
    argtemp = args[ii].split('=')
    if argtemp[0] == 't':
        filelist = argtemp[1]
        count += 1
    if argtemp[0] == 'o':
        outputfile = argtemp[1]
        count += 1

# Verifying if the mandatory arguments.
if count != 2:
    print('\nERROR - Missing or excess of arguments.\n')
    print('Syntax: '+args[0].split('\\')[-1]+ ' t=targetsfile o=outputfile')
    exit()
        
# Verifying if the targets file exist and if it is not empty.
# If all is good to go, loads the target file into a Pandas Dataframe.
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

# If all is good to go, we load the modules that will help identify the location
# of the targets and download the desired images.

ngal = 0
nhst = 0
ncomp = 0

for ii in range(targets_df.shape[0]):
    
    # This line below verifies that the important information is there. RA and Dec.
    if (str(targets_df['decsn'][ii]).split() != ['nan']) & (str(targets_df['rasn'][ii]).split() != ['nan']):
        
        # Getting the information to the correct variables.
        RAstr = str(targets_df['rasn'][ii]).split()[0]
        DECstr = str(targets_df['decsn'][ii]).split()[0]
        
        pos = SkyCoord(RAstr,DECstr,frame='fk5',unit=(u.hourangle,u.deg))
        if(chkavail==True):
            sleep(0.5)
            avail=check_avail(RAstr, DECstr, radius)
            if 'HST' in avail: nhst+=1
            if 'GALEX' in avail: ngal+=1
            if ('GALEX' in avail) & ('HST' in avail): ncomp+=1
            
            print('\nObject: '+str(targets_df['snname'][ii].split()[0])+' was observed on the following surveys: \n'+str(avail))
        print('\nGALEX: '+str(ngal)+'\nHST: '+str(nhst)+'\nBoth: '+str(ncomp))    
        
        # count, ObjIDs = check_exist_location(pos.ra.value, pos.dec.value, radius, 'GALEX', 'descript', 'project')