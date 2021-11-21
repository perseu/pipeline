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
from astroquery.sdss import SDSS
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy import units as u
from astropy.coordinates import SkyCoord



#############################################################################
# Global Variables                                                          #
#############################################################################

filelist = []
outputfile = []
count = 0
ObsIDNum = 0 # This value is differentiates Observations of PS1.
survey = ['GALEX', 'HST']
dict_bands = {'GALEX':['FUV','NUV'],'PS1':['i','g','r','y','z'],'SDSS':['u','g','r','i','z']}
descript_params = ['Background subtracted intensity map (J2000)', 
                    'Sky background image (J2000)', 
                    'Intensity map (J2000)',
                    'MIS',
                    'stack data image']

project_params = ['MIS', 'AIS']

SDSS_base_local = '/SDSS'

#############################################################################
# Control variables                                                         #
#############################################################################

chkavail = True                 # - Check if exists data from any listed survey
getall = True                   # - Gets data from all surveys, independently of
                                # coverage.
counting = True

#############################################################################
# Constants                                                                 #
#############################################################################

radius= "0.02 deg"

#############################################################################
# Functions                                                                 #
#############################################################################



def check_avail(pos, radius):
    
    avail_surveys = []
    sdss_avail = []
    out = []
    location = str(pos.ra.deg)+' '+str(pos.dec.deg)
    avail_surveys=np.unique(np.array(Observations.query_region(pos, radius=radius)['obs_collection']))
    sdss_avail= SDSS.query_region(pos, radius=radius)
    avail_surveys=list(avail_surveys)

    if ('NoneType' in str(type(sdss_avail))):
        pass
    else:
        avail_surveys.append('SDSS')

    if set(['GALEX']).issubset(set(avail_surveys)):
        out.append('GALEX')
    if set(['SDSS']).issubset(set(avail_surveys)):
        out.append('SDSS')
    if len(out) == 0:
        out='None'
    return out, location

#############################################################################

def filter_for_download(ra, dec, radius, survey):
    

    # This function takes the columns of RA and DEC to check if there are any observations 
    # of that region, and how many are there in MAST. 
    # RA and DEC must be in decimal format, NOT in DEG:MM:SS.
    # print("\nChecking...\n")
    
    results=[]
    
    location = str(ra)+' '+str(dec)
        
    results = Observations.query_region(location,radius = radius * u.arcsec)
    
    # This block of code filters the GALEX data of interest.
    if(survey == 'GALEX'):
        final = results[(results['obs_collection']==survey) & 
                        (results['dataproduct_type']=='image')]
    
        data_products_by_obs = Observations.get_product_list(final)
     	
        ObjIDs = data_products_by_obs[((data_products_by_obs['description']==descript_params[0]) | 
                                   (data_products_by_obs['description']==descript_params[1]) |
                                   (data_products_by_obs['description']==descript_params[2])) & 
                                  ((data_products_by_obs['project']==project_params[0]))] 
    
    # This block of code filters the PS1 data of interest.
    if(survey == 'SDSS'):
        final = results[(results['obs_collection']==survey)]
        data_products_by_obs = Observations.get_product_list(final)
        ObjIDs = data_products_by_obs[(data_products_by_obs['description']==descript_params[4])]
    
    # This small block verifies if any information was returned. If no information was returned,
    # then it returns -9999 as an error code.
    counts = len(ObjIDs)
    if counts == 0:
        counts=-9999
        ObjIDs=-9999
        
    # Returns the number of images found and the Object IDs for download.
    return counts, ObjIDs

#############################################################################


def strip_filename(survey,path):
    
    if(survey == 'PS1'):
        splitted = (path.split('\\')[-1]).split('.')
        band = splitted[6]
        cellx = splitted[3]
        celly = splitted[4]
        stat = splitted[7]
        return survey, band, cellx, celly, stat
    
    if(survey == 'GALEX'):
        temp = path.split('\\')
        expinfo=(temp[-1].split('-'))
        band = expinfo[-2]
        runid = (expinfo[-3]).split('_')[-1]
        einfo = (expinfo[-1]).split('.')[0]
        return einfo, runid, band
        
#############################################################################

def SDSS_download(pos, radius, datatype='image'):
    pass


#############################################################################

def filter_and_download(pos, radius, survey, mission='MIS', datatype='image'):
    
    results = []
    data_products = []
    
    if survey=='GALEX':
        results = Observations.query_region(pos, radius=radius)
        final = results[(results['obs_collection']==survey) & 
                        (results['dataproduct_type']==datatype)]
        
        data_products_by_obs = Observations.get_product_list(final)
        
        ObjIDs = data_products_by_obs[((data_products_by_obs['description']==descript_params[0]) | 
                                   (data_products_by_obs['description']==descript_params[1]) |
                                   (data_products_by_obs['description']==descript_params[2])) & 
                                  ((data_products_by_obs['project']==project_params[0]))] 
        
        data_products = Observations.download_products(ObjIDs,description=descript_params[0])
        return data_products
        
    if survey=='SDSS':
        results = SDSS.query_region(pos, radius=radius)
        # data_products = SDSS.get_images(pos) 
        pass


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
nsdss = 0
nnone = 0
ncomp = 0
location=[]

for ii in range(targets_df.shape[0]):
    
    # This line below verifies that the important information is there. RA and Dec.
    if (str(targets_df['decsn'][ii]).split() != ['nan']) & (str(targets_df['rasn'][ii]).split() != ['nan']):
        
        # Getting the information to the correct variables.
        RAstr = str(targets_df['rasn'][ii]).split()[0]
        DECstr = str(targets_df['decsn'][ii]).split()[0]
        
        pos = SkyCoord(RAstr,DECstr,frame='fk5',unit=(u.hourangle,u.deg))
        #avail,location = check_avail(pos.ra.deg, pos.dec.deg, radius)
        avail,location = check_avail(pos, radius)
        print('\nObject: '+str(targets_df['snname'][ii].split()[0])+' was observed on the following surveys: \n'+str(avail))
        
        if counting == True:
            if set(['GALEX']).issubset(set(avail)): ngal += 1
            if set(['SDSS']).issubset(set(avail)): nsdss += 1
            if set(['None']).issubset(set(avail)): nnone += 1
            if ('GALEX' in avail) & ('SDSS' in avail): ncomp += 1
            
        if 'GALEX' in avail:
            data_products = []
            downloaded_products = filter_and_download(pos, radius, 'GALEX')
            
        if 'SDSS' in avail:
            data_products = []
            downloaded_products = filter_and_download(pos, radius, 'SDSS')
            

# =============================================================================
#         avail,location = check_avail(pos.ra.deg, pos.dec.deg, radius)
#         
#         print('\nObject: '+str(targets_df['snname'][ii].split()[0])+' was observed on the following surveys: \n'+str(avail))
#         print('\nGALEX: '+str(ngal)+'\nPS1: '+str(nps1)+'\nHST: '+str(nhst)+'\nBoth: '+str(ncomp))    
#         
#         if 'PS1' in avail: 
#             nps1+=1
#             count, ObjIDs = filter_for_download(pos.ra.deg, pos.dec.deg, radius, 'PS1')
#             data_products = Observations.download_products(ObjIDs)           
#         if 'HST' in avail: 
#             nhst+=1
#         if 'GALEX' in avail: 
#             ngal+=1
#             count, ObjIDs = filter_for_download(pos.ra.deg, pos.dec.deg, radius, 'GALEX')
#             data_products = Observations.download_products(ObjIDs)
#         if ('GALEX' in avail) & (('PS1' in avail)|('HST' in avail)): 
#             ncomp+=1
# =============================================================================
        
        
        # count, ObjIDs = check_exist_location(pos.ra.value, pos.dec.value, radius, 'GALEX', 'descript', 'project')