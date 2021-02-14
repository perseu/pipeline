# -*- coding: utf-8 -*-
"""

This script performs photometric measurements on the targets provided on the
input targets list. The target list file contains information about the location
of the FITS file, the name of the object, band used in the observation. The second
input file provides the redshift of the target to calculate the the radius of the appertures for 
the photometric mask.
The script outputs a CSV file containing a line per target, each line contains 
columns with the apparent magnitudes measured per apperture/band.


Created on Sun Jan 24 20:30:57 2021

@author: Joao Aguas
"""


import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import csv
import sys
import os
from sys import exit
from os import path
#from astroquery.mast import Observations
#from astroquery.ned import Ned
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
# from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import WMAP9 as cosmo


##############################################################################
#                   Notes for this script...                                 #
##############################################################################
"""
The cosmological constants for the calculations of the distances are taken from WMAP.
GALEX image resolution = 1.5 arcsec/pixel

"""
##############################################################################
#                   Config variables, if needed...                           #
##############################################################################

pix_size = 1.5 # arcsec/pixel
mask_sizes=[5,10,15] # mask radius in kpc
targets_file = []
redshift_file = []
filesOK = True

##############################################################################
#                            The Functions                                   #
##############################################################################

def radius_pix_kpc(n_pix,arc_p_pix,z):
    """
    This function calculates the transverse distance in kpc from the distance in pixels.
    To calculate this, the correspondence between angular distance and pixel distance,
    redshift and the distance in pixels. It returns the correspondent distance in kpc.
    
    Parameters
    ----------
    n_pix : FLOAT
        Transverse distance in pixels.
    arc_p_pix : FLOAT
        Size of each pixel in arcseconds.
    z : FLOAT
        The redshift of the object.

    Returns
    -------
    x_kpc : FLOAT
        The transverse distance in kpc.

    """
    
    kpc_arcmin=cosmo.kpc_proper_per_arcmin(z).value # returns kpc/arcmin
    
    #Converting to kpc per arcsec.
    kpc_arcsec = kpc_arcmin/60
    
    x_kpc = arc_p_pix*kpc_arcsec*n_pix
    
    return x_kpc
    
##############################################################################

def radius_kpc_pix(d_kpc,arc_p_pix,z):
    """
    This function calculates the transverse distance in pixels from the distance in kpc.
    To calculate this, it takes the correspondence between angular distance and pixel distance,
    redshift and the distance in kpc. It returns the correspondent distance in pixels.
   

    Parameters
    ----------
    d_kpc : FLOAT
        Tranverse distance in kpc.
    arc_p_pix : FLOAT
        The correspondence of a pixel in arcseconds.
    z : FLOAT
        The objects redshift.

    Returns
    -------
    x_pix : FLOAT
        The transverse distance in pixels

    """
    
    kpc_arcmin=cosmo.kpc_proper_per_arcmin(z).value # returns kpc/arcmin
    
    #Converting to kpc per arcsec.
    kpc_arcsec = kpc_arcmin/60
    
    x_pix = d_kpc/(arc_p_pix*kpc_arcsec)
    
    return x_pix

##############################################################################

def load_targets(filename):
    """
    Loads the data file with the target names and respective fits files that come
    from the script batch.py.

    Parameters
    ----------
    filename : STRING
        The file containing the targets list and respective file paths.

    Returns
    -------
    df : TYPE
        This function returns a PANDAS dataframe with the data from the input file.

    """
    df = pd.read_csv(filename,header=0)
    df.columns
    
    return df

##############################################################################

def load_z_file(filename):
    """
    Loads the data file with the table of redshifts.

    Parameters
    ----------
    filename : STRING
        The file containing the targets list and respective file paths.

    Returns
    -------
    df : TYPE
        This function returns a PANDAS dataframe with the data from the input file.

    """
    df = pd.read_csv(filename,header=0)
    df.columns
    
    return df
    
##############################################################################

def load_fits(filename):
    """
    This just loads a FITS file, and return the data and the header.

    Parameters
    ----------
    filename : STRING
        path to the FITS file.

    Returns
    -------
    data : array
        The data contained in the file.
    
    header: array
        An array containing the files header information.
        
    """
    
    data=fits.open(filename)
    
    header=data[0].header
    
    return data, header

##############################################################################

def center_in_pix(data, ra, dec):
    
    
    w = wcs.WCS(data.header)
    t_loc = [[ra,dec]]
    pix_loc = w.wcs_world2pix(t_loc,0)
    x0 = int(pix_loc[0][0])
    y0 = int(pix_loc[0][1])
    
    return x0, y0
    
##############################################################################
##########                 The main!!!                            ############
##############################################################################

# Getting script inline parameters
args = sys.argv

# Debug arguments... Comment the next line when the program is production.
args = ['batch.py','t=outputfile.csv','z=zlistsample.txt']

# Parsing and interpreting the command line.
for ii in range(len(args)):
    argtemp = args[ii].split('=')
    if argtemp[0] == 't':
        targets_file = argtemp[1]
    if argtemp[0] == 'z':
        redshift_file = argtemp[1]
        
# This block of code checks the existance of the targets list and redshift list. 
# It also checks if the respective files are empty or not.
if len(args) > 1:
    if os.path.isfile(targets_file):
        if os.path.getsize(targets_file) > 0:
            targets_df = load_targets(targets_file)
        else:
            print('\n\nTarget list file '+targets_file+' it\'s empty.')
            filesOK = False            
    else:
        print('\n\nTarget list file '+targets_file+' not found.')
        filesOK=False
        
    if os.path.isfile(redshift_file):
        if os.path.getsize(redshift_file) > 0:
            redshift_df = load_z_file(redshift_file)
        else:
            print('\n\nTarget list file '+redshift_file+' it\'s empty.')
            filesOK = False            
    else:
        print('\n\nTarget list file '+redshift_file+' not found.')
        filesOK=False
        
if filesOK == False:
    exit()
    
# Creating lists of unique values of object names and bands.
object_list = targets_df.Object.unique()
band_list = targets_df.Band.unique()

for ii in range(targets_df):
    hdul, header = load_fits(targets_df['filepath'][ii])
    