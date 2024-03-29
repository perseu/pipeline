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

@author: Cepheu
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
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
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import Background2D, SExtractorBackground
from astropy.stats import SigmaClip


##############################################################################
#                   Notes for this script...                                 #
##############################################################################
"""
The cosmological constants for the calculations of the distances are taken from WMAP9.
GALEX image resolution = 1.5 arcsec/pixel
PS1 image resolution = 0.258 arcsec/pixel
"""
##############################################################################
#                   Constants                                                #
##############################################################################

zeropoint = {'fd':18.82, 'nd':20.08}

##############################################################################
#                   Config variables, if needed...                           #
##############################################################################

pix_size = 1.5 # arcsec/pixel
p_size = {'GALEX':1.5, 'PS1':0.25}
mask_sizes_kpc=[5,10,15,30] # mask radius in kpc
targets_file = []
redshift_file = []
output_file = []
filesOK = True
mask_size_pix = []
accumulator = []

# Debug activation flags
db_photo = 0 # When 1, it plots the region of interest.
silentmode = 0 # When 0, it displays on the console the calculated values.

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

def center_in_pix(header, ra, dec):
    """
    Converts RA and Dec to pixel coordinates. 

    Parameters
    ----------
    header : STRING
        FITS header information.
        
    ra : FLOAT
        Right Ascension
        
    dec : FLOAT
        Declination

    Returns
    -------
    x0 : FLOAT 
        x pixel coordinate
        
    y0 : FLOAT
        y pixel coordinate
        
    """
    
    w = wcs.WCS(header)
    t_loc = [[ra,dec]]
    pix_loc = w.wcs_world2pix(t_loc,0)
    x0 = int(pix_loc[0][0])
    y0 = int(pix_loc[0][1])
    
    return x0, y0

##############################################################################

def center_in_ra_dec(header, x0, y0):
    """
    Converts pixel coordinates to Right Ascension and Declination. 

    Parameters
    ----------
    header : STRING
        FITS header information.
        
    x0 : FLOAT
        x pixel coordinate
        
    y0 : FLOAT
        y pixel coordinate

    Returns
    -------
    ra : FLOAT 
        Right Ascension
        
    dec : FLOAT
        Declination
        
    """
    
    w = wcs.WCS(header)
    pix_loc=[[x0, y0]]
    locRADEC = w.wcs_pix2world(pix_loc,0)
    ra = locRADEC[0][0]
    dec = locRADEC[0][1]
    
    return ra, dec
    

##############################################################################

def select_images(target_df, bands, objname):
    """
    Takes the targets lists and returns a list of paths to the downloaded images. 

    Parameters
    ----------
    target_df : Dataframe
        The targets dataframe.
        
    bands : array
        The bands used on the observations
        
    objname : string
        Target name

    Returns
    -------
    img_path_list : array
        Paths for the files of the specific name objname
                
    """

    
    img_path_list = []
    
    for band in bands:
        img_path_list.append(targets_df[(targets_df['Object']==objname) & (targets_df['Band']==band)]['filepath'].values)
    
    return img_path_list

##############################################################################

def stacking(obj_cube):
    """
    This function stacks a list of images from the obj_cube and returns an
    average image.

    Parameters
    ----------
    obj_cube : List
        List of images.

    Returns
    -------
    final_list : List
        List of stacked images
                
    """
    
    final_list = []
    
    if(len(obj_cube)>1):    
        img_shape = obj_cube[0][0].data.shape
        ref_pix = [[(img_shape[0])/2,(img_shape[1])/2],[((img_shape[0])/2)+100,((img_shape[1])/2)+100]]
        ref_RA_Dec = []
        ref_vector = []
        ref_RA_Dec.append(center_in_ra_dec(obj_cube[0][0].header,(img_shape[0])/2,(img_shape[0])/2))
        ref_RA_Dec.append(center_in_ra_dec(obj_cube[0][0].header,((img_shape[0])/2)+100,((img_shape[0])/2)+100))
        ref_vector = (ref_pix[1][0]-ref_pix[0][0],ref_pix[1][1]-ref_pix[0][1])
        
        final_image = np.zeros(shape=img_shape)
        if(len(obj_cube)==1):
            return obj_cube[0][0].data
        else:
            # Alignment code goes here!!!
            for img in range(1,len(obj_cube)):
                pix_list = []
                t_vector = []
                pix_list.append(center_in_pix(obj_cube[img][0].header,ref_RA_Dec[0][0],ref_RA_Dec[0][1]))
                pix_list.append(center_in_pix(obj_cube[img][0].header,ref_RA_Dec[1][0],ref_RA_Dec[1][1]))
                t_vector = (pix_list[1][0]-pix_list[0][0],pix_list[1][1]-pix_list[0][1])
                
                # Now needs the translation and rotation part of the code.
                # I'M CONTINUING HEEERRRRREEEEEEE!!!
            
            # Computing the average.
            for img in obj_cube:
                final_image += img[0].data
                
            final_image = np.divide(final_image,len(obj_cube))
            final_list.append(final_image)
            
            return final_list                        
        
    if(len(obj_cube)==0):
        return -99

##############################################################################

def estimate_radius_pix(rr, zz, zzerr, arc_per_pix):
    """
    This function estimates the radius in pixels for the several radius in kpc
    stored in rr, using the redshift stored in zz, and taking into account the
    error in zzerr, and the correspondence between angular size and pixel size.
    
    Parameters
    ----------
    rr : Float
        Radius in kpc.
        
    zz : Float
        Redshift.
        
    zzerr : Float
        Redshift associated error.
        
    arc_per_pix : Float
        The size in arcsec for a pixel.
        
    Returns
    -------
    rr : Float
        The radius of the mask in pixels
    
    rrmin : Float
        The minimum radius taking into account the error.
    
    rrmax : Float
        The maximum radius taking into account the error.
                
    """

    
    zzvals = [zz-zzerr, zz, zz+zzerr]
    rrvals = []
    
    for ii in range(len(zzvals)):
        rrvals.append(radius_kpc_pix(rr,pix_size,zzvals[ii]))
    
    rrmax = rrvals[0]
    rr = rrvals[1]
    rrmin = rrvals[2]
    return rr , rrmin, rrmax

##############################################################################

def background_rms(data):
    """
    Computes the RMS from the image stored in the array data.

    Parameters
    ----------
    data : array
        Image array.

    Returns
    -------
    rms_median : 
        List of stacked images
                
    """
    
    sigmaClipper = SigmaClip(sigma=3, maxiters=5)
    bck_estimator = SExtractorBackground(data)
    bkg = Background2D(data,(50,50),filter_size=(3,3),sigma_clip=sigmaClipper,bkg_estimator=bck_estimator)
    print((bkg.background_median, bkg.background_rms_median))
    return bkg.background_rms_median


##############################################################################

def photo_estimate(data, background, x0, y0, radii, zp):
    """
    Computes the photometry for a mask of radius 'radii', centered at the pixel
    coordinates x0, y0, with the zeropoint zp, and taking into account the 
    background.

    Parameters
    ----------
    data : array
        Image array.
        
    background : Float
        background value.
        
    x0 : Float
        x0 pixel coordinate.
    
    y0 : Float
        y0 pixel coordinate.
        
    radii : Float
        The radius of the apperture.
        
    zp : Float
        The zeropoint.

    Returns
    -------
    mag : Float 
        The estimated magnitude.
        
    magerr : Float
        The associated error.
                
    """
    
    phot_data=[]
    mag = 0
    magerr = 0
    
    aperture = [CircularAperture((x0,y0), r) for r in radii]
    phot_data.append(aperture_photometry(data, aperture))
    
    mag = -2.5*np.log10(phot_data[0]['aperture_sum_0'][0])+zp
    mag1 = -2.5*np.log10(phot_data[0]['aperture_sum_1'][0])+zp
    mag2 = -2.5*np.log10(phot_data[0]['aperture_sum_2'][0])+zp
    
    magerr = np.abs(mag2-mag1)/2
    
    return mag, magerr

##############################################################################

def photo_measure(data, zp, x0, y0, r):
    """
    Computes the photometry for a mask of radius 'radii', centered at the pixel
    coordinates x0, y0, with the zeropoint zp, and taking into account the 
    background.

    Parameters
    ----------
    data : array
        Image array.
        
    zp : Float
        The zeropoint.
        
    x0 : Float
        x0 pixel coordinate.
    
    y0 : Float
        y0 pixel coordinate.
        
    r : Float
        The radius of the apperture.
        

    Returns
    -------
    magval : Float 
        The estimated magnitude.
        
    err : Float
        The associated error.
                
    """
    
    magval = []
    erms = []
    fval = []
    
    for rr in r:
        accum = 0
        npix = 0
        pixval = []
        roi = data[np.int(y0-rr):np.int(y0+rr),np.int(x0-rr):np.int(x0+rr)]
        for xx in range(roi.shape[0]):
            for yy in range(roi.shape[1]):
                if (np.sqrt((xx-rr/2)**2+(yy-rr/2)**2)<=rr) & (roi[yy][xx]>0):
                    pixval.append(roi[yy][xx])
                    npix +=1
                    
        avgval = np.average(pixval)
        npix=0
        pixval = []
        
        for xx in range(roi.shape[0]):
            for yy in range(roi.shape[1]):
                if (np.sqrt((xx-rr/2)**2+(yy-rr/2)**2)<=rr):
                    npix +=1
#                    if (roi[yy][xx]>0):
                    pixval.append(roi[yy][xx])
                    accum += roi[yy][xx]
#                    else:
#                        pixval.append(avgval)
#                        accum += avgval
                        
        erms.append(np.sqrt(np.multiply(npix,avgval)))
        
        fval.append(accum)
        
        magval.append(-2.5*np.log10(accum)+zp)
        
    ferr = np.abs(fval[2]-fval[1])/2
    zerr = np.abs(5/(2*fval[0]*np.log(10)))*np.abs(ferr)
    
    magerr = np.abs((-2.5*np.log10(np.average(np.array(erms)))))
    
    err = np.sqrt(magerr**2 + zerr**2)
    
    if db_photo == 1:
        plt.imshow(roi)
        plt.show()
    
    if (magval[0] == np.inf) | (magval[0] == np.NaN) | (err == np.inf) | (err == np.NaN):
        magval[0] = -99
        err = -99
        magerr = -99
        zerr = -99
    
    return magval[0], err, magerr, zerr, npix
    
    
##############################################################################
##########                 The main!!!                            ############
##############################################################################
"""
This script computes the magnitude and associated errors of a batch of objects
that are stored on a target file, taking a redshift file to compute the
distances and radius. 

Parameters
----------
t : CSV file
    Target list.
    
z : CSV file
    Redshift list.
    

Returns
-------
o : CSV Output file 
    The output with the estimated magnitudes/errors for the required appertures.
                
"""

# Getting script inline parameters
args = sys.argv

# Debug arguments... Comment the next line when the program is production.
args = ['batch.py','t=outputfile_MIS_ONLY.csv','z=list_redshift.txt', 'o=results_decomposed_error_upto30_pix.txt']

# Parsing and interpreting the command line.
for ii in range(len(args)):
    argtemp = args[ii].split('=')
    if argtemp[0] == 't':
        targets_file = argtemp[1]
    if argtemp[0] == 'z':
        redshift_file = argtemp[1]
    if argtemp[0] == 'o':
        res_file = argtemp[1]
        
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

# Discarding Objects that do not have the required information.

failed_obj = []
valid_obj = []

for obj in object_list:
    temp_z = []
    temp_t = []
    t_fail = False
    temp_z = redshift_df[redshift_df['cubename']==obj]
    temp_t = targets_df[targets_df['Object']==obj]
    if (len(temp_z) == 0) | (len(temp_t) == 0):
        t_fail = True
    else:
        if (np.array(temp_z['zhelio'])[0]==np.nan) | (np.array(temp_z['ezhelio'])[0]==np.nan):
            t_fail = True
        if (np.array(temp_t['ra'])[0]==np.nan) | (np.array(temp_t['dec'])[0]==np.nan):
            t_fail = True
            
    if t_fail == True:
        failed_obj.append(obj)
    else:
        valid_obj.append(obj)
        
object_list = valid_obj


# Starting the measurements!
# NOTE: Need to automatically generate this header depending on the number of appertures selected. Fixed header for 3 different appertures.
resaccum = [['Object','Band','mag_5kpc', 'e_5kpc', 'magerr_5kpc', 'zerr_5kpc', 'npix_5kpc','mag_10kpc', 'e_10kpc', 'magerr_10kpc', 'zerr_10kpc', 'npix_10kpc', 'mag_15kpc', 'e_15kpc', 'magerr_15kpc', 'zerr_15kpc', 'npix_15kpc', 'mag_30kpc', 'e_30kpc', 'magerr_30kpc', 'zerr_30kpc', 'npix_30kpc']]
radius_register = []


for obj in object_list:
    band_flag = 1
    mask_size_pix = []
    img_path_list=[]
    img_path_list = select_images(targets_df, band_list, obj)
    zhelio=np.array(redshift_df[redshift_df['cubename']==obj]['zhelio'])[0]    
    ezhelio=np.array(redshift_df[redshift_df['cubename']==obj]['ezhelio'])[0]
    
    # Calculating the radius in pixels that correspond to a radius of 5kpc, 10kpc and 15kpc.
    for jj in range(len(mask_sizes_kpc)):
        mask_pix=mask_sizes_kpc[jj]        
        mask_size_pix.append(estimate_radius_pix(mask_pix,zhelio,ezhelio,pix_size))
    
    for band in range(len(band_list)):
        obj_cube = []
        header_cube = []
        laccum = []    # Line accumulator. It's used to create the output line. Starts with the object name, band, and afterwards magnitude and associated error.
        laccum.append(obj)
        # laccum.append(band_list[band]) # NOTE: Need to change this band ID to the one on the photo file.
        for frame in img_path_list[band]:  # This cycle loads the filepaths to arrays depending on the object and band that is being reduced.
            tmpobj, tmphead = load_fits(frame)
            header_cube.append(tmphead)
            obj_cube.append(tmpobj)
        
        # This block of code gives the correct band name to the line. It takes it from the filename. 
        # There may be targets that only have information from a single band instead of both.
        if frame.split('-')[-2]=='fd': 
            laccum.append('FUV')
        if frame.split('-')[-2]=='nd':
            laccum.append('NUV')
        if (frame.split('-')[-2]!='nd') & (frame.split('-')[-2]!='fd'):
            band_flag = 0
        
        final_img = stacking(obj_cube)  # stacking images of the same object, and with the same band.
        if (len(np.shape(final_img)) > 0) & (band_flag == 1):
            final_img = final_img[0]
            x0, y0 = center_in_pix(header_cube[0],np.array(targets_df[targets_df['Object']==obj]['ra'])[0],np.array(targets_df[targets_df['Object']==obj]['dec'])[0])        
        
            for kk in range(len(mask_size_pix)): # In this cycle it computes the magnitudes for the desired mask sizes, and stores it first in a line accumulator laccum. At the end of the cycle the line is stored in the result accumulator resaccum.
                mag, merr, magerr, zerr, npix = photo_measure(final_img, zeropoint[frame.split('-')[-2]], x0, y0, mask_size_pix[kk])
                laccum.append(mag)
                laccum.append(merr)
                laccum.append(magerr)
                laccum.append(zerr)
                laccum.append(npix)
                if silentmode == 0:
                    print('Object='+str(obj)+',\tBand='+str(band_list[band])+',\tRedshift='+str(zhelio)+',\tRedshift Error='+str(ezhelio)+',\tRadius(in kpc)='+str(mask_sizes_kpc[kk])+',\tRadius(in pixels)='+str(mask_size_pix[kk][0])+',\tMag='+str(mag)+',\tError='+str(merr)+',\tMag_Error='+str(magerr)+',\tz_Error='+str(zerr)+',\tn_pixels='+str(npix))
            resaccum.append(laccum)
        else:
            failed_obj.append(obj) # If the stacking function does not return an image, it does not attempt to compute the magnitudes and errors, it simply adds the object to a failure list.
        
file = open(res_file, 'w+', newline ='')
with file:
    write = csv.writer(file)
    write.writerows(resaccum)
file.close()
