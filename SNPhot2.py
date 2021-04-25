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
The cosmological constants for the calculations of the distances are taken from WMAP.
GALEX image resolution = 1.5 arcsec/pixel

"""
##############################################################################
#                   Constants                                                #
##############################################################################

zeropoint = {'fd':18.82, 'nd':20.08}

##############################################################################
#                   Config variables, if needed...                           #
##############################################################################

pix_size = 1.5 # arcsec/pixel
mask_sizes_kpc=[5,10,15] # mask radius in kpc
targets_file = []
redshift_file = []
output_file = []
filesOK = True
mask_size_pix = []
accumulator = []

# Debug activation flags
db_photo = 0 # When 1, it plots the region of interest.

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
    
    
    w = wcs.WCS(header)
    t_loc = [[ra,dec]]
    pix_loc = w.wcs_world2pix(t_loc,0)
    x0 = int(pix_loc[0][0])
    y0 = int(pix_loc[0][1])
    
    return x0, y0

##############################################################################

def center_in_ra_dec(header, x0, y0):
    
    w = wcs.WCS(header)
    pix_loc=[[x0, y0]]
    locRADEC = w.wcs_pix2world(pix_loc,0)
    ra = locRADEC[0][0]
    dec = locRADEC[0][1]
    
    return ra, dec
    

##############################################################################

def select_images(target_df, bands, objname):
    
    img_path_list = []
    
    for band in bands:
        img_path_list.append(targets_df[(targets_df['Object']==objname) & (targets_df['Band']==band)]['filepath'].values)
    
    return img_path_list

##############################################################################

def stacking(obj_cube):
    
    # marking two points on the reference image, assuming the first image as the reference.
    # These points are then translated from pixel coordinates to the WCS. 
    # As all images have a WCS, it is easy to find out the location of those points 
    # on other pictures of the same field, and afterward compute the necessary
    # translation and rotation.
    
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
        return -9999

##############################################################################

def estimate_radius_pix(rr, zz, zzerr, arc_per_pix):
    
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
    # Determines the rms of the background.
    sigmaClipper = SigmaClip(sigma=3, maxiters=5)
    bck_estimator = SExtractorBackground(data)
    bkg = Background2D(data,(50,50),filter_size=(3,3),sigma_clip=sigmaClipper,bkg_estimator=bck_estimator)
    print((bkg.background_median, bkg.background_rms_median))
    return bkg.background_rms_median


##############################################################################

def photo_estimate(data, background, x0, y0, radii, zp):
    
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
    
    magval = []
    erms = []
    
    for rr in r:
        accum = 0
        npix = 0
        pixval = []
        roi = data[np.int(y0-rr):np.int(y0+rr),np.int(x0-rr):np.int(x0+rr)]
        for xx in range(roi.shape[0]):
            for yy in range(roi.shape[1]):
                if (np.sqrt((xx-rr/2)**2+(yy-rr/2)**2)<=rr) & (roi[yy][xx]>0):
                    pixval.append(roi[yy][xx])
                    accum += roi[yy][xx]
                    npix +=1
        
        erms.append(np.sqrt(np.multiply(npix**2,np.average(np.array(pixval))**2)))
        
        magval.append(-2.5*np.log10(accum)+zp)
        
    magerr = np.abs(magval[2]-magval[1])/2
    
    err = np.sqrt(magerr**2 + (-2.5*np.log10(np.average(np.array(erms))))**2)
    
    if db_photo == 1:
        plt.imshow(roi)
        plt.show()
    
    if (magval[0] == np.inf) | (magval[0] == np.NaN) | (err == np.inf) | (err == np.NaN):
        magval[0] = -9999
        err = -9999
    
    return magval[0], err
    
    
##############################################################################
##########                 The main!!!                            ############
##############################################################################

# Getting script inline parameters
args = sys.argv

# Debug arguments... Comment the next line when the program is production.
args = ['batch.py','t=outputfile.csv','z=zlistsample.txt', 'o=results.txt']

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
resaccum = []
radius_register = []

for obj in object_list:
    mask_size_pix = []
    img_path_list=[]
    img_path_list = select_images(targets_df, band_list, obj)
    zhelio=np.array(redshift_df[redshift_df['cubename']==obj]['zhelio'])[0]    # Erro no segundo objecto
    ezhelio=np.array(redshift_df[redshift_df['cubename']==obj]['ezhelio'])[0]

    for jj in range(len(mask_sizes_kpc)):
        mask_pix=mask_sizes_kpc[jj]        
        mask_size_pix.append(estimate_radius_pix(mask_pix,zhelio,ezhelio,pix_size))
    
    for band in range(len(band_list)):
        obj_cube = []
        header_cube = []
        laccum = []
        laccum.append(obj)
        laccum.append(band_list[band])
        for frame in img_path_list[band]:
            tmpobj, tmphead = load_fits(frame)
            header_cube.append(tmphead)
            obj_cube.append(tmpobj)
            
        final_img = stacking(obj_cube)
        final_img = final_img[0]
        x0, y0 = center_in_pix(header_cube[0],np.array(targets_df[targets_df['Object']==obj]['ra'])[0],np.array(targets_df[targets_df['Object']==obj]['dec'])[0])        
      #  mask = np.zeros(shape=final_img.shape)
        
        for kk in range(len(mask_size_pix)):
            mag, merr = photo_measure(final_img, zeropoint[band_list[band]], x0, y0, mask_size_pix[kk])
            laccum.append(mag)
            laccum.append(merr)
        
        resaccum.append(laccum)
        
file = open(res_file, 'w+', newline ='')
with file:
    write = csv.writer(file)
    write.writerows(resaccum)
file.close()

# for ii in range(len(targets_df)):
#     # This cycle loads the images, detects the center of the the interest area,
#     # creates the interest region for the photometry, and performs the photometric 
#     # estimates with the associated errors.
    
#     mask_size_pix = []
    
#     hdul, header = load_fits(targets_df['filepath'][ii])
#     x0, y0 = center_in_pix(header, targets_df['ra'][ii],targets_df['dec'][ii])
    
#     # Determining the background.
#     # bck = background_rms(hdul[0].data)
#     bck = 0
    
#     # This cycle estimates the radius in pixels for the interest region.
#     # It takes the cosmology of the WMAP9 to calculate the angular sizes for the 
#     # distances stored on the array mask_sizes_kpc, and converts it to the angular
#     # size in pixels.
    
#     for jj in range(len(mask_sizes_kpc)):
#         mask=mask_sizes_kpc[jj]
#         zhelio=redshift_df[redshift_df['cubename']==targets_df['Object'][ii]]['zhelio'].values[0]
#         ezhelio=redshift_df[redshift_df['cubename']==targets_df['Object'][ii]]['ezhelio'].values[0]
#         mask_size_pix.append(estimate_radius_pix(mask,zhelio,ezhelio,pix_size))

#     # This array is used to store the lines containing the calculations done with
#     # the several aperture sizes.
    
#     templine = []
#     templine.append(targets_df['Object'][ii])
#     templine.append(targets_df['Band'][ii])
    
#     for kk in range(len(mask_size_pix)):
#         print('Estimating target: '+ targets_df['Object'][ii])
#         mag, err = photo_estimate(hdul[0].data,bck,x0,y0,mask_size_pix[kk], zeropoint[targets_df['Band'][ii]])
#         templine.append(mag)
#         templine.append(err)
    
#     # The accumulator stores the lines that are formatted above,
#     # afterward this information will be properly formatted, and stored in an
#     # output file.
    
#     accumulator.append(templine)
    
#     # creating output data format.
#     # Starting with the header.
    
# outputdata = [['Object', 'Band']]
# for ll in range(len(mask_sizes_kpc)):
#     outputdata[0].append(str(mask_sizes_kpc[ll])+'kpc')
#     outputdata[0].append(str(mask_sizes_kpc[ll])+'kpc_error')
        
# for ll in range(len(accumulator)):
#     outputdata.append(accumulator[ll])