# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 02:22:32 2021

Data Visualisation script. This script takes the output data from the photometry script, and presents the relations between the important variables.

@author: Cobra
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os


#############################################################################
# Global Constants and Variables                                            #
#############################################################################

FilesOK = 0


#############################################################################
# Functions                                                                 #
#############################################################################

#############################################################################

def load_data(filename,header="infer"):
    
    dataframe = pd.read_csv(filename,header=header)
    return dataframe

#############################################################################

def clear_nan(df):
    
    df.replace(np.nan,value=-99)
    return df

#############################################################################

def list_unique(df,colname):
    
    obj_list = np.unique(np.array(df[colname]))
    return obj_list

#############################################################################

def available_bands(df):
    
    bands=list_unique(df, 'Band')
    return bands

#############################################################################

def merge_redshift(data_df, z_df):
    
    data_contents = data_df.Object.unique()
    z_contents = z_df.cubename.unique()
    
    # Adding a new z column to the data_df dataframe.
    data_df['z'] = np.zeros(len(data_df))
    
    

#############################################################################
# If this was C, then from this point on, this would be the MAIN            #
#############################################################################

# Getting script inline parameters
args = sys.argv

# Debug arguments... Comment the next line when the program is production.
args = ['batch.py','z=list_redshift.txt', 'r=results_finalMIS.txt']

# Parsing and interpreting the command line.
for ii in range(len(args)):
    argtemp = args[ii].split('=')
    if argtemp[0] == 'z':
        redshift_file = argtemp[1]
    if argtemp[0] == 'r':
        res_file = argtemp[1]
        
if os.path.isfile(redshift_file):
    if os.path.getsize(redshift_file) > 0:
        FilesOK += 1
if os.path.isfile(res_file):
    if os.path.getsize(redshift_file) > 0:
        FilesOK += 1

if FilesOK < 2:
    print('\nATTENTION: One or both argument files is missing or empty!!!\n')
    exit(0)
    
data_df = load_data(res_file)
z_df = load_data(redshift_file)

data_df = clear_nan(data_df)
z_df = clear_nan(z_df)
bands = available_bands(data_df)

