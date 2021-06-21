# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 02:22:32 2021

Data Visualisation script. This script takes the output data from the photometry script, 
and presents the relations between the important variables.

@author: João Águas
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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
    
    df.replace(np.nan,value=-99, inplace=True)
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
    
    fail_value = []
    
    data_contents = data_df.Object.unique()
    z_contents = z_df.cubename.unique()
    bands_contents = data_df.Band.unique()
    
    z = []
    ez = []
    
    # Filling the z and e_z columns with the values from the z_df dataframe
    for ii in range(len(data_df)):
        if data_df['Object'][ii] in np.array(z_df['cubename']):
            z.append(z_df[z_df.cubename==data_df.Object[ii]].zhelio.unique()[0])
            ez.append(z_df[z_df.cubename==data_df.Object[ii]].ezhelio.unique()[0])
        else:
            z.append(-99)
            ez.append(-99)

    # Adding a new z and e_z columns to the data_df dataframe.
    # the z column contains the measured redshift and the e_z contains the 
    # associated error of the redshift measurement.
    data_df['z']=z
    data_df['ez']=ez
        
    return data_df
            
#############################################################################

def preview_relations(data_df,sx,sy,hue=None):
    
    data_df[data_df[data_df.columns[2:]]<0]=np.nan
    plt.figure(figsize=(sx,sy))
    sns.pairplot(data_df, hue=hue, dropna=True)
    
    return 0

#############################################################################

def scatter_reg_1(data, x_col, y_col, hue=None, sx=20, sy=20):
    plt.figure(figsize=(sx,sy))
    #plt.title(title)
    sns.lmplot(x=x_col, y=y_col, data=data, hue=hue)
    
    return 0

#############################################################################

def multi_reg(data, x_var, y_var, hue=None, kind='reg', aspect=1, height=1):
    sns.pairplot(data=data, hue=hue, x_vars=x_var,y_vars=y_var,kind=kind,height=height,aspect=aspect)
    
    return 0

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

# Merging the redshift table with the estimates table. 
data_df=merge_redshift(data_df, z_df)

# This line creates a preview plot where we may find the relations between all columns.
# Although some relations that are presented may not have a physical meaning.
# The objective of this plot is to quickly see if there is a relation between two variables.
preview_relations(data_df,40,40,'Band')

