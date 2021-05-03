# pipeline
Data aquisition and reduction pipeline for GALEX data.

The batch_S_v2.py script is used to connect to the GALEX repository and to download the FITS files containing the data of interest.
The SNPhot2.py takes the list of data files and their location, if there is more than on photo of a certain field it preforms an average of all the photos and then performs photometry of the region for appertures that correspond to 5kpc, 10kpc and 15kpc, the correspondence for these appertures is calculated taking into account the Cosmological constant values from WMAP and the pixel size of the instrument, in this case 1.5 arcsec/pixel. Afterward it stores each set of measurements on a CSV type file. Where the first columns are the object name and the band used on the observation. The other columns are the measured magnitude and associated error for the different appertures (5kpc, 10kpc and 15kpc).

Refere to the sample files to create the input CSV files for the scripts.
