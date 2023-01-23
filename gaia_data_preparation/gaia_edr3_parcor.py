# ------------------------------------------------------------
# For applying bias correction to Gaia EDR3 parallaxes using 
# Python script associated with Lindegren et al 2020
# Gaia Early Data Release 3: Parallax bias versus magnitude, colour and position
# https://gitlab.com/icc-ub/public/gaiadr3_zeropoint
# Author: Andrew Wilson
# Institution: University of Exeter, Astrophysics Group
# Version date: 4 December 2020
# University email: aw648@exeter.ac.uk
# Personal email: andyjwilson_uk@hotmail.com
# ------------------------------------------------------------

import sys
import os

import numpy as np
import math

from astropy.io import fits
from astropy.table import Table
#from astropy import units
#from astropy.coordinates import SkyCoord

# Lindegren EDR3 parallax bias calculation
from zero_point import zpt

import datetime

# Parameters common to the suite of Python scripts
from run_parameters import *

# Hard coded parameters
#----------------------
# Limit to the probablity of YSO over not YSO
prob_limit = 0.999

# Get the user inputs
def user_inputs():

	print('--------------------------------------------------------------------')

	# Ask the user the file path
	filepath = input('Data file path: ')
	# Ask the user the filename
	filename = input('Fits filename without .fits extention: ')

	print('--------------------------------------------------------------------')

	return filepath, filename

# Read a fits file of observations with a standard set of columns
def read_observations(file_name):
	print("Reading observations file ...")
	fits_file = fits.open(file_name)
	#fits_file.info()
	fits_data = fits_file[1].data
	fits_file.close
	# Convert the data into Table format
	table_data = Table(fits_data)
	# Calculate parallax bias and corrected parallax
	table_data['parallax_bias'] = np.nan
	table_data['parallax_zero_corrected'] = np.nan
	table_data['parallax_bias'] = zpt.get_zpt(table_data['phot_g_mean_mag'], table_data['nu_eff_used_in_astrometry'], table_data['pseudocolour'], table_data['ecl_lat'], table_data['astrometric_params_solved'])
	table_data['parallax_zero_corrected'] = table_data['parallax'] - table_data['parallax_bias']
	return table_data

def output_full_fits(output_filename, data_table, source_file):
	print("Writing output fits file ...")
	# Delete file if it exists
	if os.path.isfile(output_filename):
		os.remove(output_filename)
	# Create fits header
	hdr = fits.Header()
	hdr['run_ver'] = str(run_version)
	hdr['datetime'] = str(datetime.datetime.now())
	hdr['file_obs'] = source_file
	primary_hdu = fits.PrimaryHDU(header=hdr)
	# Convert table data into HDU
	hdu_data = fits.table_to_hdu(data_table)
	# Combine header with data
	hdul = fits.HDUList([primary_hdu, hdu_data])
	# Create the output file
	hdul.writeto(output_filename)

def remove_interim_columns():
	global obs_data
	print("Removing interim processing columns ...")
	obs_data.remove_columns(['ra', 'dec', 'l', 'b', 'parallax', 'parallax_zero_corrected', 'parallax_error', 'parallax_external_error', 'phot_g_n_obs', 'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_rp_mean_mag', 'phot_rp_mean_flux_over_error', 'bp_rp', 'phot_bp_rp_excess_factor', 'astrometric_excess_noise', 'astrometric_chi2_al', 'astrometric_n_good_obs_al', 'classiiyso_scale_length', 'classiiyso_frac_non_var', 'other_scale_length', 'other_frac_non_var', 'gaia_g_var_classiiyso_l_raw', 'gaia_g_var_other_l_raw', 'gaia_g_var_classiiyso_normprob', 'gaia_g_var_other_normprob'])





#------------------------------------------------------------
# Start of master routine
#------------------------------------------------------------

print('')
print('-----------------------------------------------------------------------------')
print('------------------ Calculate Gaia G magnitude Variability -------------------')
print('-----------------------------------------------------------------------------')
print('')

# Get the parameters from the user
data_path, data_filename = user_inputs()

# Main path for the run
main_path = "/Users/andy/OneDrive/Data/catalogues/Gaia/EDR3/"

# Load Gaia EDR3 coefficient tables
zpt.load_tables()

# Read the observation data file
obs_file = data_filename + ".fits"
obs_data = read_observations(data_path + "/" + obs_file)

# Create output fits file with all data
output_filename = data_path + "/" + data_filename + "_ParCor.fits"
output_full_fits(output_filename, obs_data, obs_file)

# Remove interim processing columns not needed by Naive Bayes Classifier
#remove_interim_columns()

#------------------------------------------------------------
# End of master routine
#------------------------------------------------------------
