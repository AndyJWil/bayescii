# ------------------------------------------------------------
# For preparing WISE data for input to 
#  likelihood calculation in Naive Bayes Classifier
# Author: Andrew Wilson
# Institution: University of Exeter, Astrophysics Group
# Version date: 18 January 2023
# University email: aw648@exeter.ac.uk
# Personal email: andyjwilson_uk@hotmail.com
# ------------------------------------------------------------

import sys
import os

import numpy as np
#import math

from astropy.io import fits
from astropy.table import Table

import datetime

# Parameters common to the suite of Python scripts
from run_parameters import *


# Get the user inputs
def user_inputs():

	# Ask the user which data set to work with
	print('Choose source data:')
	for region_idx in range(0,len(region_descriptions)):
		print(' ' + region_descriptions[region_idx])
	source_data = input('Source data: ')
	# Convert source data choice into a title
	region_choice = region_list.get(source_data, "Invalid")
	if region_choice == "Invalid":
		print("Invalid choice")
		sys.exit()

	# Convert source data choice into a data file name
	data_file = wise_list.get(source_data, "Invalid")
	if data_file == "Invalid":
		print("Invalid choice")
		sys.exit()

	return data_file, region_choice

# Read a fits file of observations with a standard set of columns
def read_observations(file_name):
	print("Reading observations file ...")
	fits_file = fits.open(file_name)
	#fits_file.info()
	fits_data = fits_file[1].data
	fits_file.close
	# Convert the data into Table format
	table_data = Table(fits_data)
	# Add new columns
	table_data['wise_W1mW2'] = np.nan
	mask = np.logical_and( ~np.isnan(table_data['W1mag']) , ~np.isnan(table_data['W2mag']) )
	table_data['wise_W1mW2'][mask] = table_data['W1mag'][mask] - table_data['W2mag'][mask]
	return table_data

def output_full_fits(output_filename, data_table):
	print("Writing output fits file ...")
	# Delete file if it exists
	if os.path.isfile(output_filename):
		os.remove(output_filename)
	# Create fits header
	hdr = fits.Header()
	hdr['run_ver'] = str(run_version)
	hdr['datetime'] = str(datetime.datetime.now())
	hdr['file_obs'] = obs_file
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
	col_to_delete = []
	for col in obs_data.columns:
		if col not in ['source_id', 'l', 'b', 'parallax_zero_corrected', 'wise_W1mW2']:
			col_to_delete.append(col)
	obs_data.remove_columns(col_to_delete)



#------------------------------------------------------------
# Start of master routine
#------------------------------------------------------------

print('')
print('-----------------------------------------------------------------------------')
print('-------------------- Calculate WISE Excess Likelihoods ----------------------')
print('-----------------------------------------------------------------------------')
print('')

# Get the parameters from the user
obs_file, region_name = user_inputs()

# Main path for the run
main_path = run_filepath + region_name + "/"

# Read the observation data files
obs_data = read_observations(main_path + obs_file)

# Create output file with all columns
output_filename = main_path + "WISE_W1W2_" + region_name + "_full.fits"
output_full_fits(output_filename, obs_data)

# Remove interim processing columns not needed by Naive Bayes Classifier
remove_interim_columns()

# Create smaller output fits file for Naive Bayes Classifier
output_filename = main_path + "WISE_W1W2_" + region_name + ".fits"
output_full_fits(output_filename, obs_data)

#------------------------------------------------------------
# End of master routine
#------------------------------------------------------------