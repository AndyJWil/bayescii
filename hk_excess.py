# ------------------------------------------------------------
# For calculating Near-IR using (J-H,H-K) from 2MASS and UKIDSS
# Author: Andrew Wilson
# Institution: University of Exeter, Astrophysics Group
# Original Version date: 18 January 2023
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

import scipy.optimize as optimization

import datetime

# Parameters common to the suite of Python scripts
from run_parameters import *

# Hard coded parameters
#----------------------


# Get the user inputs
def user_inputs():

	print('--------------------------------------------------------------------')

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

	# Convert source data choice into a input data file path
	data_file = hk_file_list.get(source_data, "Invalid")
	if data_file == "Invalid":
		print("Invalid choice")
		sys.exit()
	obs_data_file = data_file

	# Convert source data choice into catalogues
	cat_choice = hk_cat_list.get(source_data, "Invalid")

	print('--------------------------------------------------------------------')

	return obs_data_file, region_choice, cat_choice

# Read a fits file of observations with a standard set of columns
def read_observations(file_name, file_type):
	print("Reading " + file_type + " file ...")
	fits_file = fits.open(file_name)
	#fits_file.info()
	fits_data = fits_file[1].data
	fits_file.close
	# Convert the data into Table format
	table_data = Table(fits_data)

	return table_data

def add_reddening():
	global main_data, red_data

	# First sort the data so combine rows for the same sources
	print("Sorting data")
	main_data.sort('source_id')
	red_data.sort('source_id')

	# Check arrays are the same size
	main_size = np.size(main_data,0)
	red_size = np.size(red_data,0)
	if main_size != red_size:
		print("Aborting as number of rows in the main data (" + str(main_size) + ") is different to the number of rows in the reddening data (" + str(red_size) + ").")
		print("Files need to match for combination to be successful.")
		sys.exit()

	# Check the source_id in each row matches the source_id in the same row of other table
	check_mask = main_data['source_id'] == red_data['source_id']
	check_sum = np.sum(check_mask)
	if check_sum != main_size:
		print("Aborting as only " + str(check_sum) + " source_ids out of " + str(main_size) + " match between the main file and the reddening file.")
		print("Files need to match for combination to be successful.")
		sys.exit()
	
	print("Merging reddening into the main data.")
	main_data['E_JmH_2MASS'] = red_data['E_JmH_2MASS']
	main_data['E_HmK_2MASS'] = red_data['E_HmK_2MASS']
	main_data['E_JmH_UKIDSS'] = red_data['E_JmH_UKIDSS']
	main_data['E_HmK_UKIDSS'] = red_data['E_HmK_UKIDSS']
	main_data['E_JmH_VVV'] = red_data['E_JmH_VVV']
	main_data['E_HmK_VVV'] = red_data['E_HmK_VVV']
	main_data['suspect_reddening'] = red_data['suspect_reddening']

	# Add a blank space
	print()

# Calculate (H-K) excess from 2MASS and UKIDSS with a standard set of columns
def hk_2mass_ukidss():
	global main_data

	print("Processing combined 2MASS and UKIDSS data ...")

	# Add results columns
	print("Calculating H-K excess ...")
	# Select which source to use for which rows
	main_data['Source'] = "No data"
	T2MASS_mask = ~np.isnan(main_data['Kmag'])
	main_data['Source'][T2MASS_mask] = "2MASS"
	# UKIDSS selection is more complex
	# Selection based on Lucas et al 2008
	# Use 2MASS photometry where J>=13.25 and Hmag>=12.75 and Kmag>=12.
	# If no 2MASS data then ensure UKIDSS photometry has same limits, though won't be as good
	# TOPCAT filter
	#(!null_KAPERMAG3 && Jmag>=13.25 && Hmag>=12.75 && Kmag>=12) || (null_Kmag && JAPERMAG3>=13.25 && HAPERMAG3>=12.75 && KAPERMAG3>=12)
	# 2MASS photometry to select unsaturaged UKIDSS data
	print("----------------------------------------------------------------------")
	print("The following error message on greater than or equal to may be ignored")
	print("----------------------------------------------------------------------")
	UKIDSS_DATA = ~np.isnan(main_data['HMKPNT'])
	Phot_2MASS_p1 = np.logical_and( main_data['Jmag']>=13.25 , main_data['Hmag']>=12.75 )
	Phot_2MASS_p2 = np.logical_and( main_data['Kmag']>=12. , UKIDSS_DATA )
	Phot_2MASS_mask = np.logical_and( Phot_2MASS_p1 , Phot_2MASS_p2 )
	print("----------------------------------------------------------------------")
	# Where no 2MASS use UKIDSS photometry to select unsaturated UKIDSS data
	No_2MASS = ~T2MASS_mask
	Phot_UKIDSS_p1 = np.logical_and( main_data['JAPERMAG3']>=13.25 , main_data['HAPERMAG3']>=12.75 )
	Phot_UKIDSS_p2 = np.logical_and( main_data['KAPERMAG3']>=12. , No_2MASS )
	Phot_UKIDSS_mask = np.logical_and( Phot_UKIDSS_p1 , Phot_UKIDSS_p2 )
	# Combined UKIDSS unsaturated selection
	UKIDSS_mask = np.logical_or( Phot_2MASS_mask , Phot_UKIDSS_mask )
	main_data['Source'][UKIDSS_mask] = "UKIDSS"
	# Deredden each source
	main_data['JmH_dered_2MASS'] = np.nan
	main_data['HmK_dered_2MASS'] = np.nan
	main_data['JmH_dered_UKIDSS'] = np.nan
	main_data['HmK_dered_UKIDSS'] = np.nan
	main_data['JmH_dered_2MASS'][T2MASS_mask] = main_data['Jmag'][T2MASS_mask] - main_data['Hmag'][T2MASS_mask] - main_data['E_JmH_2MASS'][T2MASS_mask]
	main_data['HmK_dered_2MASS'][T2MASS_mask] = main_data['Hmag'][T2MASS_mask] - main_data['Kmag'][T2MASS_mask] - main_data['E_HmK_2MASS'][T2MASS_mask]
	main_data['JmH_dered_UKIDSS'][UKIDSS_mask] = main_data['JMHPNT'][UKIDSS_mask] - main_data['E_JmH_UKIDSS'][UKIDSS_mask]
	main_data['HmK_dered_UKIDSS'][UKIDSS_mask] = main_data['HMKPNT'][UKIDSS_mask] - main_data['E_HmK_UKIDSS'][UKIDSS_mask]
	# H-K offset
	main_data['HK_offset_2MASS'] = np.nan
	main_data['HK_offset_UKIDSS'] = np.nan
	main_data['HK_offset_2MASS'][T2MASS_mask] = main_data['HmK_dered_2MASS'][T2MASS_mask] - ((T2MASS_iso_slope * main_data['JmH_dered_2MASS'][T2MASS_mask]) + T2MASS_iso_intercept)
	main_data['HK_offset_UKIDSS'][UKIDSS_mask] = main_data['HmK_dered_UKIDSS'][UKIDSS_mask] - ((UKIDSS_iso_slope * main_data['JmH_dered_UKIDSS'][UKIDSS_mask]) + UKIDSS_iso_intercept)
	main_data['HK_offset'] = np.nan
	main_data['HK_offset'][T2MASS_mask] = main_data['HK_offset_2MASS'][T2MASS_mask]
	main_data['HK_offset'][UKIDSS_mask] = main_data['HK_offset_UKIDSS'][UKIDSS_mask]

# Calculate (H-K) excess from 2MASS only with a standard set of columns
def hk_2mass():
	global main_data

	# Add results columns
	print("Calculating H-K excess for 2MASS only data...")
	# Select which source to use for which rows
	main_data['Source'] = "No data"
	T2MASS_mask = ~np.isnan(main_data['Kmag'])
	main_data['Source'][T2MASS_mask] = "2MASS"

	# Deredden each source
	main_data['JmH_dered_2MASS'] = np.nan
	main_data['HmK_dered_2MASS'] = np.nan
	main_data['JmH_dered_UKIDSS'] = np.nan
	main_data['HmK_dered_UKIDSS'] = np.nan
	main_data['JmH_dered_2MASS'][T2MASS_mask] = main_data['Jmag'][T2MASS_mask] - main_data['Hmag'][T2MASS_mask] - main_data['E_JmH_2MASS'][T2MASS_mask]
	main_data['HmK_dered_2MASS'][T2MASS_mask] = main_data['Hmag'][T2MASS_mask] - main_data['Kmag'][T2MASS_mask] - main_data['E_HmK_2MASS'][T2MASS_mask]
	# H-K offset
	main_data['HK_offset_2MASS'] = np.nan
	main_data['HK_offset_UKIDSS'] = np.nan
	main_data['HK_offset_2MASS'][T2MASS_mask] = main_data['HmK_dered_2MASS'][T2MASS_mask] - ((T2MASS_iso_slope * main_data['JmH_dered_2MASS'][T2MASS_mask]) + T2MASS_iso_intercept)
	main_data['HK_offset'] = np.nan
	main_data['HK_offset'][T2MASS_mask] = main_data['HK_offset_2MASS'][T2MASS_mask]


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
	hdr['2MASS_JH'] = "2MASS: E(B-V) to E(J-H) = " + str(T2MASS_E_BV_to_E_JH)
	hdr['2MASS_HK'] = "2MASS: E(B-V) to E(H-K) = " + str(T2MASS_E_BV_to_E_HK)
	hdr['UKIDS_JH'] = "UKIDSS: E(B-V) to E(J-H) = " + str(UKIDSS_E_BV_to_E_JH)
	hdr['UKIDS_HK'] = "UKIDSS: E(B-V) to E(H-K) = " + str(UKIDSS_E_BV_to_E_HK)
	hdr['2MASSfit'] = "2MASS: (H-K) = " + str(T2MASS_iso_slope) + " * (J-H) + " + str(T2MASS_iso_intercept)
	hdr['UKIDSfit'] = "UKDISS: (H-K) = " + str(UKIDSS_iso_slope) + " * (J-H) + " + str(UKIDSS_iso_intercept)
	primary_hdu = fits.PrimaryHDU(header=hdr)
	# Convert table data into HDU
	hdu_data = fits.table_to_hdu(data_table)
	# Combine header with data
	hdul = fits.HDUList([primary_hdu, hdu_data])
	# Create the output file
	hdul.writeto(output_filename)

def remove_interim_columns():
	global main_data
	print("Removing interim processing columns ...")
	col_to_delete = []
	for col in main_data.columns:
		if col not in ['source_id', 'l', 'b', 'parallax_zero_corrected', 'Source', 'HK_offset', 'suspect_reddening']:
			col_to_delete.append(col)
	main_data.remove_columns(col_to_delete)

#------------------------------------------------------------
# Start of master routine
#------------------------------------------------------------

print('')
print('-----------------------------------------------------------------------------')
print('-------------------------- Calculate H-K Excess -----------------------------')
print('-----------------------------------------------------------------------------')
print('')

# Get the parameters from the user
obs_file, region_name, catalogues_used = user_inputs()

# Main path for the run
main_path = run_filepath + region_name + "/"

# Read the observation data file
main_data = read_observations(main_path + obs_file, "observations")

# Read the reddening data
red_filename = "Reddening_" + region_name + ".fits"
red_data = read_observations(main_path + red_filename, "reddening")

# Add the reddening to the main data
add_reddening()

# Calculate (H-K) excess
if catalogues_used=="2MASS_UKIDSS":
	hk_2mass_ukidss()
elif catalogues_used=="2MASS":
	hk_2mass()

# Create output fits file with all data
output_filename = main_path + "HKExcess_full_" + region_name + ".fits"
output_full_fits(output_filename, main_data, obs_file)

# Remove interim processing columns not needed by Naive Bayes Classifier
remove_interim_columns()

# Create smaller output fits file for Naive Bayes Classifier
output_filename = main_path + "HKExcess_" + region_name + ".fits"
output_full_fits(output_filename, main_data, obs_file)

#------------------------------------------------------------
# End of master routine
#------------------------------------------------------------