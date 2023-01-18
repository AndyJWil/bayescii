# ------------------------------------------------------------
# For creating an output file of source classification likelihoods
# from an input file of measuremens
# for input to Naive Bayes Classifier
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

	# Ask the user which measurement set to use
	print('Choose measurement data:')
	print('1 - WISE W1-W2')
	print('2 - (H-K) Excess from 2MASS UKIDSS')
	print('3 - Isochronal Age')
	print('4 - (r - H-alpha) Excess')
	source_data = input('Source data: ')
	# Convert measurement choice into a short name
	source_list = {
			"1": "WISE_W1W2",
			"2": "HKExcess",
			"3": "IsoAge",
			"4": "HaExcess"
			}
	measurement_choice = source_list.get(source_data, "Invalid")
	if measurement_choice == "Invalid":
		print("Invalid choice")
		sys.exit()
	# Convert measurement choice into a field name
	source_list = {
			"1": "W1mW2",
			"2": "HKExcess",
			"3": "IsoAge",
			"4": "rmha_offset"
			}
	measurement_choice_field_name = source_list.get(source_data, "Invalid")

	# Ask which version of likelihoods to use
	lhood_v = input('Likelihood version: ')

	return region_choice, measurement_choice, measurement_choice_field_name, lhood_v

# Read a fits file of observations with a standard set of columns
def read_observations(file_name, meas_type):
	print("Reading observations file ...")
	fits_file = fits.open(file_name)
	#fits_file.info()
	fits_data = fits_file[1].data
	fits_file.close
	# Convert the data into Table format
	table_data = Table(fits_data)
	# Add new columns
	table_data[meas_type.lower() + '_classiiyso_l_raw'] = 1.
	table_data[meas_type.lower() + '_other_l_raw'] = 1.
	table_data[meas_type.lower() + '_classiiyso_like'] = 1.
	table_data[meas_type.lower() + '_other_like'] = 1.
	table_data[meas_type.lower() + '_classiiyso_calc'] = False
	table_data[meas_type.lower() + '_other_calc'] = False
	table_data[meas_type.lower() + '_calc'] = False
	table_data[meas_type.lower() + '_classiiyso_capped'] = False
	table_data[meas_type.lower() + '_other_capped'] = False
	return table_data

def output_full_fits(output_filename, data_table, source_file, lhood_file, lhood_ver):
	print()
	print("Writing output fits file ...")
	# Delete file if it exists
	if os.path.isfile(output_filename):
		os.remove(output_filename)
	# Create fits header
	hdr = fits.Header()
	hdr['run_ver'] = str(run_version)
	hdr['lh_ver'] = str(lhood_ver)
	hdr['datetime'] = str(datetime.datetime.now())
	hdr['data_f'] = source_file
	hdr['lhood_f'] = lhood_file
	primary_hdu = fits.PrimaryHDU(header=hdr)
	# Convert table data into HDU
	hdu_data = fits.table_to_hdu(data_table)
	# Combine header with data
	hdul = fits.HDUList([primary_hdu, hdu_data])
	# Create the output file
	hdul.writeto(output_filename)

def remove_interim_columns(file_type, meas_type):
	global obs_data
	print("Removing interim processing columns ...")
	# Remove from all
	obs_data.remove_columns(['l', 'b', 'parallax_zero_corrected', meas_type.lower() + '_classiiyso_l_raw', meas_type.lower() + '_other_l_raw', meas_type.lower() + '_classiiyso_calc', meas_type.lower() + '_other_calc'])
	# Extra columns to remove from specified features
	if file_type == "IsoAge":
		obs_data.remove_columns(['suspect_reddening', 'rmi_r_iso_mass', 'igaps_data', 'brighter_max_iso_mag', 'fainter_max_iso_mag', 'above_iso_col', 'below_iso_col', 'back_interp_colour_diff', 'back_interp_mag_diff', 'good_except_back_interp_colour', 'good_except_back_interp_mag', 'above_max_iso_red', 'above_max_iso_mass', 'below_min_iso_mass'])
	elif file_type == "HaExcess":
		obs_data.remove_columns(['suspect_reddening'])
	elif file_type == "HKExcess":
		obs_data.remove_columns(['suspect_reddening'])

# Read a csv file of likelihood data based on criteria in the observation file
def read_likelihoods(file_name):
	print("Reading likelihoods file ...")
	read_data = np.loadtxt(file_name, dtype={'names': ('Classification', 'Probability', 'Field1', 'Field1Criteria1Qualifier', 'Field1Criteria1Value', 'Field1Criteria2Qualifier', 'Field1Criteria2Value', 'Field2', 'Field2Criteria1Qualifier', 'Field2Criteria1Value', 'Field2Criteria2Qualifier', 'Field2Criteria2Value', 'Field3', 'Field3Criteria1Qualifier', 'Field3Criteria1Value', 'Field3Criteria2Qualifier', 'Field3Criteria2Value'), 'formats': ('U10', 'f8', 'U40', 'U2', 'U40', 'U2', 'U40', 'U40', 'U2', 'U40', 'U2', 'U40', 'U40', 'U2', 'U40', 'U2', 'U40')}, delimiter=',', skiprows=1, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16))
	#read_data = np.loadtxt(file_name, dtype={'names': ('Classification', 'Probability'), 'formats': ('S6', 'f8')}, delimiter=',', skiprows=2, usecols=(0, 1))
	return read_data

def bayes_likelihoods(meas_type):
	global lhood_data, obs_data
	print("Calculating likelihoods for Naive Bayes Classifier ...")

	# Iterate over the rows in the likelihood array
	for lhood_row in np.nditer(lhood_data):

		print()
		
		mask_1a = single_mask(lhood_row['Field1'], lhood_row['Field1Criteria1Qualifier'], lhood_row['Field1Criteria1Value'])

		if str(lhood_row['Field1Criteria2Qualifier']) == 'na':
			mask_1 = mask_1a
		else:
			mask_1b = single_mask(lhood_row['Field1'], lhood_row['Field1Criteria2Qualifier'], lhood_row['Field1Criteria2Value'])
			mask_1 = np.logical_and( mask_1a , mask_1b)

		if str(lhood_row['Field2']) == 'na':
			mask_2 = mask_1
		else:
			mask_2a = single_mask(lhood_row['Field2'], lhood_row['Field2Criteria1Qualifier'], lhood_row['Field2Criteria1Value'])
			if str(lhood_row['Field2Criteria2Qualifier']) == 'na':
				mask_2_temp = mask_2a
			else:
				mask_2b = single_mask(lhood_row['Field2'], lhood_row['Field2Criteria2Qualifier'], lhood_row['Field2Criteria2Value'])
				mask_2_temp = np.logical_and( mask_2a , mask_2b)
			mask_2 = np.logical_and( mask_1 , mask_2_temp)

		if str(lhood_row['Field3']) == 'na':
			mask_3 = mask_2
		else:
			mask_3a = single_mask(lhood_row['Field3'], lhood_row['Field3Criteria1Qualifier'], lhood_row['Field3Criteria1Value'])
			if str(lhood_row['Field3Criteria2Qualifier']) == 'na':
				mask_3_temp = mask_3a
			else:
				mask_3b = single_mask(lhood_row['Field3'], lhood_row['Field3Criteria2Qualifier'], lhood_row['Field3Criteria2Value'])
				mask_3_temp = np.logical_and( mask_3a , mask_3b)
			mask_3 = np.logical_and( mask_2 , mask_3_temp)

		#print(lhood_row)
		print("Count = " + str(np.sum(mask_3)))
		obs_field = meas_type.lower() + "_" + str(lhood_row['Classification']).lower()
		obs_data[obs_field + "_l_raw"][mask_3] = float(lhood_row['Probability'])
		obs_data[obs_field + "_like"][mask_3] = float(lhood_row['Probability'])
		obs_data[obs_field + "_calc"][mask_3] = True
		obs_data[meas_type.lower() + "_calc"][mask_3] = True

# Creates a mask for a single criteria
def single_mask(field_name, criteria_qualifier, criteria_value):
	global obs_data

	crit_field = str(field_name)
	crit_qual = str(criteria_qualifier)

	try:
		crit_val = float(criteria_value)
	except:
		#crit_val = str(criteria_value).ljust(40)
		crit_val = str(criteria_value)

	#print(crit_field + ": " + crit_qual + str(criteria_value))

	if crit_qual == '=':
		output_mask = obs_data[crit_field] == crit_val
	elif crit_qual == '!=':
		output_mask = obs_data[crit_field] != crit_val
	elif crit_qual == '>':
		output_mask = obs_data[crit_field] > crit_val
	elif crit_qual == '>=':
		output_mask = obs_data[crit_field] >= crit_val
	elif crit_qual == '<':
		output_mask = obs_data[crit_field] < crit_val
	elif crit_qual == '<=':
		output_mask = obs_data[crit_field] <= crit_val

	print(crit_field + ": " + crit_qual + str(criteria_value) + " = " + str(np.sum(output_mask)))

	return output_mask

# Check likelihoods calculated for both categories, reset both to one if not so as to avoid erroneous likelihood ratios
def calc_check(meas_type):
	global obs_data

	print()
	print("Checking where likelihoods are present they are for all categories ...")

	mask_calc = np.logical_and( obs_data[meas_type.lower() + '_classiiyso_calc']==True , obs_data[meas_type.lower() + '_other_calc']==True )
	print("Calculated = " + str(np.sum(mask_calc)))
	# Enforce likelihoods of unity where likelihoods are not present for both categories
	mask_not_calc = ~mask_calc
	obs_data[meas_type.lower() + '_classiiyso_like'][mask_not_calc] = 1.
	obs_data[meas_type.lower() + '_other_like'][mask_not_calc] = 1.

# Cap the ratio of the likelihoods to the upper prior ratio
def cap_lhoods(meas_type):
	global obs_data

	print()
	print("Capping the likelihoods to a maximum of the upper ratio of the priors ...")

	# Calculate the upper prior ratio
	upper_prior_ratio = prior_other / prior_classiiyso

	# Limit the likelihood ratio to the opposite of the ratio of the priors
	mask_ratio_limit = obs_data[meas_type.lower() + '_classiiyso_like'] / obs_data[meas_type.lower() + '_other_like'] > upper_prior_ratio
	print("Capped Class II YSO likelihoods = " + str(np.sum(mask_ratio_limit)))
	if np.sum(mask_ratio_limit) > 0:
		obs_data[meas_type.lower() + '_classiiyso_like'][mask_ratio_limit] = obs_data[meas_type.lower() + '_other_like'][mask_ratio_limit] * upper_prior_ratio
		obs_data[meas_type.lower() + '_classiiyso_capped'][mask_ratio_limit] = True

	# Limit the likelihood ratio to the ratio of the priors
	mask_ratio_limit = obs_data[meas_type.lower() + '_other_like'] / obs_data[meas_type.lower() + '_classiiyso_like'] > upper_prior_ratio
	print("Capped Other likelihoods = " + str(np.sum(mask_ratio_limit)))
	if np.sum(mask_ratio_limit) > 0:
		obs_data[meas_type.lower() + '_other_like'][mask_ratio_limit] = obs_data[meas_type.lower() + '_classiiyso_like'][mask_ratio_limit] * upper_prior_ratio
		obs_data[meas_type.lower() + '_other_capped'][mask_ratio_limit] = True


#------------------------------------------------------------
# Start of master routine
#------------------------------------------------------------

print('')
print('-----------------------------------------------------------------------------')
print('----------------------- Determine Source Likelihoods ------------------------')
print('-----------------------------------------------------------------------------')
print('')

# Get the parameters from the user
region_name, measurement_name, measurement_field_name, lhood_version = user_inputs()

# Main path for the run
main_path = run_filepath + region_name + "/"

# Read the observation data files
input_filename = measurement_name + "_" + region_name + ".fits"
obs_data = read_observations(main_path + input_filename, measurement_field_name)

# Read in the likelihoods file
lhood_filename = measurement_name + "_Likelihoods_v" + str(lhood_version) + ".csv"
lhood_data = read_likelihoods(main_path + lhood_filename)

# Update the likelihoods
bayes_likelihoods(measurement_field_name)

# Check likelihoods calculated for both categories, reset both to one if not so as to avoid erroneous likelihood ratios
calc_check(measurement_field_name)

# Cap the ratio of the likelihoods to the upper prior ratio
cap_lhoods(measurement_field_name)

# Create output file with all columns
output_filename = main_path + measurement_name + "_lhoods_full_" + region_name + ".fits"
output_full_fits(output_filename, obs_data, input_filename, lhood_filename, lhood_version)

# Remove interim processing columns not needed by Naive Bayes Classifier
remove_interim_columns(measurement_name, measurement_field_name)

# Create smaller output fits file for Naive Bayes Classifier
output_filename = main_path + measurement_name + "_lhoods_" + region_name + ".fits"
output_full_fits(output_filename, obs_data, input_filename, lhood_filename, lhood_version)

#------------------------------------------------------------
# End of master routine
#------------------------------------------------------------