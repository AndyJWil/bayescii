# ------------------------------------------------------------
# Naive Bayes Classifier
# Takes input source fits files of likelihoods for different
#  measures with columns for different classifications.
# Author: Andrew Wilson
# Institution: University of Exeter, Astrophysics Group
# Version date: 27 April 2021
# University email: aw648@exeter.ac.uk
# Personal email: andyjwilson_uk@hotmail.com
# ------------------------------------------------------------
#
# Code assumptions
# ----------------
#
# The code is adaptive with only the following hard coded:
# - File names, set in the user_inputs procedure
# - Categories (e.g. YSO, not YSO), set in bayes_calculation procedure
# - Priors, set in the bayes_calculation procedure
#
# For this to work properly it makes assumptions about text in the field names
# that must be in each type of field and not in other types of field.
# - The likelihood columns have _like in their names
# - Categories are named consistently and contain the following in the column names
#   YSO contains _yso
#   Not YSO contains _notyso
#
# ------------------------------------------------------------

import sys
import os

import numpy as np

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

	# Convert source data choice into a data file names
	# Coordinates
	source_files=[]
	source_files.append("Coords_" + region_choice + ".fits")
	# Reddening
	source_files.append("Reddening_" + region_choice + ".fits")
	# Gaia G variability
	source_files.append("Gaia_G_Var_calc_lhoods_" + region_choice + ".fits")
	# WISE Excess
	source_files.append("WISE_W1W2_lhoods_" + region_choice + ".fits")
	# H-K Excess
	source_files.append("HKExcess_lhoods_" + region_choice + ".fits")

	# Ask the user whether to include Isochronal Age
	print('Include Isochronal Age:')
	IsoAge_input = input('y/n: ')
	# Convert source data choice into sky bin length
	IsoAge_list = {
			"y": True,
			"Y": True,
			"n": False,
			"N": False
			}
	IsoAge_choice = IsoAge_list.get(IsoAge_input, "Invalid")
	if IsoAge_choice == "Invalid":
		print("Invalid choice")
		sys.exit()
	elif IsoAge_choice == True:
		# Isochronal Age
		source_files.append("IsoAge_lhoods_" + region_choice + ".fits")
	
	# Ask the user whether to include H-alpha Excess
	print('Include H-alpha Excess:')
	Ha_input = input('y/n: ')
	# Convert source data choice into sky bin length
	Ha_list = {
			"y": True,
			"Y": True,
			"n": False,
			"N": False
			}
	Ha_choice = Ha_list.get(Ha_input, "Invalid")
	if Ha_choice == "Invalid":
		print("Invalid choice")
		sys.exit()
	elif Ha_choice == True:
		# H-alpha Excess
		source_files.append("HaExcess_lhoods_" + region_choice + ".fits")

	# Add a blank space
	print()

	return source_files, region_choice

def load_and_check_data(source_filepath, source_files):
	global gaia_g_var_data, ir_excess_data

	for file_idx in range(0,len(source_files)):

		print("Loading file: " + source_files[file_idx])
		load_table = read_observations(source_filepath + source_files[file_idx])

		print("Sorting data")
		load_table.sort('source_id')

		# Combine the input data into a single structured array
		if file_idx==0:
			# Start off the array with the first file
			data_table = load_table
		else:
			# Merge subsequent files into the master array
			# Check arrays are the same size
			master_size = np.size(data_table,0)
			new_size = np.size(load_table,0)
			if master_size != new_size:
				print("Aborting as file " + source_files[file_idx] + " has " + str(new_size) + " rows when expecting " + str(master_size))
				sys.exit()
			# Check the source_id in each row matches the source_id in the same row of other table
			check_mask = data_table['source_id'] == load_table['source_id']
			check_sum = np.sum(check_mask)
			if check_sum != master_size:
				print("Aborting as only " + str(check_sum) + " source_ids out of " + str(master_size) + " match from the file " + source_files[file_idx])
				sys.exit()
			#print("Element 10: " + str(gaia_g_var_data[10]['source_id']) + " " + str(ir_excess_data[10]['source_id']) + " " + str(gaia_g_var_data[10]['source_id']==ir_excess_data[10]['source_id']))
			#rfn.drop_fields(temp_data, 'source_id')
			#print("Removing source_id as this should only appear once in the table")
			#temp_table = remove_field_from_structured_array(load_table, 'source_id')
			#load_table.remove_columns('source_id')
			print("Merging tables")
			#data_table = rfn.merge_arrays((data_table, temp_table), asrecarray=True, flatten=True)
			for col_name in load_table.dtype.names:
				# Do not add source_id as this should only appear once in the table
				if col_name != 'source_id':
					data_table[col_name] = load_table[col_name]
			#print(data_table.dtype.descr)
			#print(str(np.size(data_table)))
			#print(data_table[10])

	
	# Remove reddening for other filters
	for col_name in ['E_rmi', 'E_rmha', 'A_r', 'E_JmH_2MASS', 'E_HmK_2MASS', 'E_JmH_UKIDSS', 'E_HmK_UKIDSS', 'E_JmH_VVV', 'E_HmK_VVV']:
		if col_name in data_table.columns:
			data_table.remove_columns([col_name])
	
	# Add a Bayes flags
	data_table['bayes_calc'] = False
	data_table['bayes_zero_norm'] = False
	# Add counts of properties in end calculation
	data_table['properties_calc'] = 0

	# Add a blank space
	print()

	return data_table

# Read a fits file of observations with a standard set of columns
def read_observations(file_name):
	print("Reading observations file ...")
	fits_file = fits.open(file_name)
	#fits_file.info()
	fits_data = fits_file[1].data
	fits_file.close
	# Convert the data into Table format
	table_data = Table(fits_data)
	return table_data

def output_full_fits(output_filename, data_table, source_files):
	print("Writing output fits file ...")
	# Delete file if it exists
	if os.path.isfile(output_filename):
		os.remove(output_filename)
	# Create fits header
	hdr = fits.Header()
	hdr['run_ver'] = str(run_version)
	hdr['datetime'] = str(datetime.datetime.now())
	for file_idx in range(0,len(source_files)):
		source_files[file_idx]
		hdr['data_f' + str(file_idx)] = source_files[file_idx]
	primary_hdu = fits.PrimaryHDU(header=hdr)
	# Convert table data into HDU
	hdu_data = fits.table_to_hdu(data_table)
	# Combine header with data
	hdul = fits.HDUList([primary_hdu, hdu_data])
	# Create the output file
	hdul.writeto(output_filename)

def remove_interim_columns():
	global bayes_table
	print("Removing interim processing columns ...")
	for col_name in bayes_table.dtype.names:
		if '_like' in col_name or 'normalisation' in col_name or '_capped' in col_name or '_calc' in col_name or 'gaia_g_var_' in col_name or col_name == 'bayes_zero_norm':
			if 'properties_' not in col_name:
				bayes_table.remove_columns([col_name])
	# Remove observation feature value columns
	for col_name in ['wise_W1mW2', 'par_corr_r_dered', 'rmi_r_iso_log10age', 'older_max_iso_age', 'younger_min_iso_age', 'good_classiiyso_fit', 'rmha_offset', 'rmi_dered', 'Source', 'HK_offset', 'phot_g_mean_mag', 'gaia_g_fluxsd_over_flux', 'sigma_i', 'model_n_obs', 'Umg_offset', 'brighter_max_iso_mag', 'fainter_max_iso_mag']:
		if col_name in bayes_table.columns:
			bayes_table.remove_columns([col_name])

def bayes_calculation():
	global bayes_table
	print("Calculating Naive Bayes Classifier posteriors ...")

	# Set the priors
	bayes_table['prior_classiiyso'] = prior_classiiyso
	bayes_table['prior_other'] = prior_other

	# Iterate through the list to create the prior * likelihood
	print("Multiplying the priors with the likelihoods ...")
	for b_cat in bayes_cats:
		# Initiliase the priorxlikelihoods
		bayes_table['prior_x_like' + b_cat] = bayes_table['prior' + b_cat]
		for class_name in bayes_table.dtype.names:
			if (b_cat + '_like') in class_name:
				# Bayes calculation
				bayes_table['prior_x_like' + b_cat] = bayes_table['prior_x_like' + b_cat] * bayes_table[class_name]

	# Calculate the number of properties contributing to each Bayes calculation
	for class_name in bayes_table.dtype.names:
		if ('_calc') in class_name and class_name!='bayes_calc' and class_name!='properties_calc':
			bayes_table['properties_calc'] = bayes_table['properties_calc'] + bayes_table[class_name]

	# Calculate the normalisation factor
	print("Calculating the normalisation factors ...")
	bayes_table['normalisation'] = 0.
	for b_cat in bayes_cats:
		bayes_table['normalisation'] += bayes_table['prior_x_like' + b_cat]

	# Calculate the Bayes Posteriors
	print("Calculating the posteriors ...")
	# Check for zero normalisation
	zero_mask = bayes_table['normalisation'] == 0.
	bayes_table['bayes_zero_norm'][zero_mask] = True
	bayes_table['bayes_calc'][~zero_mask] = True
	for b_cat in bayes_cats:
		bayes_table['posterior' + b_cat] = bayes_table['prior' + b_cat]
		bayes_table['posterior' + b_cat][~zero_mask] = bayes_table['prior_x_like' + b_cat][~zero_mask] / bayes_table['normalisation'][~zero_mask]

	# Add a blank space
	print()



#------------------------------------------------------------
# Start of master routine
#------------------------------------------------------------

print('')
print('-----------------------------------------------------------------------------')
print('------------------------- Naive Bayes Classifier ----------------------------')
print('-----------------------------------------------------------------------------')
print('')

# Get the parameters from the user
data_files, region_name = user_inputs()

# Main path for the run
main_path = run_filepath + region_name + "/"

# Read the observation data files
bayes_table = load_and_check_data(main_path, data_files)

# Calculate the Bayes Posteriors
bayes_calculation()

# Create output file with all columns
output_filename = main_path + "Bayes_full_" + region_name + ".fits"
output_full_fits(output_filename, bayes_table, data_files)

# Remove interim processing columns not needed by Naive Bayes Classifier
remove_interim_columns()

# Create smaller output fits file for Naive Bayes Classifier
output_filename = main_path + "Bayes_" + region_name + ".fits"
output_full_fits(output_filename, bayes_table, data_files)

#------------------------------------------------------------
# End of master routine
#------------------------------------------------------------