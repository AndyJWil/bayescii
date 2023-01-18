# ------------------------------------------------------------
# For calculating Age by fitting isochrones to IGAPS CMDs
# Author: Andrew Wilson
# Institution: University of Exeter, Astrophysics Group
# Version date: 18 January 2023
# University email: aw648@exeter.ac.uk
# Personal email: andyjwilson_uk@hotmail.com
# ------------------------------------------------------------

import sys
import os

import numpy as np
import math
from scipy import interpolate

from astropy.io import fits
from astropy.table import Table
#from astropy import units
#from astropy.coordinates import SkyCoord

import datetime

# Parameters common to the suite of Python scripts
from run_parameters import *

# Hard coded parameters
#----------------------

# Start log10 age of isochrones
iso_log_age_start = 5.
# Extinction steps
ex_steps=120

# Radial Basis Function Parameters
rbf_function = 'thin_plate'
rbs_smooth = 0

# Maximum reddening avaiable from the isochrones
max_reddening = 5.999

# Convert E(B-V) reddening to index
red_to_index = 20.

# Interpolation and data limits
max_interp_age = 7.0
min_interp_age = 5.0
max_interp_mass = 2.00
min_interp_mass = 0.05
colour_interp_tolerance = 0.1
mag_interp_tolerance = 0.1


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

	# Convert source data choice into the data file
	obs_data_file = igaps_file_list.get(source_data, "Invalid")
	if obs_data_file == "Invalid":
		print("Invalid choice")
		sys.exit()

	print('--------------------------------------------------------------------')
	print('')

	return obs_data_file, region_choice

# # Calculate the distance corrected magnitude, called absolute magnitude though this can be misleading as it is not corrected for reddening
# def abs_mag(data_t):
# 	absolute_magnitude = data_t['phot_g_mean_mag'] - 5 * np.log10(data_t['r_est']) + 5
# 	return absolute_magnitude

# Read a fits file of observations with a standard set of columns
def read_observations(file_name, file_type):
	print("Reading " + file_type + " file ...")
	fits_file = fits.open(file_name)
	#fits_file.info()
	fits_data = fits_file[1].data
	fits_file.close
	# Convert the data into Table format
	table_data = Table(fits_data)
	# If the IGAPS file then add columns for derived INT AB magnitudes
	if (file_type=="IGAPS"):
		IGAPS_mask = ~np.isnan(table_data['imagAB'])
		table_data['imagAB_INT'] = np.nan
		table_data['rmi'] = np.nan
		table_data['imagAB_INT'][IGAPS_mask] = table_data['imagAB'][IGAPS_mask] - 0.06 * (table_data['rImagAB'][IGAPS_mask] - table_data['imagAB'][IGAPS_mask])
		table_data['rmi'][IGAPS_mask] = table_data['rImagAB'][IGAPS_mask] - table_data['imagAB_INT'][IGAPS_mask]
		# Add calculated columns
		table_data['par_corr_r'] = np.nan
		table_data['par_corr_r_dered'] = np.nan
		#table_data['par_corr_r_dered'] = np.nan
		table_data['par_corr_r'][IGAPS_mask] = table_data['rImagAB'][IGAPS_mask] + 5. * np.log10(table_data['parallax_zero_corrected'][IGAPS_mask]/1000.) + 5.
		table_data['rmi_r_iso_log10age'] = np.nan
		table_data['rmi_r_iso_mass'] = np.nan
#		table_data['isochrone_interpolation_result'] = "No calculation                          "
		# Add result columns
		table_data['igaps_data'] = False
		table_data['brighter_max_iso_mag'] = False
		table_data['fainter_max_iso_mag'] = False
		table_data['above_iso_col'] = False
		table_data['below_iso_col'] = False
		table_data['back_interp_colour_diff'] = np.nan
		table_data['back_interp_mag_diff'] = np.nan
		table_data['good_except_back_interp_colour'] = False
		table_data['good_except_back_interp_mag'] = False
		table_data['above_max_iso_red'] = False
		table_data['older_max_iso_age'] = False
		table_data['younger_min_iso_age'] = False
		table_data['above_max_iso_mass'] = False
		table_data['below_min_iso_mass'] = False
		table_data['good_classiiyso_fit'] = False

		# Remove unnecessary columns
		table_data.remove_columns(['Name', 'HamagAB'])

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
	hdr['iso_int'] = 'Dotter 2008 z=0.018'
	hdr['iso_atmo'] = 'BTSettl+Bell'
	hdr['iso_ext'] = 'Fitzpatrick & Massa (2007)'
	primary_hdu = fits.PrimaryHDU(header=hdr)
	# Convert table data into HDU
	hdu_data = fits.table_to_hdu(data_table)
	# Combine header with data
	hdul = fits.HDUList([primary_hdu, hdu_data])
	# Create the output file
	hdul.writeto(output_filename)

def remove_interim_columns():
	global igaps_data
	print("Removing interim processing columns ...")
	col_to_delete = []
	for col in igaps_data.columns:
		if col not in ['source_id', 'l', 'b', 'parallax_zero_corrected', 'par_corr_r_dered', 'rmi_r_iso_log10age', 'rmi_r_iso_mass', 'igaps_data', \
			'brighter_max_iso_mag', 'fainter_max_iso_mag', 'above_iso_col', 'below_iso_col', 'back_interp_colour_diff', 'back_interp_mag_diff', \
			'good_except_back_interp_colour', 'good_except_back_interp_mag', 'above_max_iso_red', 'older_max_iso_age', 'younger_min_iso_age', \
			'above_max_iso_mass', 'below_min_iso_mass', 'good_classiiyso_fit', 'suspect_reddening']:
			col_to_delete.append(col)
	igaps_data.remove_columns(col_to_delete)

def add_reddening():
	global igaps_data, red_data

	# First sort the data so combine rows for the same sources
	print("Sorting data")
	igaps_data.sort('source_id')
	red_data.sort('source_id')

	# Check arrays are the same size
	igaps_size = np.size(igaps_data,0)
	red_size = np.size(red_data,0)
	if igaps_size != red_size:
		print("Aborting as number of rows in the igaps data (" + str(igaps_size) + ") is different to the number of rows in the reddening data (" + str(red_size) + ").")
		print("Files need to match for combination to be successful.")
		sys.exit()

	# Check the source_id in each row matches the source_id in the same row of other table
	check_mask = igaps_data['source_id'] == red_data['source_id']
	check_sum = np.sum(check_mask)
	if check_sum != igaps_size:
		print("Aborting as only " + str(check_sum) + " source_ids out of " + str(igaps_size) + " match between the igaps file and the reddening file.")
		print("Files need to match for combination to be successful.")
		sys.exit()
	
	print("Merging reddening into the main data.")
	igaps_data['E_BmV'] = red_data['E_BmV']
	igaps_data['A_r'] = red_data['A_r']
	igaps_data['suspect_reddening'] = red_data['suspect_reddening']

	# Add a column for dereddened parallax corrected r
	IGAPS_mask = ~np.isnan(igaps_data['imagAB'])
	igaps_data['par_corr_r_dered'][IGAPS_mask] = igaps_data['par_corr_r'][IGAPS_mask] - igaps_data['A_r'][IGAPS_mask]

	# Add a blank space
	print()


# ------------------------------------------------------------------------------------------------------------
# Read Isochrones and Convert into Radial Basis Functions
# ------------------------------------------------------------------------------------------------------------

# Read in the Isochrones at multiple ages and reddenings for age fit
def read_multi_age_isochrones(iso_filepath):
	print("Reading isochrones for age fit ...")
	
	#iso_steps = 3250
	log_age_steps=21
	log_age_start=int(iso_log_age_start*10)
	int_count = 0
	for log_age in range(0,log_age_steps):
		log_age_number = (log_age+log_age_start)/10
		log_age_string = '0' + '{:.4f}'.format(log_age_number)
		print(' log10(age)=' + log_age_string)
		for E_BV in range(0,ex_steps):
			E_BV_string = '{:.2f}'.format(E_BV/red_to_index)
			#print(' E(B-V)=' + E_BV_string)
			# Read the data from the isochrone file into a numpy ndarray
			iso_file = iso_filepath + "r_r-i_age_5_7_Dotter08_BTSettl_Bell_E" + E_BV_string + "/r_r-i_" + log_age_string + ".iso"
			iso_data = np.loadtxt(iso_file, dtype={'names': ('colour', 'mag', 'mass', 'age', 'isgood'), 'formats': ('f8', 'f8', 'f8', 'f8', 'S4')}, skiprows=5, usecols=(0, 1, 2, 6, 8))
			# Select the good data and below 2 Msun to avoid isochrones overlapping
			mask_good = iso_data['isgood'] == b'T'
			mask_mass = iso_data['mass'] <= max_interp_mass
			mask = np.logical_and( mask_good, mask_mass )
			good_data = iso_data[mask]
			# Reduce the oversampled isochrone to 1 in every 50 pointa
			selected_data = good_data[0:len(good_data):50]
			# Add the extinction as an extra column
			new_data = np.empty(selected_data.shape, dtype=selected_data.dtype.descr + [('E_BmV', 'f8')])
			for name in selected_data.dtype.names:
				new_data[name] = selected_data[name]
			new_data['E_BmV'] = float(E_BV) / red_to_index
			# Counter used to know when to combine data
			int_count += 1
			if int_count > 1:
				combined_data = np.concatenate((combined_data, new_data), axis=0)
			else:
				combined_data = new_data
	return combined_data

# Convert isochrone CMDs by E(B-V) into bivariate B-splines to interpolate age and mass
def interp_iso_cmd(iso_array):
	print("Creating Radial Basis Functions to interpolate the isochrone CMDs ...")
	# Initialise Radial Basis Function arrays by extinction
	rbf_age_set = []
	rbf_mass_set = []
	rbf_colour_set = []
	rbf_mag_set = []
	# Loop through extinction values
	for ext_index in range(0,ex_steps):
		E_BV = float(ext_index) / red_to_index
		iso_e_bv_mask = iso_array['E_BmV'] == E_BV
		iso_e_bv = iso_array[iso_e_bv_mask]
		# Radial Basis Functions to interpolate age and mass
		# A thin plate function is used as this gives the best fit to the isochrones
		# Zero smoothing is used as this forces the function to go through the provided isochrone points
		rbf_age = interpolate.Rbf(iso_e_bv['colour'], iso_e_bv['mag'], iso_e_bv['age'],function=rbf_function,smooth=rbs_smooth)
		rbf_mass = interpolate.Rbf(iso_e_bv['colour'], iso_e_bv['mag'], iso_e_bv['mass'],function=rbf_function,smooth=rbs_smooth)
		rbf_colour = interpolate.Rbf(iso_e_bv['age'], iso_e_bv['mass'], iso_e_bv['colour'],function=rbf_function,smooth=rbs_smooth)
		rbf_mag = interpolate.Rbf(iso_e_bv['age'], iso_e_bv['mass'], iso_e_bv['mag'],function=rbf_function,smooth=rbs_smooth)
		rbf_age_set.append(rbf_age)
		rbf_mass_set.append(rbf_mass)
		rbf_colour_set.append(rbf_colour)
		rbf_mag_set.append(rbf_mag)

	return rbf_age_set, rbf_mass_set, rbf_colour_set, rbf_mag_set

# ------------------------------------------------------------------------------------------------------------
# Calculation Age
# ------------------------------------------------------------------------------------------------------------

# Estimate age of stars using splines created from isochrones
def find_best_age():
	print("Finding age of young stars ...")
	global rbf_age_arr, rbf_mass_arr, rbf_colour_arr, rbf_mag_arr, iso_yso

	# Loop through the data calculating age where there is igaps data
	for igaps_idx in range(0,len(igaps_data)):

		if igaps_idx % 10000 == 0:
			print("Processing " + str(igaps_idx) + " of " + str(len(igaps_data)))

		# Clear interpolation message
		interp_msg = ""

		# Check whether there is igaps data
		if ~np.isnan(igaps_data['rmi'][igaps_idx]):

			# igaps data
			igaps_data['igaps_data'][igaps_idx] = True

			# Reset Good YSO fit tracker
			good_classiiyso_fit = True

			# Extract the reddening
			E_BV = igaps_data['E_BmV'][igaps_idx]
			# Convert the reddening into an index for the RBF function
			ext_idx = int(E_BV * red_to_index)
			#print(str(E_BV) + '=' + str(ext_idx))
			# Only proceed with age calculation where reddening has been determined
			if E_BV >= 0.:

				# If reddening is above the maximum then do not calculate
				if E_BV > max_reddening:
					igaps_data['above_max_iso_red'][igaps_idx] = True

				else:

					# Check sources lies within the isochrone magnitude and column range to within +/-0.1
					iso_e_bv_mask = iso_yso['E_BmV'] == ext_idx / int(red_to_index)
					#print(iso_e_bv.dtype.descr)
					#print(iso_e_bv)
					#print()
					if igaps_data['par_corr_r'][igaps_idx] > np.max(iso_yso['mag'][iso_e_bv_mask]) + 0.1:
						igaps_data['fainter_max_iso_mag'][igaps_idx] = True
					elif igaps_data['par_corr_r'][igaps_idx] < np.min(iso_yso['mag'][iso_e_bv_mask]) - 0.1:
						igaps_data['brighter_max_iso_mag'][igaps_idx] = True
					elif igaps_data['rmi'][igaps_idx] > np.max(iso_yso['colour'][iso_e_bv_mask]) + 0.1:
						igaps_data['above_iso_col'][igaps_idx] = True
					elif igaps_data['rmi'][igaps_idx] < np.min(iso_yso['colour'][iso_e_bv_mask]) - 0.1:
						igaps_data['below_iso_col'][igaps_idx] = True
		
					# Source within isochrone magnitude and colour eange
					else:
		
						# Interpolate the age and mass using the Radial Basis Functions
						age_rbf = rbf_age_arr[ext_idx](igaps_data['rmi'][igaps_idx], igaps_data['par_corr_r'][igaps_idx])
						mass_rbf = rbf_mass_arr[ext_idx](igaps_data['rmi'][igaps_idx], igaps_data['par_corr_r'][igaps_idx])
						# Record the age and mass interpolated values
						igaps_data['rmi_r_iso_log10age'][igaps_idx] = age_rbf
						igaps_data['rmi_r_iso_mass'][igaps_idx] = mass_rbf

						# Check the log(age) is within the bounds of the source isochrones
						# Extrapoloated older than age of the isochrones
						if age_rbf > max_interp_age:
							igaps_data['older_max_iso_age'][igaps_idx] = True
							good_classiiyso_fit = False
						# Extrapoloated younger than age of the isochrones plus a small tolerance
						elif age_rbf < min_interp_age:
							igaps_data['younger_min_iso_age'][igaps_idx] = True
							good_classiiyso_fit = False

						# Check the mass is within the bounds of the source isochrones plus a small tolerance
						if mass_rbf > max_interp_mass:
							igaps_data['above_max_iso_mass'][igaps_idx] = True
							good_classiiyso_fit = False
						elif mass_rbf < min_interp_mass:
							igaps_data['below_min_iso_mass'][igaps_idx] = True
							good_classiiyso_fit = False

						# Back interpolate the colour and magnitude from addition Radial Basis Functions
						back_colour_rbf = rbf_colour_arr[ext_idx](age_rbf, mass_rbf)
						back_mag_rbf = rbf_mag_arr[ext_idx](age_rbf, mass_rbf)
						# Then calculate the difference between the original and back interpolated colour and magnitude
						colour_diff = abs(igaps_data['rmi'][igaps_idx] - back_colour_rbf)
						mag_diff = abs(igaps_data['par_corr_r'][igaps_idx] - back_mag_rbf)
						igaps_data['back_interp_colour_diff'][igaps_idx] = colour_diff
						igaps_data['back_interp_mag_diff'][igaps_idx] = mag_diff

						# Good on all initial checks
						if good_classiiyso_fit:

							# Check the colour and mag back interpolated differences are within tolerance
							if colour_diff >= colour_interp_tolerance:
								igaps_data['good_except_back_interp_colour'][igaps_idx] = True
								good_classiiyso_fit = False
							if mag_diff >= mag_interp_tolerance:
								igaps_data['good_except_back_interp_mag'][igaps_idx] = True
								good_classiiyso_fit = False
			
							# Still good after back calculations
							if good_classiiyso_fit:
								igaps_data['good_classiiyso_fit'][igaps_idx] = True



#------------------------------------------------------------
# Start of master routine
#------------------------------------------------------------

print('')
print('-----------------------------------------------------------------------------')
print('------------------------ Calculate Isochronal Age ---------------------------')
print('-----------------------------------------------------------------------------')
print('')
print('Start: ' + str(datetime.datetime.now()))
print('')

# Get the parameters from the user
obs_file, region_name = user_inputs()

# Main path for the run
main_path = run_filepath + region_name + "/"

# Read the observation data file
igaps_data = read_observations(main_path + obs_file, "IGAPS")

# Read the reddening data
red_filename = "Reddening_" + region_name + ".fits"
red_data = read_observations(main_path + red_filename, "reddening")

# Add the reddening to the igaps data
add_reddening()

# Read isochrone files for extinction and age
iso_yso = read_multi_age_isochrones(isoage_iso_filepath)

# Convert isochrone CMDs by E(B-V) into Radial Basis Functions to interpolate age and mass
rbf_age_arr, rbf_mass_arr, rbf_colour_arr, rbf_mag_arr = interp_iso_cmd(iso_yso)

print('')
print('Time check: ' + str(datetime.datetime.now()))
print('')

# Find the best age for each source if it is within the region defined by the YSO isochrones
find_best_age()

# Create output fits file with all data
output_filename = main_path + "IsoAge_Full_" + region_name + ".fits"
output_full_fits(output_filename, igaps_data)

# Remove interim processing columns not needed by Naive Bayes Classifier
remove_interim_columns()

# Create smaller output fits file for Naive Bayes Classifier
output_filename = main_path + "IsoAge_" + region_name + ".fits"
output_full_fits(output_filename, igaps_data)

print('')
print('Finish: ' + str(datetime.datetime.now()))
print('')

#------------------------------------------------------------
# End of master routine
#------------------------------------------------------------