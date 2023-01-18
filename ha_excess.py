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


# IPHAS
# -----
# Standard empirical isochrone with slope and intercept that is displaced in r-i
# Values from fitting Norther Galactic Plane IPHAS data within 500 pc (Bailer-Jones 2018) 
# and E(B-V)=0.05 from calculation for Isochronal Age fitting in this suite of programs
#ms_slope = 0.4492
#ms_intercept = 0.1189


# IGAPS AB with I corrected to INT
# --------------------------------
# Empirical isochrone from good quality IGAPS photometry matched to the WEAVE NGPv2 catalogue
# With E(B-V)<0.01 from Stilism matched to nearest 0.5 desgrees on the sky
# Limited by (r-i) colour to linear region 0.2<(r-i)<1.3
# Iterative fit such that data points lie within +/-0.1 in (r-Ha) of the best fit line
ms_slope = 0.4716
ms_intercept = -0.0083

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
	data_file = igaps_file_list.get(source_data, "Invalid")
	if data_file == "Invalid":
		print("Invalid choice")
		sys.exit()
	obs_data_file = data_file

	print('--------------------------------------------------------------------')

	return obs_data_file, region_choice

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
		# Add columns for derived INT AB magnitudes and colours
		IGAPS_mask = ~np.isnan(table_data['imagAB'])
		table_data['imagAB_INT'] = np.nan
		table_data['rmi'] = np.nan
		table_data['rmha'] = np.nan
		table_data['imagAB_INT'][IGAPS_mask] = table_data['imagAB'][IGAPS_mask] - 0.06 * (table_data['rImagAB'][IGAPS_mask] - table_data['imagAB'][IGAPS_mask])
		table_data['rmi'][IGAPS_mask] = table_data['rImagAB'][IGAPS_mask] - table_data['imagAB_INT'][IGAPS_mask]
		table_data['rmha'][IGAPS_mask] = table_data['rImagAB'][IGAPS_mask] - table_data['HamagAB'][IGAPS_mask]

		# Remove unnecessary columns
		table_data.remove_columns(['Name'])

	return table_data

def add_reddening():
	global igaps_data, red_data

	# First sort the data so combine rows for the same sources
	print("Sorting data")
	igaps_data.sort('source_id')
	red_data.sort('source_id')

	# Check arrays are the same size
	main_size = np.size(igaps_data,0)
	red_size = np.size(red_data,0)
	if main_size != red_size:
		print("Aborting as number of rows in the main data (" + str(main_size) + ") is different to the number of rows in the reddening data (" + str(red_size) + ").")
		print("Files need to match for combination to be successful.")
		sys.exit()

	# Check the source_id in each row matches the source_id in the same row of other table
	check_mask = igaps_data['source_id'] == red_data['source_id']
	check_sum = np.sum(check_mask)
	if check_sum != main_size:
		print("Aborting as only " + str(check_sum) + " source_ids out of " + str(main_size) + " match between the main file and the reddening file.")
		print("Files need to match for combination to be successful.")
		sys.exit()
	
	print("Merging reddening into the main data.")
	igaps_data['E_rmi'] = red_data['E_rmi']
	igaps_data['E_rmha'] = red_data['E_rmha']
	igaps_data['suspect_reddening'] = red_data['suspect_reddening']

	# Add a blank space
	print()

# Calculate (r-Ha) excess
def ha_excess():
	global igaps_data

	# Add results columns
	print("Calculating r-Ha excess...")

	# Select rows with IGAPS data
	nonan_mask = ~np.isnan(igaps_data['imagAB'])

	# Deredden each source
	igaps_data['rmi_dered'] = np.nan
	igaps_data['rmha_dered'] = np.nan
	igaps_data['rmi_dered'][nonan_mask] = igaps_data['rmi'][nonan_mask] - igaps_data['E_rmi'][nonan_mask]
	igaps_data['rmha_dered'][nonan_mask] = igaps_data['rmha'][nonan_mask] - igaps_data['E_rmha'][nonan_mask]
	# r-Ha offset
	igaps_data['rmha_offset'] = np.nan
	igaps_data['rmha_offset'][nonan_mask] = igaps_data['rmha_dered'][nonan_mask] - ((ms_slope * igaps_data['rmi_dered'][nonan_mask]) + ms_intercept)

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
	hdr['E_rI'] = "E(r-I) = " + str(E_BV_to_ArI_p0) + " + " + str(E_BV_to_ArI_p1) + "*E(B-V)" + " + "  + str(E_BV_to_ArI_p2) + "*E(B-V)^2"
	hdr['E_rHa'] = "E(r-Ha) = " + str(E_BV_to_ArHa_p0) + " + " + str(E_BV_to_ArHa_p1) + "*E(B-V)" + " + "  + str(E_BV_to_ArHa_p2) + "*E(B-V)^2"
	hdr['MSfit'] = "(r-Ha) = " + str(ms_slope) + " * (r-i) + " + str(ms_intercept)
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
		if col not in ['source_id', 'l', 'b', 'parallax_zero_corrected', 'rmi_dered', 'rmha_offset', 'suspect_reddening']:
			col_to_delete.append(col)
	igaps_data.remove_columns(col_to_delete)

#------------------------------------------------------------
# Start of master routine
#------------------------------------------------------------

print('')
print('-----------------------------------------------------------------------------')
print('------------------------ Calculate H-alpha Excess ---------------------------')
print('-----------------------------------------------------------------------------')
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

# Add the reddening to the main data
add_reddening()

# Calculate (r-Ha) excess
ha_excess()

# Create output fits file with all data
output_filename = main_path + "HaExcess_full_" + region_name + ".fits"
output_full_fits(output_filename, igaps_data, obs_file)

# Remove interim processing columns not needed by Naive Bayes Classifier
remove_interim_columns()

# Create smaller output fits file for Naive Bayes Classifier
output_filename = main_path + "HaExcess_" + region_name + ".fits"
output_full_fits(output_filename, igaps_data, obs_file)

#------------------------------------------------------------
# End of master routine
#------------------------------------------------------------