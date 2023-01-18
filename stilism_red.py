# ------------------------------------------------------------
# For calculating Interstellar Reddening to list of sources
# using Stilism
# https://stilism.obspm.fr
# Capitanio et al, A&A, 606, A65 (2017)
# ------------------------------------------------------------
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

from astropy.io import fits
from astropy.table import Table
#from astropy import units
#from astropy.coordinates import SkyCoord

import datetime

# Parameters common to the suite of Python scripts
from run_parameters import *

# Hard coded parameters
#----------------------
# Resolution of reddening on the sky
sky_res = 0.5
# Distance resolution in pc
d_res = 5


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

	# Convert source data choice into the IGAPS file
	data_file = igaps_file_list.get(source_data, "Invalid")
	if data_file == "Invalid":
		print("Invalid choice")
		sys.exit()
	igaps_filename = data_file

	print('--------------------------------------------------------------------')

	return region_choice, igaps_filename

# Read a fits file of observations with a standard set of columns
def read_observations(file_name, file_type):
	print("Reading " + file_type + " file ...")
	fits_file = fits.open(file_name)
	#fits_file.info()
	fits_data = fits_file[1].data
	fits_file.close
	# Convert the data into Table format
	table_data = Table(fits_data)

	# If the coordinates file then add reddening columns
	if (file_type=="Coords"):
		# Add reddening columns
		table_data['E_BmV'] = np.nan
		table_data['E_BmV_calc'] = False
		table_data['E_BmV_uncertainty_min'] = np.nan
		table_data['E_BmV_uncertainty_max'] = np.nan
		table_data['distance_uncertainty'] = np.nan
		table_data['E_rmi'] = np.nan
		table_data['E_rmha'] = np.nan
		table_data['A_r'] = np.nan
		table_data['E_JmH_2MASS'] = np.nan
		table_data['E_HmK_2MASS'] = np.nan
		table_data['E_JmH_UKIDSS'] = np.nan
		table_data['E_HmK_UKIDSS'] = np.nan
		table_data['E_JmH_VVV'] = np.nan
		table_data['E_HmK_VVV'] = np.nan
		table_data['rmi'] = np.nan
		table_data['rmi_dered'] = np.nan
		table_data['suspect_reddening'] = False

	# If the IGAPS file then add columns for derived INT AB magnitudes
	elif (file_type=="IGAPS"):
		# Add columns for derived INT AB magnitudes and colours
		IGAPS_mask = ~np.isnan(table_data['imagAB'])
		table_data['imagAB_INT'] = np.nan
		table_data['rmi'] = np.nan
		table_data['imagAB_INT'][IGAPS_mask] = table_data['imagAB'][IGAPS_mask] - 0.06 * (table_data['rImagAB'][IGAPS_mask] - table_data['imagAB'][IGAPS_mask])
		table_data['rmi'][IGAPS_mask] = table_data['rImagAB'][IGAPS_mask] - table_data['imagAB_INT'][IGAPS_mask]

		# Remove unnecessary columns
		table_data.remove_columns(['Name', 'HamagAB'])

	return table_data

def output_full_fits(output_filename, data_table, source_file_obs):
	print("Writing output fits file ...")
	# Delete file if it exists
	if os.path.isfile(output_filename):
		os.remove(output_filename)
	# Create fits header
	hdr = fits.Header()
	hdr['run_ver'] = str(run_version)
	hdr['datetime'] = str(datetime.datetime.now())
	hdr['file_obs'] = source_file_obs
	hdr['stilism'] = 'These reddenings come from Stilism'
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
		if col not in ['source_id', 'E_BmV', 'E_rmi', 'E_rmha', 'A_r', 'E_JmH_2MASS', 'E_HmK_2MASS', 'E_JmH_UKIDSS', 'E_HmK_UKIDSS', 'E_JmH_VVV', 'E_HmK_VVV', 'suspect_reddening']:
			col_to_delete.append(col)
	obs_data.remove_columns(col_to_delete)

# ------------------------------------------------------------------------------------------------------------
# Read Stilism Reddening Files
# ------------------------------------------------------------------------------------------------------------

# Read the reddening files
def load_readening_files(l_min, l_max, b_min, b_max, d_max, stilism_filepath):

	l_size = int((l_max-l_min)/sky_res)
	b_size = int((b_max-b_min)/sky_res)
	d_size = int(d_max/d_res)+1
	red_cube = np.zeros((l_size,b_size,d_size), dtype='f8')
	red_calc_cube = np.zeros((l_size,b_size,d_size), dtype='bool')
	red_uncert_min_cube = np.zeros((l_size,b_size,d_size), dtype='f8')
	red_uncert_max_cube = np.zeros((l_size,b_size,d_size), dtype='f8')
	d_uncert_cube = np.zeros((l_size,b_size,d_size), dtype='f8')
	d_temp = np.zeros((d_size), dtype='f8')

	# Loop through the files by Galactic longitude and latitude
	for l_idx in range(0,l_size):
		for b_idx in range(0,b_size):

			# Convert indexes to Galactic coordinates
			l_val = gal_index_to_value(l_min, l_idx)
			b_val = gal_index_to_value(b_min, b_idx)

			# Read the data file
			stilism_file = 'l' + '{:.2f}'.format(l_val) + 'b' + '{:.2f}'.format(b_val) + ".csv"
			print(' Loading Stilism file: ' + stilism_file)
			stilism_data = np.loadtxt(stilism_filepath + stilism_file, dtype={'names': ('distance', 'reddening', 'distance_uncertainty', 'reddening_uncertainty_min', 'reddening_uncertainty_max'), 'formats': ('u4', 'f8', 'f8', 'f8', 'f8')}, delimiter=',', skiprows=1, usecols=(0, 1, 2, 3, 4))

			# Select data within the maxium distance
			mask = stilism_data['distance'] <= d_max
			# Deal with files that don't go to the maximum distance
			last_distance_index = np.sum(mask)
			# By default populate the reddening with the maximum in the file
			max_red = np.max(stilism_data[mask]['reddening'])
			red_cube[l_idx,b_idx,:] = max_red
			# Update arrays with values fro the file up to the greatest distance in the file
			red_cube[l_idx,b_idx,:last_distance_index] = stilism_data[mask]['reddening']
			red_calc_cube[l_idx,b_idx,:last_distance_index] = True
			red_uncert_min_cube[l_idx,b_idx,:last_distance_index] = stilism_data[mask]['reddening_uncertainty_min']
			red_uncert_max_cube[l_idx,b_idx,:last_distance_index] = stilism_data[mask]['reddening_uncertainty_max']
			d_uncert_cube[l_idx,b_idx,:last_distance_index] = stilism_data[mask]['distance_uncertainty']

	return red_cube, red_calc_cube, red_uncert_min_cube, red_uncert_max_cube, d_uncert_cube


# Convert Galactic coordinate index into value in degrees
def gal_index_to_value(gal_min, gal_index):
	# Note adds half a resolution element to place in the centre of the box in sky coordinates
	# This is only needed when going from index to value
	# as the integer conversion from value to index rounds down
	gal_value = gal_min + (gal_index*sky_res) + (sky_res/2.)
	return gal_value

# Not used
# Convert Galactic coordinate value in degrees into index
# def gal_value_to_index(gal_min, gal_value):
# 	gal_index = int((gal_value - gal_min)/sky_res)
# 	return gal_index

# Convert distance value in pc into index
def d_value_to_index(d_value):
	d_index = int(d_value/d_res)
	return d_index

# Convert Galactic coordinate value in degrees into index
def gal_array_to_index(gal_min, gal_array):
	gal_index = ((gal_array - gal_min)/sky_res).astype(int)
	return gal_index

# Convert distance value in pc into index
def d_array_to_index(d_array):
	d_index = (d_array/d_res).astype(int)
	return d_index


# ------------------------------------------------------------------------------------------------------------
# Calculation Reddening to Sources
# ------------------------------------------------------------------------------------------------------------

# Loop through regions of sky
def add_reddening_to_sources(l_min, b_min):
	global obs_data, reddening_calc_cube, reddening_cube, reddening_min_cube, reddening_max_cube, d_uncert_cube

	print("Adding E(B-V) reddening to sources from STILISM ...")

	# Add temporary columns to index the sources to the reddening map
	#obs_data['l_idx'] = ((obs_data['l'] - l_min)/sky_res).astype(int)
	#obs_data['b_idx'] = ((obs_data['b'] - b_min)/sky_res).astype(int)
	obs_data['l_idx'] = gal_array_to_index(l_min, obs_data['l'])
	obs_data['b_idx'] = gal_array_to_index(b_min, obs_data['b'])
	obs_data['d_idx'] = d_array_to_index(1000/obs_data['parallax_zero_corrected'])

	# Add the reddening and other data
	# The reddening is for the centre of the voxel
	# This is achieved as by loading the galactic coordinates from STILISM plus half the sky resolution as the coordinates are converted to array indexes
	# And taking the average of the near and back edges in distance of the voxel
	obs_data['E_BmV'] = (reddening_cube[obs_data['l_idx'],obs_data['b_idx'],obs_data['d_idx']]+reddening_cube[obs_data['l_idx'],obs_data['b_idx'],obs_data['d_idx']+1])/2
	obs_data['E_BmV_calc'] = reddening_calc_cube[obs_data['l_idx'],obs_data['b_idx'],obs_data['d_idx']]
	obs_data['E_BmV_uncertainty_min'] = reddening_min_cube[obs_data['l_idx'],obs_data['b_idx'],obs_data['d_idx']]
	obs_data['E_BmV_uncertainty_max'] = reddening_max_cube[obs_data['l_idx'],obs_data['b_idx'],obs_data['d_idx']]
	obs_data['distance_uncertainty'] = d_uncert_cube[obs_data['l_idx'],obs_data['b_idx'],obs_data['d_idx']]

	# Calculate reddening values for other filters
	print("Calculating reddening for IGAPS, 2MASS, UKIDSS and VVV")
	nonan_mask = ~np.isnan(obs_data['E_BmV'])
	# Reddening for IGAPS
	obs_data['E_rmi'][nonan_mask] = E_BV_to_ArI_p0 + E_BV_to_ArI_p1 * obs_data['E_BmV'][nonan_mask] + E_BV_to_ArI_p2 * obs_data['E_BmV'][nonan_mask] * obs_data['E_BmV'][nonan_mask]
	obs_data['E_rmha'][nonan_mask] = E_BV_to_ArHa_p0 + E_BV_to_ArHa_p1 * obs_data['E_BmV'][nonan_mask] + E_BV_to_ArHa_p2 * obs_data['E_BmV'][nonan_mask] * obs_data['E_BmV'][nonan_mask]
	obs_data['A_r'][nonan_mask] = obs_data['E_BmV'][nonan_mask] * E_BmV_to_A_r
	# Reddening for (H-K)
	obs_data['E_JmH_2MASS'][nonan_mask] = T2MASS_E_BV_to_E_JH * obs_data['E_BmV'][nonan_mask]
	obs_data['E_HmK_2MASS'][nonan_mask] = T2MASS_E_BV_to_E_HK * obs_data['E_BmV'][nonan_mask]
	obs_data['E_JmH_UKIDSS'][nonan_mask] = UKIDSS_E_BV_to_E_JH * obs_data['E_BmV'][nonan_mask]
	obs_data['E_HmK_UKIDSS'][nonan_mask] = UKIDSS_E_BV_to_E_HK * obs_data['E_BmV'][nonan_mask]
	obs_data['E_JmH_VVV'][nonan_mask] = UKIDSS_E_BV_to_E_JH * obs_data['E_BmV'][nonan_mask]
	obs_data['E_HmK_VVV'][nonan_mask] = UKIDSS_E_BV_to_E_HK * obs_data['E_BmV'][nonan_mask]


# ------------------------------------------------------------------------------------------------------------
# Identify Suspect Reddening Values
# ------------------------------------------------------------------------------------------------------------

def Identify_suspect_reddening():
	global obs_data, igaps_data

	print("")
	print("Identifying suspect reddening values")
	print("------------------------------------")

	# Check arrays are the same size
	red_size = np.size(obs_data,0)
	igaps_size = np.size(igaps_data,0)
	if red_size != igaps_size:
		print("Aborting as file reddening file has " + str(red_size) + " rows and IGAPS file has " + str(igaps_size))
		sys.exit()
	
	# Merge IGAPS rmi into the reddening table
	print("Sorting reddening data")
	obs_data.sort('source_id')
	print("Sorting IGAPS data")
	igaps_data.sort('source_id')
	# Check the source_id in each row matches the source_id in the same row of other table
	check_mask = obs_data['source_id'] == igaps_data['source_id']
	check_sum = np.sum(check_mask)
	if check_sum != red_size:
		print("Aborting as only " + str(check_sum) + " source_ids out of " + str(red_size) + " match between the reddening file and IGAPS")
		sys.exit()
	print("Adding (r-i) to the reddening table")
	obs_data['rmi'] = igaps_data['rmi']


	# Select rows with IGAPS data
	nonan_mask = ~np.isnan(obs_data['rmi'])
	obs_data['rmi_dered'][nonan_mask] = obs_data['rmi'][nonan_mask] - obs_data['E_rmi'][nonan_mask]

	# Identify suspect reddening values
	print("Identifying suspect reddening as (r-i)<=-0.3 where IGAPS data exists")
	suspect_mask = np.logical_and( obs_data['rmi_dered'] <= -0.3 , ~np.isnan(obs_data['rmi_dered']) )
	obs_data['suspect_reddening'][suspect_mask] = True

	# Add a blank space
	print()

#------------------------------------------------------------
# Start of master routine
#------------------------------------------------------------

print('')
print('-----------------------------------------------------------------------------')
print('--------------------------- Calculate Reddening -----------------------------')
print('-----------------------------------------------------------------------------')
print('')

# Get the parameters from the user
region_name, igaps_file = user_inputs()

# Main path for the run
main_path = run_filepath + region_name + "/"


# Reddening determination
# -----------------------

# Import filename
obs_file = "Coords_" + region_name + ".fits"

# Read the coordinates data file
obs_data = read_observations(main_path + obs_file, "Coords")

# Determine the limits of the data for use with loops
l_min = math.floor(np.min(obs_data['l']))
l_max = math.ceil(np.max(obs_data['l']))
b_min = math.floor(np.min(obs_data['b']))
b_max = math.ceil(np.max(obs_data['b']))
d_max = math.ceil(np.max(1000/obs_data['parallax_zero_corrected'])/d_res)*d_res

# Read isochrone files for to fit the extinction
reddening_cube, reddening_calc_cube, reddening_min_cube, reddening_max_cube, d_uncert_cube = load_readening_files(l_min, l_max, b_min, b_max, d_max, stilism_filepath)

# Update the observations with reddening
add_reddening_to_sources(l_min, b_min)


# Identify suspect reddening values
# ---------------------------------

# Read the IGAPS data file
igaps_data = read_observations(main_path + igaps_file, "IGAPS")

# Identify suspect reddening values
Identify_suspect_reddening()

# Create output files
# -------------------

# Create output fits file with all data
output_filename = main_path + "Reddening_full_" + region_name + ".fits"
output_full_fits(output_filename, obs_data, obs_file)

# Remove interim processing columns not needed by Naive Bayes Classifier
remove_interim_columns()

# Create smaller output fits file for Naive Bayes Classifier
output_filename = main_path + "Reddening_" + region_name + ".fits"
output_full_fits(output_filename, obs_data, obs_file)

#print(reddening_cube)
#print(reddening_cube[0,0,0])
#print(reddening_cube[1,1,1])
#print(str(l_max) + ' ' + str(b_max) + ' ' + str(d_max))
#print(str(gal_value_to_index(l_min,l_max)) + ' ' + str(gal_value_to_index(b_min,b_max)) + ' ' + str(d_value_to_index(d_max)))
#print(reddening_cube[gal_value_to_index(l_min,l_max),gal_value_to_index(b_min,b_max),d_value_to_index(d_max)])

#------------------------------------------------------------
# End of master routine
#------------------------------------------------------------