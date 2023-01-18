# ------------------------------------------------------------
# For calculating variability in G mag using Gaia DR2
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
from scipy import stats
from scipy.optimize import minimize

from astropy.io import fits
from astropy.table import Table
#from astropy import units
#from astropy.coordinates import SkyCoord

import datetime

# Custom code for calculating the variabiilty function pdf
# Written by Prof Tim Naylor with some edits by Andrew Wilson
import dP_dr_fast

# Parameters common to the suite of Python scripts
from run_parameters import *

# Hard coded parameters
#----------------------
# None

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

	# Convert source data choice into a data file path
	data_file = gaia_g_var_file_list.get(source_data, "Invalid")
	if data_file == "Invalid":
		print("Invalid choice")
		sys.exit()

	# Ask which version of likelihoods to use
	lhood_v = input('Likelihood version: ')
	# Convert to parameter file
	par_file = 'GaiaGVariabilityParameters_V' + str(lhood_v) + '.csv'


	print('--------------------------------------------------------------------')

	return data_file, region_choice, par_file

# Read a fits file of observations with a standard set of columns
def read_observations(file_name):
	print("Reading observations file ...")
	fits_file = fits.open(file_name)
	#fits_file.info()
	fits_data = fits_file[1].data
	fits_file.close
	# Convert the data into Table format
	table_data = Table(fits_data)
	# Calculated column needed in both calibration and observation files
#	table_data['gaia_g_fluxsd'] = (table_data['phot_g_n_obs']**0.5) * table_data['phot_g_mean_flux_error']
	table_data['gaia_g_fluxsd_over_flux'] = (table_data['phot_g_n_obs']**0.5) * (table_data['phot_g_mean_flux_error']/table_data['phot_g_mean_flux_corrected'])
#	table_data['gaia_g_fluxsd_over_sqrtflux'] = (table_data['phot_g_n_obs']**0.5) * (table_data['phot_g_mean_flux_error']/(table_data['phot_g_mean_flux']**0.5))
#	table_data['gaia_g_fmsv'] = (table_data['phot_g_n_obs'] * table_data['phot_g_mean_flux_error']**2.) / table_data['phot_g_mean_flux']**2.
	# TOPCAT
	# VARmag - y-axis
	# 2.5*log10((gaia_g_var_fluxsd_over_sqrtflux/sqrt(phot_g_mean_flux))+1)
	# d(VARratio)/d(VARmag) - y-axix
	# (ln(10)/2.5)*(gaia_g_var_fluxsd_over_sqrtflux+sqrt(phot_g_mean_flux))
	# Non-YSO 15<=G<16
	# 0.8/(0.087401*sqrt(2*PI))*exp(-pow(x-0.767382,2)/(2*pow(0.087401,2)))+0.2*(140*2.5*log10(E))/sqrt(11836)*pow(sqrt(11836)/(sqrt(pow(x,2)-pow(0.767386,2))+sqrt(11836)),(140*2.5*log10(E))+1)
	# Double exponential function
	# 0.8/(0.087401*sqrt(2*PI))*exp(-pow(x-0.767382,2)/(2*pow(0.087401,2)))+0.19*(100*2.5*log10(E))/sqrt(11836)*pow(sqrt(11836)/(sqrt(pow(x,2)-pow(0.767386,2))+sqrt(11836)),(100*2.5*log10(E))+1)+0.01*(18*2.5*log10(E))/sqrt(11836)*pow(sqrt(11836)/(sqrt(pow(x,2)-pow(0.767386,2))+sqrt(11836)),(18*2.5*log10(E))+1)
	# Add results columns
	table_data['sigma_i'] = np.nan
	# Effective number of observations
	table_data['model_n_obs'] = np.nan
	mask_G_less_16 = table_data['phot_g_mean_mag_corrected'] < 16
	table_data['model_n_obs'][mask_G_less_16] = ( ( table_data['phot_g_mean_mag_corrected'][mask_G_less_16] * np.power(table_data['phot_g_n_obs'][mask_G_less_16],0.93) - 13.0 * np.power(table_data['phot_g_n_obs'][mask_G_less_16],0.93) ) / 3. ).astype(int)
	table_data['model_n_obs'][~mask_G_less_16] = np.power(table_data['phot_g_n_obs'][~mask_G_less_16],0.93).astype(int)
	# Defaults for other values
	table_data['classiiyso_scale_length'] = np.nan
	table_data['classiiyso_frac_non_var'] = np.nan
	table_data['other_scale_length'] = np.nan
	table_data['other_frac_non_var'] = np.nan
	table_data['gaia_g_var_classiiyso_l_raw'] = np.nan
	table_data['gaia_g_var_other_l_raw'] = np.nan
	table_data['gaia_g_var_classiiyso_like'] = 1.
	table_data['gaia_g_var_other_like'] = 1.
	table_data['gaia_g_var_classiiyso_capped'] = False
	table_data['gaia_g_var_other_capped'] = False
	table_data['gaia_g_var_calc'] = False
	return table_data

def output_full_fits(output_filename, data_table, source_file, para_file):
	print("Writing output fits file ...")
	# Delete file if it exists
	if os.path.isfile(output_filename):
		os.remove(output_filename)
	# Create fits header
	hdr = fits.Header()
	hdr['run_ver'] = str(run_version)
	hdr['datetime'] = str(datetime.datetime.now())
	hdr['file_obs'] = source_file
	hdr['para_f'] = para_file
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
		if col not in ['source_id', 'phot_g_mean_mag_corrected', 'gaia_g_fluxsd_over_flux', 'sigma_i', 'model_n_obs', 'gaia_g_var_classiiyso_like', 'gaia_g_var_other_like', 'gaia_g_var_classiiyso_capped', 'gaia_g_var_other_capped', 'gaia_g_var_calc']:
			col_to_delete.append(col)
	obs_data.remove_columns(col_to_delete)

# Read a csv file of likelihood data based on criteria in the observation file
def read_Gvar_parameters(file_name):
	print("Reading likelihoods file ...")
	read_data = np.loadtxt(file_name, dtype={'names': ('start_mag', 'end_mag', 'inst_noise', 'other_scale', 'classiiyso_scale', 'other_frac_non_var', 'classiiyso_frac_non_var'), 'formats': ('f4', 'f4', 'f8', 'f8', 'f8', 'f4', 'f4')}, delimiter=',', skiprows=1, usecols=(0, 1, 2, 3, 4, 5, 6))
	return read_data


#------------------------------------------------------------
# Variability calculations
#------------------------------------------------------------

# Calculate the category likelihoods based on variability
def calc_lhoods():
	global obs_data, para_data

	# Use the Gaussian approximation to the chi-squared distribution
	use_gaussian_approximation = True

	# Calculate how many sources will have probability calculations
	#total_mask = np.logical_and(np.min(para_data['start_mag'])<=obs_data['phot_g_mean_mag'] , obs_data['phot_g_mean_mag']<np.max(para_data['end_mag']))
	total_count = np.size(obs_data,0)
	count_idx = 0

	# Calculate the upper prior ratio
	upper_prior_ratio = prior_other / prior_classiiyso

	print('Calculating likelihoods for ' + str(total_count) + ' sources ...')

	# Loop through the observations in the bin
	for idx, obs_row in enumerate(obs_data):

		# Select magnitude bin
		mag_mask = np.logical_and(para_data['start_mag']<=obs_row['phot_g_mean_mag_corrected'] , obs_row['phot_g_mean_mag_corrected']<para_data['end_mag'])
		
		# Check observation data within the magnitude range
		# And enforce number of G-band observations >= 70
		if np.sum(mag_mask)==1 and obs_row['phot_g_n_obs']>=70:
		
			# Indicate calculation performed
			obs_data['gaia_g_var_calc'][idx] = True

			# Add model parameters to the observation file
			obs_data['sigma_i'][idx] = para_data['inst_noise'][mag_mask]
			obs_data['other_scale_length'][idx] = para_data['other_scale'][mag_mask]
			obs_data['other_frac_non_var'][idx] = para_data['other_frac_non_var'][mag_mask]
			obs_data['classiiyso_scale_length'][idx] = para_data['classiiyso_scale'][mag_mask]
			obs_data['classiiyso_frac_non_var'][idx] = para_data['classiiyso_frac_non_var'][mag_mask]

			# Give a count every 10,0000 observations
			if count_idx > 0 and count_idx % 10000 == 0:
				print('  Processed ' + str(count_idx) + ' of ' + str(total_count))

			# Calculate the non-YSO likelihood
			Other_lhood = dP_dr_fast.variability_pdf(obs_row['gaia_g_fluxsd_over_flux'], obs_row['sigma_i'], obs_row['other_scale_length'], 
							obs_row['other_frac_non_var'], obs_row['model_n_obs'], use_gaussian_approximation)

			# Calculate the YSO likelihood
			ClassIIYSO_lhood = dP_dr_fast.variability_pdf(obs_row['gaia_g_fluxsd_over_flux'], obs_row['sigma_i'], obs_row['classiiyso_scale_length'], 
							obs_row['classiiyso_frac_non_var'], obs_row['model_n_obs'], use_gaussian_approximation)

			# Clean likelihoods for problem cases of
			#  Zero likelihood
			#  Unphysical ratio between YSO and non-YSO likelihoods

			# Deal with rare case of zero likelihoods
			if ClassIIYSO_lhood == 0.0 or Other_lhood == 0.0:
				if ClassIIYSO_lhood == 0.0 and Other_lhood == 0.0:
					ClassIIYSO_lhood_fin = 1.
					Other_lhood_fin = 1.
				elif ClassIIYSO_lhood == 0.0:
					ClassIIYSO_lhood_fin = 0.
					Other_lhood_fin = 1.
				else:
					ClassIIYSO_lhood_fin = 1.
					Other_lhood_fin = 0.

			else:

				# Enforce power limit on Class II YSOs
				if ClassIIYSO_lhood/Other_lhood > upper_prior_ratio:
					ClassIIYSO_lhood_fin = Other_lhood * upper_prior_ratio
					obs_data['gaia_g_var_classiiyso_capped'][idx] = True
				else:
					ClassIIYSO_lhood_fin = ClassIIYSO_lhood

				# Enforce power limit on Not YSO
				if Other_lhood/ClassIIYSO_lhood > upper_prior_ratio:
					Other_lhood_fin = ClassIIYSO_lhood * upper_prior_ratio
					obs_data['gaia_g_var_other_capped'][idx] = True
				else:
					Other_lhood_fin = Other_lhood


			# Update data table
			obs_data['gaia_g_var_other_l_raw'][idx] = Other_lhood
			obs_data['gaia_g_var_classiiyso_l_raw'][idx] = ClassIIYSO_lhood
			obs_data['gaia_g_var_other_like'][idx] = Other_lhood_fin
			obs_data['gaia_g_var_classiiyso_like'][idx] = ClassIIYSO_lhood_fin

		count_idx += 1



#------------------------------------------------------------
# Start of master routine
#------------------------------------------------------------

print('')
print('-----------------------------------------------------------------------------')
print('------------------ Calculate Gaia G magnitude Variability -------------------')
print('-----------------------------------------------------------------------------')
print('')

# Get the parameters from the user
obs_file, region_name, param_file = user_inputs()

# Main path for the run
main_path = run_filepath + region_name + "/"

# Read the observation data file
obs_data = read_observations(main_path + obs_file)

# Read the parameter data
para_data = read_Gvar_parameters(main_path + param_file)

# Calculate the category likelihoods based on variability
calc_lhoods()

# Create output fits file with all data
output_filename = main_path + "Gaia_G_Var_calc_full_" + region_name + ".fits"
output_full_fits(output_filename, obs_data, obs_file, param_file)

# Remove interim processing columns not needed by Naive Bayes Classifier
remove_interim_columns()

# Create smaller output fits file for Naive Bayes Classifier
output_filename = main_path + "Gaia_G_Var_calc_lhoods_" + region_name + ".fits"
output_full_fits(output_filename, obs_data, obs_file, param_file)

#------------------------------------------------------------
# End of master routine
#------------------------------------------------------------