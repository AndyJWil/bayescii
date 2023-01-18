# --------------------------------------------------------------
# For converting the full output file into a publishable version
# Author: Andrew Wilson
# Institution: University of Exeter, Astrophysics Group
# Version date: 18 Jan 2023
# University email: aw648@exeter.ac.uk
# Personal email: andyjwilson_uk@hotmail.com
# --------------------------------------------------------------

import sys
import os

import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy import units

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

	# Ask the user which sky area to use
	print('What do you wish to do?')
	print('E - Export field names and descriptions from fits file')
	print('U - Create published FITS file with updated columns names and field descriptions')
	print('V - Create VizieR data file')
	print('R - Create file format section of VizieR read me file')
	user_process_type = input('Process: ')
	# Convert source data choice into sky bin length
	process_list = {
			"E": "E",
			"e": "E",
			"U": "U",
			"u": "U",
			"V": "V",
			"v": "V",
			"r": "R",
			"R": "R"
			}
	selected_process_type = process_list.get(user_process_type, "Invalid")
	if selected_process_type == "Invalid":
		print("Invalid choice")
		sys.exit()
	
	return region_choice, selected_process_type

# Read a fits file of data with a standard set of columns
def read_fits(file_name):
	print("Reading data file ...")
	fits_file = fits.open(file_name)
	#fits_file.info()
	#table_hdr = fits_file[1].header
	fits_data = fits_file[1].data
	fits_file.close
	# Convert the data into Table format
	table_data = Table(fits_data)

	return table_data

# Read a fits file of data with a standard set of columns
def read_fits_header(file_name):
	print("Reading data file ...")
	fits_file = fits.open(file_name)
	table_hdr = fits_file[1].header
	#print(table_hdr.comments['TTYPE1'])
	#print(fits_file[1].header.comments['TTYPE1'])
	fits_file.close

	return table_hdr

# Read a fits file of observations with a standard set of columns
def read_field_csv(file_name):
	print("Reading csv file field names and descriptions ...")
	read_data = np.loadtxt(file_name, dtype={'names': ('Old name', 'VizieR name', 'FITS data type', 'Dictionary data type', 'Numpy data type', 'VizieR data type', 'VizieR data range', 'VizieR data nullable', 'Unit', 'Symbol', 'FITS comment', 'VizieR description', 'Long description'), 'formats': ('U80', 'U80', 'U2', 'U40', 'U40', 'U40', 'U20', 'U1', 'U10', 'U20', 'U80', 'U80', 'U80')}, delimiter=',', skiprows=1, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
	#read_data = np.loadtxt(file_name, dtype={'names': ('Classification', 'Probability'), 'formats': ('S6', 'f8')}, delimiter=',', skiprows=2, usecols=(0, 1))
	return read_data

# Function to output a 2D array from 2 lists as a file
def export_csv(fname, out_data, header=''):
	print("Writing output file ...")
	# Write file
	with open(fname, 'w') as out_file:
		# If header is provided then use it as the first row
		if header != '':
			#out_file.write(','.join(header) + '\n');
			out_file.write(header + '\n');

		# Loop through each row in the data array
		for out_row in out_data:
			# Create output row
			out_row = ','.join(out_row) + '\n'
			# Output the row
			out_file.write(out_row);
	out_file.closed

# Function to output an AstroPy fits table in VizieR format one row and field at a time
def export_vizier_data(fname, out_data, updated_headers):
	print("Preparing VizieR data file ...")
	num_fields = updated_headers.shape[0]
	first_field_name = updated_headers[0]['VizieR name']
	num_rows = out_data[first_field_name].shape[0]
	row_length = 0
	# Loop through the columns to work out row length (width of the data file)
	for int_index in range(0, num_fields):
		if int_index > 0:
			row_length += 1
		if updated_headers[int_index]['VizieR data type'][0] == 'F':
			decimal_pos = updated_headers[int_index]['VizieR data type'].find('.')
			field_len = int(updated_headers[int_index]['VizieR data type'][1:decimal_pos])
			row_length += field_len
		elif updated_headers[int_index]['VizieR data type'][0] == 'I':
			if updated_headers[int_index]['Numpy data type'] == 'bool':
				row_length += 1
			else:
				field_len = int(updated_headers[int_index]['VizieR data type'][1:])
				row_length += field_len
		elif updated_headers[int_index]['VizieR data type'][0] == 'A':
			field_len = int(updated_headers[int_index]['VizieR data type'][1:])
			row_length += field_len
		
	print(row_length)
	# Create blank array for the file export
	#  1 row per data row with a width to match the new output file
	formatted_row = np.empty((num_rows), dtype='U'+str(row_length))
	
	# Loop through the data columns with sub loops to add the column to each row of the output array
	for int_index in range(0, num_fields):

		# Record the VizieR name for this data column
		field_name = updated_headers[int_index]['VizieR name']
		data_field = out_data[field_name]
		print(field_name)
		#print(data_field.shape)
		#print(data_field.dtype)
		print(updated_headers[int_index]['VizieR data type'])

		# Add a space between data columns
		if int_index > 0:
			for int_row in range(0,data_field.shape[0]):
				#print(formatted_row[int_row])
				formatted_row[int_row] += ' '

		# If statements to deal with different data types
		# F = Floating point
		if updated_headers[int_index]['VizieR data type'][0] == 'F':
			# Extract the two portions of the variable format F6.3 to field_dec=6 and field_int=3
			decimal_pos = updated_headers[int_index]['VizieR data type'].find('.')
			field_int = int(updated_headers[int_index]['VizieR data type'][1:decimal_pos])
			field_dec = int(updated_headers[int_index]['VizieR data type'][decimal_pos+1:])
			int_pad = field_int-field_dec-1
			# Loop through the data rows
			for int_row in range(0,data_field.shape[0]):
				# Test for null
				if np.isnan(data_field[int_row]):
					formatted_row[int_row] += ' ' * field_int
				else:
					formatted_row[int_row] += np.format_float_positional(data_field[int_row], unique=False, precision=field_dec, trim='k', pad_left=int_pad)
		# I = Integer
		elif updated_headers[int_index]['VizieR data type'][0] == 'I':
			# Check for Boolean
			if updated_headers[int_index]['Numpy data type'] == 'bool':
				# Loop through the data rows
				for int_row in range(0,data_field.shape[0]):
					if data_field[int_row] == True:
						formatted_row[int_row] += '1'
					else:
						formatted_row[int_row] += '0'
			else:
				field_len = int(updated_headers[int_index]['VizieR data type'][1:])
				# Loop through the data rows
				for int_row in range(0,data_field.shape[0]):
					str_field = str(data_field[int_row])
					int_pad = field_len - len(str_field)
					extra_spaces = ' ' * int_pad
					formatted_row[int_row] += extra_spaces + str_field
		# A = Character/Text
		elif updated_headers[int_index]['VizieR data type'][0] == 'A':
			int_pad = int(updated_headers[int_index]['VizieR data type'][1:])
			# Loop through the data rows
			for int_row in range(0,data_field.shape[0]):
				str_field = str(data_field[int_row])
				extra_spaces = int_pad - len(str_field)
				formatted_row[int_row] += str_field.rjust(int_pad, ' ')
		else:
			# Data type not found
			print(updated_headers[int_index]['VizieR name'] + ' - Data type not found')


		print('---')

	print('Output rows:')
	print(formatted_row.shape)
	np.savetxt(fname, formatted_row, fmt='%s', delimiter=' ')

# Create the VizieR file using a simple method
# Not used as didn't work out how to deal with NULLs
# def export_vizier_simple(fname, out_data, updated_headers):
# 	print("Preparing VizieR data file ...")
# 	num_fields = updated_headers.shape[0]
# 	first_field_name = updated_headers[0]['VizieR name']
# 	fmt_str = ''
# 	fmt_arr = np.empty(num_fields, dtype='U10')
# 	# Loop through the columns
# 	for int_index in range(0, num_fields):
# 		if int_index > 0:
# 			fmt_str += ' '
# 		if updated_headers[int_index]['VizieR data type'][0] == 'F':
# 			decimal_pos = updated_headers[int_index]['VizieR data type'].find('.')
# 			field_int = int(updated_headers[int_index]['VizieR data type'][1:decimal_pos])
# 			fmt_str += "%"+updated_headers[int_index]['VizieR data type'][1:]+"f"
# 			fmt_arr[int_index] = '%'+updated_headers[int_index]['VizieR data type'][1:]+'f'
# 		elif updated_headers[int_index]['VizieR data type'][0] == 'I':
# 			if updated_headers[int_index]['Numpy data type'] == 'bool':
# 				fmt_str += "%1d"
# 				fmt_arr[int_index] = '%1d'
# 			else:
# 				int_pad = int(updated_headers[int_index]['VizieR data type'][1:])
# 				fmt_str += "%"+str(int_pad)+"d"
# 				fmt_arr[int_index] = '%'+str(int_pad)+'d'
# 		elif updated_headers[int_index]['VizieR data type'][0] == 'A':
# 			int_pad = int(updated_headers[int_index]['VizieR data type'][1:])
# 			fmt_str += "%"+str(int_pad)+"s"
# 			fmt_arr[int_index] = '%'+str(int_pad)+'s'
# 	fmt_str += ""
# 	print(fmt_str)
# 	np.savetxt(fname, out_data, fmt=fmt_str, delimiter=' ')

# Create the file format section of the VizieR read me file
def export_vizier_readme(fname, updated_headers):
	print("Preparing VizieR read me file ...")

	# Initialise variables
	num_fields = updated_headers.shape[0]
	first_field_name = updated_headers[0]['VizieR name']
	row_length = 0
	vizier_readme = np.empty((num_fields), dtype='U110')
	max_label_length = len(max(updated_headers['VizieR name'], key=len))
	
	# First loop through the columns - Creates ReadMe but without wrapping lines that exceed 80 characters
	for int_index in range(0, num_fields):
		# Add gap of a single space between fields
		if int_index > 0:
			row_length += 1
		# Add the start column of the field
		pad_spaces = (3 - len(str(row_length+1))) * ' '
		vizier_readme[int_index] += pad_spaces + str(row_length+1)
		# If statement to deal with different data types
		# F = Floating point
		if updated_headers[int_index]['VizieR data type'][0] == 'F':
			decimal_pos = updated_headers[int_index]['VizieR data type'].find('.')
			field_len = int(updated_headers[int_index]['VizieR data type'][1:decimal_pos])
			row_length += field_len
		# I = Integer
		elif updated_headers[int_index]['VizieR data type'][0] == 'I':
			if updated_headers[int_index]['Numpy data type'] == 'bool':
				row_length += 1
			else:
				field_len = int(updated_headers[int_index]['VizieR data type'][1:])
				row_length += field_len
		# A = Character/Text
		elif updated_headers[int_index]['VizieR data type'][0] == 'A':
			field_len = int(updated_headers[int_index]['VizieR data type'][1:])
			row_length += field_len
		# Add the end column of the field
		pad_spaces = (3 - len(str(row_length))) * ' '
		vizier_readme[int_index] += '-' + pad_spaces + str(row_length)
		# Add the data type
		vizier_readme[int_index] += ' '
		pad = 6 - len(updated_headers[int_index]['VizieR data type'])
		vizier_readme[int_index] += updated_headers[int_index]['VizieR data type'] + (pad*' ')
		# Add the units
		vizier_readme[int_index] += ' '
		if updated_headers[int_index]['Unit'] == '':
			vizier_unit = '---   '
		else:
			pad = 6 - len(updated_headers[int_index]['Unit'])
			vizier_unit = updated_headers[int_index]['Unit'] + (pad * ' ')
		vizier_readme[int_index] += vizier_unit
		# Add the label
		vizier_readme[int_index] += ' '
		pad = max_label_length - len(updated_headers[int_index]['VizieR name'])
		vizier_readme[int_index] += updated_headers[int_index]['VizieR name'] + (pad*' ')
		# Add explanation -
		#  Limits + Nullable + Description
		vizier_readme[int_index] += ' '
		vizier_explanation = ''
		if updated_headers[int_index]['VizieR data range'] != '':
			vizier_explanation = updated_headers[int_index]['VizieR data range']
		if updated_headers[int_index]['VizieR data nullable'] != '':
			vizier_explanation += updated_headers[int_index]['VizieR data nullable']
		if vizier_explanation != '':
			vizier_explanation += ' '
		vizier_explanation += updated_headers[int_index]['VizieR description']
		vizier_readme[int_index] += vizier_explanation

	# Second loop through the columns - Creates second row to wrap description onto a new line where the column width exceeds 80 characters
	vizier_readme_wrap = np.empty((num_fields), dtype='U80')
	int_wrap = 0
	for int_index in range(0, num_fields):
		if len(vizier_readme[int_index])>80:
			vizier_readme_wrap = np.append(vizier_readme_wrap, [''], axis=0)
			last_space = vizier_readme[int_index][:80].rfind(' ')
			temp_left = vizier_readme[int_index][:last_space]
			temp_right = vizier_readme[int_index][last_space:]
			vizier_readme_wrap[int_wrap] = temp_left
			int_wrap += 1
			vizier_readme_wrap[int_wrap] = (' ' * 48) + temp_right
		else:
			vizier_readme_wrap[int_wrap] = vizier_readme[int_index]
		int_wrap += 1

	# Create the ReadMe data table dictionary
	np.savetxt(fname, vizier_readme_wrap, fmt='%s', delimiter=' ')

def extract_fields(field_count):
	global obs_header, output_table

	# Loop through the fields
	for int_index in range(0, field_count):
		output_table[int_index][0] = obs_header['TTYPE' + str(int_index + 1)]
		output_table[int_index][1] = obs_header['TFORM' + str(int_index + 1)]
		output_table[int_index][2] = obs_header.comments['TTYPE' + str(int_index + 1)]

# Create a new table to change the fits columns
def update_fits_columns(input_file_name, header_descriptions):
	global fits_data, updated_headers

	print("Updating fits columns ...")

	# Create new table with the first columns
	#first_cols = ['source_id', 'dr2_source_id', 'random_index', 'ra', 'dec', 'l', 'b', 'parallax_zero_corrected']
	first_cols = ['source_id']
	new_fits = fits_data[first_cols]

	# Initial fields
	new_fits['source_id'].name = 'GaiaEDR3SourceId'
	#new_fits['dr2_source_id'] = fits_data['dr2_source_id']
	#new_fits['random_index'] = fits_data['random_index']
	new_fits['RAdeg'] = fits_data['ra'].astype('float64')
	new_fits['DEdeg'] = fits_data['dec'].astype('float64')
	new_fits['GLON'] = fits_data['l'].astype('float64')
	new_fits['GLAT'] = fits_data['b'].astype('float64')
	new_fits['plx'] = fits_data['parallax_zero_corrected'].astype('float64')
	new_fits['EBmV'] = fits_data['E_BmV'].astype('float32')
	new_fits['SuspectReddening'] = fits_data['suspect_reddening'].astype('bool')
	#new_fits['suspect_reddening'].format('I1')
	new_fits['GmagCorr'] = fits_data['phot_g_mean_mag_corrected'].astype('float32')
	new_fits['rmagParCorrDered'] = fits_data['par_corr_r_dered'].astype('float32')
	new_fits['rmiDered'] = fits_data['rmi_dered'].astype('float32')

	# Bayes calculations
	new_fits['FeaturesCalc'] = fits_data['properties_calc'].astype('int16')
	new_fits['PriorCII'] = fits_data['prior_classiiyso'].astype('float32')
	new_fits['PriorOther'] = fits_data['prior_other'].astype('float32')
	new_fits['PosteriorCII'] = fits_data['posterior_classiiyso'].astype('float64')
	new_fits['PosteriorOther'] = fits_data['posterior_other'].astype('float64')

	# Variability feature
	new_fits['GFracStdDevFlux'] = fits_data['gaia_g_fluxsd_over_flux'].astype('float64')
	#new_fits['sigma_i'] = fits_data['sigma_i']
	#new_fits['model_n_obs'] = fits_data['model_n_obs']
	new_fits['GvarCIILike'] = fits_data['gaia_g_var_classiiyso_like'].astype('float64')
	new_fits['GvarOtherLike'] = fits_data['gaia_g_var_other_like'].astype('float64')
	new_fits['GvarLikeRatio'] = fits_data['gaia_g_var_classiiyso_like'].astype('float64') / fits_data['gaia_g_var_other_like'].astype('float64')
	new_fits['GvarCalc'] = fits_data['gaia_g_var_calc'].astype('bool')
	new_fits['GvarCIICapped'] = fits_data['gaia_g_var_classiiyso_capped'].astype('bool')
	new_fits['GvarOtherCapped'] = fits_data['gaia_g_var_other_capped'].astype('bool')
	# WISE feature
	new_fits['W1mW2'] = fits_data['wise_W1mW2'].astype('float32')
	new_fits['W1mW2CIILike'] = fits_data['w1mw2_classiiyso_like'].astype('float64')
	new_fits['W1mW2OtherLike'] = fits_data['w1mw2_other_like'].astype('float64')
	new_fits['W1mW2LikeRatio'] = fits_data['w1mw2_classiiyso_like'].astype('float64') / fits_data['w1mw2_other_like'].astype('float64')
	new_fits['W1mW2Calc'] = fits_data['w1mw2_calc'].astype('bool')
	new_fits['W1mW2CIICapped'] = fits_data['w1mw2_classiiyso_capped'].astype('bool')
	new_fits['W1mW2OtherCapped'] = fits_data['w1mw2_other_capped'].astype('bool')
	# (H-K) excess feature
	new_fits['HmKExSource'] = fits_data['Source'].astype('U7')
	new_fits['HmKEx'] = fits_data['HK_offset'].astype('float32')
	new_fits['HmKExCIILike'] = fits_data['hkexcess_classiiyso_like'].astype('float64')
	new_fits['HmKExOtherLike'] = fits_data['hkexcess_other_like'].astype('float64')
	new_fits['HmKExLikeRatio'] = fits_data['hkexcess_classiiyso_like'].astype('float64') / fits_data['hkexcess_other_like'].astype('float64')
	new_fits['HmKExCalc'] = fits_data['hkexcess_calc'].astype('bool')
	new_fits['HmKExCIICapped'] = fits_data['hkexcess_classiiyso_capped'].astype('bool')
	new_fits['HmKExOtherCapped'] = fits_data['hkexcess_other_capped'].astype('bool')
	# (r-Ha) feature
	new_fits['rmHaEx'] = fits_data['rmha_offset'].astype('float32')
	new_fits['rmHaExCIILike'] = fits_data['rmha_offset_classiiyso_like'].astype('float64')
	new_fits['rmHaExOtherLike'] = fits_data['rmha_offset_other_like'].astype('float64')
	new_fits['rmHaExLikeRatio'] = fits_data['rmha_offset_classiiyso_like'].astype('float64') / fits_data['rmha_offset_other_like'].astype('float64')
	new_fits['rmHaExCalc'] = fits_data['rmha_offset_calc'].astype('bool')
	new_fits['rmHaExCIICapped'] = fits_data['rmha_offset_classiiyso_capped'].astype('bool')
	new_fits['rmHaExOtherCapped'] = fits_data['rmha_offset_other_capped'].astype('bool')
	# Isochronal Age feature
	new_fits['log10Age'] = fits_data['rmi_r_iso_log10age'].astype('float32')
	new_fits['OlderMaxIsoAge'] = fits_data['older_max_iso_age'].astype('bool')
	new_fits['YoungerMinIsoAge'] = fits_data['younger_min_iso_age'].astype('bool')
	new_fits['GoodClassIIYSOFit'] = fits_data['good_classiiyso_fit'].astype('bool')
	new_fits['IsoAgeCIILike'] = fits_data['isoage_classiiyso_like'].astype('float64')
	new_fits['IsoAgeOtherLike'] = fits_data['isoage_other_like'].astype('float64')
	new_fits['IsoAgeLikeRatio'] = fits_data['isoage_classiiyso_like'].astype('float64') / fits_data['isoage_other_like'].astype('float64')
	new_fits['IsoAgeCalc'] = fits_data['isoage_calc'].astype('bool')
	new_fits['IsoAgeCIICapped'] = fits_data['isoage_classiiyso_capped'].astype('bool')
	new_fits['IsoAgeOtherCapped'] = fits_data['isoage_other_capped'].astype('bool')

	return new_fits

# Set the DR2 ids to Null where they do not exist
def dr2_null():
	global replacement_fits

	print("Setting DR2 ids to zero where they do not exist")

	mask = (replacement_fits['dr2_source_id'] == -9223372036854775808)
	replacement_fits['dr2_source_id'][mask] = 0


# Create the new fits file and add the field comments
def output_full_fits(output_table, output_filename, source_file, updated_headers):
	print("Writing output fits file ...")
	# Delete file if it exists
	if os.path.isfile(output_filename):
		os.remove(output_filename)
	# Create fits header
	hdr = fits.Header()
	hdr['run_ver'] = str(run_version)
	hdr['datetime'] = str(datetime.datetime.now())
	hdr['data_f'] = source_file
	primary_hdu = fits.PrimaryHDU(header=hdr)
	# Convert table data into HDU
	hdu_data = fits.table_to_hdu(output_table)
	# Combine header with data
	hdul = fits.HDUList([primary_hdu, hdu_data])

	# Add comments to the fits fields
	num_fields = hdul[1].header['TFIELDS']
	all_fields_match = True
	for int_index in range(0, num_fields):
		#print(fits_file[1].header['TTYPE' + str(int_index+1)] + '=' + header_descriptions[int_index][0])
		if hdul[1].header['TTYPE' + str(int_index+1)] == updated_headers[int_index]['VizieR name']:
			#print(updated_headers[int_index][0] + ': ' + updated_headers[int_index][1])
			hdul[1].header.comments['TTYPE' + str(int_index+1)] = updated_headers[int_index]['FITS comment']
			hdul[1].header['TCOMM' + str(int_index+1)] = updated_headers[int_index]['FITS comment']
		else:
			print('Fail: ' + hdul[1].header['TTYPE' + str(int_index+1)] + '<>' + updated_headers[int_index]['VizieR name'])
			all_fields_match = False

	# Create the output file
	hdul.writeto(output_filename)



#------------------------------------------------------------
# Start of master routine
#------------------------------------------------------------

print('')
print('-----------------------------------------------------------------------------')
print('------------- Read/Add Description to Fields in a Fits File -----------------')
print('-----------------------------------------------------------------------------')
print('')

# Get the parameters from the user
region_name, process_type = user_inputs()

# Main path for the run
main_path = run_filepath + region_name + "/"
input_filename = main_path + "Bayes_full_" + region_name + ".fits"

if process_type == 'E':

	# Read the observation data file
	obs_header = read_fits_header(input_filename)

	# Create an empty array for the field names and descriptions
	num_fields = obs_header['TFIELDS']
	#num_fields = np.size(obs_data.colnames)

	#num_fields = np.size(obs_data,1)
	output_table = [["" for x in range(0,4)] for y in range(0,num_fields)]
	#print(num_fields)

	# Extract the field names and descriptions
	extract_fields(num_fields)
	#print(output_table)

	output_filename = main_path + "Bayes_" + region_name + "_Desc" + ".csv"
	file_header = "Name,FITS data type,FITS comment,Long description"
	export_csv(output_filename, output_table, file_header)

elif process_type == 'U':

	# Read the updated csv file of headers and descriptions
	header_filename = main_path + "Bayes_" + region_name + "_DescUpdated" + ".csv"
	updated_headers = read_field_csv(header_filename)

	# Read in the fits file
	fits_data = read_fits(input_filename)

	# Update the columns
	replacement_fits = update_fits_columns(input_filename, updated_headers)

	# Update DR2 ids to null where the nulls have been replaced by default value
	#dr2_null()

	# Create the new fits file and add the field comments
	output_filename = main_path + "Bayes_full_" + region_name + "_pub.fits"
	output_full_fits(replacement_fits, output_filename, input_filename, updated_headers)

	# Create sample file with first ten rows
	fits_ten = replacement_fits[:10]
	output_filename = main_path + "Bayes_full_" + region_name + "_pub_first10.fits"
	output_full_fits(fits_ten, output_filename, input_filename, updated_headers)

	# Create sample file with first ten rows of the P(CII)>0.5 set
	mask = fits_data['posterior_classiiyso']>0.5
	fits_cii = replacement_fits[mask]
	fits_cii_ten = fits_cii[:10]
	output_filename = main_path + "Bayes_full_" + region_name + "_pub_PCII05_first10.fits"
	output_full_fits(fits_cii_ten, output_filename, input_filename, updated_headers)

elif process_type == 'V':

	# Read the updated csv file of headers and descriptions
	header_filename = main_path + "Bayes_" + region_name + "_DescUpdated" + ".csv"
	updated_headers = read_field_csv(header_filename)

	# Read in the published fits file
	input_filename = main_path + "Bayes_full_" + region_name + "_pub.fits"
	fits_data = read_fits(input_filename)

	# Export the VizieR data file
	output_filename = main_path + "catalog.dat"
	export_vizier_data(output_filename, fits_data, updated_headers)

elif process_type == 'R':

	# Read the updated csv file of headers and descriptions
	header_filename = main_path + "Bayes_" + region_name + "_DescUpdated" + ".csv"
	updated_headers = read_field_csv(header_filename)

	# Export the VizieR data file
	output_filename = main_path + "Bayes_full_" + region_name + "_readme"
	export_vizier_readme(output_filename, updated_headers)

else:

	print('Nothing done! Please check your selections.')

#------------------------------------------------------------
# End of master routine
#------------------------------------------------------------