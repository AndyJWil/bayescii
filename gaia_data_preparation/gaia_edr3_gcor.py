# ------------------------------------------------------------
# For applying corrections to Gaia EDR3 G magnitude and fluxes
# Python script associated with Riello et al 2020
# Gaia Early Data Release 3: Photometric content and validation
# https://github.com/agabrown/gaiaedr3-6p-gband-correction
# Author: Andrew Wilson
# Institution: University of Exeter, Astrophysics Group
# Version date: 7 December 2020
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
	table_data['phot_g_mean_mag_corrected'] = np.nan
	table_data['phot_g_mean_flux_corrected'] = np.nan
	gmag_corrected, gflux_corrected = correct_gband(table_data['bp_rp'], table_data['astrometric_params_solved'], table_data['phot_g_mean_mag'], table_data['phot_g_mean_flux'])
	table_data['phot_g_mean_mag_corrected'] = gmag_corrected
	table_data['phot_g_mean_flux_corrected'] = gflux_corrected
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


# Function from Riello et al 2020 Jupter notebook
# Taken from github 07 Dec 2020
# https://github.com/agabrown/gaiaedr3-6p-gband-correction
def correct_gband(bp_rp, astrometric_params_solved, phot_g_mean_mag, phot_g_mean_flux):
    """
    Correct the G-band fluxes and magnitudes for the input list of Gaia EDR3 data.
    
    Parameters
    ----------
    
    bp_rp: float, numpy.ndarray
        The (BP-RP) colour listed in the Gaia EDR3 archive.
    astrometric_params_solved: int, numpy.ndarray
        The astrometric solution type listed in the Gaia EDR3 archive.
    phot_g_mean_mag: float, numpy.ndarray
        The G-band magnitude as listed in the Gaia EDR3 archive.
    phot_g_mean_flux: float, numpy.ndarray
        The G-band flux as listed in the Gaia EDR3 archive.
        
    Returns
    -------
    
    The corrected G-band magnitudes and fluxes. The corrections are only applied to
    sources with a 6-parameter astrometric solution fainter than G=13, for which a
    (BP-RP) colour is available.
    
    Example
    -------
    
    gmag_corr, gflux_corr = correct_gband(bp_rp, astrometric_params_solved, phot_g_mean_mag, phot_g_mean_flux)
    """

    if np.isscalar(bp_rp) or np.isscalar(astrometric_params_solved) or np.isscalar(phot_g_mean_mag) \
                    or np.isscalar(phot_g_mean_flux):
        bp_rp = np.float64(bp_rp)
        astrometric_params_solved = np.int64(astrometric_params_solved)
        phot_g_mean_mag = np.float64(phot_g_mean_mag)
        phot_g_mean_flux = np.float64(phot_g_mean_flux)
    
    if not (bp_rp.shape == astrometric_params_solved.shape == phot_g_mean_mag.shape == phot_g_mean_flux.shape):
        raise ValueError('Function parameters must be of the same shape!')
    
    do_not_correct = np.isnan(bp_rp) | (phot_g_mean_mag<=13) | (astrometric_params_solved != 95)
    bright_correct = np.logical_not(do_not_correct) & (phot_g_mean_mag>13) & (phot_g_mean_mag<=16)
    faint_correct = np.logical_not(do_not_correct) & (phot_g_mean_mag>16)
    bp_rp_c = np.clip(bp_rp, 0.25, 3.0)
    
    correction_factor = np.ones_like(phot_g_mean_mag)
    correction_factor[faint_correct] = 1.00525 - 0.02323*bp_rp_c[faint_correct] + \
        0.01740*np.power(bp_rp_c[faint_correct],2) - 0.00253*np.power(bp_rp_c[faint_correct],3)
    correction_factor[bright_correct] = 1.00876 - 0.02540*bp_rp_c[bright_correct] + \
        0.01747*np.power(bp_rp_c[bright_correct],2) - 0.00277*np.power(bp_rp_c[bright_correct],3)
    
    gmag_corrected = phot_g_mean_mag - 2.5*np.log10(correction_factor)
    gflux_corrected = phot_g_mean_flux * correction_factor
    
    return gmag_corrected, gflux_corrected


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

# Read the observation data file
obs_file = data_filename + ".fits"
obs_data = read_observations(data_path + "/" + obs_file)

# Create output fits file with all data
output_filename = data_path + "/" + data_filename + "_MagCor.fits"
output_full_fits(output_filename, obs_data, obs_file)

# Remove interim processing columns not needed by Naive Bayes Classifier
#remove_interim_columns()

#------------------------------------------------------------
# End of master routine
#------------------------------------------------------------
