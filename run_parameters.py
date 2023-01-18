# ------------------------------------------------------------
# File of run paramters for the Naive Bayes Classifier
# Author: Andrew Wilson
# Institution: University of Exeter, Astrophysics Group
# Version date: 18 Jan 2023 (file paths and region options updated 4 Nov 2022)
# University email: aw648@exeter.ac.uk
# Personal email: andyjwilson_uk@hotmail.com
# ------------------------------------------------------------

run_version = "8.28.00"

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# File Paths
# ----------

# Master run path where results files are saved
run_filepath = "/Users/Shared/BayesCIIYSO/"
# Input data file paths
isoage_iso_filepath = "/Users/Shared/BayesCIIYSO/Isochrones/"
stilism_filepath = "/Users/Shared/BayesCIIYSO/STILISM/"

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# Region lists
# ------------
region_descriptions=['1 - Northern Galactic Plane Gaia EDR3: 20<l<220 ; -4<b<4', 
	'2 - 4 sq deg test region Gaia EDR3: 60<l<62 ; -1<b<1']
# This is how the regions are called in all file names
region_list = {
	"1": "NGPv3",
	"2": "TestEDR3"
	}

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# Source File Lists
# -----------------

# IGAPS data file for H-alpha excess and Isochronal Age
igaps_file_list = {
	"1": "IGAPS_GoodPhot_riHaAB_NGPv3_v828.fits",
	"2": "IGAPS_GoodPhot_riHaAB_TestEDR3_v828.fits"
	}

# Gaia G Variability source files
gaia_g_var_file_list = {
	"1": "GaiaEDR3_NGPv3_RP18_RUWE14_ParCor_PoverE3_P05_r18_Classifier.fits",
	"2": "GaiaEDR3_TestEDR3_RP18_RUWE14_ParCor_PoverE3_P05_r18_Classifier.fits"
	}

# (H-K) Excess source files
hk_file_list = {
	"1": "2MASS_UKIDSS_Comb_NGPv3.fits",
	"2": "2MASS_UKIDSS_Comb_TestEDR3.fits"
	}

# (H-K) catalogues
hk_cat_list = {
	"1": "2MASS_UKIDSS",
	"2": "2MASS_UKIDSS"
	}

# WISE source files
wise_list = {
	"1": "ALLWISE_GaiaEDR3_NGPv3_NoPhot_MatchProbGt50.fits",
	"2": "ALLWISE_GaiaEDR3_TestEDR3_NoPhot_MatchProbGt50.fits"
	}

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# Bayes parameters
# ----------------
# Set the priors
prior_classiiyso = 0.001
prior_other = 0.999

# Create an python list of the catgories
bayes_cats = ['_classiiyso', '_other']

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# Reddening photometry conversions
# --------------------------------

# Set the conversion factors from E(B-V) to other colours
#E_BV_to_E_ri = 0.71
#E_BV_to_E_rHa = 0.22
# E(B-V) to A_r-A_I from Tim 10/10/2020
E_BV_to_ArI_p0 = 0.0047801677
E_BV_to_ArI_p1 = 0.73329043
E_BV_to_ArI_p2 = -0.0096184695
# E(B-V) to A_r-A_Ha for hot stars at 10,000K from Tim 10/10/2020
E_BV_to_ArHa_p0 = 0.0064770393
E_BV_to_ArHa_p1 = 0.24482922
E_BV_to_ArHa_p2 = -0.022416774

# Convert E(B-V) to A(r)
E_BmV_to_A_r = 2.38

# Set the conversion factors from E(B-V) to other colours
#J-H: 0.2963 +/- 0.0002
#H-K: 0.1668 +/- 0.0004
T2MASS_E_BV_to_E_JH = 0.2963
T2MASS_E_BV_to_E_HK = 0.1668
#J-H: 0.2747 +/- 0.0006
#H-K: 0.1815 +/- 0.0005
UKIDSS_E_BV_to_E_JH = 0.2747
UKIDSS_E_BV_to_E_HK = 0.1815

# Fit to linear region of 1 Myr Baraffe 2015 BT-Settl isochrone at E(B-V)=0.0
# for masses 0.13Msol<mass<0.35Msol
T2MASS_iso_slope = -0.3733
T2MASS_iso_intercept = 0.5188
UKIDSS_iso_slope = -0.4835
UKIDSS_iso_intercept = 0.6003

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------