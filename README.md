# bayescii
Naive Bayes Classifier for indentifying Class II YSOs

# ------------------------------------------------------------
# Naive Bayes Classifier for indentifying Class II YSOs
# Full details in Wilson et al (2023) in MNRAS
# Author: Andrew Wilson
# Institution: University of Exeter, Astrophysics Group
# Date: 18 January 2022
# University email: aw648@exeter.ac.uk
# Personal email: andyjwilson_uk@hotmail.com
# ------------------------------------------------------------

This Read Me file is a brief introduction to the Python scripts used to generate Naive Bayes Classifier results for Class II YSOs.
It has been written for version 8.28 of the code and will only be updated for major code releases.

The Python scripts have been written in Python 3.7 on a Mac. They may work in other Python 3.x versions and on Windows. It is intended to be run from the command line and all the Python scripts should be saved in the same folder.

The code can be run by updating the paths in the run_parameters.py file and obtaining the input files. The run paramters file also contains common variables such as the Bayes priors, categories and the list of available regions on which the code can run. Add to or amend the regions by extending region_descriptions and igaps_file_list variables, and the file list variables.

The run_version is important as this is saved in the fits file headers and is used to identify the run. It is suggested the version numebr is incremented as appropriate with a minor/major update to code, data or likelihoods.

The Gaia source_id field is used as the index to data in all files, and to combine results at the end. All files should contain a full list of the source ids, even if not all source ids contain data for the feature. For example the WISE input file should contain blank entries in the WISE columns where no data exists for the source. This approach allows the code to run quickly, without needing to combine data by searching the rows of different feature files.

Run Order
---------

The individual feature code can be run independently, followed by the likelihoods and bayes scripts. Below is the suggested run order for a full end-to-end run.

* The version of the data, code and parameters is 8.28. This will need to be entered into some of the Python scripts when prompted.

1 - Obtain the input files
	Feature data:
		Gaia - Save in dedicated folder
		IGAPS or VPHAS - Save in dedicated folder
		WISE - Save in dedicated folder
		2MASS / UKIDSS / VVV - Manually as appropriate - Save in run path
	Other data
		Coordinates - Any columns appearing in this file will be added to the final results files - Save in the run path
		Likelihood csv files - Save in the run path
		Gaia G Variability parameters csv file (includes likelihoods) - Save in the run path
		Stilism files for 0.5 degree to a side - https://stilism.obspm.fr - Save in a dedicated folder
		Isochrones from the Cluster Isochrone server or the ark iso FORTRAN code (or from Tim Naylor) (https://www.astro.ex.ac.uk/people/timn/isochrones/) - Save in a dedicated folder

2 - Update the parameters file with the run version and filepaths.
		run_parameters.py
		You will need to create folder(s) for the run results. These start with the run_filepath variable followed by subfolders that match the region_list Python list.

3 - Run the reddening.
		stilism_red.py

4 - Run the feature scripts.
		gaia_g_var.py
		ha_excess.py
		hk_excess.py
		isoage.py
		wise_w1w2.py

5 - Run the likelihood script once for each feature/likelihood change.
	It is advisable to only have the current likelihood files in the run folder to avoid any risk of accidentally generating results using old likelihoods.
	Note this script does not generate the Gaia G-band variability likelihoods as they are incorported in the gaia_g_var script run in the previous step.
		lhoods.py

6 - Run the Bayes script.
		bayes.py

7 - Optional step that updates the file layout to match that published in MNRAS, CDS and other locations.
	First run option U to create the publilshed fits file.
	Second run option V to create the VizieR format text file.
		bayes_pub.py - Converts the columns to those in the published catalogue.


Other Python Scripts
--------------------

gaia_edr3_gcor.py
Corrects the Gaia EDR3 G-band fluxes and magnitudes. Adapated copy of the Python code from the Gaia EDR3 website.

gaia_edr3_parcor.py
Updates the Gaia EDR3 parallaxes for the biases found by Lindegren et al (2021). Adapated copy of the Python code from the Gaia EDR3 website.

dP_dr_fast.py
Script called automatically in the calculation of Gaia G variability.

var_model_run.py
Not directly used but provided as useful for exploring the Gaia G variability model.