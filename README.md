# Naive Bayes Classifier for indentifying Class II YSOs
## Full details in Wilson et al (2023) in MNRAS
## Author: Andrew Wilson
## Institution: University of Exeter, Astrophysics Group
## Date: 18 January 2022

This Read Me file is a brief introduction to the Python scripts used to generate Naive Bayes Classifier results for Class II YSOs.
It has been written for version 8.28 of the code and will only be updated for major code releases.

The Python scripts have been written in Python 3.7 on a Mac. They may work in other Python 3.x versions and on Windows. It is intended to be run from the command line and all the Python scripts should be saved in the same folder.

The code can be run by updating the paths in the run_parameters.py file and obtaining the input files. The run parameters file also contains common variables such as the Bayes priors, categories and the list of available regions on which the code can run. Add to or amend the regions by extending region_descriptions and igaps_file_list variables, and the file list variables.

The run_version is important as this is saved in the fits file headers and is used to identify the run. It is suggested the version number is incremented as appropriate with a minor/major update to code, data or likelihoods.

The Gaia source_id field is used as the index to data in all files, and to combine results at the end. All files should contain a full list of the source ids, even if not all source ids contain data for the feature. For example, the WISE input file should contain blank entries in the WISE columns where no data exists for the source. This approach allows the code to run quickly, without needing to combine data by searching the rows of different feature files.

Data Files
----------

The following data files are needed to run the code. Provided as links to a zip.

Isochrones - A zip of the files and folders of the isochrones (see Wilson et al 2023 for details) in the correct structure to be used by the classifier. Unzipped there are 3.2 GB of files. There are 120 folders at different reddenings, each containing a set of 21 isochrones at different ages.
https://1drv.ms/u/s!AoS2nZUeme4DmIVIGSGAx8Q_eVYxIA?e=RNynTW

Stilism - A zip of the Stilism reddening files (Lallement et al 2014; Capitanio et al 2017; Lallement et al 2018) used for the Wilson et al (2021) catalogue covering the footprint in Galactic coordinates 20<l<220 and |b|<4. Unzipped there are 17,282 files totalling 214 MB.
https://1drv.ms/u/s!AoS2nZUeme4DmIVHbPDzhM4gz4hvxw?e=qhWzB9

TestEDR3 - A folder of test input files in the correct format to run the classifier. These also demonstrate the required columns and names to run the classifier. Unzipped there are 12 files totalling 8.9 MB.
https://1drv.ms/u/s!AoS2nZUeme4DmIVGTRoNnlhqMJBhyg?e=it3LBv

Run Order
---------

The individual feature code can be run independently, followed by the likelihoods and bayes scripts. Below is the suggested run order for a full end-to-end run.

**The version of the data, code and parameters is 8.28. This will need to be entered into some of the Python scripts when prompted.**

1. Obtain the input files and save them in your chosen run folder, except for the Stilism files and isochrones which should be saved in dedicated folders. These folders are specified in the run_parameters.py file as described in the next step.  
   - Feature data:
     - Gaia
     - IGAPS
     - WISE
     - Combined 2MASS and UKIDSS
   - Other data:
     - Coordinates - Any columns appearing in this file will be added to the final results files
     - Likelihood csv files
     - Gaia G Variability parameters csv file (includes likelihoods)
     - Stilism files for 0.5 degree to a side (link above to download and save in dedicated folder)
     - Isochrones (link above to download and save in dedicated folder)

2. Update the parameters file with the run version and filepaths.  
You will need to create a folder for the run results. This starts with the run_filepath variable followed by a subfolder name that matches the region_list Python list.
   - **run_parameters.py**

3. Run the reddening.
   - stilism_red.py

4. Run the feature scripts.
   - gaia_g_var.py
   - ha_excess.py
   - hk_excess.py
   - isoage.py
   - wise_w1w2.py

5. Run the likelihood script once for each feature/likelihood change.  
It is advisable to only have the current likelihood files in the run folder to avoid any risk of accidentally generating results using old likelihoods.  
Note this script does not generate the Gaia G-band variability likelihoods as they are incorporated in the gaia_g_var script run in the previous step.
   - lhoods.py

6. Run the Bayes script.
   - bayes.py

7. Optional step that updates the file layout to match that published in MNRAS, CDS and other locations.  
First run option U to create the published fits file.  
Second run option V to create the VizieR format text file.
   - bayes_pub.py (Converts the columns to those in the published catalogue.)


Other Python Scripts
--------------------

**gaia_data_preparation/gaia_edr3_gcor.py**  
Corrects the Gaia EDR3 G-band fluxes and magnitudes. Adapted copy of the Python code from the Gaia EDR3 website.

**gaia_data_preparation/gaia_edr3_parcor.py**  
Updates the Gaia EDR3 parallaxes for the biases found by Lindegren et al (2021). Adapted copy of the Python code from the Gaia EDR3 website.

**dP_dr_fast.py**  
Script called automatically in the calculation of Gaia G variability.

**var_model_run.py**  
Not directly used but provided as useful for exploring the Gaia G variability model.

Other Files
-----------

**stilism_example/Stilism_Code_Example.txt**  
An example of the command line code to download the Stilism map in Galactic coordinates using the Stilism_Coords_Example.csv file.

**stilism_example/Stilism_Coords_Example**  
Example coordinates file for use with Stilism_Script_Example.txt.
