# ------------------------------------------------------------------
# For calculating variability in G mag using Gaia DR2
# Author: Andrew Wilson created this by modifying original code written by Prof Tim Naylor
# Institution: University of Exeter, Astrophysics Group
# Version date: 15 September 2020
# University email: aw648@exeter.ac.uk
# Personal email: andyjwilson_uk@hotmail.com
# ------------------------------------------------------------------
# This code calculates the run of first quartile values of the PDF as a
# function of the median of the PDF.
# ------------------------------------------------------------------


import numpy as np
import dP_dr_fast

def locate_median(sigma_i, eFold, nObs, frac, approx, write_file=False, filepath='', mag=''):

  # The range of the output pdf is harder.  If we want it to integrate to one, lets try and get about
  # 99.9% of the variables.  There are two parts to consider.  Let's start with the chi-squared disribution
  # in the case that there are no intrinsic variables.  Approximate it as a Gaussian and then go out to
  # 3.3 sigma, which gets 99.9 percent.
  gauss_sig = sigma_i/np.sqrt(np.float32(2*(nObs-1)))
  upper = sigma_i + 3.3*gauss_sig
  # But maybe the tail is dominated by the intrinsic variables, so also go out to 7 times b
  # (gets 99.9% of the integral).
  upper_best = max(upper, 7.0*eFold)

  # And likewise we can do a lower limit based on the chi-squared distrubution.
  lower = max (0.0, sigma_i-3.3*gauss_sig)

  # The width of the bins in the output pdf, set to 5% of the width of the chi-squared distribution.
  binwidth=0.05*gauss_sig

  # So now we can calculate the number of points.
  npts = 1 + np.int((upper_best-lower)/binwidth)

  # Define the output array, which is the PDF of the root-mean-square values.
  dP_dr_pdf=np.zeros(shape=(npts), dtype=float)
  # And its x-axis.
  ipts = np.arange(0, npts, dtype=float)
  r_val = lower + binwidth*ipts

  # Now calculate the PDF at each value of the RMS.
  for ipts in range (0, npts):
    dP_dr_pdf[ipts]=dP_dr_fast.variability_pdf(r_val[ipts], sigma_i, eFold, frac, nObs, approx)

  if write_file:
    filename = 'dP_dr_mag' + mag + '_sig' + str(round(sigma_i,4)) + '_sca' + str(round(eFold,3)) + '_n' + str(nObs).zfill(3) + '_frac' + str(round(frac,3)) + '.dat'
    file = open(filepath+filename,'w')
    file.write('# The PDF of the RMS for a model which has an instrumental sigma of %g \n' % sigma_i)
    file.write(
        '# %i Gaia observations, an e-folding width of %g and a zero variability fraction of %g \n'
               % (nObs, eFold, frac))
    file.write('# r dP_dr \n')

    for ipts in range (0, npts):
      file.write('%14.8f' % r_val[ipts])
      file.write('%14.8f' % dP_dr_pdf[ipts])
      file.write('\n')

    file.close()

  # Now find the quartiles.
  median_pdf = dP_dr_fast.pdf_median(dP_dr_pdf, r_val, binwidth)

  return median_pdf
