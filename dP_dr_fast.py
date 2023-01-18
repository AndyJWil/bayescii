# ------------------------------------------------------------------
# For calculating variability in G mag using Gaia DR2
# Author: Prof Tim Naylor with some updates by Andrew Wilson
# Institution: University of Exeter, Astrophysics Group
# Version date: 2 March 2021
# University email: aw648@exeter.ac.uk
# Personal email: andyjwilson_uk@hotmail.com
# ------------------------------------------------------------------
# This is the library of routines which calculate the predicted PDFs
# for the distribution of RMS values from Gaia.
# ------------------------------------------------------------------

################################################################################

import numpy as np
from scipy.special import gamma
from scipy.stats import chi2

################################################################################

def variability_pdf(rms, sigma_i, expScale, frac, n, approx):

  # Calculate the value of df_dr at rms given an instrumental noise sigma_i
  # a number of observations n 
  # and parameters for the model of the astrophysical variability b and frac.
  # Use the approximation of the chi-squared function for large n if approx is true.

  # First question.  Are we so far out in the chi-squared tail that the instrumental
  # variability does not matter?  Approximate it as a Gaussion and go out 50 times
  # the width.   

  long_approx=True
  
  if rms > sigma_i+(10.0*sigma_i/np.sqrt(np.float32(2*(n-1)))) and long_approx:
      
    integ_norm = (1.0-frac)*intrinsic(expScale, np.sqrt(rms*rms-sigma_i*sigma_i))*rms/np.sqrt(rms*rms-sigma_i*sigma_i)

  else:

    # Define the step in the intrinsic variability.  Here assume expScale defines a
    # scale length, and let's make 20 samples of it.  This seems to give an
    # answer to a percent or so for an exponential.
    n_samples = 20
    binScale = expScale/np.float(n_samples)
    # This seems to cause ringing, If the scale of the chi-squared function is much tighter.
    binWidth = min(binScale, 3.0*sigma_i/np.sqrt(np.float32(2*(n-1))))

    # Go out to four times the scale length of the exponential, this gets
    # 99.5% of the integral
    #n_sigma_a=np.int(4*n_samples) - Error corrected 19 Feb 2021. Appears not to make any difference as worst case it integrates over more points
    n_sigma_a = np.int(4.0*expScale/binWidth)
        
    j_sigma_a = np.arange(0, n_sigma_a, dtype=float)

    # So we pick a group of stars with a given variability.
    # Note the binWidth cancels so the maximum is simply
    #  sigma_a = 4.0*exp_Scale
    sigma_a = j_sigma_a*binWidth
    # We are integrating, so mulitply by the width.
    intr = intrinsic(expScale,sigma_a)*binWidth

    # Now consider what fraction of those stars will have a variability of
    # sigma_a.  Normally this would be given by the chi-squared function, but
    # we have not divided the x-axis by the expected sigma^2, but by the
    # instrumental variability squared.  So we have to move the peak of the
    # chi-squared distribution right by a factor of sigma_i^2+sigma_a^2.
    # Hence...  

    sigma2_astro = sigma_i*sigma_i + sigma_a*sigma_a
    if approx:
      conv_astro_temp = gauss(rms*rms*np.float32(n)/sigma2_astro, np.float32(n-1),
                      np.sqrt(2.0*np.float32(n-1)))
    else:
      conv_astro_temp = chi2.pdf(rms*rms*np.float32(n)/sigma2_astro, np.float32(n-1))
    # Apply part of the normalisation of the chi-squared/Gaussian.
    conv_astro_norm = conv_astro_temp/sigma2_astro

    # And finally multiply this by the fraction of stars at this variability,
    # and add it to the running total for this value of the RMS.
    integ_var = np.sum(conv_astro_norm*intr)


    # Now add the effect of a delta function.
    sigma2_delta = sigma_i*sigma_i
    if approx:
      integ_delta_temp = gauss(rms*rms*np.float32(n)/sigma2_delta,
                                         np.float32(n-1), np.sqrt(2.0*np.float32(n-1)))
    else:
      integ_delta_temp = chi2.pdf(rms*rms*np.float32(n)/sigma2_delta, np.float32(n-1))
    # Apply part of the normalisation of the chi-squared/Gaussian.
    integ_delta = integ_delta_temp/sigma2_delta
    
    # Combine the integrals
    integ_comb = (1.0-frac)*integ_var + frac*integ_delta

    # Then apply the rest of the normalisation of the chi-squared/Gaussian, and
    # apply the factor to change from mean-square to RMS.
    integ_norm = integ_comb*2.0*rms*np.float32(n)

  return integ_norm


################################################################################

# The distribution of intrinsic variabilities. Must be normalised, and must not be
# called for f < 0.

def intrinsic(b, f):
    return (np.exp(-f/b))/b

################################################################################

# A Gaussian distribution.

def gauss(x, mu=0.0, sigma=1.0):
    # Original function for x replaced as float not working with x as an array - Andy Wilson 16 Dec 2020
    #x = float(x - mu) / sigma
    exp_power = (x - mu) / sigma
    return np.exp(-exp_power*exp_power/2.0) / (sigma * np.sqrt(2.0*np.pi) )

################################################################################

# Get the first two quartiles.

def quartiles(df_dr_pdf, r, binwidth):

  point25=-1.0
  point5=-1.0

  sum=0.0

  for ipts in range (0, len(df_dr_pdf)):

    sum+=df_dr_pdf[ipts]*binwidth
    if sum > 0.1 and point25 < 0.0:
      point25=r[ipts]
    if sum > 0.5 and point5 < 0.0:
      point5=r[ipts]

  print('Integral (which may not cover the whole rage) is ', sum)

  return point25, point5

################################################################################

# Get the median

def pdf_median(df_dr_pdf, r, binwidth):

  point5=-1.0

  sum=0.0

  for ipts in range (0, len(df_dr_pdf)):

    sum+=df_dr_pdf[ipts]*binwidth
    if sum > 0.5 and point5 < 0.0:
      point5=r[ipts]

  return point5

################################################################################
