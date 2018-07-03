#####################################################################
# This sample.R routine shows how to calculate synthetic SSP spectra
#####################################################################

#---------------------------------
# Define functions
#---------------------------------
# For a description of the functions 
#  and parameters, see manual.pdf
source('stpars.R')
source('pfant12.R')
source('SSP_model.R')

#---------------------------------
# Set parameters
#---------------------------------
n_ms = 4          # number of desired main sequence stars
n_rg = 3          # number of desired red giant stars
feh = -0.5         # [Fe/H]
afe = 0.33         # [alpha/Fe]
age = 12.5        # age (Gyr)
lmin = 4600      # lower lambda
lmax = 5600       # upper lambda
imf = 'Kroupa'    # IMF


# Select pairs of Teff, logg from Dartmouth isochrones
  stpars(n_ms, n_rg, feh, afe, age, fig = F)

# Calculate the synthetic stellar spectra for pairs of Teff, logg
#   (this step takes a few minutes)
  pfant12(feh, afe, lmin, lmax, age, n_ms = n_ms, n_rg = n_rg)

# Add the synthetic stellar spectra according to an IMF
  ssp.model(feh, afe, age, imf, fwhm = 2.5, n_ms = n_ms, n_rg = n_rg, lmin = lmin, lmax = lmax)

# Plot the spectra
  file_ssp <- set.ssp.filename(feh, afe, age = age, imf = imf)
  t = read.table(file_ssp)
  xx1 = median(t$V2[t$V1 > lmin & t$V1 < lmax])
  plot(t$V1, t$V2/xx1, col = 'blue', type = 'l', lty = 1, lwd = 2.5,
       xlab = expression(paste(lambda, ' (', ring(A), ')')), 
       ylab = 'Flux', ylim = c(0, 1.5))

  # t = read.table('Miles_Ku1.30Zm0.35T06.dat')
  # xx2 = median(t$V2[t$V1 > lmin & t$V1 < lmax])
  # lines(t$V1, t$V2/xx2, type = 'l', col ='red', lwd = 0.8)

#------------------------------------------------------------------------------------------------
# Changing abundances ratios
#  The abundance ratios [O, Mg, Si, Ca, Ti/Fe] are set to [a/Fe], unless you specify other value
#  Example [Ca/Fe] = 0.4: the two lines below will calculate synthetic spectra with 
#    [O, Si, Mg, Ti/Fe] = [a/Fe] = 0.10 and
#    [Ca/Fe] = 0.4
#------------------------------------------------------------------------------------------------
  # pfant12(feh, afe, lmin, lmax, age, CaFe = 0.4, n_ms = n_ms, n_rg = n_rg)
  # ssp.model(feh, afe, age, imf, fwhm = 2.5, CaFe = 0.4, n_ms = n_ms, n_rg = n_rg, lmin = lmin, lmax = lmax)

# Overplot the new spectrum
  # file_ssp <- set.ssp.filename(feh, afe, age = age, CaFe = 0.4, imf = imf)
  # t = read.table(file_ssp)
  # xx3 = median(t$V2[t$V1 > lmin & t$V1 < lmax])
  # lines(t$V1, t$V2/xx3, col = 'black', type = 'l', lty = 2, lwd = 2.5)

#------------------------------------------------------------------------------------------------
# Changing the IMF
# The options are Salpeter, Kroupa and Unimodal 
#------------------------------------------------------------------------------------------------
# imf   = 'Unimodal'
# slope = 3.3
#   ssp.model(feh, afe, age, imf, fwhm = 2.5, slope = slope, CaFe = 0.4, n_ms = n_ms,
#           n_rg = n_rg, lmin = lmin, lmax = lmax)

# Overplot the new spectrum
  # file_ssp <- set.ssp.filename(feh, afe, age = age, CaFe = 0.4, imf = imf, slope = slope)
  # t = read.table(file_ssp)
  # xx4 = median(t$V2[t$V1 > lmin & t$V1 < lmax])
  # lines(t$V1, t$V2/xx4, col = 'green', type = 'l', lty = 2, lwd = 1.5)

# legend(4600, 0.68, c('[Ca/Fe] = [a/Fe] = 0.1, Kroupa',
#                     '[Ca/Fe] = 0.4, Kroupa', '[Ca/Fe] = 0.4, Uni slope = 3.3',
#                     'MIUSCAT'),
#        lty = c(1, 2, 2, 1), col = c('blue', 'black', 'green', 'red'),
#        lwd = c(2, 2, 1, 1), bty = 'n')
