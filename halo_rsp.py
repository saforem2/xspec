# Python script to load a synphot FITS spectrum and display it

# import needed extensions
from numpy import *
# old version should be
# import pyfits
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt # plotting package

# read in the transmission data
ener, c2, c1 = transpose(genfromtxt('C1C2eff.dat'))
# remove points below ethreshold
ethreshold = 0.1  # keV
q = where(ener >= ethreshold)
ener = ener[q]
c1 = c1[q]
c2 = c2[q]
# replace with evenly spaced energies
t = max(ener) + ener[1] - ener[0]
ener = arange(min(ener), t, (t-min(ener))/len(ener))
# should add silicon dead layer
area = 0.25 # detector area in cm^2

plt.ion() # do plots in interactive mode
# plot the spectrum
fig1 = plt.figure(1)
plt.clf()
plt.xlabel('Energy (keV)')
plt.ylabel('Response (cm^2)')
# plt.xscale('log')
plt.plot(ener, area*c1, '-b')
plt.plot(ener, area*c2, '-r')
plt.show() # display the plot
plt.savefig('effective_area_energy.pdf', format='pdf')

# translate into energy bin lower and upper bounds and average response
n = len(ener)
emin = float32(ener[0:n-1])

