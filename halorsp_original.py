# Python program to load a synphot FITS spectrum and display it

# import needed extensions
from numpy import *
import pyfits
import matplotlib.pyplot as plt # plotting package

# read in the transmission data
ener, c2, c1 = transpose(genfromtxt('C1C2eff.dat'))
# remove points below ethreshold
ethreshold = 0.2 # keV
q = where(ener >= 0.2)
ener = ener[q]
c1 = c1[q]
c2 = c2[q]
# replace with evenly spaced energies
t = max(ener)+ener[1]-ener[0]
ener = arange(min(ener), t, (t-min(ener))/len(ener))
# should add silicon dead layer
area = 0.25 # detector area in cm^2

plt.ion() # do plots in interactive mode
# plot the spectrum
plt.figure(1)
plt.clf()
plt.xlabel('Energy (keV)')
plt.ylabel('Response (cm^2)')
# plt.xscale('log')
plt.plot(ener, area*c1 , '-b')
plt.plot(ener, area*c2 , '-r')
plt.show() # display the plot

# translate into energy bin lower and upper bounds and average response
n = len(ener)
emin = float32(ener[0:n-1])
emax = float32(ener[1:n])
resp = float32(area*0.5*(c1[0:n-1]+c1[1:n]))
plt.plot(0.5*(emin+emax), resp, '-k')


# make .arf file
# create FITS table with ARF data
col1 = pyfits.Column(name='E_MIN', format='1E', unit='keV', array=emin)
col2 = pyfits.Column(name='E_MAX', format='1E', unit='keV', array=emax)
col3 = pyfits.Column(name='SPECRESP', format='1E', unit='cm**2', array=resp)
# create a ColDefs (column-definitions) object for all columns:
cols = pyfits.ColDefs([col1, col2, col3])
# create a new binary table HDU object
arfhdu = pyfits.new_table(cols)
arfhdu.header.update('EXTNAME', 'SPECRESP', 'Extension type')
arfhdu.header.update('TELESCOP', 'HaloSat', 'Mission name')
arfhdu.header.update('INSTRUME', 'SDD1', 'Instrument')
#arfhdu.header.update('CHANTYPE', 'PHA', 'Channels from detector electronics')
#arfhdu.header.update('DETCHANS', len(emin), 'total number of raw detector PHA channels')
arfhdu.header.update('HDUCLASS', 'OGIP', 'file format is OGIP standard')
arfhdu.header.update('HDUCLAS1', 'RESPONSE', 'extension contains response data')
arfhdu.header.update('HDUCLAS2', 'SPECRESP', 'extension contains a response matrix')
arfhdu.header.update('HDUVERS', '1.1.0', 'version of the file format')
# create primary FITS header
prihdr = pyfits.Header()
prihdu = pyfits.PrimaryHDU(header=prihdr)
# create HDUList with primary HDU and the data table extension
thdulist = pyfits.HDUList([prihdu, arfhdu])
# write to file
thdulist.writeto('halo_arf.fits', clobber=True)

# make .rmf file
# create EBOUNDS extension
chan = arange(0, len(emin), dtype='int16')
col1 = pyfits.Column(name='CHANNEL', format='I', array=chan)
col2 = pyfits.Column(name='E_MIN', format='1E', array=emin)
col3 = pyfits.Column(name='E_MAX', format='1E', array=emax)
# create a ColDefs (column-definitions) object for all columns:
cols = pyfits.ColDefs([col1, col2, col3])
# create a new binary table HDU object
eboundshdu = pyfits.new_table(cols)
eboundshdu.name = 'EBOUNDS'
eboundshdu.header.update('EXTNAME', 'EBOUNDS', 'Extension type')
eboundshdu.header.update('TELESCOP', 'HaloSat', 'Mission name')
eboundshdu.header.update('INSTRUME', 'SDD1', 'Instrument')
eboundshdu.header.update('CHANTYPE', 'PHA', 'Channels from detector electronics')
eboundshdu.header.update('DETCHANS', len(emin), 'total number of raw detector PHA channels')
eboundshdu.header.update('HDUCLASS', 'OGIP', 'file format is OGIP standard')
eboundshdu.header.update('HDUCLAS1', 'RESPONSE', 'extension contains response data')
eboundshdu.header.update('HDUCLAS2', 'EBOUNDS', 'extension contains a response matrix')
eboundshdu.header.update('HDUVERS', '1.3.0', 'version of the file format')
# create MATRIX extension
n = len(emin)
de = emin[1]-emin[0] # assumes uniform bins
ener = 0.5*(emax+emin)
enoise = 0.090/2.35 # electronic noise (rms) in keV
rwidth = enoise*3.0 # half width of response in keV
rn = int(ceil(rwidth/de)) # width of response in channels = 2*rn+1
ngrp = 1+zeros(n, dtype=int16)
fchan = zeros(n)
nchan = (2*rn+1)+zeros(n, dtype=int16)
matr = zeros((n, 2*rn+1), dtype=float32)
des = de*arange(-rn, rn+de)
for i in range(n):
  e0 = max([ener[i]+min(des), ener[0]])
  f = int((e0-ener[0])/de)
  fchan[i] = f
  sigma = enoise # add Poission term
  norm = 1/sqrt(2*pi*sigma**2)
  r = de*norm*exp(-(des**2/(2*sigma**2)))
  matr[i,:] = r[:] # should integrate over de then divide by de
'''
  if i == 100:
    print ener[i], e0, f, ener[f]
    print ener[i]+des
    print ener[f+arange(0,2*rn+1,1,dtype=int32)]
    print r
    print sum(r),
    plt.clf()
    plt.xlabel('Energy (keV)')
    plt.ylabel('Response (cm^2)')
    plt.plot(ener[f+arange(0,2*rn+1,1,dtype=int32)], r, '-b')
    plt.show() # display the plot
'''

col1 = pyfits.Column(name='ENERG_LO', format='1E', array=emin)
col2 = pyfits.Column(name='ENERG_HI', format='1E', array=emax)
col3 = pyfits.Column(name='N_GRP', format='I', array=ngrp)
col4 = pyfits.Column(name='F_CHAN', format='I', array=fchan)
col5 = pyfits.Column(name='N_CHAN', format='I', array=nchan)
col6 = pyfits.Column(name='MATRIX', format=str(2*rn+1)+'E', array=matr)
# create a ColDefs (column-definitions) object for all columns:
cols = pyfits.ColDefs([col1, col2, col3, col4, col5, col6])
# create a new binary table HDU object
matrixhdu = pyfits.new_table(cols)
matrixhdu.name = 'MATRIX'
matrixhdu.header.update('EXTNAME', 'SPECRESP', 'Extension type')
matrixhdu.header.update('TELESCOP', 'HaloSat', 'Mission name')
matrixhdu.header.update('INSTRUME', 'SDD1', 'Instrument')
matrixhdu.header.update('CHANTYPE', 'PHA', 'Channels from detector electronics')
matrixhdu.header.update('DETCHANS', len(emin), 'total number of raw detector PHA channels')
matrixhdu.header.update('HDUCLASS', 'OGIP', 'file format is OGIP standard')
matrixhdu.header.update('HDUCLAS1', 'RESPONSE', 'extension contains response data')
matrixhdu.header.update('HDUCLAS2', 'RSP_MATRIX', 'extension contains a response matrix')
matrixhdu.header.update('HDUVERS', '1.3.0', 'version of the file format')
matrixhdu.header.update('TLMIN4', 1, 'first channel in the response')
matrixhdu.header.update('HDUCLAS3', 'REDIST', 'response for photon redistribution only')
# create primary FITS header
prihdr = pyfits.Header()
prihdu = pyfits.PrimaryHDU(header=prihdr)
# create HDUList with primary HDU and the data table extension
thdulist = pyfits.HDUList([prihdu, eboundshdu, matrixhdu])
# write to file
thdulist.writeto('halo_rmf.fits', clobber=True)


