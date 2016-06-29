# Python program to load a synphot FITS spectrum and display it

# import needed extensions
from numpy import *
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt # plotting package
from simcrab import fitCrab
from xspec import *
import os, time

# read in the transmission data
ener, c2, c1 = transpose(genfromtxt('C1C2eff.dat'))
# remove points below ethreshold
ethreshold = 0.1 # keV
q = where(ener >= ethreshold)
ener = ener[q]
c1 = c1[q]
c2 = c2[q]
# replace with evenly spaced energies
t = max(ener)+ener[1]-ener[0]
ener = arange(min(ener), t, (t-min(ener))/len(ener))
# should add silicon dead layer
area = 0.25 # detector area in cm^2

'''
#plt.ion() # do plots in interactive mode
# plot the spectrum
fig1 = plt.figure(1)
plt.clf()
plt.xlabel('Energy (keV)')
plt.ylabel('Response (cm^2)')
# plt.xscale('log')
plt.plot(ener, area*c1 , '-b')
plt.plot(ener, area*c2 , '-r')
#plt.show() # display the plot
plt.savefig('effective_area_energy.pdf',format='pdf')
'''
# translate into energy bin lower and upper bounds and average response
n = len(ener)
emin = float32(ener[0:n-1])
emax = float32(ener[1:n])
resp = float32(area*0.5*(c1[0:n-1]+c1[1:n]))
'''
fig2 = plt.figure()
plt.plot(0.5*(emin+emax), resp, '-k')
#plt.show()
plt.savefig('e_min_e_max_response.pdf',format='pdf')
'''

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

# createMatrix returns the total counts with error bounds
# for a given energy resolution, e_res
def createMatrix(e_res):
    n = len(emin)
    de = emin[1]-emin[0] # assumes uniform bins
    ener = 0.5*(emax+emin)
    #enoise = 0.090/2.35 # electronic noise (rms) in keV
    enoise = e_res/2.35
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
    return fitCrab()

# create evenly spaced values of the energy resolution
# from 50 eV to 150 eV in steps of 5 eV
ener_res = arange(1E-3,500E-3,5E-3)
# initialize empty array to store [total counts, lower bound, upper bound]
tot_count_err = []

# Loop over all energy resolutions
# and append NEW values of [total counts, lower bound, upper bound]
for i in range(len(ener_res)):
    vals = createMatrix(ener_res[i])
    norm_val = vals[0]
    err_m = list(vals[1])[0]
    err_p = list(vals[1])[1]
    tot_count_err.append([norm_val,err_m,err_p])
    #Plot.setRebin(3,12)
    #Plot.add = 1
    Plot.xAxis = "keV"
    #Plot("data model")
#    time.sleep(1)

# convert the list (tot_count_err) to numpy array (x)
x = array(tot_count_err)
# define error bars
yerr_m = x[:,0] - x[:,1]
yerr_p = x[:,2] - x[:,0]

# scale energy resolution from keV ---> eV for plotting
ener_res_plot = 1000*ener_res
fig3 = plt.figure()
plt.errorbar(ener_res_plot, x[:,0], yerr = (yerr_m,yerr_p), marker='o')
plt.grid(True)
plt.xlabel('$\Delta E$ (eV)')
plt.ylabel('counts')
plt.title('Total Counts vs. Energy Resolution $(\Delta E)$')
plt.savefig('counts_vs_ener_res.pdf',type='pdf')
plt.show()

# clean up the unnecessary .fak data-files created from fakeit command
filelist = [ f for f in os.listdir(".") if f.endswith(".fak") ]
for f in filelist:
    os.remove(f)
