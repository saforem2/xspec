#!/usr/bin/python2

from xspec import *
import matplotlib.pyplot as plt

'''
NOTE: many of the plotting commands below have been commented out
in order to speed up runtime
'''

def fitData():

    # Xspec Settings
    Xset.abund = "angr"
    Xset.cosmo = "70 0 0.73"
    Xset.xsect = "bcmc"

    norm_parameter = 0.25
    m1 = Model("TBabs(vapec + powerlaw + gaussian) + vapec", setPars=(0.018,0.225,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0.278,1.52,.376,.597,0.0503,norm_parameter,0.099,1,1,1,1,1,1,1,1,1,1,1,1,0,0.66))

    for i in range(1,39):
        m1(i).frozen = True
        if i == 20 or i == 21 or i == 22:
            m1(i).frozen = False

    fs1 = FakeitSettings(response="halo_rmf.fits", arf="halo_arf.fits",exposure=10000.)
    AllData.fakeit(1,[fs1])
    #AllData.ignore("**-0.1")
    data = AllData(1)
    Fit.statMethod = "cstat"
    Fit.query = "yes"
    Fit.perform()
    Fit.error("22")

    Plot.device = "/xs"
    Plot.setRebin(3.0,12)
    Plot.add = 1
    Plot.xAxis = "keV"
    Plot("data ratio")

    '''
    To access the value of a parameter (k) defined in the model (m1)
    par = m1(k) = [val, delta, min, bot, top, max]
    then, the value of interest is
    par[0] = parameter_value
    '''
    # store 'norm' value, the total counts for our gaussian model
    norm_par = m1(22)
    norm_val = norm_par.values[0]
    err_pm = list(norm_par.error)
    return_vals = [norm_val, err_pm]


    #plt.plot(data.noticed,data.values,'r.',data.noticed,m1.folded(1),ms=0.6)
    #plt.xlabel('Channels')
    #plt.ylabel('counts/cm$^{2}$/sec/chan')
    #plt.errorbar(data.values,m1.folded(1),yerr=yerr).v
    #plt.xlim(0,200)
    #plt.savefig('sim3.pdf',type='pdf')
    #plt.show()

    return return_vals
