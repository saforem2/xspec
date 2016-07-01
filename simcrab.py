#!/usr/bin/python2

from xspec import *
import matplotlib.pyplot as plt

Plot.device = "/xs"
global z

def fitCrab():
    # Xspec Settings
    Xset.abund = "angr"
    Xset.cosmo = "70 0 0.73"
    Xset.xsect = "bcmc"
    global m1
    # PyXspec operations:
    m1 = Model("wabs(pow + gaussian)", setPars=(0.018, 1.52, 0.3, .7, 0.1,.25))
    sigma_par = m1.gaussian.Sigma
    sigma_par.values = "0.1, ,75E-3,75E-3 , , "
    #m1 = Model("wabs(gaussian)", setPars=(0.018, 0.7, 0.1, 0.25))
    #m1 = Model("wabs(pow + gaussian)", setPars=(0.018))
    m1(1).frozen = True
    m1(2).frozen = True
    #m1(3).frozen = True
    m1(4).frozen = True
    m1(6).frozen = True

   # or, can also set parameters by
    # m1.setPars(a,b,c,...)

    fs1 = FakeitSettings(response="halo_rmf.fits",arf="halo_arf.fits",exposure=10000.)
    AllData.fakeit(1,[fs1])
    data = AllData(1)
    #data.ignore("**-0.35 1.75-**")
    Fit.statMethod = "cstat"
    Fit.query = "yes"
    Fit.perform()
    Plot.setRebin(3.0,12)
    Plot.add = 1
    Plot.xAxis = "keV"
    Plot("data model")
    # store 'sigma' value, the width of the gaussian in spectrum
    sig_par = m1(5)
    sig_val = sig_par.values[0]
    err_pm = list(sig_par.error)
    return_vals = [sig_val, err_pm,AllModels.eqwidth(3)]

    return return_vals
