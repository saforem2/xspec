#!/usr/bin/python2

from xspec import *
import matplotlib.pyplot as plt


def fitCrab():
    # Xspec Settings
    Xset.abund = "angr"
    Xset.cosmo = "70 0 0.73"
    Xset.xsect = "bcmc"

    # PyXspec operations:
    m1 = Model("wabs(pow + gaussian)", setPars=(0.018, 1.52, 0.376, 0.7,0.1,0.15))
    #m1 = Model("wabs(pow + gaussian)", setPars=(0.018))
    m1(1).frozen = True
    m1(4).frozen = True
    m1(5).frozen = True

   # or, can also set parameters by
    # m1.setPars(a,b,c,...)

    fs1 = FakeitSettings(response="halo_rmf.fits",arf="halo_arf.fits",exposure=10000.)
    AllData.fakeit(1,[fs1])
    data = AllData(1)
    Fit.statMethod = "cstat"
    Fit.query = "yes"
    Fit.perform()

    Plot.device = "/xs"
    Plot.setRebin(3.0,12)
    Plot.add = 1
    Plot.xAxis = "keV"
    Plot("data residuals")

    # store 'sigma' value, the width of the gaussian in spectrum
    sig_par = m1(5)
    sig_val = sig_par.values[0]
    err_pm = list(sig_par.error)
    return_vals = [sig_val, err_pm]

    return return_vals
