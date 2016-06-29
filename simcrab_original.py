#!/usr/bin/python2

from xspec import *
import matplotlib.pyplot as plt

# Xspec Settings
Xset.abund = "angr"
Xset.cosmo = "70 0 0.73"
Xset.xsect = "bcmc"




# PyXspec operations:
m1 = Model("wabs(pow)", setPars=(0.018, 1, 1))

# or, can also set parameters by
# m1.setPars(a,b,c,...)

fs1 = FakeitSettings(response="halo_rmf.fits",arf="halo_arf.fits",exposure=10000.)
AllData.fakeit(1,[fs1])
data = AllData(1)
Fit.statMethod = "cstat"
Fit.perform()


Plot.device = "/xs"
#Plot.yAxis = "Normalized counts s$^{-1}$ keV$^{-1}$"
Plot.setRebin(3.0,12)
Plot.xAxis = "keV"
Plot("data")
