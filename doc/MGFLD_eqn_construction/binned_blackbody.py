#!/usr/bin/env python
# matplotlib-based plotting script
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import numpy as np
import matplotlib

# update font sizes
matplotlib.rcParams.update({'font.size': 24})  # tick labels
fs = 26                                        # other text

# declare some constants
h = 6.6260693e-27;         # Planck's constant [ergs*s]
kb = 1.3806504e-16;        # Boltzmann's constant [ergs/K]
c = 2.99792458e10;         # speed of light [cm/s]
ev2erg = 1.60217653e-12;   # conversion constant from eV to ergs

# declare plot/bin bounds
nuL = 1.0
nuR = 100.0

# blackbody spectrum function
def blackbody_sed(hnu):
    """ Accepts argument in eV """
    nu = hnu*ev2erg/h
    sigma = 8.0*np.pi*h*(nu/c)**3/(exp(h*nu/kb/1.e5)-1.0)*1e11
    return sigma

# set plot frequencies
N = 1000
hnu = np.linspace(nuL,nuR,N)

# set true spectrum
blackbody = np.zeros(N)
for i in arange(N):
    blackbody[i] = blackbody_sed(hnu[i])

# plot true spectrum
fig = figure()
#plot(hnu, blackbody)
fill_between(hnu, 0, blackbody)
xlabel(r'$h\nu$', fontsize=fs)
ylabel('SED', fontsize=fs)
title('$T=10^5$ Blackbody', fontsize=fs)
grid()
savefig('blackbody.pdf')


### 5-bin approximation ###

# set bins for approximation
Nbins = 5
hnus = np.linspace(nuL,nuR,Nbins+1)

# set bin values by filling each with Simpson's rule from true
# spectrum, and dividing by width to get average value
bin_vals = np.zeros(Nbins)
for i in arange(Nbins):
    hnu_l = hnus[i]
    hnu_r = hnus[i+1]
    hnu_c = 0.5*(hnu_l+hnu_r)
    bin_vals[i] = 1.0/6.0*(1.0*blackbody_sed(hnu_l) +
                           4.0*blackbody_sed(hnu_c) +
                           1.0*blackbody_sed(hnu_r))


# binned blackbody spectrum function
def binned_sed(bins, bsed, hnu):
    """ Accepts argument: 
        bins - frequency bounds for each bin
        bsed - SED values for each bin
        hnu  - evaluation frequency (in eV) """
    ibin = 0
    for i in arange(len(bins)-1):
        if ((hnu >= bins[i]) and (hnu < bins[i+1])):
            ibin = i
            break
    sigma = bsed[ibin]
    return sigma


# set binned spectrum
binned_blackbody = np.zeros(N)
for i in arange(N):
    binned_blackbody[i] = binned_sed(hnus, bin_vals, hnu[i])

# plot binned spectrum
fig = figure()
#plot(hnu, binned_blackbody)
fill_between(hnu, 0, binned_blackbody)
xlabel(r'$h\nu$', fontsize=fs)
ylabel('SED', fontsize=fs)
title('5-Bin Approximation', fontsize=fs)
grid()
savefig('blackbody-5bin.pdf')




### 10-bin approximation ###

# set bins for approximation
Nbins = 10
hnus = np.linspace(nuL,nuR,Nbins+1)

# set bin values by filling each with Simpson's rule from true
# spectrum, and dividing by width to get average value
bin_vals = np.zeros(Nbins)
for i in arange(Nbins):
    hnu_l = hnus[i]
    hnu_r = hnus[i+1]
    hnu_c = 0.5*(hnu_l+hnu_r)
    bin_vals[i] = 1.0/6.0*(1.0*blackbody_sed(hnu_l) +
                           4.0*blackbody_sed(hnu_c) +
                           1.0*blackbody_sed(hnu_r))


# binned blackbody spectrum function
def binned_sed(bins, bsed, hnu):
    """ Accepts argument: 
        bins - frequency bounds for each bin
        bsed - SED values for each bin
        hnu  - evaluation frequency (in eV) """
    ibin = 0
    for i in arange(len(bins)-1):
        if ((hnu >= bins[i]) and (hnu < bins[i+1])):
            ibin = i
            break
    sigma = bsed[ibin]
    return sigma


# set binned spectrum
binned_blackbody = np.zeros(N)
for i in arange(N):
    binned_blackbody[i] = binned_sed(hnus, bin_vals, hnu[i])

# plot binned spectrum
fig = figure()
#plot(hnu, binned_blackbody)
fill_between(hnu, 0, binned_blackbody)
xlabel(r'$h\nu$', fontsize=fs)
ylabel('SED', fontsize=fs)
title('10-Bin Approximation', fontsize=fs)
grid()
savefig('blackbody-10bin.pdf')

