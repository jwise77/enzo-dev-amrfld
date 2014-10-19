#!/usr/bin/env python
# Python script to compute the ionization energy for each 
# group within a multigroup simulation.
#
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import numpy as np

# disable overflow warning messages (happens inside exp)
np.seterr(over='ignore')

# relevant physical constants
pi = 3.14159265358979323846
hplanck = 6.6260693e-27
ev2erg = 1.60217653e-12
clight = 2.99792458e10
kboltz = 1.3806504e-16
nu0_HI = 13.6 / hplanck * ev2erg    # ionization threshold of HI (in Hz)
nu0_HeI = 24.6 / hplanck * ev2erg   # ionization threshold of HeI (in Hz)
nu0_HeII = 54.4 / hplanck * ev2erg  # ionization threshold of HeII (in Hz)


###### USER-MODIFIABLE DATA ######

# define the radiation energy groups (in eV)
energy_bands = ( (13.6, 24.6), (24.6, 54.4), (54.4, 100.0) )

# set the total source luminosity (in photon/s)
NGammaDot = 1e50

# define the radiation spectrum object
class SED:
    """ An SED object must have four member functions:

          monochromatic() -- returns true or false, depending on 
              whether the SED is monochromatic

          lower_bound() -- if the SED cuts off below a given frequency, 
              return that cutoff frequency here.  If there is no lower 
              bound, return a negative value.  If the SED is monochromatic, 
              return that frequency here.  The return value should be given 
              in eV, e.g. if the SED exists from the ionization threshold 
              of Hydrogen, this should be 13.6.

          upper_bound() -- if the SED cuts off above a given frequency, 
              return that cutoff frequency here.  If there is no upper 
              bound, return a negative value.  If the SED is monochromatic, 
              this is unused.  The value should be given in eV.

          value(hnu) -- the main SED function: implements the non-
              monochromatic SED, and should return the SED value at a given 
              frequency (the input is in units of eV).  If the SED is 
              monochromatic, this function should return a value of 1.

        This example SED implements a blackbody spectrum. """
    def __init__(self, T):
        self.T = T
    def monochromatic(self):
        return False
    def lower_bound(self):
        return 0.0
    def upper_bound(self):
        return -1.0
    def value(self, hnu):
        nu = hnu*ev2erg/hplanck   # convert frequency to Hz
        return (8.0*pi*hplanck*(nu/clight)**3/(np.exp(hplanck*nu/kboltz/self.T)-1.0))


# create SED object for T=1e5 blackbody spectrum to use in subsequent integration;
# MUST name the object "source_sed"
source_sed = SED(1e5) 



###### NON-MODIFIABLE DATA ######

# define the radiation spectrum moment object
class SED_moment:
    """ An SED_moment object is built off of a SED object, and is used for 
        computing the relevant integrals that determine the integrals that 
        define each radiation group emission energy. """
    def __init__(self, sed):
        # ensure that the sed passed in meets qualifications
        try:
            tst = sed.monochromatic()
        except:
            raise Exception('non-functional SED argument: monochromatic() failure')
        try:
            tst = sed.lower_bound()
        except:
            raise Exception('non-functional SED argument: lower_bound() failure')
        try:
            tst = sed.upper_bound()
        except:
            raise Exception('non-functional SED argument: upper_bound() failure')
        try:
            tst = sed.value(13.6)
        except:
            raise Exception('non-functional SED argument: value() failure')
        self.baseSED = sed
    def monochromatic(self):
        return (self.baseSED.monochromatic())
    def lower_bound(self):
        return (self.baseSED.lower_bound())
    def upper_bound(self):
        return (self.baseSED.upper_bound())
    def value(self, hnu):
        return (self.baseSED.value(hnu)*hnu*ev2erg/hplanck)

# create SED moment object based off of user-defined source_se
source_sed_moment = SED_moment(source_sed)



###################
# define the numerical integration utility functions

def SED_composite_integral(sed, a, b, n):
    """ O(h^16) accurate composite Gaussian quadrature of a given integrand
        over a definite interval [a,b], using n subintervals.
    
        Arguments:
           sed -- SED class object
           a, b -- end points of definite interval [a,b]; Require: a<b and a>=0
                   These should be given in units of eV.
           n -- number of subintervals to use in Gaussian quadrature.  """

    # nodes/weights defining quadrature method
    hwid = (b-a)/n      # subinterval width
    x = ( -0.18343464249564980493, 0.18343464249564980493, -0.52553240991632898581,
          0.52553240991632898581, -0.79666647741362673959, 0.79666647741362673959,
          -0.96028985649753623168, 0.96028985649753623168 )
    w = ( 0.36268378337836198296, 0.36268378337836198296, 0.31370664587788728733,
          0.31370664587788728733, 0.22238103445337447054, 0.22238103445337447054,
          0.10122853629037625915, 0.10122853629037625915 )

    F = 0.0                              # initialize result
    for i in range(n):
        xmid = a + (i+0.5)*hwid          # subinterval midpoint
        for j in range(len(x)):          # loop over quadrature nodes
            z = xmid + 0.5*hwid*x[j]
            F += w[j]*sed.value(z)       # add contribution to result

    return (0.5*hwid*F)                  # return final result, scaled by hwid/2
 

def SED_composite_integral2(sed, a, n):
    """ O(h^16) accurate composite Gaussian quadrature of a given integrand
        over an indefinite interval [a,infinity], that is remapped to [0,1], 
        using n subintervals
    
        Arguments:
           sed -- SED class object
           a -- lower limit of integration, [a,infinity]; Require: a>=0.
                Should be given in units of eV.
           n -- number of subintervals to use in Gaussian quadrature.  """

    # remapped interval
    a2 = 0.0
    b2 = 1.0

    # nodes/weights defining quadrature method
    hwid = (b2-a2)/n         # subinterval width
    nodes = 8                # num quadrature nodes
    x = ( -0.18343464249564980493, 0.18343464249564980493, -0.52553240991632898581,
          0.52553240991632898581, -0.79666647741362673959, 0.79666647741362673959,
          -0.96028985649753623168, 0.96028985649753623168 )
    w = ( 0.36268378337836198296, 0.36268378337836198296, 0.31370664587788728733,
          0.31370664587788728733, 0.22238103445337447054, 0.22238103445337447054,
          0.10122853629037625915, 0.10122853629037625915 )

    F = 0.0                                # initialize result
    for i in range(n):
        xmid = a2 + (i+0.5)*hwid           # subinterval midpoint
        for j in range(len(x)):            # loop over quadrature nodes
            z = xmid + 0.5*hwid*x[j]
            F += w[j]*sed.value(a/z)/z/z   # add contribution to result

    return (a*0.5*hwid*F)                  # return final result, scaled by a*hwid/2


def SED_integral(sed, a, b, convertHz):
    """ main integration function: calls either SED_composite_integral() 
        or SED_composite_integral2(), depending on the nature of the 
        integration interval.
    
        Arguments:
           sed -- SED class object
           a, b -- end points of definite interval [a,b]; ; Require: a<b and a>=0
                   or for an indefinite interval [a,infinity], set b<0.
                   These should be given in units of eV.
           convertHz -- flag denoting whether the output from the integral
                        (performed over hnu in eV) should be converted to Hz. """


    # initialize result
    R = 0.0

    # return zero on an illegal interval
    if ((a > b) and (a > 0) and (b > 0)):
        print "illegal interval: a =", a, ", b =", b, ", returning zero"
        return R

    # if the sed is monochromatic:
    # return sed value if the sed frequency lies inside the interval, otherwise return 0
    if (sed.monochromatic()):
        if (sed.lower_bound() <= 0):
            R = 0
        elif ((a > 0) and (sed.lower_bound() < a)):
            R = 0
        elif ((b > 0) and (sed.lower_bound() > b)):
            R = 0
        else:
            R = sed.value(sed.lower_bound())
        return R

    # set local variables for integration bounds
    nu_L = max(max(a, 0), sed.lower_bound())
    nu_R = b

    b_infinite = False
    if ((b <= 0) and (sed.upper_bound() <= 0)):
        b_infinite = True

    if ((not b_infinite) and (b > 0) and (sed.upper_bound() > 0)):
        nu_R = min(sed.upper_bound(), b)
    if ((not b_infinite) and (b > 0) and (sed.upper_bound() <= 0)):
        nu_R = b
    if ((not b_infinite) and (b <= 0) and (sed.upper_bound() > 0)):
        nu_R = sed.upper_bound()

    # if the integral has now disappeared, return 0
    if ((not b_infinite) and (nu_L >= nu_R)):
        return R

    # set the number of subintervals and compute approximation
    N = 5000
    if (b_infinite):   # improper integral (infinite upper bound)
        R = SED_composite_integral2(sed, nu_L, N)
    else:              # definite integral
        R = SED_composite_integral(sed, nu_L, nu_R, N)

    # R is computed based on integration in eV; scale result if Hz was desired
    if (convertHz):
        R *= ev2erg/hplanck

    return R



###################
# perform the integration and output results

# compute the total SED integral over all ionizing frequencies
total_integral = SED_integral(source_sed, 13.6, -1.0, True)

# loop over radiation bins
nbins = len(energy_bands)
for i in range(nbins):

    # set this bin
    freq_bin = energy_bands[i]
    
    # skip monochromatic radiation bands
    if (freq_bin[1] <= freq_bin[0]):
        print ' Group', i,'energy = 0 (monochromatic, so cannot integrate)'
        continue

    # computing source contributions to this radiation field
    SourceGroupEnergy = (NGammaDot * hplanck / total_integral * 
                         SED_integral(source_sed_moment, freq_bin[0], freq_bin[1], True))

    print ' Group',i,'energy =', "{0:.16e}".format(SourceGroupEnergy)
 


# end of script
