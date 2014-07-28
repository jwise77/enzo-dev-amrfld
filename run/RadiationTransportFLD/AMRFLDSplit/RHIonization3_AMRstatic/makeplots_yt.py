#!/usr/bin/env 
# yt-based plotting script for Iliev et al. test #3
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
from yt.mods import *
from os import *
import matplotlib.colorbar as cb


# set the graphics output type
#pictype = 'pdf'
pictype = 'png'


# Define some derived fields
#   Neutral Hydrogen fraction (log plot)
def _xHI(field, data):
    return (data["HI_Density"]/(data["HII_Density"]+data["Density"]))

add_field("xHI", take_log=True, function=_xHI, 
          display_name="Neutral\; Fraction")

#   Ionized Hydrogen fraction (log plot)
def _xHII(field, data):
    return (data["HII_Density"]/(data["HI_Density"]+data["HII_Density"]))

add_field("xHII", take_log=True, function=_xHII, 
          display_name="Ionized\; Fraction")

#   Temperature (log plot)
def _logT(field, data):
    tmp = (data["Temperature"])
    return tmp

add_field("logT", take_log=True, function=_logT, 
          display_name="Temperature", units=r"\rm{K}")

#   Pressure (log plot)
def _logP(field, data):
    mp = 1.67262171e-24
    kb = 1.3806504e-16
    tmp = (data["Temperature"] * (data["HI_Density"] + 2.0*data["HII_Density"]))
    return (tmp*kb/mp)

add_field("logP", take_log=True, function=_logP, 
          display_name="Pressure", units=r"\rm{g} \rm{cm}^{-1} \rm{s}^{-2}")




# loop over snapshots, loading values and times
tsteps = (1, 5, 10, 50)
for it in range(len(tsteps)):
    tstep = tsteps[it]
    
    # load relevant information
    sdump = repr(tstep).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    pf = load(pfile)
    t = pf.current_time * pf["TimeUnits"]
    
    # set time label
    if (tstep == tsteps[0]):
        Myr = '01'
        linestyle = ':'
        style = 'r:'
        lab = '1 Myr'
    elif (tstep == tsteps[1]):
        Myr = '05'
        linestyle = '-.'
        style = 'g-.'
        lab = '5 Myr'
    elif (tstep == tsteps[2]):
        Myr = '10'
        linestyle = '--'
        style = 'b--'
        lab = '10 Myr'
    else:
        Myr = '50'
        linestyle = '-'
        style = 'k-'
        lab = '50 Myr'
        
    # determine domain "center" for plotting
    xL = pf["DomainLeftEdge"][0]
    xR = pf["DomainRightEdge"][0]
    xC = 0.5*(xL + xR)
    yL = pf["DomainLeftEdge"][1]
    yR = pf["DomainRightEdge"][1]
    yC = 0.5*(yL + yR)
    zL = pf["DomainLeftEdge"][2]
    zR = pf["DomainRightEdge"][2]
    zC = 0.5*(zL + zR)

    # begin multi-plot plot collection
#    fig, axes, colorbars = get_multi_plot(2, 2, colorbar='vertical', bw=4)
    fig, axes, colorbars = get_multi_plot(2, 2, colorbar=None, bw=4)
    pc = PlotCollection(pf, [xC,yC,zC])
        
    # add slices to plot
    p = pc.add_slice("xHI", 'z', figure=fig, axes=axes[0][0], use_colorbar=False)
    p = pc.add_slice("logP", 'z', figure=fig, axes=axes[0][1], use_colorbar=False)
    p = pc.add_slice("Density", 'z', figure=fig, axes=axes[1][0], use_colorbar=False)
    p = pc.add_slice("logT", 'z', figure=fig, axes=axes[1][1], use_colorbar=False)
    for p, cax in zip(pc.plots, colorbars):
        cbar = cb.Colorbar(cax, p.image, orientation='vertical')
        p.colorbar = cbar
        p._autoset_label()

    fig.savefig('slices_' + Myr + 'Myr' % pf)

    # generate profile plots along line (x,0,0)
    #    extract data along ray through center of box
    r = pf.h.ray([xL, yC, zC], [xR, yC, zC])
    dProf = r["Density"]
    TProf = r["logT"]
    pProf = r["logP"]
    xProf = r["xHI"]
    rvals = r["x"]

    #    sort results by x (since that isn't quite working correctly from yt)
    ptype = [('r', float), ('x', float), ('p', float), ('T', float), ('d', float)]
    pdata = np.zeros(rvals.shape, dtype=ptype);
    nr = (rvals.shape)[0]
    for ir in range(0,nr):
        pdata[ir] = (rvals[ir], xProf[ir], pProf[ir], TProf[ir], dProf[ir])
    pdata = np.sort(pdata, order='r')
            
    # interpolate results into output arrays
    xvals = np.linspace(xL, xR, 200)
    xHIProfile = np.interp(xvals, pdata['r'], pdata['x'])
    PresProfile = np.interp(xvals, pdata['r'], pdata['p'])
    TempProfile = np.interp(xvals, pdata['r'], pdata['T'])
    DensProfile = np.interp(xvals, pdata['r'], pdata['d'])

    #    generate stored profile plots
    figure(2)
    subplot(221)
    semilogy(xvals,DensProfile,style,hold=True)
    ylabel('Density (g/cm$^3$)')
    ax = gca()
    ax.yaxis.set_ticks_position("left")
    ax.yaxis.set_label_position("left")
#    axis([ 0.0, 1.0, 1e-28, 1e-24 ])
    subplot(222)
    semilogy(xvals,TempProfile,style,hold=True,label=lab)
    ylabel('Temperature (K)')
    ax = gca()
    ax.yaxis.set_ticks_position("right")
    ax.yaxis.set_label_position("right")
#    axis(( 0.0, 1.0, 1e2, 1e5 ))
    subplot(223)
    semilogy(xvals,xHIProfile,style,hold=True)
    ylabel('Neutral Fraction')
    ax = gca()
    ax.yaxis.set_ticks_position("left")
    ax.yaxis.set_label_position("left")
#    axis(( 0.0, 1.0, 1e-5, 2 ))
    xlabel('x/L$_{box}$')
    subplot(224)
    semilogy(xvals,PresProfile,style,hold=True)
    ylabel('Pressure (g cm$^{-1}$ s$^{-2}$)')
    ax = gca()
    ax.yaxis.set_ticks_position("right")
    ax.yaxis.set_label_position("right")
#    axis(( 0.0, 1.0, 1e-16, 1e-12 ))
    xlabel('x/L$_{box}$')

    

# save accumulated profile plot
figure(2)
subplot(222)
legend( loc=3 )
savefig('profiles.' + pictype)
