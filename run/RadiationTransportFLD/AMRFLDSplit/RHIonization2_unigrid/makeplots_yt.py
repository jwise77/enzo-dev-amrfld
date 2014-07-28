#!/usr/bin/env python
# yt-based plotting script for Iliev et al. test #2
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
from yt.mods import *
from os import *

# set the total number of snapshots
te = 50

# set the graphics output type
#pictype = 'pdf'
pictype = 'png'

# set some constants
Ngammadot = 5.0e48     # ionization source strength [photons/sec]
aHII = 2.52e-13        # recombination rate coefficient
mp = 1.67262171e-24    # proton mass [g]
Myr = 3.15576e13       # duration of a Megayear [sec]
nH = 1.0e-3            # input hydrogen number density [cm^(-3)]
trec = 1.0/(aHII*nH)   # recombination time [sec]
rs0 = (3.0*Ngammadot/4/pi/aHII/nH/nH)**(1.0/3.0)   # Stromgren radius


# Define some derived fields
#   Neutral Hydrogen fraction (log plot)
def _xHI(field, data):
    return (data["HI_Density"]/data["Density"])
add_field("xHI", take_log=True, function=_xHI, 
          display_name="Neutral\; Fraction")

#   Ionized Hydrogen fraction (log plot)
def _xHII(field, data):
    return (data["HII_Density"]/data["Density"])
add_field("xHII", take_log=True, function=_xHII, 
          display_name="Ionized\; Fraction")

#   Radiation energy density (log plot)
def _logE0(field, data):
    return (data["Radiation0"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE0", take_log=True, function=_logE0, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Group\; 0", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Radiation energy density (log plot)
def _logE1(field, data):
    return (data["Radiation1"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE1", take_log=True, function=_logE1, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Group\; 1", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Radiation energy density (log plot)
def _logE2(field, data):
    return (data["Radiation2"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE2", take_log=True, function=_logE2, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Group\; 2", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Radiation energy density (log plot)
def _logE3(field, data):
    return (data["Radiation3"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE3", take_log=True, function=_logE3, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Group\; 3", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Radiation energy density (log plot)
def _logE4(field, data):
    return (data["Radiation4"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE4", take_log=True, function=_logE4, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Group\; 4", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Radiation energy density (log plot)
def _logE5(field, data):
    return (data["Radiation5"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE5", take_log=True, function=_logE5, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Group\; 5", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Radiation energy density (log plot)
def _logE6(field, data):
    return (data["Radiation6"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE6", take_log=True, function=_logE6, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Group\; 6", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Radiation energy density (log plot)
def _logE7(field, data):
    return (data["Radiation7"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE7", take_log=True, function=_logE7, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Group\; 7", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Radiation energy density (log plot)
def _logE8(field, data):
    return (data["Radiation8"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE8", take_log=True, function=_logE8, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Group\; 8", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Radiation energy density (log plot)
def _logE9(field, data):
    return (data["Radiation9"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE9", take_log=True, function=_logE9, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Group\; 9", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Temperature (log plot)
def _logT(field, data):
    mp = 1.67262171e-24
    kb = 1.3806504e-16
    gamma = 5.0/3.0
    tmp = (data["TotalEnergy"] * data["Density"] / 
           (2.0*data["Density"] - data["HI_Density"]))
    return ((gamma - 1.0)*mp/kb * tmp)
add_field("logT", take_log=True, function=_logT, 
          display_name="Temperature", units=r"\rm{K}")

#   Radius from domain center
def _radius(field, data):
    return (np.sqrt(data["x"]*data["x"] + data["y"]*data["y"] +
                    data["z"]*data["z"]))
add_field("radius", take_log=False, function=_radius, 
          display_name="radius", units=r"{r/L_{box}}")



# initialize time-history outputs
#    row 1: time (t)
#    row 2: computed i-front radius
rdata = zeros( (2, te+1), dtype=float);


# loop over snapshots, loading values and times
for tstep in range(0,te+1):
    
    # load relevant information
    sdump = repr(tstep).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    pf = load(pfile)
    t = pf.current_time * pf["TimeUnits"]

    # determine if simulation was run with source in center or corner
    spherical = (2.0**(pf.domain_left_edge[0]/pf.domain_right_edge[0]+1.0) 
               * 2.0**(pf.domain_left_edge[1]/pf.domain_right_edge[1]+1.0) 
               * 2.0**(pf.domain_left_edge[2]/pf.domain_right_edge[2]+1.0))

    # compute I-front radius (assuming spherical)
    sp = pf.h.sphere([0.0, 0.0, 0.0], 1.0)
    HIIvolume = (sp["xHII"]*sp["CellVolumeCode"]*pf["cm"]**3).sum()*spherical
    radius = (3.0/4.0*HIIvolume/pi)**(1.0/3.0)
    
    # store data
    rdata[0][tstep] = t/trec
    rdata[1][tstep] = radius
    
    # generate 2D plots at certain times
    if (tstep == 1) or (tstep == 5) or (tstep == 10) or (tstep == 20) or (tstep == 50):
        
        # set time label
        if (tstep == 1):
            Myr = '010'
        elif (tstep == 5):
            Myr = '050'
        elif (tstep == 10):
            Myr = '100'
        elif (tstep == 20):
            Myr = '200'
        else:
            Myr = '500'
        
        # determine domain "center" for plotting
        xC = 0.5*(pf["DomainLeftEdge"][0] + pf["DomainRightEdge"][0])
        yC = 0.5*(pf["DomainLeftEdge"][1] + pf["DomainRightEdge"][1])
        zC = 0.5*(pf["DomainLeftEdge"][2] + pf["DomainRightEdge"][2])

        # begin plot collection
        pc = PlotCollection(pf, [xC,yC,0.0])
        
        # xHI slice through z=0
        p = pc.add_slice("xHI",'z')
        p.modify["title"]('HI fraction, t =' + Myr + ' Myr')

        p = pc.add_slice("xHII",'z')
        p.modify["title"]('HII fraction, t =' + Myr + ' Myr')

        numradfields = 0
        try:
            p = pc.add_slice("logE0",'z')
            p.modify["title"]('radiation density 0, t =' + Myr + ' Myr')
            numradfields += 1
        except:
            pass

        try:
            p = pc.add_slice("logE1",'z')
            p.modify["title"]('radiation density 1, t =' + Myr + ' Myr')
            numradfields += 1
        except:
            pass

        try:
            p = pc.add_slice("logE2",'z')
            p.modify["title"]('radiation density 2, t =' + Myr + ' Myr')
            numradfields += 1
        except:
            pass

        try:
            p = pc.add_slice("logE3",'z')
            p.modify["title"]('radiation density 3, t =' + Myr + ' Myr')
            numradfields += 1
        except:
            pass

        try:
            p = pc.add_slice("logE4",'z')
            p.modify["title"]('radiation density 4, t =' + Myr + ' Myr')
            numradfields += 1
        except:
            pass

        try:
            p = pc.add_slice("logE5",'z')
            p.modify["title"]('radiation density 5, t =' + Myr + ' Myr')
            numradfields += 1
        except:
            pass

        try:
            p = pc.add_slice("logE6",'z')
            p.modify["title"]('radiation density 6, t =' + Myr + ' Myr')
            numradfields += 1
        except:
            pass

        try:
            p = pc.add_slice("logE7",'z')
            p.modify["title"]('radiation density 7, t =' + Myr + ' Myr')
            numradfields += 1
        except:
            pass

        try:
            p = pc.add_slice("logE8",'z')
            p.modify["title"]('radiation density 8, t =' + Myr + ' Myr')
            numradfields += 1
        except:
            pass

        try:
            p = pc.add_slice("logE9",'z')
            p.modify["title"]('radiation density 9, t =' + Myr + ' Myr')
            numradfields += 1
        except:
            pass

        p = pc.add_slice("logT",'z')
        p.modify["title"]('temperature, t =' + Myr + ' Myr')

        pc.save(Myr + 'Myr', format=pictype)

        # rename generated files
        if (numradfields > 0):
            f1 = Myr + 'Myr_Slice_z_logE0.' + pictype
            f2 = 'E0contour_' + Myr + 'Myr.' + pictype
            rename(f1,f2)
        if (numradfields > 1):
            f1 = Myr + 'Myr_Slice_z_logE1.' + pictype
            f2 = 'E1contour_' + Myr + 'Myr.' + pictype
            rename(f1,f2)
        if (numradfields > 2):
            f1 = Myr + 'Myr_Slice_z_logE2.' + pictype
            f2 = 'E2contour_' + Myr + 'Myr.' + pictype
            rename(f1,f2)
        if (numradfields > 3):
            f1 = Myr + 'Myr_Slice_z_logE3.' + pictype
            f2 = 'E3contour_' + Myr + 'Myr.' + pictype
            rename(f1,f2)
        if (numradfields > 4):
            f1 = Myr + 'Myr_Slice_z_logE4.' + pictype
            f2 = 'E4contour_' + Myr + 'Myr.' + pictype
            rename(f1,f2)
        if (numradfields > 5):
            f1 = Myr + 'Myr_Slice_z_logE5.' + pictype
            f2 = 'E5contour_' + Myr + 'Myr.' + pictype
            rename(f1,f2)
        if (numradfields > 6):
            f1 = Myr + 'Myr_Slice_z_logE6.' + pictype
            f2 = 'E6contour_' + Myr + 'Myr.' + pictype
            rename(f1,f2)
        if (numradfields > 7):
            f1 = Myr + 'Myr_Slice_z_logE7.' + pictype
            f2 = 'E7contour_' + Myr + 'Myr.' + pictype
            rename(f1,f2)
        if (numradfields > 8):
            f1 = Myr + 'Myr_Slice_z_logE8.' + pictype
            f2 = 'E8contour_' + Myr + 'Myr.' + pictype
            rename(f1,f2)
        if (numradfields > 9):
            f1 = Myr + 'Myr_Slice_z_logE9.' + pictype
            f2 = 'E9contour_' + Myr + 'Myr.' + pictype
            rename(f1,f2)
        f1 = Myr + 'Myr_Slice_z_logT.' + pictype
        f2 = 'TempContour_' + Myr + 'Myr.' + pictype
        rename(f1,f2)
        f1 = Myr + 'Myr_Slice_z_xHI.' + pictype
        f2 = 'HIcontour_' + Myr + 'Myr.' + pictype
        rename(f1,f2)
        f1 = Myr + 'Myr_Slice_z_xHII.' + pictype
        f2 = 'HIIcontour_' + Myr + 'Myr.' + pictype
        rename(f1,f2)

        # generate profile plots by averaging results from multiple rays
#         rays = np.array( ((1.0,0.0,0.0), (0.0,1.0,0.0), (0.0,0.0,1.0), 
#                           (-1.0,0.0,0.0), (0.0,-1.0,0.0), (0.0,0.0,-1.0), 
#                           (1.0,1.0,0.0), (1.0,0.0,1.0), (0.0,1.0,1.0), 
#                           (-1.0,1.0,0.0), (-1.0,0.0,1.0), (0.0,-1.0,1.0), 
#                           (1.0,-1.0,0.0), (1.0,0.0,-1.0), (0.0,1.0,-1.0), 
#                           (-1.0,-1.0,0.0), (-1.0,0.0,-1.0), (0.0,-1.0,-1.0), 
#                           (1.0,1.0,1.0), (-1.0,1.0,1.0), (1.0,-1.0,1.0), 
#                           (1.0,1.0,-1.0), (-1.0,-1.0,1.0), (-1.0,1.0,-1.0), 
#                           (1.0,-1.0,-1.0), (-1.0,-1.0,-1.0)) )
#         nrays = 26
        rays = np.array( ((1.0,0.0,0.0), (0.0,1.0,0.0), (0.0,0.0,1.0), 
                          (1.0,1.0,0.0), (1.0,0.0,1.0), (0.0,1.0,1.0), 
                          (1.0,1.0,1.0)) )
        nrays = 7
        nradii = 200
        rvals = np.linspace(0,1*pf["cm"]/rs0,nradii)
        E0profile  = np.zeros(rvals.shape, dtype=float)
        E1profile  = np.zeros(rvals.shape, dtype=float)
        E2profile  = np.zeros(rvals.shape, dtype=float)
        E3profile  = np.zeros(rvals.shape, dtype=float)
        E4profile  = np.zeros(rvals.shape, dtype=float)
        E5profile  = np.zeros(rvals.shape, dtype=float)
        E6profile  = np.zeros(rvals.shape, dtype=float)
        E7profile  = np.zeros(rvals.shape, dtype=float)
        E8profile  = np.zeros(rvals.shape, dtype=float)
        E9profile  = np.zeros(rvals.shape, dtype=float)
        HIprofile  = np.zeros(rvals.shape, dtype=float)
        HIIprofile = np.zeros(rvals.shape, dtype=float)
        Tprofile   = np.zeros(rvals.shape, dtype=float)

        # generate 1D profiles from ray emanating out from box center to corner
        for iray in range(0,nrays):
            rvec = rays[iray,:]
            rvec = rvec / sqrt(rvec[0]**2 + rvec[1]**2 + rvec[2]**2)
            r = pf.h.ray([0.0, 0.0, 0.0], rvec)
            ptype = [('r', float), ('xHI', float), ('xHII', float), ('T', float),
                     ('E0', float), ('E1', float), ('E2', float), ('E3', float), 
                     ('E4', float), ('E5', float), ('E6', float), ('E7', float), 
                     ('E8', float), ('E9', float)]
            HIprof  = log10(r["xHI"])
            HIIprof = log10(r["xHII"])
            Tprof   = log10(r["logT"])
            E0prof = E1prof = E2prof = E3prof = E4prof = 0.0*HIprof
            E5prof = E6prof = E7prof = E8prof = E9prof = 0.0*HIprof
            if (numradfields > 0):
                E0prof = log10(r["logE0"])
            if (numradfields > 1):
                E1prof = log10(r["logE1"])
            if (numradfields > 2):
                E2prof = log10(r["logE2"])
            if (numradfields > 3):
                E3prof = log10(r["logE3"])
            if (numradfields > 4):
                E4prof = log10(r["logE4"])
            if (numradfields > 5):
                E5prof = log10(r["logE5"])
            if (numradfields > 6):
                E6prof = log10(r["logE6"])
            if (numradfields > 7):
                E7prof = log10(r["logE7"])
            if (numradfields > 8):
                E8prof = log10(r["logE8"])
            if (numradfields > 9):
                E9prof = log10(r["logE9"])
            Hradii  = r["radius"]
        
            # sort results by radius (since that isn't quite working correctly from yt)
            pdata = np.zeros(Hradii.shape, dtype=ptype);
            nrad = (Hradii.shape)[0]
            for irad in range(0,nrad):
                E0p = E1p = E2p = E3p = E4p = E5p = E6p = E7p = E8p = E9p = 0
                if (numradfields > 0):
                    E0p = E0prof[irad]
                if (numradfields > 1):
                    E1p = E1prof[irad]
                if (numradfields > 2):
                    E2p = E2prof[irad]
                if (numradfields > 3):
                    E3p = E3prof[irad]
                if (numradfields > 4):
                    E4p = E4prof[irad]
                if (numradfields > 5):
                    E5p = E5prof[irad]
                if (numradfields > 6):
                    E6p = E6prof[irad]
                if (numradfields > 7):
                    E7p = E7prof[irad]
                if (numradfields > 8):
                    E8p = E8prof[irad]
                if (numradfields > 9):
                    E9p = E9prof[irad]
                pdata[irad] = (Hradii[irad], HIprof[irad], HIIprof[irad], 
                               Tprof[irad], E0p, E1p, E2p, E3p, E4p, E5p, 
                               E6p, E7p, E8p, E9p)
            pdata = np.sort(pdata, order='r')
            
            # interpolate results into output arrays
            tmp = np.interp(rvals, pdata['r'], pdata['xHI'])
            HIprofile += tmp
            tmp = np.interp(rvals, pdata['r'], pdata['xHII'])
            HIIprofile += tmp
            tmp = np.interp(rvals, pdata['r'], pdata['T'])
            Tprofile += tmp
            if (numradfields > 0):
                tmp = np.interp(rvals, pdata['r'], pdata['E0'])
                E0profile += tmp
            if (numradfields > 1):
                tmp = np.interp(rvals, pdata['r'], pdata['E1'])
                E1profile += tmp
            if (numradfields > 2):
                tmp = np.interp(rvals, pdata['r'], pdata['E2'])
                E2profile += tmp
            if (numradfields > 3):
                tmp = np.interp(rvals, pdata['r'], pdata['E3'])
                E3profile += tmp
            if (numradfields > 4):
                tmp = np.interp(rvals, pdata['r'], pdata['E4'])
                E4profile += tmp
            if (numradfields > 5):
                tmp = np.interp(rvals, pdata['r'], pdata['E5'])
                E5profile += tmp
            if (numradfields > 6):
                tmp = np.interp(rvals, pdata['r'], pdata['E6'])
                E6profile += tmp
            if (numradfields > 7):
                tmp = np.interp(rvals, pdata['r'], pdata['E7'])
                E7profile += tmp
            if (numradfields > 8):
                tmp = np.interp(rvals, pdata['r'], pdata['E8'])
                E8profile += tmp
            if (numradfields > 9):
                tmp = np.interp(rvals, pdata['r'], pdata['E9'])
                E9profile += tmp
        HIprofile  /= nrays
        HIIprofile /= nrays
        Tprofile   /= nrays
        if (numradfields > 0):
            E0profile /= nrays
        if (numradfields > 1):
            E1profile /= nrays
        if (numradfields > 2):
            E2profile /= nrays
        if (numradfields > 3):
            E3profile /= nrays
        if (numradfields > 4):
            E4profile /= nrays
        if (numradfields > 5):
            E5profile /= nrays
        if (numradfields > 6):
            E6profile /= nrays
        if (numradfields > 7):
            E7profile /= nrays
        if (numradfields > 8):
            E8profile /= nrays
        if (numradfields > 9):
            E9profile /= nrays
        
        # chemistry profiles
        figure()
        plot(rvals,HIprofile,'b-',rvals,HIIprofile,'r--')
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(xHI), log(xHII)')
        title('HI, HII Profiles, t =' + Myr + ' Myr')
        legend( ('xHI','xHII'), 'lower right' )
        axis([ 0.0, 1.2, -7.0, 1.0 ])
        savefig('profiles_' + Myr + 'Myr.' + pictype)
        
        # Temperature profile
        figure()
        plot(rvals,Tprofile)
        grid()
        xlabel('$r/r_S$')
        ylabel('log(T) [K]')
        title('Temperature Profile, t =' + Myr + ' Myr')
        axis([ 0.0, 1.2, 3.5, 4.6 ])
        savefig('TempProfile_' + Myr + 'Myr.' + pictype)
        
        # radiation profiles
        if (numradfields > 0):
            figure()
            plot(rvals, E0profile); 
            labels = ['E0']
            if (numradfields > 1):
                plot(rvals, E1profile); 
                labels.append('E1')
            if (numradfields > 2):
                plot(rvals, E2profile); 
                labels.append('E2')
            if (numradfields > 3):
                plot(rvals, E3profile); 
                labels.append('E3')
            if (numradfields > 4):
                plot(rvals, E4profile); 
                labels.append('E4')
            if (numradfields > 5):
                plot(rvals, E5profile); 
                labels.append('E5')
            if (numradfields > 6):
                plot(rvals, E6profile); 
                labels.append('E6')
            if (numradfields > 7):
                plot(rvals, E7profile); 
                labels.append('E7')
            if (numradfields > 8):
                plot(rvals, E8profile); 
                labels.append('E8')
            if (numradfields > 9):
                plot(rvals, E9profile); 
                labels.append('E9')
            grid()
            xlabel('$r/r_S$')
            ylabel('log(E)')
            title('Radiation Profiles, t =' + Myr + ' Myr')
            legend( labels, 'lower right' )
            #axis([ 0.0, 1.2, 3.5, 4.6 ])
            savefig('EProfiles_' + Myr + 'Myr.' + pictype)


# I-front radius
figure()
plot(rdata[0],rdata[1]/rs0,'b-')
xlabel('$t/t_{rec}$')
ylabel('$r_I/r_S$')
title('Propagation of HII Region')
axis([ 0.0, 4.25, 0.0, 1.2 ])
grid()
savefig('rad_vs_time.' + pictype)

