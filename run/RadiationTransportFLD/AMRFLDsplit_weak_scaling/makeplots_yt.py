# matplotlib-based plotting script for Iliev et al. test #2
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
from yt.mods import *
from os import *

# set the total number of snapshots
te = 10

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
def _logE(field, data):
    return (data["Grey_Radiation_Energy"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE", take_log=True, function=_logE, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Energy\; Density", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Temperature (log plot)
def _logT(field, data):
    mp = 1.67262171e-24
    kb = 1.3806504e-16
    gamma = 5.0/3.0
    tmp = (data["Total_Energy"] * data["Density"] / 
           (2.0*data["Density"] - data["HI_Density"]))
    return ((gamma - 1.0)*mp/kb * tmp)
add_field("logT", take_log=True, function=_logT, 
          display_name="Temperature", units=r"\rm{K}")

#   Radius from domain center
def _radius(field, data):
    return (np.sqrt(data["x"]*data["x"] + data["y"]*data["y"] +
                    data["z"]*data["z"]))
def _convertradius(data):
    return (data.convert("cm"))
add_field("radius", take_log=False, function=_radius, 
          convert_function=_convertradius, 
          display_name="radius", units=r"\rm{cm}")



# loop over snapshots, loading values and times
for tstep in range(0,te+1):
    
    # load relevant information
    sdump = repr(tstep).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    pf = load(pfile)
    t = pf.current_time * pf["TimeUnits"]

    # generate 2D plots at certain times
    if (tstep == 0) or (tstep == 1) or (tstep == 5) or (tstep == 10) or (tstep == 20) or (tstep == 50):
        
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

        p = pc.add_slice("logE",'z')
        p.modify["title"]('radiation density, t =' + Myr + ' Myr')

        p = pc.add_slice("logT",'z')
        p.modify["title"]('temperature, t =' + Myr + ' Myr')

        pc.save(Myr + 'Myr', format=pictype)

        # rename generated files
        f1 = Myr + 'Myr_Slice_z_logE.' + pictype
        f2 = 'Econtour_' + Myr + 'Myr.' + pictype
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



