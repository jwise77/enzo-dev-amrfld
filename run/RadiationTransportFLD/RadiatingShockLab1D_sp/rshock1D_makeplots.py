# matplotlib-based plotting script for radiating shock tests
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import numpy as np

# set the graphics output type
pictype = '.png'

# store exact solution values
xi = [8.000838354813387e-02, 8.284041833623391e-02, 8.394370939772287e-02, 
      8.501057335585425e-02, 8.604153195060078e-02, 8.700271757035795e-02,
      8.751663034854217e-02, 8.800393249934589e-02, 8.850802709734139e-02,
      8.900280804043151e-02, 8.920192262961063e-02, 8.940133858887539e-02,
      8.960058691879980e-02, 8.980358834109224e-02, 8.990275612051976e-02,
      9.000000000000000e-02, 9.000000000000000e-02, 9.050181717389578e-02,
      9.100363434779157e-02, 9.200726869558316e-02, 9.301090304337475e-02,
      9.401453739116633e-02, 9.501817173895791e-02, 9.602180608674950e-02,
      9.702544043454109e-02, 9.802907478233266e-02, 9.903270913012424e-02,
      1.000363434779158e-01]
x_ex = [8.700271757035795e-02,  8.721289728521497e-02,  8.741110537269545e-02,  
        8.761411040468112e-02,  8.780312958806377e-02,  8.800393249934589e-02,  
        8.820939716865782e-02,  8.840635840794454e-02,  8.860325618672420e-02,  
        8.880510806722709e-02,  8.900280804043151e-02,  8.920192262961063e-02,  
        8.940133858887539e-02,  8.960058691879980e-02,  8.980358834109224e-02,  
        9.000000000000000e-02,  9.000000000000000e-02,  9.022302985506479e-02,  
        9.044605971012959e-02,  9.061333210142819e-02,  9.083636195649299e-02,  
        9.100363434779157e-02,  9.122666420285638e-02,  9.144969405792117e-02,  
        9.161696644921977e-02,  9.183999630428456e-02,  9.200726869558316e-02,  
        9.223029855064796e-02,  9.245332840571276e-02,  9.262060079701136e-02,  
        9.284363065207615e-02,  9.301090304337475e-02]
vi = [-6.283244113400758e+02, -7.661597371058387e+03, -2.029904324144572e+04, 
      -5.208449113051507e+04, -1.295032018305447e+05, -3.029091395730755e+05, 
      -4.773022575921432e+05, -7.348903072922676e+05, -1.149074215131261e+06, 
      -1.783124926830709e+06, -2.128526220776094e+06, -2.541962412231995e+06, 
      -3.036037600331375e+06, -3.640022928812999e+06, -3.978601944258832e+06, 
      -4.342581582537412e+06, -1.513357745471116e+07, -1.757990948196923e+07, 
      -1.859085440588360e+07, -1.927046311292141e+07, -1.942626081456771e+07, 
      -1.946302254503200e+07, -1.947175433815990e+07, -1.947383159167752e+07, 
      -1.947432594221196e+07, -1.947444360160807e+07, -1.947447160488609e+07, 
      -1.947447826937688e+07] 
v_ex = [-3.029091395730755e+05, -3.648035828119376e+05, -4.347404425190039e+05, 
        -5.203212093986827e+05, -6.151265166072322e+05, -7.348903072922676e+05, 
        -8.816856670425900e+05, -1.049961367295480e+06, -1.250410737341433e+06, 
        -1.495860346651081e+06, -1.783124926830709e+06, -2.128526220776094e+06, 
        -2.541962412231995e+06, -3.036037600331375e+06, -3.640022928812999e+06, 
        -4.342581582537412e+06, -1.513357745471122e+07, -1.651366078418055e+07, 
        -1.740674580823335e+07, -1.788069204316191e+07, -1.833938124041248e+07, 
        -1.859085440588360e+07, -1.883926562006382e+07, -1.901639407487639e+07, 
        -1.911545895091638e+07, -1.921467104252037e+07, -1.927046311292141e+07, 
        -1.932655628985928e+07, -1.936715248529808e+07, -1.939007608727717e+07, 
        -1.941319186813242e+07, -1.942626081456771e+07]
di = [1.000018151049645e+00,  1.000221373374236e+00,  1.000586732624816e+00,  
      1.001506857838764e+00,  1.003755071296794e+00,  1.008827530598369e+00,  
      1.013980836750794e+00,  1.021689593771764e+00,  1.034333542580215e+00,  
      1.054307392556115e+00,  1.065516246670119e+00,  1.079250362294035e+00,  
      1.096134900308262e+00,  1.117507169805125e+00,  1.129856504895149e+00,  
      1.143440384852932e+00,  1.776739179147769e+00,  2.031857908173327e+00,  
      2.160029080434648e+00,  2.255683864151188e+00,  2.278818299637701e+00,  
      2.284346420618191e+00,  2.285663426185862e+00,  2.285976959403660e+00,  
      2.286051587577088e+00,  2.286069350399400e+00,  2.286073578043371e+00,  
      2.286074584181251e+00]
rho_ex = [1.008827530598369e+00,  1.010650500550850e+00,  1.012718287092141e+00,  
          1.015260149830882e+00,  1.018090932238415e+00,  1.021689593771764e+00,  
          1.026135352363984e+00,  1.031279487452762e+00,  1.037474876629775e+00,  
          1.045163241572978e+00,  1.054307392556115e+00,  1.065516246670119e+00,  
          1.079250362294035e+00,  1.096134900308262e+00,  1.117507169805125e+00,  
          1.143440384852932e+00,  1.776739179147775e+00,  1.912186079245962e+00,  
          2.011414111823347e+00,  2.068374059474956e+00,  2.126658865315985e+00,  
          2.160029080434648e+00,  2.194037366425768e+00,  2.218948268817144e+00,  
          2.233128736909006e+00,  2.247513096313387e+00,  2.255683864151188e+00,  
          2.263958839834421e+00,  2.269985641510433e+00,  2.273403014083470e+00,  
          2.276859471059020e+00,  2.278818299637701e+00]
egasi = [2.696542574809257e+14,  2.699642853963563e+14,  2.705209739927604e+14,  
         2.719190314454886e+14,  2.753115710530059e+14,  2.828452022975568e+14,  
         2.903308914170840e+14,  3.012210206428349e+14,  3.183149130955234e+14,  
         3.434875880463936e+14,  3.566938035092925e+14,  3.720313963864993e+14,  
         3.896889383682254e+14,  4.102811291577671e+14,  4.213465812917462e+14,  
         4.328591502673416e+14,  5.939404791217856e+14,  5.818598091904346e+14,  
         5.716248784856595e+14,  5.630212190809336e+14,  5.608536552518507e+14,  
         5.603315842279227e+14,  5.602069840557834e+14,  5.601773085716539e+14,  
         5.601702444112664e+14,  5.601685629763826e+14,  5.601681627838805e+14,  
         5.601680675418579e+14]
egas_ex = [2.828452022975568e+14,  2.855124005862494e+14,  2.885123612229592e+14,  
           2.921634402273271e+14,  2.961824494793456e+14,  3.012210206428349e+14,  
           3.073382340808936e+14,  3.142711681779321e+14,  3.224189609804636e+14,  
           3.322319934795784e+14,  3.434875880463936e+14,  3.566938035092925e+14,  
           3.720313963864993e+14,  3.896889383682254e+14,  4.102811291577671e+14,  
           4.328591502673416e+14,  5.939404791217856e+14,  5.893325698206294e+14,  
           5.833053901915996e+14,  5.791350352139968e+14,  5.744572827543389e+14,  
           5.716248784856595e+14,  5.686407142512102e+14,  5.663998192836241e+14,  
           5.651054887522295e+14,  5.637797312028029e+14,  5.630212190809336e+14,  
           5.622492005500327e+14,  5.616845821803565e+14,  5.613635732059496e+14,  
           5.610382769336721e+14,  5.608536552518507e+14]
eradi = [3.001398564464470e+10,  3.063487268084298e+10,  3.174986183701081e+10,  
         3.455070778699154e+10,  4.135141619234581e+10,  5.647543730315511e+10,  
         7.153548121759041e+10,  9.350764194538100e+10,  1.281653507008099e+11,  
         1.796401053595543e+11,  2.068842100993477e+11,  2.387564956298731e+11,  
         2.757822927153420e+11,  3.194398115048245e+11,  3.431167802708561e+11,  
         3.679068119138071e+11,  3.679068119138071e+11,  4.640395830628364e+11,  
         5.120246018991776e+11,  5.471361045577606e+11,  5.555192773137209e+11,  
         5.575158340548044e+11,  5.579911054684722e+11,  5.581042292712521e+11,  
         5.581311541205513e+11,  5.581375626357207e+11,  5.581390878830764e+11,  
         5.581394508595256e+11]
erad_ex = [5.647543730315511e+10,  6.183760018908191e+10,  6.787378035190774e+10,  
           7.522748823046490e+10,  8.333202229695720e+10,  9.350764194538100e+10,  
           1.058853091159983e+11,  1.199467038841325e+11,  1.365200077395836e+11,  
           1.565541472957868e+11,  1.796401053595543e+11,  2.068842100993477e+11,  
           2.387564956298731e+11,  2.757822927153420e+11,  3.194398115048245e+11,  
           3.679068119138071e+11,  3.679068119138071e+11,  4.187194729125132e+11,  
           4.563151055718904e+11,  4.777975688812859e+11,  4.996212470928763e+11,  
           5.120246018991776e+11,  5.245867919848723e+11,  5.337345553458995e+11,  
           5.389205130508717e+11,  5.441647032246688e+11,  5.471361045577606e+11,  
           5.501398106802782e+11,  5.523238825166902e+11,  5.535609654633843e+11,  
           5.548111981993607e+11,  5.555192773137209e+11]

# set some constants
kb_div_everg = 1.3806505/1.60217646 * 1.0e-4
CvRL         = 2.218056e12
Cv           = CvRL * kb_div_everg
gamma        = 5.0/3.0

# problem information
Length   =  0.1
V_inflow =  1.94744804e-2
tstop    =  1.73325

# shock velocity
machno = 2.0
Vshock = machno * sqrt( gamma*(gamma-1)*CvRL*121.6 ) * 1e-9

# define some helpful functions
def get_params(file):
    """Returns dUnit, tUnit, lUnit, vUnit from a given parameter file"""
    import shlex
    f = open(file)
    for line in f:
        text = shlex.split(line)
        if ("DensityUnits" in text):
            dUnit = float(text[len(text)-1])
        elif ("TimeUnits" in text):
            tUnit = float(text[len(text)-1])
        elif ("LengthUnits" in text):
            lUnit = float(text[len(text)-1])
    vUnit = lUnit/tUnit
    return [dUnit, tUnit, lUnit, vUnit]

def load_vals(tdump):
    """Returns Eg, etot, ke from a given data dump"""
    import h5py
    import numpy as np
    sdump = repr(tdump).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    hfile = pfile + '.cpu0000'
    dUnit, tUnit, lUnit, vUnit = get_params(pfile)
    f = h5py.File(hfile,'r')
    vx = f.get('/Grid00000001/x-velocity')
    vy = f.get('/Grid00000001/y-velocity')
    vz = f.get('/Grid00000001/z-velocity')
    Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
    etot = f.get('/Grid00000001/TotalEnergy')
    vx = np.multiply(vx,vUnit)
    vy = np.multiply(vy,vUnit)
    vz = np.multiply(vz,vUnit)
    Eg = np.multiply(Eg,dUnit*vUnit*vUnit)
    etot = np.multiply(etot,vUnit*vUnit)
    ke = np.multiply(np.multiply(vx,vx) + np.multiply(vy,vy) + np.multiply(vz,vz),0.5)
    return [Eg, etot, ke]

# load relevant data snapshots
r_e, t_e, k_e = load_vals(100)

# generate figure with solution data (Figure 6 from paper)
#    reference solution data
x_sol = np.subtract(xi,0.09)
tgas_sol = np.multiply(egasi, kb_div_everg/121.6/Cv)
trad_sol = np.multiply(kb_div_everg/121.6,np.sqrt(np.sqrt(np.multiply(eradi,1.0e15/7.56))))
#    internal energies
i_e = t_e - k_e
#    compute gas temperature from internal energy, in units of 121.6 eV
Tgas = i_e * ( kb_div_everg / 121.6 / Cv )
Trad = ( kb_div_everg / 121.6 ) * sqrt( sqrt( 1.0e15 * r_e / 7.56 ) )
#    generate x-coordinates of Enzo data
len = t_e.size
dx  = Length/len
x = linspace(dx/2,Length-dx/2,len)
#    shock location
Xshock = Length - Vshock * tstop
#    transform x into frame comoving with the INFLOWING material.
x_cmi = np.subtract(x, V_inflow*tstop)
#    transform x_cmi into frame comoving with the SHOCK
x_cms = np.subtract(x_cmi, Xshock)


# plot computed solutions overlaid with analytical solution data
plot(x_cms, Tgas, 'b-', x_cms, Trad, 'r--')
plot(x_sol, tgas_sol, 'bo', x_sol, trad_sol, 'rs')
xlabel('Distance From Shock (cm)'); 
ylabel('T');
title('Analytical and Computed Solutions')
axis([-0.01, 0.01, 1.0, 2.25])
grid()
savefig('rshock_sol' + pictype)
