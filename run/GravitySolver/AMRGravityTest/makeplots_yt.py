import os.path
import numpy as na
import matplotlib.pyplot as plt

# tolerance to pass test
tol = 0.12

# load data, based on serial vs parallel test
Data = []
if (os.path.isfile("TestGravityCheckResults.out")):
    Data = na.loadtxt("TestGravityCheckResults.out")
else:
    for iproc in range(0,10000):
        sproc = repr(iproc).zfill(4)
        fname = 'TestGravityCheckResults.out' + sproc
        if (os.path.isfile(fname)):
            nprocs = iproc
            DataP = na.loadtxt(fname)
            if (iproc == 0):
                Data = DataP
            else:
                Data = na.append(Data, DataP, axis=0)
        else:
            break

radius = Data[:,0]
ForceTangentialComputed = Data[:,1]
ForceRadialComputed = Data[:,2]
ForceRadialTrue = Data[:,3]

Error = (ForceRadialComputed-ForceRadialTrue)/ForceRadialTrue
indices = na.where((radius > 1.0) & (radius < 8.0))
rmsError = na.std(Error[indices])
print "rms error = "+str(rmsError)
if (rmsError < tol):
    print "pass"
else:
    print "fail"
    

# Plot the computed radial force again the r-2 profile 
#  (which should be equal except for outer part where the periodic 
#   boundary conditions are important, and inner part, where grid
#   softening is important).
# Also plot tangential component of force, which should be zero
#  (or at least small compared to radial component).

plt.loglog(radius, ForceRadialComputed, label='Frad', ls='None', marker='+')
plt.loglog(radius, ForceRadialTrue, label='Frad_true', ls='None', marker='.')
plt.loglog(radius, ForceTangentialComputed, label='Ftang', ls='None', marker='o')
plt.xlabel('r (cells)')
plt.ylabel('Force')
plt.axis([0.1,20.0,1.0e-7,1.0e0])
plt.legend()
plt.savefig('GravityTest.png')
