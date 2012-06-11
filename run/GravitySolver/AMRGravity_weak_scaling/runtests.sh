#!/bin/bash

cd p01
mpiexec_intel -n 1 ./enzo.exe -d GravityTest.enzo &> output.txt
export FIELD=Grav_Potential
plotgrids.sh
mv frames potential_frames
export FIELD=Dark_Matter_Density
plotgrids.sh
mv frames dm_density_frames
cd ..

cd p02
mpiexec_intel -n 2 ./enzo.exe -d GravityTest.enzo &> output.txt
export FIELD=Grav_Potential
plotgrids.sh
mv frames potential_frames
export FIELD=Dark_Matter_Density
plotgrids.sh
mv frames dm_density_frames
cd ..

cd p04
mpiexec_intel -n 4 ./enzo.exe -d GravityTest.enzo &> output.txt
export FIELD=Grav_Potential
plotgrids.sh
mv frames potential_frames
export FIELD=Dark_Matter_Density
plotgrids.sh
mv frames dm_density_frames
cd ..

cd p06
mpiexec_intel -n 6 ./enzo.exe -d GravityTest.enzo &> output.txt
export FIELD=Grav_Potential
plotgrids.sh
mv frames potential_frames
export FIELD=Dark_Matter_Density
plotgrids.sh
mv frames dm_density_frames
cd ..

cd p08
mpiexec_intel -n 8 ./enzo.exe -d GravityTest.enzo &> output.txt
export FIELD=Grav_Potential
plotgrids.sh
mv frames potential_frames
export FIELD=Dark_Matter_Density
plotgrids.sh
mv frames dm_density_frames
cd ..

# cd p10
# mpiexec_intel -n 10 ./enzo.exe -d GravityTest.enzo &> output.txt
# export FIELD=Grav_Potential
# plotgrids.sh
# mv frames potential_frames
# export FIELD=Dark_Matter_Density
# plotgrids.sh
# mv frames dm_density_frames
# cd ..

# cd p12
# mpiexec_intel -n 12 ./enzo.exe -d GravityTest.enzo &> output.txt
# export FIELD=Grav_Potential
# plotgrids.sh
# mv frames potential_frames
# export FIELD=Dark_Matter_Density
# plotgrids.sh
# mv frames dm_density_frames
# cd ..

# cd p14
# mpiexec_intel -n 14 ./enzo.exe -d GravityTest.enzo &> output.txt
# export FIELD=Grav_Potential
# plotgrids.sh
# mv frames potential_frames
# export FIELD=Dark_Matter_Density
# plotgrids.sh
# mv frames dm_density_frames
# cd ..

# cd p16
# mpiexec_intel -n 16 ./enzo.exe -d GravityTest.enzo &> output.txt
# export FIELD=Grav_Potential
# plotgrids.sh
# mv frames potential_frames
# export FIELD=Dark_Matter_Density
# plotgrids.sh
# mv frames dm_density_frames
# cd ..

