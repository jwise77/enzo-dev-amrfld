#!/bin/bash

# determine mpiexec name and number of procs
nargs=$#
if [ "$nargs" -lt 2 ]
then
    echo " "
    echo "Script utilizes command-line arguments:"
    echo "   <mpiexec> <numprocs>"
    echo "where"
    echo "   <mpiexec> is the name of the mpiexec executable to use"
    echo "   <numprocs> is the number of processors to use"
    echo " "
    echo "No arguments given, so using defaults:  \"mpiexec_intel 1\""
    echo " "
    MPIEXEC=mpiexec_intel
    NP=1
else
    MPIEXEC=$1
    NP=$2
fi

# now set up yt for remaining tests
source /usr/local/yt_dev/bin/activate

# array of tests to run
tests=( \
    TurnerStoneEquil1 \
    TurnerStoneEquil2 \
    RadiationStreamX0 \
    RadiationStream1D \
    RadiatingShockLab \
    RadiatingShockLab1D \
    RHIonization1 \
    RHIonization2 \
    CosmoIonization_q5z4 \
    CosmoIonization_q05z4 \
    CosmoIonization_q5z10 \
    CosmoIonization_q05z10 \
)

echo "    "
echo "Date:"
date
echo "Running tests with $NP processes:"
echo "    "

# run the tests
for i in "${tests[@]}"
do
  echo "Running test: $i"
  pushd $i > /dev/null
  ln -fs ../runtest.sh
  ./runtest.sh $MPIEXEC $NP
  popd > /dev/null
done

echo "Finished!!"
