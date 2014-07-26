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
    RHIonization1_enzochem \
    RHIonization2_enzochem \
    CosmoIonization_q5z4_enzochem \
    CosmoIonization_q05z4_enzochem \
    CosmoIonization_q5z10_enzochem \
    CosmoIonization_q05z10_enzochem \
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
