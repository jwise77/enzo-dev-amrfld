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
# tests=( \
#     RHIonization1_unigrid \
#     RHIonization2_unigrid \
#     RHIonization3_unigrid \
#     RHIonization5_unigrid \
#     RHIonization6_unigrid \
#     RHIonization7_unigrid \
#     CosmoIonization_q5z4_unigrid \
#     CosmoIonization_q05z4_unigrid \
#     CosmoIonization_q5z10_unigrid \
#     CosmoIonization_q05z10_unigrid \
#     RHIonization1_AMRstatic \
#     RHIonization2_AMRstatic \
#     RHIonization3_AMRstatic \
#     RHIonization5_AMRstatic \
#     RHIonization6_AMRstatic \
#     RHIonization7_AMRstatic \
#     CosmoIonization_q5z4_AMRstatic \
#     CosmoIonization_q05z4_AMRstatic \
#     CosmoIonization_q5z10_AMRstatic \
#     CosmoIonization_q05z10_AMRstatic \
#     RHIonization1_AMRdynamic \
#     RHIonization2_AMRdynamic \
#     RHIonization3_AMRdynamic \
#     RHIonization5_AMRdynamic \
#     RHIonization6_AMRdynamic \
#     RHIonization7_AMRdynamic \
#     CosmoIonization_q5z4_AMRdynamic \
#     CosmoIonization_q05z4_AMRdynamic \
#     CosmoIonization_q5z10_AMRdynamic \
#     CosmoIonization_q05z10_AMRdynamic \
# )
tests=( \
    RHIonization1_unigrid \
    RHIonization2_unigrid \
    RHIonization3_unigrid \
    RHIonization5_unigrid \
    RHIonization6_unigrid \
    RHIonization7_unigrid \
    CosmoIonization_q5z4_unigrid \
    CosmoIonization_q05z4_unigrid \
    CosmoIonization_q5z10_unigrid \
    CosmoIonization_q05z10_unigrid \
    RHIonization1_AMRstatic \
    RHIonization2_AMRstatic \
    CosmoIonization_q5z10_AMRstatic \
    CosmoIonization_q05z10_AMRstatic \
    RHIonization1_AMRdynamic \
    RHIonization2_AMRdynamic \
    CosmoIonization_q5z10_AMRdynamic \
    CosmoIonization_q05z10_AMRdynamic \
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
