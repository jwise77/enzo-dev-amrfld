#!/bin/bash

# array of tests to clean up
tests=( \
    TurnerStoneEquil1_sp \
    TurnerStoneEquil2_sp \
    RadiationStreamX0_sp \
    RadiatingShockLab_sp \
    RadiatingShockLab1D_sp \
    RHIonization1_sp \
    RHIonization2_sp \
    CosmoIonization_q5z4_sp \
    CosmoIonization_q05z4_sp \
    CosmoIonization_q5z10_sp \
    CosmoIonization_q05z10_sp\
)

echo "  Cleaning up regression test results:"

for i in "${tests[@]}"
do
  pushd $i > /dev/null
  ln -fs ../cleanup.sh
  ./cleanup.sh
  popd > /dev/null
done

echo "  Finished!"
