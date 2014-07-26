#!/bin/bash

# array of tests to clean up
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

echo "  Cleaning up regression test results:"

for i in "${tests[@]}"
do
  pushd $i > /dev/null
  ln -fs ../cleanup.sh
  ./cleanup.sh
  popd > /dev/null
done

echo "  Finished!"
