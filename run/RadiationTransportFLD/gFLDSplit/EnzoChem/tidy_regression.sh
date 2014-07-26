#!/bin/bash

# array of tests to clean up
tests=( \
    RHIonization1_enzochem \
    RHIonization2_enzochem \
    CosmoIonization_q5z4_enzochem \
    CosmoIonization_q05z4_enzochem \
    CosmoIonization_q5z10_enzochem \
    CosmoIonization_q05z10_enzochem \
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
