#!/bin/bash

# array of tests to clean up
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

echo "  Cleaning up regression test results:"

for i in "${tests[@]}"
do
  pushd $i > /dev/null
  ln -fs ../cleanup.sh
  ./cleanup.sh
  popd > /dev/null
done

echo "  Finished!"
