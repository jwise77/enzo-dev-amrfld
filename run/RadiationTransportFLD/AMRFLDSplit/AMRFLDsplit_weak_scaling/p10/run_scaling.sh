#!/bin/bash

NUMPROCS=10


sed '/StopTime/ d' RHIonization2_scaling.enzo > tmp.enzo
echo "StopTime = 2.0e-8" >> tmp.enzo

mpiexec_intel -n $NUMPROCS ./enzo.exe -d tmp.enzo &> output_init.txt

sed 's:StaticHierarchy     = 0:StaticHierarchy     = 1:' DD0001/data0001 > tmp
grep StopTime RHIonization2_scaling.enzo > stoptime.txt
sed '/StopTime/ d' tmp > tmp.enzo
\rm tmp
cat stoptime.txt >> tmp.enzo
\rm stoptime.txt
mv tmp.enzo DD0001/data0001

mpiexec_intel -n $NUMPROCS ./enzo.exe -d -r DD0001/data0001 &> output.txt

./plotgrids.sh
grep Wallclock output.txt > stats.txt
grep cumulative output.txt | tail -n 1 >> stats.txt
grep CommunicationTime output.txt >> stats.txt

