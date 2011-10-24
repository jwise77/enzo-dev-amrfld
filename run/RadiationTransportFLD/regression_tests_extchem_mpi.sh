#!/bin/sh
source /usr/local/yt_dev/bin/activate

echo "    "
echo "Date:"
date
echo "    "

echo "    "
echo "Running Iliev et al. extchem Test 1"
cd RHIonization1_extchem
ln -fs ../../../src/enzo/enzo.exe enzo
mpiexec -n 4 ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Iliev et al. extchem Test 2"
cd RHIonization2_extchem
ln -fs ../../../src/enzo/enzo.exe enzo
mpiexec -n 4 ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Shapiro & Giroux q0=0.5 z0=4 extchem Test"
cd CosmoIonization_q5z4_extchem
ln -fs ../../../src/enzo/enzo.exe enzo
mpiexec -n 4 ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Shapiro & Giroux q0=0.05 z0=4 extchem Test"
cd CosmoIonization_q05z4_extchem
ln -fs ../../../src/enzo/enzo.exe enzo
mpiexec -n 4 ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Shapiro & Giroux q0=0.5 z0=10 extchem Test"
cd CosmoIonization_q5z10_extchem
ln -fs ../../../src/enzo/enzo.exe enzo
mpiexec -n 4 ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Shapiro & Giroux q0=0.05 z0=10 extchem Test"
cd CosmoIonization_q05z10_extchem
ln -fs ../../../src/enzo/enzo.exe enzo
mpiexec -n 4 ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Finished!!"
