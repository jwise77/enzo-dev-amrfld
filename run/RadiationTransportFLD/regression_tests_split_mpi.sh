#!/bin/sh
yt_activate

echo "    "
echo "Date:"
date
echo "    "

echo "    "
echo "Running Iliev et al. Split Test 1"
cd RHIonization1_sp
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
echo "Running Iliev et al. Split Test 2"
cd RHIonization2_sp
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
echo "Running Shapiro & Giroux q0=0.5 z0=4 Split Test"
cd CosmoIonization_q5z4_sp
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
echo "Running Shapiro & Giroux q0=0.05 z0=4 Split Test"
cd CosmoIonization_q05z4_sp
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
echo "Running Shapiro & Giroux q0=0.5 z0=10 Split Test"
cd CosmoIonization_q5z10_sp
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
echo "Running Shapiro & Giroux q0=0.05 z0=10 Split Test"
cd CosmoIonization_q05z10_sp
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
