#!/bin/sh

echo "    "
echo "Date:"
date
echo "    "

echo "    "
echo "Running Turner & Stone 1 Split Test"
cd TurnerStoneEquil1_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Turner & Stone 2 Split Test"
cd TurnerStoneEquil2_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiation Stream X0 Split Test"
cd RadiationStreamX0_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiation Stream Y1 Split Test"
cd RadiationStreamY1_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiation Stream Z0 Split Test"
cd RadiationStreamZ0_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiation Stream 1D Split Test"
cd RadiationStream1D_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiating Shock Lab Split Test"
cd RadiatingShockLab_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiating Shock Lab 1D Split Test"
cd RadiatingShockLab1D_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Iliev et al. Split Test 1"
cd RHIonization1_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Iliev et al. Split Test 2"
cd RHIonization2_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Shapiro & Giroux q0=0.5 z0=4 Split Test"
cd CosmoIonization_q5z4_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Shapiro & Giroux q0=0.05 z0=4 Split Test"
cd CosmoIonization_q05z4_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Shapiro & Giroux q0=0.5 z0=10 Split Test"
cd CosmoIonization_q5z10_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Shapiro & Giroux q0=0.05 z0=10 Split Test"
cd CosmoIonization_q05z10_sp
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Finished!!"
