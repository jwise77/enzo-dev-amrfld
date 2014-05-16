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
    echo "No arguments given, so using defaults:  \"mpiexec 1\""
    echo " "
    MPIEXEC=mpiexec
    NP=1
else
    MPIEXEC=$1
    NP=$2
fi

echo "    "
echo "Date:"
date
echo "Running tests with $NP processes:"
echo "    "

echo "    "
echo "Running Turner & Stone 1 Split Test"
cd TurnerStoneEquil1_sp
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots.py 1> /dev/null 2> /dev/null
python ./*check.py 1> PASS_FAIL.txt 2> /dev/null
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Turner & Stone 2 Split Test"
cd TurnerStoneEquil2_sp
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
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
$MPIEXEC -n $NP ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "


# echo "    "
# echo "Running Radiation Stream 1D Split Test"
# cd RadiationStream1D_sp
# ln -fs ../../../src/enzo/enzo.exe enzo
# $MPIEXEC -n $NP ./enzo -d *.enzo &> output.txt 
# grep Wallclock output.txt
# grep "Successful run" output.txt
# python ./*makeplots.py &> /dev/null
# python ./*check.py &> PASS_FAIL.txt
# echo "error checking result:"
# cat PASS_FAIL.txt
# cd ../
# echo "    "

echo "    "
echo "Running Radiating Shock Lab Split Test"
cd RadiatingShockLab_sp
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
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
$MPIEXEC -n $NP ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "


# now set up yt for remaining tests
source /usr/local/yt_dev/bin/activate


echo "    "
echo "Running Ionization Split Test 1"
cd RHIonization1_sp
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Ionization Split Test 2"
cd RHIonization2_sp
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Cosmological Ionization q0=0.5 z0=4 Split Test"
cd CosmoIonization_q5z4_sp
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Cosmological Ionization q0=0.05 z0=4 Split Test"
cd CosmoIonization_q05z4_sp
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Cosmological Ionization q0=0.5 z0=10 Split Test"
cd CosmoIonization_q5z10_sp
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Cosmological Ionization q0=0.05 z0=10 Split Test"
cd CosmoIonization_q05z10_sp
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Finished!!"
