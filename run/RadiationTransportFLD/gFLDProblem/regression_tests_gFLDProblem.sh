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
echo "Running Turner & Stone 1 Test"
cd TurnerStoneEquil1
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Turner & Stone 2 Test"
cd TurnerStoneEquil2
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiation Stream X0 Test"
cd RadiationStreamX0
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "


echo "    "
echo "Running Radiation Stream 1D Test"
cd RadiationStream1D
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiating Shock Lab Test"
cd RadiatingShockLab
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots.py &> /dev/null
python ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiating Shock Lab 1D Test"
cd RadiatingShockLab1D
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
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
echo "Running Ionization Test 1"
cd RHIonization1
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Ionization Test 2"
cd RHIonization2
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Cosmological Ionization q0=0.5 z0=4 Test"
cd CosmoIonization_q5z4
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Cosmological Ionization q0=0.05 z0=4 Test"
cd CosmoIonization_q05z4
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Cosmological Ionization q0=0.5 z0=10 Test"
cd CosmoIonization_q5z10
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Cosmological Ionization q0=0.05 z0=10 Test"
cd CosmoIonization_q05z10
ln -fs ../../../src/enzo/enzo.exe enzo
$MPIEXEC -n $NP ./enzo -d *.enzo 1> output.txt 2> error.txt
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
