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

source /usr/local/yt_dev/bin/activate

echo "    "
echo "Date:"
date
echo "Running tests with $NP processes:"
echo "    "

echo "    "
echo "Running Ionization enzochem Test 1"
cd RHIonization1_enzochem
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
echo "Running Ionization enzochem Test 2"
cd RHIonization2_enzochem
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
echo "Running Cosmological Ionization q0=0.5 z0=4 enzochem Test"
cd CosmoIonization_q5z4_enzochem
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
echo "Running Cosmological Ionization q0=0.05 z0=4 enzochem Test"
cd CosmoIonization_q05z4_enzochem
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
echo "Running Cosmological Ionization q0=0.5 z0=10 enzochem Test"
cd CosmoIonization_q5z10_enzochem
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
echo "Running Cosmological Ionization q0=0.05 z0=10 enzochem Test"
cd CosmoIonization_q05z10_enzochem
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
