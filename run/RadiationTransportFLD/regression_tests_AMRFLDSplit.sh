#!/bin/bash
source /usr/local/yt_dev/bin/activate

echo "    "
echo "Date:"
date
echo "    "

echo "    "
echo "Running Ionization static AMRFLDSplit Test 1"
cd RHIonization1_AMRstatic
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Ionization static AMRFLDSplit Test 2"
cd RHIonization2_AMRstatic
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Ionization dynamic AMRFLDSplit Test 1"
cd RHIonization1_AMRdynamic
ln -fs ../../../src/enzo/enzo.exe enzo
mpiexec -n 4 ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
./plotgrids.sh
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Ionization dynamic AMRFLDSplit Test 2"
cd RHIonization2_AMRdynamic
ln -fs ../../../src/enzo/enzo.exe enzo
mpiexec -n 4 ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
./plotgrids.sh
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Cosmological Ionization q0=0.5 z0=4 dynamic AMRFLDSplit Test"
cd CosmoIonization_q5z4_AMRdynamic
ln -fs ../../../src/enzo/enzo.exe enzo
mpiexec -n 4 ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
./plotgrids.sh
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Cosmological Ionization q0=0.05 z0=4 dynamic AMRFLDSplit Test"
cd CosmoIonization_q05z4_AMRdynamic
ln -fs ../../../src/enzo/enzo.exe enzo
mpiexec -n 4 ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
./plotgrids.sh
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Cosmological Ionization q0=0.5 z0=10 dynamic AMRFLDSplit Test"
cd CosmoIonization_q5z10_AMRdynamic
ln -fs ../../../src/enzo/enzo.exe enzo
mpiexec -n 4 ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
./plotgrids.sh
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Cosmological Ionization q0=0.05 z0=10 dynamic AMRFLDSplit Test"
cd CosmoIonization_q05z10_AMRdynamic
ln -fs ../../../src/enzo/enzo.exe enzo
mpiexec -n 4 ./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep "Successful run" output.txt
./plotgrids.sh
python ./*makeplots_yt.py &> /dev/null
python ./*check_yt.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Finished!!"
