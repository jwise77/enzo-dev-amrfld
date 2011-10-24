#!/bin/bash
yt_activate

echo "  "
echo "cleaning up old tests"
cd FSRadWave/nx128; \rm -rf DD* En* Ev* IO* Out* *.out Run* *.txt perf* *.pdf; cd ../..
cd FSRadPoint/nx128; \rm -rf DD* En* Ev* IO* Out* *.out Run* *.txt perf* *.pdf; cd ../..
cd FSRadWave/nx256; \rm -rf DD* En* Ev* IO* Out* *.out Run* *.txt perf* *.pdf; cd ../..
cd FSRadPoint/nx256; \rm -rf DD* En* Ev* IO* Out* *.out Run* *.txt perf* *.pdf; cd ../..


echo "running nx=128 test"
cd FSRadWave/nx128; mpiexec -n 8 ./enzo -d *.enzo &> output.txt; 
python *.py &> /dev/null; cd ../..
cd FSRadPoint/nx128; mpiexec -n 8 ./enzo -d *.enzo &> output.txt; 
python *.py &> /dev/null; cd ../..

echo "running nx=256 test"
cd FSRadWave/nx256; mpiexec -n 8 ./enzo -d *.enzo &> output.txt; 
python *.py &> /dev/null; cd ../..
cd FSRadPoint/nx256; mpiexec -n 8 ./enzo -d *.enzo &> output.txt; 
python *.py &> /dev/null; cd ../..

echo "  "
echo "finished!"
echo "  "
