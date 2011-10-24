#!/bin/bash
source /usr/local/yt_dev/bin/activate

echo "  "
echo "cleaning up old tests"
cd FSRadWave/nx16; \rm -rf DD* En* Ev* IO* Out* *.out Run* *.txt perf* *.pdf; cd ../..
cd FSRadPoint/nx16; \rm -rf DD* En* Ev* IO* Out* *.out Run* *.txt perf* *.pdf; cd ../..

cd FSRadWave/nx32; \rm -rf DD* En* Ev* IO* Out* *.out Run* *.txt perf* *.pdf; cd ../..
cd FSRadPoint/nx32; \rm -rf DD* En* Ev* IO* Out* *.out Run* *.txt perf* *.pdf; cd ../..

cd FSRadWave/nx64; \rm -rf DD* En* Ev* IO* Out* *.out Run* *.txt perf* *.pdf; cd ../..
cd FSRadPoint/nx64; \rm -rf DD* En* Ev* IO* Out* *.out Run* *.txt perf* *.pdf; cd ../..

# cd FSRadWave/nx128; \rm -rf DD* En* Ev* IO* Out* *.out Run* *.txt perf* *.pdf; cd ../..
# cd FSRadPoint/nx128; \rm -rf DD* En* Ev* IO* Out* *.out Run* *.txt perf* *.pdf; cd ../..


echo "  "
echo "running nx=16 tests"
cd FSRadWave/nx16; ./enzo -d *.enzo &> output.txt; 
python *.py &> /dev/null; cd ../..
cd FSRadPoint/nx16; ./enzo -d *.enzo &> output.txt; 
python *.py &> /dev/null; cd ../..

echo "running nx=32 tests"
cd FSRadWave/nx32; ./enzo -d *.enzo &> output.txt; 
python *.py &> /dev/null; cd ../..
cd FSRadPoint/nx32; ./enzo -d *.enzo &> output.txt; 
python *.py &> /dev/null; cd ../..

echo "running nx=64 tests"
cd FSRadWave/nx64; ./enzo -d *.enzo &> output.txt; 
python *.py &> /dev/null; cd ../..
cd FSRadPoint/nx64; ./enzo -d *.enzo &> output.txt; 
python *.py &> /dev/null; cd ../..

# echo "running nx=128 tests"
# cd FSRadWave/nx128; ./enzo -d *.enzo &> output.txt; 
# python *.py &> /dev/null; cd ../..
# cd FSRadPoint/nx128; ./enzo -d *.enzo &> output.txt; 
# python *.py &> /dev/null; cd ../..

echo "  "
echo "finished!"
echo "  "
