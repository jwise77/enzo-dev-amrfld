#!/bin/sh

echo "    "
echo "    "
echo "Turner & Stone 1 Test:"
cd TurnerStoneEquil1
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Turner & Stone 1 Split Test"
cd TurnerStoneEquil1_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Turner & Stone 2 Test"
cd TurnerStoneEquil2
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Turner & Stone 2 Split Test"
cd TurnerStoneEquil2_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream X0 Test"
cd RadiationStreamX0
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream X0 Split Test"
cd RadiationStreamX0_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream Y1 Test"
cd RadiationStreamY1
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream Y1 Split Test"
cd RadiationStreamY1_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream Z0 Test"
cd RadiationStreamZ0
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream Z0 Split Test"
cd RadiationStreamZ0_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream 1D Test"
cd RadiationStream1D
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream 1D Split Test"
cd RadiationStream1D_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiating Shock Lab Test"
cd RadiatingShockLab
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiating Shock Lab Split Test"
cd RadiatingShockLab_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiating Shock Lab 1D Test"
cd RadiatingShockLab1D
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiating Shock Lab 1D Split Test"
cd RadiatingShockLab1D_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Iliev et al. Test 1"
cd RHIonization1
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Iliev et al. Split Test 1"
cd RHIonization1_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Iliev et al. Test 2"
cd RHIonization2
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Iliev et al. Split Test 2"
cd RHIonization2_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.5 z0=4 Test"
cd CosmoIonization_q5z4
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.5 z0=4 Split Test"
cd CosmoIonization_q5z4_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.05 z0=4 Test"
cd CosmoIonization_q05z4
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.05 z0=4 Split Test"
cd CosmoIonization_q05z4_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.5 z0=10 Test"
cd CosmoIonization_q5z10
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.5 z0=10 Split Test"
cd CosmoIonization_q5z10_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.05 z0=10 Test"
cd CosmoIonization_q05z10
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.05 z0=10 Split Test"
cd CosmoIonization_q05z10_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Finished!!"
