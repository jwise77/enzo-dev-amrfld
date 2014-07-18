#!/bin/bash
source /usr/local/yt_dev/bin/activate

for i in {0..9}
do
    if [ -f DD000$i/data000$i ]
    then
	yt plot -p -f Density -a 2 --show-grids DD000$i/data000$i
	yt plot -p -f HI_Density -a 2 --show-grids DD000$i/data000$i
	yt plot -p -f Temperature -a 2 --show-grids DD000$i/data000$i
#	yt plot -p -f TotalEnergy -a 2 --show-grids DD000$i/data000$i
    else
	break
    fi
done

for i in {10..99}
do
    if [ -f DD00$i/data00$i ]
    then
	yt plot -p -f Density -a 2 --show-grids DD00$i/data00$i
	yt plot -p -f HI_Density -a 2 --show-grids DD00$i/data00$i
	yt plot -p -f Temperature -a 2 --show-grids DD00$i/data00$i
#	yt plot -p -f TotalEnergy -a 2 --show-grids DD00$i/data00$i
    else
	break
    fi
done

for i in {100..999}
do
    if [ -f DD0$i/data0$i ]
    then
	yt plot -p -f Density -a 2 --show-grids DD0$i/data0$i
	yt plot -p -f HI_Density -a 2 --show-grids DD0$i/data0$i
	yt plot -p -f Temperature -a 2 --show-grids DD0$i/data0$i
#	yt plot -p -f TotalEnergy -a 2 --show-grids DD0$i/data0$i
    else
	break
    fi
done

for i in {1000..9999}
do
    if [ -f DD$i/data$i ]
    then
	yt plot -p -f Density -a 2 --show-grids DD$i/data$i
	yt plot -p -f HI_Density -a 2 --show-grids DD$i/data$i
	yt plot -p -f Temperature -a 2 --show-grids DD$i/data$i
#	yt plot -p -f TotalEnergy -a 2 --show-grids DD$i/data$i
    else
	break
    fi
done

cd frames
avconv -r 1 -i data%04d_Projection_z_Temperature.png Temperature.avi
avconv -r 1 -i data%04d_Projection_z_Density.png Density.avi
avconv -r 1 -i data%04d_Projection_z_HI_Density.png HI_Density.avi
#avconv -r 1 -i data%04d_Projection_z_TotalEnergy.png TotalEnergy.avi
mkdir Temperature Density HI_Density #TotalEnergy
#mv *_TotalEnergy.png TotalEnergy
mv *_Temperature.png Temperature
mv *_HI_Density.png HI_Density
mv *_Density.png Density


