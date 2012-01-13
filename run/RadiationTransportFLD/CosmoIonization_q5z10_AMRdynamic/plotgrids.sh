#!/bin/bash
source /usr/local/yt/bin/activate

for i in {0..9}
do
    if [ -f DD000$i/data000$i ]
    then
	yt plot -f HII_Density -a 2 --show-grids DD000$i/data000$i
    else
	break
    fi
done

for i in {10..99}
do
    if [ -f DD00$i/data00$i ]
    then
	yt plot -f HII_Density -a 2 --show-grids DD00$i/data00$i
    else
	break
    fi
done

for i in {100..999}
do
    if [ -f DD0$i/data0$i ]
    then
	yt plot -f HII_Density -a 2 --show-grids DD0$i/data0$i
    else
	break
    fi
done

for i in {1000..9999}
do
    if [ -f DD$i/data$i ]
    then
	yt plot -f HII_Density -a 2 --show-grids DD$i/data$i
    else
	break
    fi
done
