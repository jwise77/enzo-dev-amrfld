#!/bin/bash
###############################################################
#
# YT-based plotting script to quickly generate projections 
# of solution fields.
# 
# By default, this script will generate projections of the 
# Grav_Potential field orthogonal to the x-axis.  
#
# The field may be changed through setting the "FIELD" 
# environment variable to a corresponding field name in the
# HDF5 output files.
#
# The axis may be changed through setting the "AXIS" 
# environment variable, with allowed values of {0,1,2} 
# corresponding to the {x,y,z} axes, respectively.
#
# Daniel R. Reynolds
# June 2012
#
###############################################################

source /usr/local/yt/bin/activate

# set plotting field if undefined
if [ $FIELD ]
then
    SKIP=1
else
    FIELD=Grav_Potential
fi


# set plotting axis
if [ $AXIS ]
then
    SKIP=1
else
    AXIS=0
fi


# generate plots
for i in {0..9}
do
    if [ -f DD000$i/data000$i ]
    then
	yt plot -f $FIELD -a $AXIS --show-grids DD000$i/data000$i
    else
	break
    fi
done

for i in {10..99}
do
    if [ -f DD00$i/data00$i ]
    then
	yt plot -f $FIELD -a $AXIS --show-grids DD00$i/data00$i
    else
	break
    fi
done

for i in {100..999}
do
    if [ -f DD0$i/data0$i ]
    then
	yt plot -f $FIELD -a $AXIS --show-grids DD0$i/data0$i
    else
	break
    fi
done

for i in {1000..9999}
do
    if [ -f DD$i/data$i ]
    then
	yt plot -f $FIELD -a $AXIS --show-grids DD$i/data$i
    else
	break
    fi
done
