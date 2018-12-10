#!/bin/bash

inputfile=input_amgtest.param
inputconfigfile=../../include/transport3d1d.cpp
run="./M3D1D $inputfile"
#------------------------------------------
sed -i -r "s/\/\/#define CSC_INTERFACE/#define CSC_INTERFACE/" $inputconfigfile 
sed -i -r "s/\/\/#define AMG_STAND_ALONE/#define AMG_STAND_ALONE/" $inputconfigfile 
#------------------------------------------

echo "========================================"
echo "First simulation amg stand-alone"
echo "========================================"
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[60,60,60]\';/" $inputfile
resu_folder=60sub_amgstalone
printfile=60sub_amgstalone_noidump.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
#------------------------------------------
sed -i -r "s/#define CSC_INTERFACE/\/\/#define CSC_INTERFACE/" $inputconfigfile 
sed -i -r "s/#define AMG_STAND_ALONE/\/\/#define AMG_STAND_ALONE/" $inputconfigfile 
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
sed -i -r "s/\/\/#define CSC_INTERFACE/#define CSC_INTERFACE/" $inputconfigfile 
sed -i -r "s/\/\/#define AMG_ACCELERATED/#define AMG_ACCELERATED/" $inputconfigfile 
#------------------------------------------

echo "========================================"
echo "First simulation amg accelerated"
echo "========================================"
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[60,60,60]\';/" $inputfile
resu_folder=60sub_amgaccel
printfile=60sub_amgacce_noidump.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

#------------------------------------------
sed -i -r "s/#define CSC_INTERFACE/\/\/#define CSC_INTERFACE/" $inputconfigfile 
sed -i -r "s/#define AMG_ACCELERATED/\/\/#define AMG_ACCELERATED/" $inputconfigfile 
#------------------------------------------


