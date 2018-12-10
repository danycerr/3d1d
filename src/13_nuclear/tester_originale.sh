#!/bin/bash

inputfile=input_stationary_nano_linf_merge.param
inputconfigfile=../../include/transport3d1d.cpp
run="./M3D1D $inputfile"
#------------------------------------------
sed -i -r "s/\/\/#define CSC_INTERFACE/#define CSC_INTERFACE/" $inputconfigfile 
sed -i -r "s/\/\/#define AMG_STAND_ALONE/#define AMG_STAND_ALONE/" $inputconfigfile 
#------------------------------------------

echo "========================================"
echo "First simulation amg stand-alone"
echo "========================================"
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[20,20,20]\';/" $inputfile
resu_folder=20sub_amgstalone
printfile=20sub_amgstalone.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
echo "========================================"
echo "Second simulation amg stand-alone"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[40,40,40]\';/" $inputfile
resu_folder=40sub_amgstalone
printfile=40sub_amgstalone.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

echo "========================================"
echo "Third simulation amg-stand alone"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[60,60,60]\';/" $inputfile
resu_folder=60sub_amgstalone
printfile=60sub_amgstalone.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

echo "========================================"
echo "Forth simulation amg-stand alone"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[80,80,80]\';/" $inputfile
resu_folder=80sub_amgstalone
printfile=80sub_amgstalone.log
make
$run >  $printfile
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
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[20,20,20]\';/" $inputfile
resu_folder=20sub_amgaccel
printfile=20sub_amgaccel.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
echo "========================================"
echo "Second simulation amg accelerated"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[40,40,40]\';/" $inputfile
resu_folder=40sub_amgaccel
printfile=40sub_amgaccel.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
echo "========================================"
echo "Third simulation amg accelerated"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[60,60,60]\';/" $inputfile
resu_folder=60sub_amgaccel
printfile=60sub_amgaccel.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
echo "========================================"
echo "Forth simulation amg accelerated"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[80,80,80]\';/" $inputfile
resu_folder=80sub_amgaccel
printfile=80sub_amgaccel.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

#------------------------------------------
sed -i -r "s/#define CSC_INTERFACE/\/\/#define CSC_INTERFACE/" $inputconfigfile 
sed -i -r "s/#define AMG_ACCELERATED/\/\/#define AMG_ACCELERATED/" $inputconfigfile 
#------------------------------------------

#------------------------------------------------------------------------------------
sed -i -r "s/\/\/#define CSC_INTERFACE/#define CSC_INTERFACE/" $inputconfigfile 
sed -i -r "s/\/\/#define DIRECT_SOLVER/#define DIRECT_SOLVER/" $inputconfigfile 
#------------------------------------------

echo "========================================"
echo "First simulation Direct solver"
echo "========================================"
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[20,20,20]\';/" $inputfile
resu_folder=20sub_direct
printfile=20sub_direct.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
echo "========================================"
echo "Second simulation Direct Solver"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[40,40,40]\';/" $inputfile
resu_folder=40sub_direct
printfile=40sub_direct.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
echo "========================================"
echo "Third simulation Direct Solver"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[60,60,60]\';/" $inputfile
resu_folder=60sub_direct
printfile=60sub_direct.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
echo "========================================"
echo "Forth simulation Direct Solver"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[80,80,80]\';/" $inputfile
resu_folder=80sub_direct
printfile=80sub_direct.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

#------------------------------------------
sed -i -r "s/#define CSC_INTERFACE/\/\/#define CSC_INTERFACE/" $inputconfigfile 
sed -i -r "s/#define DIRECT_SOLVER/\/\/#define DIRECT_SOLVER/" $inputconfigfile 
#------------------------------------------
