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
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[250,250,250]\';/" $inputfile
resu_folder=250sub_amgstalone
printfile=250sub_amgstalone.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
echo "========================================"
echo "Second simulation amg stand-alone"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[320,320,320]\';/" $inputfile
resu_folder=320sub_amgstalone
printfile=320sub_amgstalone.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

echo "========================================"
echo "Third simulation amg-stand alone"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[300,300,300]\';/" $inputfile
resu_folder=300sub_amgstalone
printfile=300sub_amgstalone.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

echo "========================================"
echo "Forth simulation amg-stand alone"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[380,380,380]\';/" $inputfile
resu_folder=380sub_amgstalone
printfile=380sub_amgstalone.log
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
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[250,250,250]\';/" $inputfile
resu_folder=250sub_amgaccel
printfile=250sub_amgaccel.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
echo "========================================"
echo "Second simulation amg accelerated"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[300,300,300]\';/" $inputfile
resu_folder=300sub_amgaccel
printfile=300sub_amgaccel.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
echo "========================================"
echo "Third simulation amg accelerated"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[320,320,320]\';/" $inputfile
resu_folder=320sub_amgaccel
printfile=320sub_amgaccel.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
echo "========================================"
echo "Forth simulation amg accelerated"
echo "========================================"
#grid sostitution
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[380,380,380]\';/" $inputfile
resu_folder=380sub_amgaccel
printfile=380sub_amgaccel.log
make
$run >  $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

#------------------------------------------
sed -i -r "s/#define CSC_INTERFACE/\/\/#define CSC_INTERFACE/" $inputconfigfile 
sed -i -r "s/#define AMG_ACCELERATED/\/\/#define AMG_ACCELERATED/" $inputconfigfile 
#------------------------------------------


