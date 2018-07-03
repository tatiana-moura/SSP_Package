#!/bin/sh

echo "Compiling DARTMOUTH codes"
# codes from http://stellar.dartmouth.edu/models/programs.html
cd ./DARTMOUTH
gfortran iso_interp_feh.f -o iso_interp_feh
gfortran isolf_split.f -o isolf_split
cd ..

echo "Compiling innewmarcs"
cd atm_models
gfortran innewmarcs2_grid_dwarfs2.f -o innewmarcs2_grid_dwarfs2
gfortran innewmarcs2_grid_giants_g_gt_1.f -o innewmarcs2_grid_giants_g_gt_1
gfortran innewmarcs2_grid_giants_g_lt_1.f -o innewmarcs2_grid_giants_g_lt_1
cd ..

echo "Compiling pfantgrade"
rm -f *.o
gfortran -c -g pfantgrade.f
gfortran -c -g pkapgeralgrade.f
gfortran -o pfantgrade pfantgrade.o pkapgeralgrade.o

echo "Compiling nulbadgrade"
gfortran -c nulbadgrade.f
gfortran -o nulbadgrade nulbadgrade.o

echo "Compiling hydro2"
gfortran -c -g hydro2.f
gfortran -c -g absor.f
gfortran -c -g calhy.f
gfortran -o hydro2 hydro2.o absor.o calhy.o 


if [ ! -d "SSP_Spectra" ]; then
  mkdir SSP_Spectra
  echo "Folder SSP_Spectra created"
fi

if [ ! -d "Stellar_Spectra" ]; then
  mkdir Stellar_Spectra
  echo "Folder Stellar_Spectra created"
fi

if [ ! -d "Stellar_pars" ]; then
  mkdir Stellar_pars
  echo "Folder Stellar_pars created"
fi


