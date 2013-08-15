#!/bin/bash

set -e

# Wrapper for Autodock Vina this shell script will import a receptor a ligand and the center corrdiates 
# and convert the output to a .mol2 format. 
# Version 1.2 - Rob Quick - Added file clean up to streamline output
# Version 1.3 - Mats Rynge - Write output files to a sub directory based on the subwf id

if [ $# -ne 6 ]; then 
    echo "Usage: vina  "
    exit 1
fi

chmod a+x vina

# run vina
./vina --cpu 1 --receptor $2.pdbqt --ligand $3.pdbqt --center_x $4 --center_y $5 --center_z $6 --size_x 21 --size_y 21 --size_z 21 --out $2-$3.pdbqt

# run pdbqt to mol2 conversion
./pdbqt2mol2.py $2-$3.pdbqt $3.mol2 $2-$3.mol2

