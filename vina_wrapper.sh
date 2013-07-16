#!/bin/bash

set -e

# Wrapper for Autodock Vina this shell script will import a receptor a ligand and the center corrdiates 
# and convert the output to a .mol2 format. Afterward it will remove the unnecessary files so only the
# .mol2 is returned
# Version 1.2 - Rob Quick - Added file clean up to streamline output

if [ $# -ne 5 ]; then 
    echo "Usage: vina  "
    exit 1
fi

chmod a+x vina

#Run vina
./vina --cpu 1 --receptor $1.pdbqt --ligand $2.pdbqt --center_x $3 --center_y $4 --center_z $5 --size_x 21 --size_y 21 --size_z 21 --out $1-$2.pdbqt

#Run pdbqt to mol2 conversion
./pdbqt2mol2.py $1-$2.pdbqt $2.mol2 $1-$2.mol2

#Get rid of unnecessary files
rm $1-$2.pdbqt
rm $1-$2.pdbqt-new

echo $HOSTNAME

