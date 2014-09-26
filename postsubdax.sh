#!/bin/bash

set -e

RUN_DIR=$1
RUN_ID=$2
SUBWF_ID=$3

START_DIR=`pwd`

echo "Creating combined dat file..."
cd $RUN_DIR/scratch/$RUN_ID/$SUBWF_ID/
find . -type f -name '*.dat' -exec cat {} \; >ALL.fet 

cd ..
find $SUBWF_ID -name \*.mol2 -print | grep -v ALL.mol2 >$SUBWF_ID/ALL.mol2

echo "Creating output tarball..."
cd $RUN_DIR/scratch/$RUN_ID
mkdir -p $RUN_DIR/outputs-tars
tar czf $RUN_DIR/outputs-tars/$SUBWF_ID.tar.gz $SUBWF_ID

exit 0

