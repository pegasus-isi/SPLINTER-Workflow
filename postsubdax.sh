#!/bin/bash

set -e

RUN_DIR=$1
RUN_ID=$2
SUBWF_ID=$3

START_DIR=`pwd`

echo "Creating combined dat file..."
mkdir -p $RUN_DIR/scratch/$RUN_ID/$SUBWF_ID
cd $RUN_DIR/scratch/$RUN_ID/$SUBWF_ID/

rm -f $SUBWF_ID.idx
rm -f $SUBWF_ID.fet
rm -f $SUBWF_ID.mol2
for DAT in `find . -type f -name '*.dat'`; do
    BASE=`echo $DAT | sed 's/.dat$//' | sed 's/^\.\///'`
    echo $BASE >>$SUBWF_ID.idx
    cat $BASE.dat >>$SUBWF_ID.fet
    cat $BASE.mol2 >>$SUBWF_ID.mol2
done

mkdir -p $RUN_DIR/outputs-tars
mv $SUBWF_ID.idx $SUBWF_ID.fet $SUBWF_ID.mol2 $RUN_DIR/outputs-tars/

exit 0

