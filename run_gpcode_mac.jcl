#!/bin/bash

cdir=`pwd`
homedir=$cdir
cd $homedir
mkdir CDOM
cd CDOM
pwd

datestamp=`date +%C%y%m%d`
timestamp=`date +%H%M%S`

datetime=${datestamp}_${timestamp}

echo "date_time  $datetime"

mkdir   testGP_${datestamp}_${timestamp}
cd testGP_${datestamp}_${timestamp}

echo "datetime ${datestamp}_${timestamp} "

cp ../../main.x  .
cp ../../GPGA_cntl_vars.in  .
cp ../../CDOM.data  .

mkdir Trees
mkdir pts
mkdir pts/Trees

/opt/local/bin/mpirun-openmpi-mp -np 16 ./main.x
#mpirun -np 32 ./GPCODE_Tree

exit 0
