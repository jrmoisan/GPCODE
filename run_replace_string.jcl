#!/bin/bash



for file in *.f90
do


./replace_string.py $file

mv UC_${file} $file


done


exit
