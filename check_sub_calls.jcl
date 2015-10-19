#!/bin/bash

#set -vx

rm call_list.txt


for file in *.f90
do

#echo " "
#echo "file = $file "
name=${file%.f90}

#echo "name = $name "

#grep -c -i "call *${name}"  *.f90 | grep -v ':0'  >> call_list.txt
grep  -i "^ *call *${name}"  *.f90   >> call_list.txt
grep  -i "^ *use *${name}"  *.f90   >> call_list.txt

done

#echo "  "
#echo "call_list.txt  "
#cat   call_list.txt
#echo "  "
echo "  "
echo "------------------------------------------------------------"
echo "the following not found in CALL or USE statements           "
echo "------------------------------------------------------------"
echo "  "

for fil in *f90
do

name2=${fil%.f90}

#echo "name2 = $name2 "

result=`grep -i "^.*:.* $name2"  call_list.txt`

if [[ $result == '' ]]
then
echo " $fil     $result "
fi

done
echo "  "

#rm -f call_list.txt

exit

