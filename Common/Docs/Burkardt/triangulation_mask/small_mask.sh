#! /bin/bash
#
g++ -c -Wall small_mask.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ ~/libcpp/triangulation_mask.o small_mask.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm small_mask.o
#
chmod ugo+x a.out
mv a.out triangulation_mask
./triangulation_mask small > small_mask_output.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm triangulation_mask
#
echo "Normal end of execution."
