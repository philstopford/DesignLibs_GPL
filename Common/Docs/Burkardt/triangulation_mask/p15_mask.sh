#! /bin/bash
#
g++ -c -Wall p15_mask.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ ~/libcpp/triangulation_mask.o p15_mask.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm p15_mask.o
#
chmod ugo+x a.out
mv a.out triangulation_mask
./triangulation_mask p15 > p15_mask_output.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm triangulation_mask
#
echo "Normal end of execution."
