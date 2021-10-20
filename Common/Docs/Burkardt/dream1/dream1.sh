#! /bin/bash
#
cp dream1.hpp /$HOME/include
#
g++ -c -Wall -I/$HOME/include dream1.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv dream1.o ~/libcpp/dream1.o
#
echo "Normal end of execution."
