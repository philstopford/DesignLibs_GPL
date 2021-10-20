#! /bin/bash
#
g++ -c -Wall -I/$HOME/include dream1_test.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ dream1_test.o \
  /$HOME/libcpp/dream1.o \
  /$HOME/libcpp/rnglib.o -lm
#
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm dream1_test.o
#
mv a.out dream1_test
./dream1_test > dream1_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm dream1_test
#
echo "Normal end of execution."
