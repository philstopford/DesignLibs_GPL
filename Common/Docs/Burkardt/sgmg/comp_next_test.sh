#! /bin/bash
#
g++ -c -Wall -I/$HOME/include comp_next_test.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ comp_next_test.o /$HOME/libcpp/sgmg.o \
                    /$HOME/libcpp/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm comp_next_test.o
#
mv a.out comp_next_test
./comp_next_test > comp_next_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm comp_next_test
#
echo "Normal end of execution."
