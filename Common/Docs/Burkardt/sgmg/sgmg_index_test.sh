#! /bin/bash
#
g++ -c -Wall -I/$HOME/include sgmg_index_test.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ sgmg_index_test.o /$HOME/libcpp/sgmg.o /$HOME/libcpp/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm sgmg_index_test.o
#
mv a.out sgmg_index_test
./sgmg_index_test > sgmg_index_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sgmg_index_test
#
echo "Normal end of execution."
