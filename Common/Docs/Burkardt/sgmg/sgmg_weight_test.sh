#! /bin/bash
#
g++ -c -Wall -I/$HOME/include sgmg_weight_test.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ sgmg_weight_test.o /$HOME/libcpp/sgmg.o /$HOME/libcpp/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm sgmg_weight_test.o
#
mv a.out sgmg_weight_test
./sgmg_weight_test > sgmg_weight_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sgmg_weight_test
#
echo "Normal end of execution."
