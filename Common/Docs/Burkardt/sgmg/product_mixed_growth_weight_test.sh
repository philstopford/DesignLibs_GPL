#! /bin/bash
#
g++ -c -Wall -I/$HOME/include product_mixed_growth_weight_test.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ product_mixed_growth_weight_test.o /$HOME/libcpp/sgmg.o /$HOME/libcpp/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm product_mixed_growth_weight_test.o
#
mv a.out product_mixed_growth_weight_test
./product_mixed_growth_weight_test > product_mixed_growth_weight_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm product_mixed_growth_weight_test
#
echo "Normal end of execution."
