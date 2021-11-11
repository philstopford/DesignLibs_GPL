#! /bin/bash
#
g++ -c -Wall -I/$HOME/include sgmg_cc_sl.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ sgmg_cc_sl.o /$HOME/libcpp/sgmg.o /$HOME/libcpp/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm sgmg_cc_sl.o
#
mv a.out sgmg_cc_sl
./sgmg_cc_sl > sgmg_cc_sl_output.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sgmg_cc_sl
#
echo "Normal end of execution."

