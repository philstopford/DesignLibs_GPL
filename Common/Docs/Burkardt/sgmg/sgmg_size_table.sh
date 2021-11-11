#! /bin/bash
#
g++ -c -Wall -I/$HOME/include sgmg_size_table.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ sgmg_size_table.o /$HOME/libcpp/sgmg.o /$HOME/libcpp/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm sgmg_size_table.o
#
mv a.out sgmg_size_table
./sgmg_size_table > sgmg_size_table_output.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sgmg_size_table
#
echo "Normal end of execution."
