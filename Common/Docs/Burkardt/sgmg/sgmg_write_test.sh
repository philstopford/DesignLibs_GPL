#! /bin/bash
#
g++ -c -Wall -I/$HOME/include sgmg_write_test.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ sgmg_write_test.o /$HOME/libcpp/sgmg.o /$HOME/libcpp/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm sgmg_write_test.o
#
mv a.out sgmg_write_test
./sgmg_write_test > sgmg_write_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sgmg_write_test
#
#  Move sparse grid files to dataset directory.
#
mv *_n.txt ../../datasets/sgmg
mv *_p.txt ../../datasets/sgmg
mv *_r.txt ../../datasets/sgmg
mv *_w.txt ../../datasets/sgmg
mv *_x.txt ../../datasets/sgmg
#
echo "Normal end of execution."
