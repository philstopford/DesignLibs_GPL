#! /bin/bash
#
g++ -c -Wall -I/$HOME/include tetrahedron_monte_carlo_test.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ tetrahedron_monte_carlo_test.o /$HOME/libcpp/tetrahedron_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm tetrahedron_monte_carlo_test.o
#
mv a.out tetrahedron_monte_carlo_test
./tetrahedron_monte_carlo_test > tetrahedron_monte_carlo_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm tetrahedron_monte_carlo_test
#
echo "Normal end of execution."
