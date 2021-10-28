#! /bin/bash
#
cp tetrahedron_monte_carlo.hpp /$HOME/include
#
g++ -c -Wall -I /$HOME/include tetrahedron_monte_carlo.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv tetrahedron_monte_carlo.o ~/libcpp/tetrahedron_monte_carlo.o
#
echo "Normal end of execution."
