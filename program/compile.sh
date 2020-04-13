#!/bin/bash

rm -rf build
mkdir build
mpic++ -std=c++11 -O3 -o build/lagrangian_distance lagrangian_distance.cpp ../hungarian_algorithm/include/core.hpp ../hungarian_algorithm/include/cost.hpp ../hungarian_algorithm/include/distance.hpp ../hungarian_algorithm/source/core.cpp ../hungarian_algorithm/source/cost.cpp ../hungarian_algorithm/source/distance.cpp ../util/include/general.hpp ../util/include/io.hpp ../util/source/general.cpp ../util/source/io.cpp
