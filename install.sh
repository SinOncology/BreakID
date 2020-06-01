#!/bin/bash
cd thirdparty/googletest/googletest-1.8.0/
cmake .
make clean
make -j 32
cd ../../../
cd thirdparty/samtools/samtools-1.3.1/
make clean
make -j 32
cd ../../../


mkdir -p build
cd build && rm -rf *
cmake ..
cmake ..
cd ../
 ./generate_installDIR.sh
cd build
make -j  16
make install
cd ../

