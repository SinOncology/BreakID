#!/bin/bash
bash ./generate_installDIR.sh
# compile samtools
cd thirdparty/samtools/samtools-1.3.1/
make -j 16
cd ../../../


mkdir -p build
cd build && rm -rf *
cmake ..

make -j  16
make install
cd ..



