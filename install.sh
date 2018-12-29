#!/bin/bash
cd thirdparty/googletest/googletest-1.8.0/
cmake .
make -j 16
cd ../../../
cd thirdparty/samtools/samtools-1.3.1/
make
cd ../../../


mkdir -p build
cd build && rm -rf *
cmake ..
cmake ..
cd ../
bash ./make_ref_data.sh
cd build
make -j  16
make install
cd ../

echo "generate version.txt file"

branch=`git rev-parse --abbrev-ref HEAD`
shaID=`git log|head -1|awk '{print $2}'`

echo "#BreakID_branch: ${branch}" >version.txt
echo "#BreakID_shaID: ${shaID}" >>version.txt
cat version.txt

