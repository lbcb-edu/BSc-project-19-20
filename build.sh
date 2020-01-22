mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Relese ..
make

cd test
ctest