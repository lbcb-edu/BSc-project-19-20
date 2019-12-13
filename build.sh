mkdir build
cd build

cmake cmake -DCMAKE_BUILD_TYPE=Release ..
make

cd test
ctest