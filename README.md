# Normal installation:

mkdir build
cd build
cmake .. -D OMP=ON -D XDR=ON
make -j
cd ../

# Debug mode:
mkdir build_debug
cd build_debug
cmake -D CMAKE_BUILD_TYPE=DEBUG ..
make -j
cd ../
