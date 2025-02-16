#/bin/bash

# Start by clearing out the 'cache'
rm -rf build
mkdir build
cd build

# Create the solution
cmake -G "Unix Makefiles" ..

# Now build the solution
make -j4


