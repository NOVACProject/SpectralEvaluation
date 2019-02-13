REM Start by clearing out the 'cache'
rmdir build
mkdir build
cd build

REM Create the solution
cmake -G "Visual Studio 15 2017 Win64" ..

