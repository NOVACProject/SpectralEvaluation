REM Start by clearing out the 'cache'
rmdir build
mkdir build
cd build

REM Create the solution
cmake -G "Visual Studio 17 2022" -A x64 ..

