mkdir build
cd build
cmake -GNinja -DCMAKE_BUILD_TYPE=Release ..
ninja
ninja install

#qmake polyview.pro -spec win32-msvc
#qmake polyview.pro -spec win32-g++
#make
