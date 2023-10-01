#!/bin/bash

qmake                                \
  QMAKE_CXXFLAGS="-I$PREFIX/include" \
  QMAKE_CC=$CC_FOR_BUILD             \
  QMAKE_CXX=$CXX_FOR_BUILD           \
  QMAKE_LINK=$CXX_FOR_BUILD          \
  polyview.pro
make -j${CPU_COUNT}
make install INSTALL_ROOT=$PREFIX
