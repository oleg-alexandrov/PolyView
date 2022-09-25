#!/bin/bash

qmake polyview.pro

make -j${CPU_COUNT}

make install INSTALL_ROOT=$PREFIX
