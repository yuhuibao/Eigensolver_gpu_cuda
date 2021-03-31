#!/bin/bash

cd lib_eigsolve
# make clean
make
cd ../test_driver
make clean
make test_dsytrd
./test_dsytrd
cd ..