#!/bin/sh

SRC=`dirname $0`/..
echo "Repertoire de source : ${SRC}"

export CC=icc CXX=icpc FC=ifort
ccmake -DOPT_NEW_GLOBAL_METHOD:bool=ON -DOPT_MPI:bool=ON -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo ${SRC}

