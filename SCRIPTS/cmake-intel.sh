#!/bin/sh

SRC=`dirname $0`/..
echo "Repertoire de source : ${SRC}"

export CC=icc FC=ifort F77=ifort CXX=icpc
ccmake -DOPT_MPI:bool=ON -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo ${SRC}

