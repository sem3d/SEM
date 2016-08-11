#!/bin/bash

DIR=`dirname $0`
DTD=${DIR}/Xdmf.dtd

xmllint --xinclude --output /dev/null --dtdvalid ${DTD} $*
