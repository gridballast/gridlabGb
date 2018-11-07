#! /bin/sh

prefix=$(cd $(dirname "$1") && pwd -P) 
#prefix="/home/justin/omf/omf/solvers/gldv990Linux"

export GLPATH="${prefix}"
#echo $GLPATH

"${prefix}/gridlabd.bin" "$@"
