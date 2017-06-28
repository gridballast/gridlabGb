#!/bin/bash

rm -rf local_gd/*
cp -r /tmp/temp/* local_gd/

sed -i '254 a DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"' local_gd/bin/gridlabd
sed -i '255 a prefix=`dirname $DIR`' local_gd/bin/gridlabd
sed -i '254d' local_gd/bin/gridlabd

sed -i '$d' local_gd/lib/gridlabd/glxengine.la
sed -i '$ a DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"' local_gd/lib/gridlabd/glxengine.la
sed -i '$ a libdir=`dirname $DIR`' local_gd/lib/gridlabd/glxengine.la

