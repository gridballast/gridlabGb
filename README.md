# utilitySimGB

Utility scale simulator for the GridBallast load control devices.

# Installation and Walkthrough

macOS binaries are available on the [releases page](https://github.com/gridballast/utilitySimGB/releases).

Details about how to configure and use the simulator are available in the [Jupyter notebook walkthrough](https://github.com/gridballast/utilitySimGB/blob/master/controller_usage_demonstration.ipynb).

# Build Instructions

```sh
git clone https://github.com/gridballast/utilitySimGB
cd workdir
git checkout feature/730
autoreconf -isf
./configure --prefix=$PWD/install --enable-silent-rules
make
make install
cp ezgridlab.sh install/gridlabd.sh
chmod +x install/gridlabd.sh
install/gridlabd.sh --version
```
