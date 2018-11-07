# utilitySimGB

Source code for the utility scale simulator for the GridBallast load control devices.

# Installation and Walkthrough

macOS binaries are available on the releases page.

Details about how to configure and use the simulator are available in the Jupyter notebook walkthrough.

# Build Instructions

```sh
git clone https://github.com/gridballast/utilitySimGB
cd workdir
git checkout feature/730
autoreconf -isf
./configure --prefix=$PWD/install --enable-silent-rules
make
make install
cp gridlabd.sh install
chmod +x install/gridlabd.sh
install/gridlabd.sh --version
```
