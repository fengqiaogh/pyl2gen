# install ocssw

```bash
./install_ocssw --install_dir=$HOME/ocssw --tag V2024.8 –-src --goci --common
```

```bash
apt install cmake gdb git gcc g++ gfortran tcsh bison flex zlib1g-dev libx11-dev pkg-config build-essential cmake libpthread-stubs0-dev
```

## install proj
```bash
wget https://download.osgeo.org/proj/proj-data-1.19.tar.gz
mkdir -p proj_data
tar -xvzf proj-data-1.19.tar.gz -C proj_data
```

```bash
vim ~/.bashrc

#ocssw config
export OCSSWROOT=$HOME/ocssw
source $OCSSWROOT/ocssw_src/OCSSW_bash.env
export PROJ_DATA=$HOME/proj_data
export OCSSW_DEBUG=1
```


```bash
cd $OCSSWROOT/opt/src
./BuildIt.py
```
```bash
cd $OCSSWROOT/ocssw_src
mkdir build 
cd build
cmake ..
make
```