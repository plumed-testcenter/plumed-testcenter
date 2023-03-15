#!/bin/bash

set -e
set -x

plumed_mode="static"
suffix=""

for opt
do
case "$opt" in
  (mode=*) mode="${opt#mode=}" ;;
  (suffix=*) suffix="${opt#suffix=}" ;;
  (*) echo "unknown option $opt" ; exit 1 ;;
esac
done

# Cloning the lammps repository
repo=https://github.com/lammps/lammps.git
echo "cloning repoisitory $repo"
git clone $repo

# Finding the latest stable version of lammps to build
cd lammps
version=$(git tag --sort=-creatordate | grep stable | head -n 1)
echo "installing latest stable lammps $version"
git checkout $version

# Building lammps using cmake
mkdir build
cd build
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HOME/opt/lib/pkgconfig
cmake -D PKG_MANYBODY=yes -D PKG_KSPACE=yes -D PKG_MOLECULE=yes -D PKG_RIGID=yes -D PKG_USER-PLUMED=yes -D PLUMED_MODE=$plumed_mode ../cmake
make
if [ ! -f ./lmp ] ; then
     exit 1
fi
cp ./lmp $HOME/opt/bin/lammps$suffix
