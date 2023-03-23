#!/bin/bash

set -e
set -x

mode="static"
suffix=""
basedir=`pwd`

# Copying arguments from input
# Mode is how plumed is linked static/runtime
# Suffix is which version of plumed to use stable/master
for opt
do
case "$opt" in
  (mode=*) mode="${opt#mode=}" ;;
  (suffix=*) suffix="${opt#suffix=}" ;;
  (*) echo "unknown option $opt" ; exit 1 ;;
esac
done

# Cloning the espresso repository
repo=https://gitlab.com/QEF/q-e.git
echo "cloning repoisitory $repo"
git clone https://gitlab.com/QEF/q-e.git q-e$suffix

# Lets build quantum expression
cd q-e$suffix
# We build the interface with QE 7.0 because we are using the patch 
git checkout qe-7.0
./configure --prefix="$HOME/opt"
plumed$suffix patch --engine qespresso-7.0 -p
make pw
make install

# Checking that lammps has been built correctly
if [ ! -f $home/opt/bin/pw.x ] ; then
     echo install_plumed$suffix: broken >> $basedir/tests/lammps/info.yml
else 
     echo install_plumed$suffix: working >> $basedir/tests/lammps/info.yml
fi
cp ./pw.x $HOME/opt/bin/pw$suffix
