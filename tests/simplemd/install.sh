#!/bin/bash

set -x

# Copying arguments from input
# Mode is how plumed is linked static/runtime
# Suffix is which version of plumed to use stable/master
mode="static"
suffix=""
basedir=`pwd`

# Reading arguments
for opt
do
case "$opt" in
  (mode=*) mode="${opt#mode=}" ;;
  (suffix=*) suffix="${opt#suffix=}" ;;
  (*) echo "unknown option $opt" ; exit 1 ;;
esac
done

# Some dummy code as we don't need to do anything to build simplemd 
# if plumed has been built
echo Hello simplemd world

# Check the required version of plumed exists
if ! command -v plumed$suffix &> /dev/null ; then
   echo install_plumed$suffix: broken >> $basedir/tests/simplemd/info.yml
else
   echo install_plumed$suffix: working >> $basedir/tests/simplemd/info.yml
fi
