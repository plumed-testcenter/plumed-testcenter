#!/bin/bash

suffix=`plumed info --version`
suffix=`echo _v$suffix`
set -x

for opt
do
case "$opt" in
  (suffix=*) suffix="${opt#suffix=}" ;;
  (code=*)  code="${opt#code=}" ;;
  (*) echo "unknown option $opt" ; exit 1 ;;
esac
done

executible=`grep executible tests/$code/info.yml | sed -e s/"executible: "//`

if [ ! -f $HOME/opt/bin/$executible ] ; then
     echo install_plumed$suffix: broken >> tests/$code/info.yml
else
     echo install_plumed$suffix: working >> tests/$code/info.yml
     # Ensures we do not overwrite master version of PLUMED when testing simplemd
     if [ ! -f $HOME/opt/bin/$executible$suffix ] ; then
          cp $HOME/opt/bin/$executible $HOME/opt/bin/$executible$suffix
     fi
fi
