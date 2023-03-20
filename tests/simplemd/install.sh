#!/bin/bash

set -e
set -x

mode="static"
suffix=""

for opt
do
case "$opt" in
  (mode=*) mode="${opt#mode=}" ;;
  (suffix=*) suffix="${opt#suffix=}" ;;
  (*) echo "unknown option $opt" ; exit 1 ;;
esac
done

echo Hello simplemd world

if [ "$suffix" == "_master" ] ; then
   echo install_plumed$suffix: broken >> tests/simplemd/info.yml
else 
   echo install_plumed$suffix: working >> tests/simplemd/info.yml
fi
