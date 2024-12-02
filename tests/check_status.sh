#!/bin/bash
# formatted with shfmt_v3.6.0

suffix=$(plumed info --version)
suffix=_v$suffix
set -x

for opt; do
     case "$opt" in
     suffix=*) suffix="${opt#suffix=}" ;;
     code=*) code="${opt#code=}" ;;
     *)
          echo "unknown option $opt"
          exit 1
          ;;
     esac
done

executible=$(grep executible "tests/$code/info.yml" | sed -e 's/executible: //')

echo "Looking for '$HOME/opt/bin/$executible' or  for '$HOME/opt/bin/$executible$suffix'"
#thest that at least the non suffixed version exists
if [[ -x $HOME/opt/bin/$executible || -x $HOME/opt/bin/$executible$suffix ]]; then
     echo "install_plumed$suffix: working" >>"tests/$code/info.yml"
     # Ensures we do not overwrite master version of PLUMED when testing simplemd
     if [[ ! -x $HOME/opt/bin/$executible$suffix ]]; then
          #this should be set up in the install script
          cp "$HOME/opt/bin/$executible" "$HOME/opt/bin/$executible$suffix"
     fi
else
     echo "install_plumed$suffix: broken" >>"tests/$code/info.yml"
fi
