#!/bin/bash
# formatted with shfmt_v3.6.0

suffix=$(plumed info --version)
suffix=_v$suffix

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

info_yml=tests/$code/info.yml

executible=$(grep executible "$info_yml" | sed -e 's/executible: //')
executible=$HOME/opt/bin/$executible
executible_suffixed=$executible$suffix

if [[ -x $executible_suffixed ]]; then
     echo "Found the executable $executible_suffixed"
     echo "install_plumed$suffix: working" >>"$info_yml"

     exit 0
fi

if [[ -x $executible ]]; then
     echo "Found the executable $executible"
     echo "install_plumed$suffix: working" >>"$info_yml"

     # no need for the extra check, because of the if above
     #if [[ ! -x $executible_suffixed ]]; then
     echo "Creating $executible_suffixed"
     #this should be set up in the install script
     cp "$executible" "$executible_suffixed"

     exit 0
fi

echo "install_plumed$suffix: broken" >>"$info_yml"
