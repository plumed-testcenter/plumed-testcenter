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
executible_suffixed=${executible}$suffix

echo -n "Looking for executable ${executible_suffixed} or ${executible}..."

#now I need the suffix with no underscores
suffix=${suffix/_/}
if [[ -x $executible ]] || [[ -x $executible_suffixed ]]; then

     python updateYaml.py "$info_yml" "install_plumed" "${suffix}" working
     if [[ -x $executible_suffixed ]]; then
          echo "found $executible_suffixed"
          # the install script should have done the homework (see below)
          exit 0
     else
          echo "found $executible"
          # the install script did not do the homework
          # this may be risky because means that the install script for _master may override the installed things from _stable
          # In case of shared libraries this may be a problem if the plumed patch is invasive
          # Or in the case the wrapper has changed between versions and it is not statically linked in the executable

          # no need for the extra check, because of the if above
          echo "Copying $executible to $executible_suffixed"
          echo "(this should be set up in the install script)"
          cp "$executible" "$executible_suffixed"

          exit 0
     fi
else
     echo "not found"

fi
echo "Something is wrong with the installation of the patched $code with plumed$suffix"

python updateYaml.py "$info_yml" "install_plumed" "${suffix}" broken
