# shellcheck disable=SC2154,2148
# Cloning the gromacs repository
# formatted with shfmt_v3.6.0
# Clone direclty into the wanted gromacs version
git clone --depth 1 --branch v2024.2 https://gitlab.com/gromacs/gromacs.git "gromacs$suffix"

# setting up the versionsuffix as in runtests scripts
versionSuffix

exeSuffix=v$("plumed$suffix" info --version)

if suffix="_master"; then
   exeSuffix=""
fi

cd "gromacs$suffix" || exit 1

# Patch with PLUMED
"plumed$suffix" patch --engine gromacs-2024.2 -p --mode "$mode"

prefix=$HOME/opt/gromacs$suffix

# Run cmake
mkdir build
cd build || exit 1
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX="$prefix" -DGMX_ENABLE_CCACHE=ON
# -j is atomatic number of concurrent cores
cmake --build . -j4
cmake --install .

if [ -x "${prefix}/bin/gmx" ]; then
   # Write a script to execute gromacs calculations
   cat <<EOF >"$HOME/opt/bin/gromacs$exeSuffix"
#!/bin/bash
mygmx=\$HOME/opt/gromacs$suffix/bin/gmx
"\$mygmx" grompp -p topol.top -c conf.gro -f md.mdp
"\$mygmx" mdrun -nt 1 -plumed plumed.dat
echo 5 | "\$mygmx" energy
EOF
   chmod u+x "$HOME/opt/bin/gromacs$exeSuffix"
fi
