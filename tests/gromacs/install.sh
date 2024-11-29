# shellcheck disable=SC2154
# Cloning the gromacs repository
# formatted with shfmt_v3.6.0
git clone --depth 1 --branch v2024.2 https://gitlab.com/gromacs/gromacs.git

# Checkout the correct version of gromacs
cd gromacs || exit 1

# Patch with PLUMED
"plumed$suffix" patch --engine gromacs-2024.2 -p --mode "$mode"

# Run cmake
mkdir build
cd build || exit 1
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX="$HOME/opt/gromacs" -DGMX_ENABLE_CCACHE=ON
# -j is atomatic number of concurrent cores
cmake --build . -j
cmake --install .

if [ -x "$HOME/opt/gromacs/bin/gmx" ]; then
   # Write a script to execute gromacs calculations
   cat <<EOF >"$HOME/opt/bin/gromacs"
#!/bin/bash
"\$HOME/opt/gromacs/bin/gmx" grompp -p topol.top -c conf.gro -f md.mdp
"\$HOME/opt/gromacs/bin/gmx" mdrun -nt 1 -plumed plumed.dat
echo 5 | "\$HOME/opt/gromacs/bin/gmx" energy
EOF
   chmod u+x "$HOME/opt/bin/gromacs"
fi
