# Cloning the gromacs repository
# shellcheck disable=SC2154,2148
git clone --depth 1 https://gitlab.com/gromacs/gromacs.git --branch release-2025 --single-branch "gromacs$suffix"
# Checkout the correct version of gromacs
cd "gromacs$suffix" || {
   echo "Something failed with cloning gromacs into gromacs$suffix" >&2
   # we want to build the site regardless
   exit 0
}

# look mum, no patch!

prefix=$HOME/opt/gromacs-vanilla$suffix

# Run cmake
mkdir build
cd build || {
   echo "cannot cd to build" >&2
   # we want to build the site regardless
   exit 0
}

# this more or less works with thread_mpi IF you keep ONLY one thread
cmake .. -DGMX_USE_PLUMED=ON -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX="$prefix" -DGMX_ENABLE_CCACHE=ON  -Wno-dev
#I want plumed on
cmake . -DGMX_USE_PLUMED=ON

# Now make
cmake --build . -j11
cmake --install .

# plumedKernel is defined in the header that will be prepended to this by the CI
if [ -x "${prefix}/bin/gmx" ]; then
   # Write a script to execute gromacs calculations
   executible=$HOME/opt/bin/gromacs$exeSuffix
   cat <<EOF >"$executible"
#!/bin/bash
export PLUMED_KERNEL=$plumedKernel
mygmx=$prefix/bin/gmx
"\$mygmx" grompp -p topol.top -c conf.gro -f md.mdp
"\$mygmx" mdrun -nt 1 -plumed plumed.dat
echo 5 | "\$mygmx" energy
EOF
   chmod u+x "$executible"
else
   echo "${prefix}/bin/gmx" have not been compiled or is not executable
fi
