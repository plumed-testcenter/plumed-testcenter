# Cloning the gromacs repository
git clone https://gitlab.com/gromacs/gromacs.git

# Checkout the correct version of gromacs
cd gromacs
git checkout v2023.2

# Patch with PLUMED
plumed$suffix patch --engine gromacs-2023.2 -p --mode $mode

# Run cmake
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=$HOME/opt/gromacs 

# Now make
make
make install

if [ -f $HOME/opt/gromacs/bin/gmx ] ; then
   # Write a script to execute gromacs calculations
   echo "#!/bin/bash" > $HOME/opt/bin/gromacs
   echo $HOME/opt/gromacs/bin/gmx grompp -p topol.top -c conf.gro -f md.mdp >> $HOME/opt/bin/gromacs
   echo $HOME/opt/gromacs/bin/gmx mdrun -nt 1 -plumed plumed.dat >> $HOME/opt/bin/gromacs
   chmod u+x $HOME/opt/bin/gromacs
fi
