# Cloning the lammps repository
repo=https://github.com/gtribello/lammps.git
# https://github.com/lammps/lammps.git
echo "cloning repoisitory $repo"
git clone $repo lammps$suffix

# Finding the latest stable version of lammps to build
cd lammps$suffix
#version=$(git tag --sort=-creatordate | grep stable | head -n 1)
version="fix-plumed"
echo "installing latest stable lammps $version"
git checkout $version

# Building lammps using cmake
mkdir build
cd build
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HOME/opt/lib/pkgconfig
cmake -D PKG_MANYBODY=yes -D PKG_KSPACE=yes -D PKG_MOLECULE=yes -D PKG_RIGID=yes -D PKG_PLUMED=yes -D PLUMED_MODE=$mode ../cmake
make

# Checking that lammps has been built correctly
if [ ! -f ./lmp ] ; then
     echo install_plumed$suffix: broken >> $basedir/tests/lammps/info.yml
else 
     echo install_plumed$suffix: working >> $basedir/tests/lammps/info.yml
     cp ./lmp $HOME/opt/bin/lammps$suffix
fi
