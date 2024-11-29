# shellcheck disable=SC2154,2148
# formatted with shfmt_v3.6.0 
# Cloning the lammps repository
#repo=https://github.com/gtribello/lammps.git
repo=https://github.com/lammps/lammps.git
sourcedir=lammps$suffix
git clone --depth 1 "$repo" "$sourcedir"

# Finding the latest stable version of lammps to build
cd "$sourcedir" || exit 1
#TODO: go back to stable lammps
#version=$(git tag --sort=-creatordate | grep stable | head -n 1)
#version="fix-plumed-cmake"
#echo "installing latest stable lammps $version"
#git checkout $version
echo "installing latest lammps version"

# Building lammps using cmake
mkdir -p build
cd build || exit 1
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HOME/opt/lib/pkgconfig
export CMAKE_C_COMPILER_LAUNCHER=ccache
export CMAKE_CXX_COMPILER_LAUNCHER=ccache
cmake -D CMAKE_CXX_COMPILER=mpic++ -D PKG_MANYBODY=yes -D PKG_KSPACE=yes -D PKG_MOLECULE=yes -D PKG_RIGID=yes -D PKG_PLUMED=yes -D PLUMED_MODE=$mode -D PLUMED_SUFFIX=$suffix ../cmake
cmake --build . -j4

if [[ -x lmp ]]; then
   cp lmp "$HOME/opt/bin/lammps$exeSuffix"
else
   echo "lmp is not executable or does not exist"
   # we want to build the site regardless
   # exit 1
fi
