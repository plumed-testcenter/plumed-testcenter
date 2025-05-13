# Clone the dlpoly repository 
repo=https://gitlab.com/ccp5/dl-poly.git
sourcedir=dl-poly$suffix
echo cloning dlpoly repository
git clone --depth 1 "$repo" "$sourcedir"
echo cloned dlpoly repository

# Add plumed to PKG_CONFIG_PATH
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HOME/opt/lib/pkgconfig

# Change to dl-poly directory
echo building dlpoly
cd dl-poly$suffix
mkdir build-mpi
cd build-mpi
echo created dlpoly directory
FFLAGS="-O3 -fallow-argument-mismatch" cmake .. -DCMAKE_BUILD_TYPE=Release -DWITH_PLUMED=ON -DINTERNAL_PLUMED=OFF &> cmake.log
echo ran cmake
cat cmake.log
echo Making code
make &> make.log
echo built dlpoly
cat make.log

if [[ -x bin/DLPOLY.Z ]]; then
   cp bin/DLPOLY.Z "$HOME/opt/bin/DLPOLY.Z$exeSuffix"
else
   echo "DLPOLY.Z is not executable or does not exist"
fi
