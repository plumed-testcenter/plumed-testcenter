# Clone the dlpoly repository 
repo=https://gitlab.com/ccp5/dl-poly.git
sourcedir=dl-poly$suffix
git clone --depth 1 "$repo" "$sourcedir"

# Add plumed to PKG_CONFIG_PATH
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HOME/opt/lib/pkgconfig

# Change to dl-poly directory
cd dl-poly$suffix
mkdir build-mpi
cd build-mpi
FFLAGS="-O3 -fallow-argument-mismatch" cmake .. -DCMAKE_BUILD_TYPE=Release -DWITH_PLUMED=ON -DINTERNAL_PLUMED=OFF &> cmake.log
cat cmake.log
make &> make.log
cat make.log

if [[ -x bin/DLPOLY.Z ]]; then
   cp bin/DLPOLY.Z "$HOME/opt/bin/DLPOLY.Z$exeSuffix"
else
   echo "DLPOLY.Z is not executable or does not exist"
fi
