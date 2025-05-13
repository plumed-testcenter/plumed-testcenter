# Clone the dlpoly repository 
repo=https://gitlab.com/ccp5/dl-poly.git
sourcedir=dl-poly$suffix
echo cloning dlpoly repository
git clone --depth 1 "$repo" "$sourcedir"
echo cloned dlpoly repository

# Change to dl-poly directory
echo building dlpoly
cd dl-poly$suffix
mkdir build-mpi
echo created dlpoly directory
cmake ../ -DCMAKE_BUILD_TYPE=Release -DWITH_PLUMED=ON -DINTERNAL_PLUMED=OFF
echo ran cmake
make
echo built dlpoly

if [[ -x bin/DLPOLY.Z ]]; then
   cp bin/DLPOLY.Z "$HOME/opt/bin/DLPOLY.Z$exeSuffix"
else
   echo "DLPOLY.Z is not executable or does not exist"
fi
