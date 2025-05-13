# Clone the dlpoly repository 
repo=https://gitlab.com/ccp5/dl-poly.git
sourcedir=dl-poly$suffix
git clone --depth 1 "$repo" "$sourcedir"

# Change to dl-poly directory
cd dl-poly$suffix
mkdir build-mpi
cmake ../ -DCMAKE_BUILD_TYPE=Release -DWITH_PLUMED=ON -DINTERNAL_PLUMED=OFF

if [[ -x bin/DLPOLY.Z ]]; then
   cp bin/DLPOLY.Z "$HOME/opt/bin/DLPOLY.Z$exeSuffix"
else
   echo "DLPOLY.Z is not executable or does not exist"
fi
