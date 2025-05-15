#This script will be embedded in one that has the relevant variables set:
#shellcheck disable=2154,2148
#formatted with shfmt_v3.6.0

# Clone the dlpoly repository
repo=https://gitlab.com/ccp5/dl-poly.git
sourcedir=dl-poly$suffix
git clone --depth 1 "$repo" "$sourcedir"

# Add plumed to PKG_CONFIG_PATH
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HOME/opt/lib/pkgconfig

# Change to dl-poly directory
cd "dl-poly$suffix" || {
   echo "Something failed with cloning dl-poly into dl-poly$suffix" >&2
   exit 1
}

if [[ -n ${suffix} ]]; then
   #this "patches" the pkg-config search in cmake
   #here I am assuming that ${suffix} is appended also to plumed.pc (like plumed_master.pc)
   sed -i CMakeLists.txt \
      -e "s/pkg_check_modules(PLUMED REQUIRED plumed>=\${PLUMED_VERSION})/pkg_check_modules(PLUMED REQUIRED plumed${suffix}>=\${PLUMED_VERSION})/"
fi

mkdir build-mpi
cd build-mpi || {
   echo "cannot cd to build-mpi" >&2
   exit 1
}

FFLAGS="-O3 -fallow-argument-mismatch" cmake .. -DCMAKE_BUILD_TYPE=Release -DWITH_PLUMED=ON -DINTERNAL_PLUMED=OFF &>cmake.log
cat cmake.log
make &>make.log
cat make.log

if [[ -x bin/DLPOLY.Z ]]; then
   cp bin/DLPOLY.Z "$HOME/opt/bin/DLPOLY.Z$exeSuffix"
else
   echo "DLPOLY.Z is not executable or does not exist"
fi
