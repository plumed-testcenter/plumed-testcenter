# Cloning the espresso repository
repo=https://gitlab.com/QEF/q-e.git
echo "cloning repoisitory $repo"
git clone https://gitlab.com/QEF/q-e.git q-e$suffix

# Lets build quantum espresso
cd q-e$suffix
if [[ $mode = "static" ]]; then
  extraLIBS="LD_LIBS='-lmpi_cxx'"
fi
# We build the interface with QE 7.0 because we are using the patch
git checkout qe-7.0
./configure --prefix="$HOME/opt" $extraLIBS
plumed$suffix patch --engine qespresso-7.0 -p --mode $mode
make pw
make install
