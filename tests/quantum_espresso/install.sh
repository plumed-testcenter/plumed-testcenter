# Cloning the espresso repository
repo=https://gitlab.com/QEF/q-e.git
echo "cloning repoisitory $repo"
git clone https://gitlab.com/QEF/q-e.git q-e$suffix

# Lets build quantum espresso
cd q-e$suffix
# We build the interface with QE 7.0 because we are using the patch 
git checkout qe-7.0
./configure --prefix="$HOME/opt"
plumed$suffix patch --engine qespresso-7.0 -p --mode $mode
make pw
make install

