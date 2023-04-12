# Cloning the espresso repository
repo=https://gitlab.com/QEF/q-e.git
echo "cloning repoisitory $repo"
git clone https://gitlab.com/QEF/q-e.git q-e$suffix

# Lets build quantum expression
cd q-e$suffix
# We build the interface with QE 7.0 because we are using the patch 
git checkout qe-7.0
./configure --prefix="$HOME/opt"
plumed$suffix patch --engine qespresso-7.0 -p
make pw
make install

# Checking that quantum espresso has been built correctly
if [ ! -f $HOME/opt/bin/pw.x ] ; then
     echo install_plumed$suffix: broken >> $basedir/tests/quantum_espresso/info.yml
else 
     echo install_plumed$suffix: working >> $basedir/tests/quantum_espresso/info.yml
     cp $HOME/opt/bin/pw.x $HOME/opt/bin/pw$suffix
fi
