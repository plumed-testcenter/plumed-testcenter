# Cloning the i-pi repository
git clone https://github.com/i-pi/i-pi.git

# Build the fortran drivers
cd i-pi/drivers/f90
make
cd ../../../

# Copy i-pi and i-pi driver to $HOME/opt/bin
cp bin/i-pi $HOME/opt/bin/i-pi-code
cp bin/i-pi-driver $HOME/opt/bin/

# Make a script to run i-pi
echo "#!/bin/bash" > $HOME/opt/bin/i-pi
echo "PYTHONPATH=$HOME/opt/lib/plumed$suffix/python" >> $HOME/opt/bin/i-pi
echo "$HOME/opt/bin/i-pi-code install.xml & sleep 5; $HOME/opt/bin/i-pi-driver -m sg -h localhost -o 15 -p 31415" >> $HOME/opt/bin/i-pi
chmod u+x $HOME/opt/bin/i-pi
