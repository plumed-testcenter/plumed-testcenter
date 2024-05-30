# Cloning the i-pi repository
git clone https://github.com/i-pi/i-pi.git

# Build the fortran drivers
cd i-pi/drivers/f90
make
cd ../../../../

# Copy i-pi to $HOME/opt
cp -pr i-pi $HOME/opt

if [ ! -f $HOME/opt/i-pi/bin/i-pi ] ; then
   echo could not find i-pi
else 
   echo installed i-pi in $HOME/opt/i-pi/bin/i-pi
fi

if [ ! -f $HOME/opt/i-pi/bin/i-pi-driver ] ; then
   echo could not find i-pi-driver
else 
   echo installed i-pi-driver in $HOME/opt/i-pi/bin/i-pi-driver
fi


# Make a script to run i-pi
echo "#!/bin/bash" > $HOME/opt/bin/i-pi
echo "PYTHONPATH=$HOME/opt/lib/plumed$suffix/python" >> $HOME/opt/bin/i-pi
echo "$HOME/opt/i-pi/bin/i-pi install.xml & sleep 5; $HOME/opt/i-pi/bin/i-pi-driver -m sg -h localhost -o 15 -p 31415" >> $HOME/opt/bin/i-pi
chmod u+x $HOME/opt/bin/i-pi
