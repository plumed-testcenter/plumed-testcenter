# shellcheck disable=SC2154,2148
# Cloning the i-pi repository
git clone --depth 1 https://github.com/i-pi/i-pi.git

# Build the fortran drivers
(
   cd i-pi/drivers/f90 || echo "failed to cd to i-pi/drivers/f90"
   make
)

# Copy i-pi to $HOME/opt
cp -pr i-pi "$HOME/opt"

if [ -d "$HOME/opt/lib/plumed$suffix/python" ]; then
   echo FOUND PYTHON DIRECTORY FOR RUNNING PLUMED
else
   echo DID NOT FIND PYTHON DIRECTORY FOR RUNNING PLUMED
fi

executable=$HOME/opt/bin/i-pi
# Make a script to run i-pi
cat <<EOF >"$executable"
#!/bin/bash
export PYTHONPATH=$HOME/opt/lib/plumed$suffix/python
$HOME/opt/i-pi/bin/i-pi input.xml & sleep 5
$HOME/opt/i-pi/bin/i-pi-driver -m sg -h localhost -o 15 -p 31415
EOF
chmod u+x "$executable"
