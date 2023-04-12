How to contribute tests to the PLUMED-TESTCENTER
------------------------------------------------

You can add your own MD code to the PLUMED-TESTCENTER. If you do so this website will then test the interface between PLUMED and your MD code once a month.
The results from these tests will then be displayed on [the dashboard](browse.md).  To build your tests you will need to collect together the input files that will
allow you to generate a trajectory with your MD code and PLUMED.  These input files should be gathered together and added to the plumed-testcenter repository in a 
directory called:

```
tests/<your code name>/input
```

You will then need to write three further files in `tests/<your code name>`

* `install.sh` - A shell script that downloads and installs a version of your MD code patched with lammps
* `info.yml` - A yaml file with basic information on your MD code including what tests should be performed
* `mdcode.py` - A python class that tells the testcenter how to run your MD code and how to get various outputs from you code to compare with output from PLUMED.

You are well advised to look the versions of these files that have been written for testing other codes.  In the sections that follow we will try to provide some rationale 
that explains the content of the three files described above.

# Writing an `install.sh` file

The `install.sh` that you will write to download and install your code must be written in bash and will have a structure similar to the one shown below:

```bash
#!/bin/bash

set -x

mode="static"
suffix=""
basedir=`pwd`

# Copying arguments from input
# Mode is how plumed is linked static/runtime
# Suffix is which version of plumed to use stable/master
for opt
do
case "$opt" in
  (mode=*) mode="${opt#mode=}" ;;
  (suffix=*) suffix="${opt#suffix=}" ;;
  (*) echo "unknown option $opt" ; exit 1 ;;
esac
done

# Add code here to download your code from some repository 

# Add code here to build your code so that final executible is in $HOME/opt/bin

# Check that your executible has been generated
if [ ! -f $HOME/opt/bin/<your executible> ] ; then
     echo install_plumed$suffix: broken >> $basedir/tests/<your code name>/info.yml
else 
     echo install_plumed$suffix: working >> $basedir/tests/<your code name>/info.yml
     cp $HOME/opt/bin/<your executible> $HOME/opt/bin/<your executible>$suffix
fi
```

The first part of the script above is common to all the tested MD codes.  This part essentially gets the name of the plumed library that we are linking with (`plumed$suffix`) and the 
mode that we are using to link PLUMED.  You will need to use this information when you build your MD code.

The final part of the code above is also common to all tested MD codes.  This part of the code tells the python script (`build.py`) that constructs [the dashboard page](browse.md) whether or not your code compiled 
sucessfully so that the appropriate badges can be put on the dashboard.  Notice that this information is transferred by adding lines to the `info.yml` file for your code.
