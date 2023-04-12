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

# Running your code from python

As explained above the `mdcode.py` file that you must contains a class with functions that serve two purposes:

1. There are functions that allow the testcenter codes to run your MD code from a python script and thus test PLUMED.
2. There are functions that recover various quantities from the output files that your code generates for comparison with the output with PLUMED.

In this section we will cover the functions that perform thee first set of operations; namely, the functions `setParams` and `runMD`.  
For the tests of `simplemd` (the MD code that is part of PLUMED) these two functions look like this:

```python
class mdcode :
  def setParams( self ) :
       params = {
         "temperature": 1.0,
         "tstep": 0.005,
         "friction": 1.0,
         "pressure": 1.0,
         "pfriction": 4
       }
       return params

  def runMD( self, mdparams ) :
       # Prepare a string that contains the input for simplemd
       inp = "inputfile input.xyz\n"
       inp = inp + "outputfile output.xyz\n"
       inp = inp + "temperature " + str(mdparams["temperature"]) + "\n"
       inp = inp + "tstep " + str(mdparams["tstep"]) + "\n"
       inp = inp + "friction " + str(mdparams["friction"]) + "\n"
       inp = inp + "forcecutoff 2.5\n"
       inp = inp + "listcutoff  3.0\n"
       inp = inp + "nstep " + str(mdparams["nsteps"]) + "\n"
       inp = inp + "nconfig 1 trajectory.xyz\n"
       inp = inp + "nstat   1 energies.dat\n"
       of=open("in","w+")
       of.write(inp)
       of.close()
       # Code to deal with restraint 
       if "restraint" in mdparams and mdparams["restraint"]>0 :
          f = open("plumed.dat","w+")
          f.write("dd: DISTANCE ATOMS=1,2 \nRESTRAINT ARG=dd KAPPA=2000 AT=" + str(mdparams["restraint"]) + "\nPRINT ARG=dd FILE=colvar FMT=%8.4f\n")
          f.close()
       # Work out the name of the plumed executable 
       executible = mdparams["executible"]
       cmd = [executible, "simplemd"]
       # Now run the calculation using subprocess
       with open("stdout","w") as stdout:
        with open("stderr","w") as stderr:
           plumed_out = subprocess.run(cmd, text=True, input=inp, stdout=stdout, stderr=stderr )
       return plumed_out.returncode 
```

To run the tests on your MD code the testcenter codes needs to be able to adjust the following parameters.  It does this by passing a dictionary called 
`mdparams` to the function `runMD`.   This function then creates an input file from the information in `mdparams` and runs the calculation with your MD code. 
The key parameters that are passed within the `mdparams` are:

* __ensemble__ - either npt or nvt.  If this variable is npt then the calculation should be run at constant pressure.  If it is nvt then the calculation should be run at constant volume.  If there is not a barostat in your code you can safely ignore this variable and run everything with nvt (as is done for simplemd in the example above).  If you set the tests correctly (see next section) you will not be testing any features that require you to run at constant pressure.
* __temperature__ - the temperature 
* __friction__ - the thermostat friction
* __pressure__ - the pressure  
* __pfriction__ - the barstat frction
* __tstep__ - the integration timestep
* __nsteps__ - the number of steps of MD to perform
* __restraint__ - this parameter is used to test that PLUMED is applying forces on atoms correctly.  If the value passed in this variable is positive then this should be used as the $d_0$ parameter of a harmonic restraint with the form $\frac{1}{2}\kappa(d_{12} - d_0)^2$, where $d_{12}$ is the distance between the first two atoms in your system and $\kappa = 2000 kJ mol$^{-1}$ nm$^{-2}$.  In the example above with simplemd this restraint is applied by PLUMED (because simplemd has no functionality of its own to apply a restraint of this form).  In most other MD codes you should be able to apply the harmonic restraint within the MD code.  The test here then determines whether the time series of distances that is returned when the restraint is applied by PLUMED is the same as the time series that is returned when the restraint is applied within the MD code.
* __executible__ - the name of the MD codes executible.  This is necessary as we test the interface between each code, the latest stable version of PLUMED and the master version of PLUMED.  Two versions of each MD code (with different names) are thus compiled and tested.  The name of the executible that is to be tested is controlled by the underlying code plumed testcenter code.

Notice that although some of these variables (e.g. nsteps) are set by the underlying plumed testcenter code, there are others that must be given sensible initial values.  This process of giving sensible initial values to variables is done by `setParams`.  Parameters here should be set in the units of the MD code (and not in PLUMED units).  Importantly, however, __the pressure must be set equal to 1 bar in whatever internal units your MD code employs__ as we assume that the pressure has been set to 1 bar when we test the virial.
