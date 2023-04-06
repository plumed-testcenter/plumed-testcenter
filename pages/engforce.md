Check force on energy
---------------------

It is common practise to use the potential energy as a collective energy. Some MD codes thus pass the potential energy to PLUMED and
PLUMED can then apply forces on this collective variable.  We test that any forces that PLUMED applies on the potential energy are 
correctly passed back to the MD code by doing the following test.  We first run a short simulation at $T$ K with a timestep of $\tau$ ps.
During the course of this simulation we monitor the potential energy using the following PLUMED input:

```plumed 
e: ENERGY
PRINT ARG=e FILE=energy1
```

We then run a second simulation (starting from identical conditions) at a temperature of $T\alpha^2$ and with a timestep of $\tau/\alpha$.
The thermostat and barostat relaxation times are similarly divided by $\alpha$.  In the tests that are run on this website we set $\alpha=1.1$.

For this second MD run the following PLUMED input file is used:

```plumed
e: ENERGY
# slope is such that 
PRINT ARG=e FILE=energy2
# slope should be (alpha-1)=0.1
RESTRAINT AT=0.0 ARG=e SLOPE=0.1
```

When forces are passed correctly the time series of the energies from these two calculations should be identical. An NPT version of this calculation 
is performed as well as an NVT calculation if the virial is passed to PLUMED.

# Trajectories

# Results

The table below contains the PLUMED outputs from the two calculations described above.
If the PLUMED interface is working correctly these two sets of numbers should be identical.
