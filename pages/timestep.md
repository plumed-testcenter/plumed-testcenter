Timestep is passed correctly
----------------------------

PLUMED must receive the timestep from an MD code in order to correctly print the times at which the CV took particular values in COLVAR files. 
To test that the timestep is passed correctly we run a short trajectory and output the time after each sstep using the following command:

```plumed
t1: TIME
PRINT ARG=t1 FILE=colvar 
```

# Trajectory

# Results

The table below includes some of the results from the calculation.  The columns contain:

1. The timestep that was obtained from the MD code, $x_{md}$.
2. The timestep that was obtained from PLUMED, $x_{pl}$.
3. The tolerances that were used when comparing these quantities, $\delta$.
4. The values of $100\frac{\vert x_{md} - x_{pl}\vert }{ \delta}$.

If the PLUMED interface is working correctly the first two sets of numbers should be identical and the final column should be filled with zeros.
