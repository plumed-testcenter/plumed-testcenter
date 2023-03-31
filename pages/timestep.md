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

The table below contains times at which each trajectory step should have ended based on the input timestep and the times 
that were output by the command above.  If the PLUMED interface is working correctly these two sets of numbers should be identical.
