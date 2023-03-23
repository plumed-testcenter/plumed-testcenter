Timestep is passed correctly
----------------------------

PLUMED must receive the timestep from an MD code in order to correctly print the times at which the CV took particular values in COLVAR files. 
When testing MD codes are working with PLUMED we must, therfore, check the timestep that PLUMED is using against the timestep in the MD code.

We can produce a file that contains the times from PLUMED using the following input.  

```plumed
t1: TIME
PRINT ARG=t1 FILE=colvar 
```

The python function `checkTimestep` that you can read in the file `test.py` reads the first column of this file (which contains the times) and compares the difference
in time between each pair of outputs with the timestep that is recovered from the MD code using the function `getTimestep`.  Notice that times output in the by PLUMED output file above
are in ps.  If you are writing tests you should thus ensure that timestep that is passed using `getTimestep` is in ps.
