Energy is passed correctly
--------------------------

It is common practise to use the potential energy as a collective energy.  Some MD codes thus pass the potential energy to PLUMED. 
To check that this quantity has been passed correctly we can output the passed energy from PLUMED using the following input.  

```plumed
e: ENERGY 
PRINT ARG=e FILE=colvar
```

We can then also output the energy from the MD code and check this matches the value output by PLUMED.  We ran a short trajectory to 
test that the energy is passed correctly.

# Trajectory

# Results

The table below contains the energies that were output by PLUMED and the energies that were ouptut by the MD code in its first two columns.  
The third column is the difference between these two numbers expressed as a percentage of $tolerance.  If the PLUMED interface is 
working correctly the first two columns of numbers should be identical and the third should be zero.
