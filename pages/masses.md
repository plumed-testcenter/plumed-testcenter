Masses are passed correctly
---------------------------

PLUMED must receive the masses from an MD code in order to calculate many CVs correctly.
To test these masses are passed correctly to PLUMED we run a short trajectory and output the masses of all the atoms that 
are passed to PLUMED using the following command: 

```plumed
DUMPMASSCHARGE FILE=mq_plumed
```

# Trajectory

# Results 

The table below contains some of the masses that were output by the above command and the masses of the corresponding atoms 
that were output by the MD code.  If the PLUMED interface is working correctly these two sets of numbers should be identical.

