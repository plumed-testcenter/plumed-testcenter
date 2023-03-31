Charges are passed correctly
---------------------------

PLUMED must receive the charges from an MD code in order to calculate some CVs correctly.
To test these charges are passed correctly to PLUMED we run a short trajectory and output the charges of all the atoms that 
are passed to PLUMED using the following command: 

```plumed
DUMPMASSCHARGE FILE=mq_plumed
```

# Trajectory

# Results 

The table below contains some of the charges that were output by the above command and the charges of the corresponding atoms 
that were output by the MD code.  If the PLUMED interface is working correctly these two sets of numbers should be identical.

