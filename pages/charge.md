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
that were output by the MD code in the first two columns.  The third column is the difference between the numbers in the first two columns expressed as a percentage of 
0.001.  If the PLUMED interface is working correctly the first two columns of numbers should be identical and the third column should be zero.

