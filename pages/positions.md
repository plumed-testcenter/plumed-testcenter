Positions are passed correctly
------------------------------

PLUMED must receive the positions from an MD code in order to calculate CVs correctly.  
To test these positions are passed correctly to PLUMED we run a short trajectory and output the positions of all the atoms 
that are passed to PLUMED using the following command: 

```plumed
DUMPATOMS ATOMS=@mdatoms PRECISION=4 FILE=plumed.xyz
```

# Trajectory

# Results

The table below contains some of the positions that were output by the above command and the positions of the corresponding atoms 
that were output by the MD code.  If the PLUMED interface is working correctly these two sets of numbers should be identical.
