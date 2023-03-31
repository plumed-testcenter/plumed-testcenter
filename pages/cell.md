Cell vectors are passed correctly
---------------------------------

PLUMED must receive the cell vectors from the MD code in order to calculate CVs correctly.  
To test these cell vectors are passed correctly to PLUMED we run a short trajectory and output the celll vectors 
that are passed to PLUMED using the following command: 

```plumed
c: CELL 
```

# Trajectory

# Results

The table below contains the cell vectors that were output by the above command and the cell vectors 
that were output by the MD code.  If the PLUMED interface is working correctly these two sets of numbers should be identical.
