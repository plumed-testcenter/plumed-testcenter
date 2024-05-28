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
that were output by the MD code in the first two columns.  The third column is the difference between the two sets of cell vectors reported as a percentage of 
0.001.  If the PLUMED interface is working correctly the first two sets of numbers should be identical and the final column should be filled with zeros.
