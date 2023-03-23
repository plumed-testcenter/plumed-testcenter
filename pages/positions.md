Positions are passed correctly
------------------------------

PLUMED must receive the positions from an MD code in order to calculate CVs correctly.  The first test that is performed 
on many MD codes is to check the positions that are output by the MD code against those that are output by PLUMED.  

We can output the positions from PLUMED using the following input.  This command will output the positions and the cell parameters.

```plumed
DUMPATOMS ATOMS=@mdatoms FILE=plumed.xyz
```

The python functions `checkPositions` and `checkCell` that you can read in the file `test.py` then check the positions and cell parameters that are output by PLUMED against those that are output 
by the MD code.  

Notice that the command above outputs positions in nm and that many MD codes will output positions in a different unit.  If you are writing tests
you should thus ensure that positions output by the MD code are converted into nanometers in the `getPositions` functions that you write.  Similarly, cell parameters returned from the function 
`getCellParameters` should be given in nm and in a  three by three matrix.  In other words all unit appropriate conversions must be done _before_ checks against the PLUMED positions and cell parameters.
