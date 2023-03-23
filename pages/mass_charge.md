Masses and charges are passed correctly
---------------------------------------

PLUMED should receive the masses and charges from an MD code in order to calculate centers of mass and dipoles correctly.  If we are interfacing PLUMED with an 
MD code we should thus test that the masses and charges that are passed to PLUMED are the same as those in the MD code 

We can output the masses and charges from PLUMED using the following input. 

```plumed
DUMPMASSCHARGE FILE=mq_plumed
```

The python functions `checkMasses` and `checkCharges` that you can read in the file `test.py` then check the masses and charges that are output by PLUMED against those that are output 
by the MD code.  
