Energy is passed correctly
--------------------------

It is common practise to use the potential energy as a collective energy.  Some MD codes thus pass the potential energy to PLUMED. 
To check that this quantity has been passed correctly we can output the passed energy from PLUMED using the following input.  

```plumed
e: ENERGY 
PRINT ARG=e FILE=colvar
```

The python function `checkEnergy` that you can read in the file `test.py` then checks whether the potential energy that PLUMED has received is the same 
as the potential energy in the MD code.  

Notice that the command above outputs energies in kJ/mol and that many MD codes will output energies in a different unit.  If you are writing tests
you should thus ensure that energies output by the MD code are converted into kJ/mol in the `getEnergy` functions that you write.  
In other words all unit appropriate conversions must be done _before_ checks against the PLUMED energy.
