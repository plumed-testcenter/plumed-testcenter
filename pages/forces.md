Checking PLUMED can add a restraint
-----------------------------------

To check that the forces are being correctly passed from PLUMED to the MD code we test that the result from a PLUMED
calculation where the following plumed input is used produces results that are equivalent to the result that is obtained
when the restraint is applied within the MD code.

```plumed
dd: DISTANCE ATOMS=1,2 
RESTRAINT ARG=dd KAPPA=2000 AT=0.6
PRINT ARG=dd FILE=plumed_restraint FMT=%8.4f
```

To check that the two ways of applying the restraint are equivalent we use the following PLUMED input for the calculation
where the MD code applies the restraint:

```plumed
dd: DISTANCE ATOMS=1,2 
PRINT ARG=dd FILE=mdcode_restraint FMT=%8.4f
```

The test passes if the time series for the distance in the file `plumed_restraint` is the same as the time series in the file
`mdcode_restraint`.
