Checking PLUMED can add a restraint
-----------------------------------

To check that the forces are being correctly passed from PLUMED to the MD code we test that the result from a PLUMED
calculation where the following plumed input is used produces results that are equivalent to the result that is obtained
when the restraint is applied within the MD code.

```plumed
dd: DISTANCE ATOMS=1,2 
RESTRAINT ARG=dd KAPPA=2000 AT=0.6
PRINT ARG=dd FILE=plumed_restraint
```

To check that the two ways of applying the restraint are equivalent we use the following PLUMED input for the calculation
where the MD code applies the restraint:

```plumed
dd: DISTANCE ATOMS=1,2 
PRINT ARG=dd FILE=mdcode_restraint
```

# Trajectories

# Results

The table below includes some of the results from the calculation.  The columns contain:

1. The time series of distance values that were obtained when the restraint was applied by the MD code, $x_{md}$.
2. The time series of distance values that were obtained when the restriant was applied by PLUMED, $x_{pl}$.
3. The tolerances that were used when comparing these quanitites, $\delta$.
4. The values of $100\frac{|x_{md} - x_{pl}| }{ \delta}$.

If the PLUMED interface is working correctly the first two sets of numbers should be identical and the final column should be filled with zeros.
