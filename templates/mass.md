Masses are passed correctly
---------------------------

PLUMED must receive the masses from an MD code in order to calculate many CVs correctly.
To test that masses are passed correctly to PLUMED we run a short trajectory and output the masses of all the atoms that 
are passed to PLUMED using the following command: 

```plumed
DUMPMASSCHARGE FILE=mq_plumed
```

# Trajectory

# Results 

The table below includes some of the results from the calculation.  The columns contain:

1. The masses that were obtained from the MD code, $x_{{md}}$.
2. The masses that were obtained from PLUMED, $x_{{pl}}$.
3. The tolerances that were used when comparing these quantities, $\delta$.
4. The values of $100\frac{{\vert x_{{md}} - x_{{pl}}\vert }}{{\delta}}$.

If the PLUMED interface is working correctly the first two sets of numbers should be identical and the final column should be filled with zeros.

