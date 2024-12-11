Check virial contribution
-------------------------

If you are running simulations at constant pressure then the virial forces cause the cell parameters 
to change with time.  Any CVs calculated by PLUMED contribute to these virial forces and PLUMED must,
therefore, have a mechanism to pass virial forces back to the MD code. 

To debug this mechanism we run a constant pressure simulation at 1 bar using the MD code.  During this simulation
we use the following PLUMED input to monitor the cell volume:

```plumed
v: VOLUME
PRINT ARG=v FILE=volume
```

We then run a second constant pressure MD simulation at a pressure of 1001 bar and the input above.

If the virial has been implemented correctly within PLUMED the following PLUMED restraint will apply a negative pressure of 1000bar, which should compensate the fact that the
second calculation was run at higher pressure.  We thus run a third MD calculation with the following input file:

```plumed
v: VOLUME 
# slope should be just 10 times the Avogadro constant:
RESTRAINT AT=0.0 ARG=v SLOPE=-60.2214129
PRINT ARG=v FILE=volume2
```

The time series for the volumes that are output by the files `volume` and `volume2` above should thus be close to identical. 

# Trajectories

{Trajectories}

# Results

{Results}

The table below includes some of the results from the calculation.  The columns contain:

1. The time series for the volume that was obtained from the simulation in that was performed at 1 bar, $x_{{md}}$.
2. The time series for the volume that was obtained from the simulation that was performed at 1001 bar and in which PLUMED applied a restraint on the volume, $x_{{pl}}$.
3. The absolute value of the difference between the time series of volumes that were obtained from the simulations running at 1001 bar and 1 bar, $\vert x_{{md}}'-x_{{md}}\vert$.  No PLUMED restraints were applied in either of these simulations.
4. The values of $100\frac{{\vert x_{{md}} - x_{{pl}}\vert }}{{ \vert x_{{md}}'-x_{{md}} \vert }}$.

If the PLUMED interface is working correctly the first two sets of numbers should be identical and the final column should be filled with zeros.
