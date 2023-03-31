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

We then run a second constant pressure MD simulation at a pressure of 1001 bar and with the following PLUMED restraint 
applied:

```plumed
v: VOLUME 
# slope should be just 10 times the Avogadro constant:
RESTRAINT AT=0.0 ARG=v SLOPE=-60.2214129
PRINT ARG=v FILE=volume2
```

The PLUMED restraint in this second calculation applies a negative pressure of 1000bar, which should compensate the fact that the 
second calculation was run at higher pressure.  The time series for the volumes that are output by the files `volume` and `volume2`
above should thus be close to identicial. 

# Trajectories

# Results

The table below contains the PLUMED outputs from the two calculations described above.
If the PLUMED interface is working correctly these two sets of numbers should be identical.
