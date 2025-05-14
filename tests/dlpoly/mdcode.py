import numpy as np
import MDAnalysis as mda
import subprocess

class mdcode :
   def setParams( self ) :
       params = {
         "temperature": 300,
         "tstep": 0.001,
         "relaxtime": 0.5,
         "pressure": 0.0,
         "prelaxtime": 1 
       }
       return params
 
   def runMD( self, mdparams ) :
       ensemble_stuff = """
ensemble                      nvt
ensemble_method               hoover
ensemble_thermostat_coupling  {mdparams["relaxtime"} ps
"""
       if mdparams["ensemble"]=="npt" :
          ensemble_stuff = """
ensemble                      npt
ensemble_method               hoover
ensemble_thermostat_coupling  {mdparams["relaxtime"} ps
ensemble_barostat_coupling    {mdparams["pelaxtime"} ps
"""
       inp = """
title                 PLUMED test calc
io_file_config        CONFIG
io_file_field         FIELD
io_file_statis        STATIS
io_file_revive        REVIVE
io_file_revcon        REVCON
io_file_tabvdw        TABLE
temperature           {mdparams["temperature"]} K
plumed                on
plumed_input          plumed.dat
print_frequency       1 steps
stats_frequency       1 steps
vdw_cutoff            8.0 ang
padding               0.25 ang
cutoff                8.0 ang
coul_method           spme
spme_precision        1e-06
time_run              {mdparams["nsteps"]} steps
time_equilibration    0 steps
time_job              3600.0 s
time_close            100.0 s
timestep              {mdparams["tstep"]} ps
pressure_hydrostatic  0.0 katm
rescale_frequency     10 steps
{ensemble_stuff}
"""
       cf = open("CONTROL","w+")
       cf.write(inp)
       cf.close()
       # Work out the name of the dlpoly executable
       executible = mdparams["executible"]
       # Now run the calculation using subprocess
       with open("stdout","w") as stdout:
        with open("stderr","w") as stderr:
          out = subprocess.run([executible], text=True, input=inp, stdout=stdout, stderr=stderr )
       return out.returncode

   def getTimestep( self ) :
       return 0.001

   def getNumberOfAtoms( self, rundir ) :
       natoms, traj = [], mda.coordinates.DLPoly.HistoryReader( rundir + "/HISTORY" )
       for frame in traj.trajectory : natoms.append( frame.positions.shape[0] ) 
       return natoms 

   def getPositions( self, rundir ) :
       first, traj = True, mda.coordinates.DLPoly.HistoryReader( rundir + "/HISTORY" )
       for frame in traj.trajectory :
          if first : pos, first = 0.1*frame.positions.copy(), False
          else : pos = np.concatenate( (pos, 0.1*frame.positions), axis=0 ) 
       return pos

   def getCell( self, rundir ) :
       cell = np.zeros([nframes,9])
       return cell

   def getMasses( self, rundir ) :
       masses = np.zeros(4*64)
       for i in range(64) : 
           masses[4*i+0] = 107.868
           masses[4*i+1] = 126.905
           masses[4*i+2] = 107.868
           masses[4*i+3] = 126.905  
       return masses

   def getCharges( self, rundir ) :
       charges = np.zeros(4*64)
       for i in range(64) : 
           charges[4*i+0] = 0.6 
           charges[4*i+1] = -0.6 
           charges[4*i+2] = 0.6 
           charges[4*i+3] = -0.6  
       return charges

   def getEnergy( self, rundir ) :
       return 1
