import numpy as np
import MDAnalysis as mda
import subprocess

class mdcode :
   def getName( self ) :
       return "simplemd"

   def setParams( self ) :
       params = {
         "temperature": 1.0,
         "tstep": 0.005,
         "friction": 1.0
       }
       return params
  
   def runMD( self, mdparams ) :
       # Prepare a string that contains the input for simplemd
       inp = "inputfile input.xyz\n"
       inp = inp + "outputfile output.xyz\n"
       inp = inp + "temperature " + str(mdparams["temperature"]) + "\n"
       inp = inp + "tstep " + str(mdparams["tstep"]) + "\n"
       inp = inp + "friction " + str(mdparams["friction"]) + "\n"
       inp = inp + "forcecutoff 2.5\n"
       inp = inp + "listcutoff  3.0\n"
       inp = inp + "nstep " + str(mdparams["nsteps"]) + "\n"
       inp = inp + "nconfig 1 trajectory.xyz\n"
       inp = inp + "nstat   1 energies.dat\n"
       of=open("in","w+")
       of.write(inp)
       of.close()
       # Work out the name of the plumed executable 
       executible = "plumed_" + mdparams["version"] 
       if mdparams["stable_version"] : executible = "plumed"
       cmd = [executible, "simplemd"]
       # Now run the calculation using subprocess
       with open("stdout","w") as stdout:
        with open("stderr","w") as stderr:
           plumed_out = subprocess.run(cmd, text=True, input=inp, stdout=stdout, stderr=stderr )
       return plumed_out.returncode

   def getTimestep( self ) :
       return 0.005

   def getNumberOfAtoms( self, rundir ) :
       natoms, traj = [], mda.coordinates.XYZ.XYZReader( rundir + "/trajectory.xyz") 
       for frame in traj.trajectory : natoms.append( frame.positions.shape[0] )
       return natoms 

   def getPositions( self, rundir ) :
       first, traj = True, mda.coordinates.XYZ.XYZReader( rundir + "/trajectory.xyz") 
       for frame in traj.trajectory :
          if first : pos, first = frame.positions, False
          else : pos = np.concatenate( (pos, frame.positions), axis=0 )
       return pos

   def getMasses( self, rundir ) :
       natoms = self.getNumberOfAtoms( rundir )
       return np.ones( natoms[0] )

   def getCharges( self, rundir ) :
       natoms = self.getNumberOfAtoms( rundir )
       return np.zeros( natoms[0] )
