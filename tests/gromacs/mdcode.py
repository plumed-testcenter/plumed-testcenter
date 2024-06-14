import os
import time
import numpy as np
import MDAnalysis as mda
import subprocess

class mdcode :
   def __init__( self ) :
       AngToNm = 0.1

   def setParams( self ) :
       params = {
         "temperature": 300,
         "tstep": 0.002,
         "relaxtime": 2.0,
         "pressure": 1,
         "prelaxtime": 2.0
       }
       return params

   def runMD( self, mdparams ) : 
       npt_stuff = ""
       if mdparams["ensemble"]=="npt" : 
          npt_stuff = f"""
pcoupl                   = Parrinello-Rahman
tau_p                    = {mdparams["prelaxtime"]}
compressibility          = 4.46e-5
ref_p                    = {mdparams["pressure"]}
"""

       inp = f"""
integrator               = md        
dt                       = {mdparams["tstep"]}     
nsteps                   = {mdparams["nsteps"]} 

nstenergy                = 1
nstlog                   = 1
nstxout-compressed       = 1
compressed-x-precision   = 100000

continuation             = yes
constraint-algorithm     = lincs
constraints              = h-bonds

cutoff-scheme            = Verlet

coulombtype              = PME
rcoulomb                 = 1.0

vdwtype                  = Cut-off
rvdw                     = 1.0
DispCorr                 = EnerPres

tcoupl                   = v-rescale
tc-grps                  = System
tau-t                    = {mdparams["relaxtime"]}
ref-t                    = {mdparams["temperature"]}

{npt_stuff}
"""
       of = open("md.mdp", "w+")
       of.write(inp)
       of.close()
       # Work out the script that will run ipi for us
       executible = mdparams["executible"] 
       # Now run the calculation using subprocess
       with open("stdout","w") as stdout:
        with open("stderr","w") as stderr:
           out = subprocess.run([executible], text=True, input=inp, stdout=stdout, stderr=stderr )
       return out.returncode

   def getTimestep( self ) :
       return 0.002

   def getNumberOfAtoms( self, rundir ) :
       natoms = [] 
       with mda.coordinates.XTC.XTCFile( rundir + "/traj_comp.xtc") as xtc : 
         for frame in xtc : natoms.append( xtc.n_atoms )
       return natoms
       
   def getPositions( self, rundir ) :
       fnum = 0
       with mda.coordinates.XTC.XTCFile( rundir + "/traj_comp.xtc") as xtc :
         for frame in xtc :
             if fnum==0 : pos = frame.x
             else : pos = np.concatenate( (pos, frame.x), axis=0 )
             fnum = fnum + 1
       return pos

   def getCell( self, rundir ) :
       first, traj = True, mda.coordinates.XTC.XTCReader( rundir + "/traj_comp.xtc") 
       for frame in traj : 
           dim = frame.dimensions
           if first : first, cell = False, np.array( [[dim[0]/10,0,0,0,dim[1]/10,0,0,0,dim[2]/10]] )
           else :
              box = np.array( [[dim[0]/10,0,0,0,dim[1]/10,0,0,0,dim[2]/10]] )
              cell = np.concatenate( (cell, box), axis=0 )  
       return cell

   def getMasses( self, rundir ) :
       natoms = self.getNumberOfAtoms( rundir )
       nmols, masses, k = int(natoms[0] / 4), np.zeros(natoms[0]), 0
       for i in range(nmols) : 
           masses[k+0] = 16.00000
           masses[k+1] = 1.00800
           masses[k+2] = 1.00800
           masses[k+3] = 0.00000
       return masses

   def getCharges( self, rundir ) :
       natoms = self.getNumberOfAtoms( rundir )
       nmols, charges, k = int(natoms[0] / 4), np.zeros(natoms[0]), 0
       for i in range(nmols) : 
           charges[k+0] = 0.0
           charges[k+1] = 0.52422
           charges[k+2] = 0.52422
           charges[k+3] = -1.04844
       return charges
  
   def getEnergy( self, rundir ) :
       eng = np.loadtxt( rundir + "/energy.xvg", comments=["#","@"] )
       return eng[:,1]
