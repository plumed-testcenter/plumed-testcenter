import os
import sys
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
          # The box mass of Parrinello-Rahman goes as compressibility/tau_p^2: the energy tests shorten tau_p
          # together with the timestep, so the compressibility follows tau_p^2 to keep the same box mass
          compressibility = 4.46e-5 * (mdparams["prelaxtime"] / self.setParams()["prelaxtime"])**2
          npt_stuff = f"""
pcoupl                   = Parrinello-Rahman
tau_p                    = {mdparams["prelaxtime"]}
compressibility          = {compressibility}
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
ld-seed                  = 4321

{npt_stuff}
"""
       # The energy tests compare runs at temperatures T and alpha*T, which must start from velocities scaled by
       # sqrt(alpha), and have the same thermostat noise (ld-seed)
       t0 = self.setParams()["temperature"]
       if mdparams["temperature"] != t0 :
          self.scaleVelocities( np.sqrt( mdparams["temperature"] / t0 ) )
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

   def scaleVelocities( self, factor ) :
       # conf.gro, with the velocities multiplied by factor, written with more decimals so that the scaling is
       # exact (GROMACS takes the precision from the distance between the decimal points of the positions)
       with open("conf.gro") as f : lines = f.readlines()
       natoms = int(lines[1])
       for i in range(2, 2 + natoms) :
          line = lines[i]
          x = [float(line[20 + 8*k:28 + 8*k]) for k in range(3)]
          v = [factor * float(line[44 + 8*k:52 + 8*k]) for k in range(3)]
          lines[i] = line[:20] + "".join(f"{a:12.7f}" for a in x) + "".join(f"{a:12.8f}" for a in v) + "\n"
       with open("conf.gro", "w") as f : f.writelines(lines)

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
           masses[4*k+0] = 16.00000
           masses[4*k+1] = 1.00800
           masses[4*k+2] = 1.00800
           masses[4*k+3] = 0.00000
           k = k + 1
       return masses

   def getCharges( self, rundir ) :
       natoms = self.getNumberOfAtoms( rundir )
       nmols, charges, k = int(natoms[0] / 4), np.zeros(natoms[0]), 0
       for i in range(nmols) : 
           charges[4*k+0] = 0.0
           charges[4*k+1] = 0.52422
           charges[4*k+2] = 0.52422
           charges[4*k+3] = -1.04844
           k = k + 1
       return charges
  
   def getEnergy( self, rundir ) :
       eng = np.loadtxt( rundir + "/energy.xvg", comments=["#","@"] )
       return eng[:,1]

   def runRegtests( self, executible, rundir ) :
       # The regression tests in regtests/: executible runs gmx_mpi with its arguments, with the PLUMED kernel
       # of the version being tested
       sys.path.insert( 0, os.path.join( os.path.dirname( os.path.abspath(__file__) ), "regtests" ) )
       import regtest
       return regtest.run( executible, None, rundir )

