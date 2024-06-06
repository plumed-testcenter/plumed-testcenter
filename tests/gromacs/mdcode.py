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
         "temperature": 25,
         "tstep": 1,
         "relaxtime": 25,
         "pressure": 1,
         "prelaxtime": 250
       }
       return params

   def runMD( self, mdparams ) : 
       # Work out the script that will run ipi for us
       executible = mdparams["executible"] 
       # Now run the calculation using subprocess
       with open("stdout","w") as stdout:
        with open("stderr","w") as stderr:
           out = subprocess.run([executible], text=True, input=inp, stdout=stdout, stderr=stderr )
       return out.returncode

   def getTimestep( self ) :
       return 0.001  # Convert timestep in fs to ps

   def getNumberOfAtoms( self, rundir ) :
       natoms = 0
       return natoms
       
   def getPositions( self, rundir ) :
       pos=0
       return pos

   def getCell( self, rundir ) :
       cell=0
       return cell

   def getMasses( self, rundir ) :
       return 0

   def getCharges( self, rundir ) :
       return 0
  
   def getEnergy( self, rundir ) :
       return 0
