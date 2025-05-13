import numpy as np
import MDAnalysis as mda
import subprocess

class mdcode :
   def setParams( self ) :
       params = {
         "temperature": 1.0,
         "tstep": 0.005,
         "relaxtime": 1.0,
         "pressure": 1.0,
         "prelaxtime": 4
       }
       return params
 
   def runMD( self, mdparams ) :
       return 1

   def getTimestep( self ) :
       return 0.005

   def getNumberOfAtoms( self, rundir ) :
       natoms = 1 
       return natoms 

   def getPositions( self, rundir ) :
       pos = np.zeros([1,3])
       return pos

   def getCell( self, rundir ) :
       cell = np.zeros([nframes,9])
       return cell

   def getMasses( self, rundir ) :
       return np.ones( natoms[0] )

   def getCharges( self, rundir ) :
       return np.zeros( natoms[0] )

   def getEnergy( self, rundir ) :
       return 1
