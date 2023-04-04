import numpy as np
import subprocess

class mdcode :
   def getName( self ) :
       return "quantum_expresso"

   def setParams( self ) :
       params = {
         "temperature": 1.0,
         "tstep": 20,
         "friction": 1.0
       }
       return params

   def runMD( self, mdparams ) : 
       inp = " &control \n"
       inp = inp + "    calculation='md' \n"
       inp = inp + "    pseudo_dir='./' \n"
       inp = inp + "    dt=" + str(mdparams["tstep"]) + ", \n"
       inp = inp + "    nstep=" + str(mdparams["nsteps"]) + " \n"
       inp = inp + " / \n"
       inp = inp + " &system \n"
       inp = inp + "    ibrav= 2, celldm(1)=10.18, nat=  2, ntyp= 1, \n"
       inp = inp + "    ecutwfc = 8.0, nosym=.true. \n"
       inp = inp + " / \n"
       inp = inp + " &electrons \n"
       inp = inp + "    conv_thr =  1.0e-8 \n"
       inp = inp + "    mixing_beta = 0.7 \n"
       inp = inp + " / \n"
       inp = inp + " &ions \n"
       inp = inp + " / \n"
       inp = inp + "ATOMIC_SPECIES \n"
       inp = inp + " Si  28.086  Si.pz-vbc.UPF \n"
       inp = inp + "ATOMIC_POSITIONS {alat} \n"
       inp = inp + " Si -0.123 -0.123 -0.123 \n"
       inp = inp + " Si  0.123  0.123  0.123 \n"
       inp = inp + "K_POINTS {automatic} \n"
       inp = inp + " 1 1 1 0 0 0 \n"
       of = open("md.in","w+")
       of.write(inp)
       of.close()
       # Work out the name of the espresso executable
       executible = "pw_" + mdparams["version"]
       if mdparams["stable_version"] : executible = "pw"
       # Now run the calculation using subprocess
       with open("stdout","w") as stdout:
        with open("stderr","w") as stderr:
           out = subprocess.run([executible, "-plumed"], text=True, input=inp, stdout=stdout, stderr=stderr )
       return out.returncode

   def getTimestep( self ) :
       return 0.005*(6.62607015E-34/2*np.pi/4.3597447222071E-18)*1E12

   def getNumberOfAtoms( self, rundir ) :
       natoms = []
       return natoms
       
   def getPositions( self, rundir ) :
       pos = []
       return pos

   def getMasses( self, rundir ) :
       return []

   def getCharges( self, rundir ) :
       return [] 
