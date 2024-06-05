import os
import numpy as np
import xml.etree.ElementTree as ET
import subprocess

class mdcode :
   def __init__( self ) :
       self.bohrToNm = 0.0529177249
       self.HaToKJ = 2*1312.75  # I think the output in the xml file is in Hartrees and not Rydbergs

   def setParams( self ) :
       params = {
         "temperature": 1.0,
         "tstep": 20,
         "relaxtime": 10,
         "pressure": 0.001,
         "prelaxtime": 4.0
       }
       return params

   def runMD( self, mdparams ) :
       
       calculation='md'
       cell_dynamics=""
       if mdparams["ensemble"]=="nvt" :
           calculation='md'
       elif mdparams["ensemble"]=="npt" :
           calculation='vc-md'
           cell_dynamics=f"""
 &cell
     cell_dynamics = 'w',
     press = {mdparams["pressure"]}
 /"""
       inp=f""" &control
    calculation='{calculation}'
    pseudo_dir='./'
    dt={mdparams["tstep"]},
    nstep={mdparams["nsteps"]}
 /
 &system
    ibrav= 2, celldm(1)=10.18, nat=  2, ntyp= 1,
    ecutwfc = 8.0, nosym=.true.
    { "" if "restraint" in mdparams and mdparams["restraint"]>0 else ""}
 /
 &electrons 
    conv_thr =  1.0e-8 
    mixing_beta = 0.7 
 / 
 &ions 
    ion_temperature = 'berendsen' 
    tempw ={mdparams["temperature"]}
    nraise = 10
 / {cell_dynamics}
 ATOMIC_SPECIES
  Si  28.086  Si.pz-vbc.UPF
 ATOMIC_POSITIONS {{alat}}
  Si -0.123 -0.123 -0.123
  Si  0.123  0.123  0.123
 K_POINTS {{automatic}}
  1 1 1 0 0 0
"""
       of = open("md.in","w+")
       of.write(inp)
       of.close()
       # Make sure the PLUMED libarary is available in the LD_LIBRARY_PATH
       os.environ["LD_LIBRARY_PATH"] = os.environ["LD_LIBRARY_PATH"] + ":$(HOME)/opt/lib/" 
       # Work out the name of the espresso executable
       executible = mdparams["executible"] 
       # Now run the calculation using subprocess
       with open("stdout","w") as stdout:
        with open("stderr","w") as stderr:
           out = subprocess.run([executible, "-plumed"], text=True, input=inp, stdout=stdout, stderr=stderr )
       return out.returncode

   def getTimestep( self ) :
       return 20*4.8378E-5

   def getNumberOfAtoms( self, rundir ) :
       natoms, root = [], ET.parse( rundir + "/pwscf.xml")
       for step in root.findall("step") :
           adict = step.find("atomic_structure").attrib
           natoms.append( int(adict['nat']) )
       return natoms
       
   def getPositions( self, rundir ) :
       natoms = self.getNumberOfAtoms( rundir )
       pos = np.zeros( [sum(natoms), 3] )

       n, root = 0, ET.parse( rundir + "/pwscf.xml")
       for step in root.findall("step") :
           struct = step.find("atomic_structure")
           apos = struct.find("atomic_positions")
           allatoms = apos.findall("atom")
           for atom in allatoms : 
               strpos = atom.text.split()
               pos[n][0], pos[n][1], pos[n][2] = self.bohrToNm*float(strpos[0]), self.bohrToNm*float(strpos[1]), self.bohrToNm*float(strpos[2])
               n = n + 1 
       return pos

   def getCell( self, rundir ) :
       natoms = self.getNumberOfAtoms( rundir )
       cell = np.zeros( [len(natoms), 9] )
       n, root = 0 , ET.parse( rundir + "/pwscf.xml")
       for step in root.findall("step") :
           struct = step.find("atomic_structure")
           cellf = struct.find("cell")
           astr = cellf.find("a1").text.split()
           cell[n][0], cell[n][1], cell[n][2] = self.bohrToNm*float(astr[0]), self.bohrToNm*float(astr[1]), self.bohrToNm*float(astr[2]) 
           astr = cellf.find("a2").text.split()
           cell[n][3], cell[n][4], cell[n][5] = self.bohrToNm*float(astr[0]), self.bohrToNm*float(astr[1]), self.bohrToNm*float(astr[2])
           astr = cellf.find("a3").text.split()
           cell[n][6], cell[n][7], cell[n][8] = self.bohrToNm*float(astr[0]), self.bohrToNm*float(astr[1]), self.bohrToNm*float(astr[2])
           n = n + 1
       return cell

   def getMasses( self, rundir ) :
       root = ET.parse( rundir + "/pwscf.xml")
       species, inpt = {}, root.find("input")
       for spec in inpt.find("atomic_species").findall("species") : 
           species[spec.attrib["name"]] = float( spec.find("mass").text )
 
       masses = []
       for atom in inpt.find("atomic_structure").find("atomic_positions").findall("atom") : 
           masses.append( species[ atom.attrib["name"] ] )
       return masses

   def getCharges( self, rundir ) :
       natoms = self.getNumberOfAtoms( rundir )
       return np.zeros( natoms[0] )
  
   def getEnergy( self, rundir ) :
       energies, root = [], ET.parse( rundir + "/pwscf.xml")
       for step in root.findall("step") :
           ebrac = step.find("total_energy")
           energies.append( self.HaToKJ*float(ebrac.find("etot").text) )
       return energies
   