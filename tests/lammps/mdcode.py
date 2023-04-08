import numpy as np
import MDAnalysis as mda
import subprocess

class mdcode :
   def getName( self ) :
       return "lammps"

   def setParams( self ) :
       params = {
         "temperature": 275,
         "tstep": 2,
         "friction": 100
       }
       return params
  
   def getExecutibleName( self ) :
       return "lammps"
 
   def runMD( self, mdparams ) :
       inp = "units           real\n"
       inp = inp + "atom_style      full\n"
       inp = inp + "pair_style      lj/charmm/coul/long 8.0 10.0 10.0\n"
       inp = inp + "bond_style      harmonic\n"
       inp = inp + "angle_style     charmm\n"
       inp = inp + "dihedral_style  charmm\n"
       inp = inp + "improper_style  harmonic\n"
       inp = inp + "kspace_style    pppm 0.0001\n"
       inp = inp + "read_data       data.peptide\n"
       inp = inp + "neighbor        2.0 bin\n"
       inp = inp + "neigh_modify    delay 5\n"
       inp = inp + "timestep        " + str(mdparams["tstep"]) + "\n"
       inp = inp + "group           peptide type <= 12\n"
       inp = inp + "group           one id 2 4 5 6\n"
       inp = inp + "group           two id 80 82 83 84\n"
       inp = inp + "group           ref id 37\n"
       inp = inp + "group           colvar union one two ref\n"
       inp = inp + "fix             2 all plumed plumedfile plumed.dat outfile p.log\n"
       inp = inp + "fix             1 all nvt temp  " + str(mdparams["temperature"]) + " " + str(mdparams["temperature"]) + " " + str(mdparams["friction"]) + " tchain 1\n"
       inp = inp + "fix             2a ref setforce 0.0 0.0 0.0\n"
       # Code to deal with restraint 
       if "restraint" in mdparams and mdparams["restraint"]>0 : inp = inp + "fix 6 all restrain bond 1 2 1000.0 1000.0 " + str(10*mdparams["restraint"]) + "\n"
       inp = inp + "fix             4 all shake 0.0001 10 100 b 4 6 8 10 12 14 18 a 31\n"
       inp = inp + "thermo_style    custom step temp etotal pe ke epair ebond f_2\n"
       inp = inp + "thermo          10\n"
       inp = inp + "dump            dd all xyz 1 lammps.xyz\n"
       inp = inp + "variable        step equal step\n"
       inp = inp + "variable        pe equal pe\n"
       inp = inp + "fix             5 all print 1 \"$(v_step) $(v_pe)\" file lammps_energy\n"
       inp = inp + "dump            mq all custom 200 mq_lammps id mass q\n"
       inp = inp + "run             " + str(mdparams["nsteps"]) + "\n"
       of = open("in.peptide-plumed","w+")
       of.write(inp)
       of.close()
       # Work out the name of the lammps executable
       executible = mdparams["executible"] 
       # Now run the calculation using subprocess
       with open("stdout","w") as stdout:
        with open("stderr","w") as stderr:
          out = subprocess.run([executible], text=True, input=inp, stdout=stdout, stderr=stderr )
       return out.returncode

   def getTimestep( self ) :
       return 0.002

   def getNumberOfAtoms( self, rundir ) :
       natoms, traj = [], mda.coordinates.XYZ.XYZReader( rundir + "/lammps.xyz")
       for frame in traj.trajectory : natoms.append( frame.positions.shape[0] )
       return natoms

   def getPositions( self, rundir ) :
       first, traj = True, mda.coordinates.XYZ.XYZReader( rundir + "/lammps.xyz")
       for frame in traj.trajectory :
          if first : pos, first = frame.positions / 10, False
          else : pos = np.concatenate( (pos, frame.positions / 10), axis=0 )
       return pos

   def getCell( self, rundir ) :
       return []

   def getMasses( self, rundir ) :
       return [] 
       #np.loadtxt( rundir + "/mq_lammps")[:,2]

   def getCharges( self, rundir ) :
       return [] 
       # np.loadtxt( rundir + "/mq_lammps")[:,3]

   def getEnergy( self, rundir ) :
       return [] 
