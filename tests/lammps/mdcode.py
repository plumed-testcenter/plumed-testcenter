import numpy as np
import MDAnalysis as mda
import subprocess

class mdcode :
   def getName( self ) :
       return "lammps"

   def setParams( self ) :
       params = {
         "temperature": 275,
         "itemperature": 275,
         "tstep": 2,
         "friction": 100,
         "pressure": 0.987,
         "pfriction": 400
       }
       return params
  
   def getExecutibleName( self ) :
       return "lammps"
 
   def runMD( self, mdparams ) :
       inp = "units           real\n"
       inp = inp + "variable        seed equal 74581\n"
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
       if mdparams["ensemble"]=="nvt" : inp = inp + "fix             1 all nvt temp  " + str(mdparams["temperature"]) + " " + str(mdparams["temperature"]) + " " + str(mdparams["friction"]) + " tchain 1\n"
       elif mdparams["ensemble"]=="npt" : inp = inp + "fix           1 all npt temp  " + str(mdparams["temperature"]) + " " + str(mdparams["temperature"]) + " " + str(mdparams["friction"]) + " iso " str(mdparams["pressure"]) + " " + str(mdparams["pressure"]) + " " + str(mdparams["pfriction"]) + " tchain 1 \n"
       inp = inp + "fix             2a ref setforce 0.0 0.0 0.0\n"
       # Code to deal with restraint 
       if "restraint" in mdparams and mdparams["restraint"]>0 : inp = inp + "fix 6 all restrain bond 1 2 10.0 10.0 " + str(10*mdparams["restraint"]) + "\n"
       inp = inp + "fix             4 all shake 0.0001 10 100 b 4 6 8 10 12 14 18 a 31\n"
       inp = inp + "thermo_style    custom step temp etotal pe ke epair ebond f_2\n"
       inp = inp + "thermo          10\n"
       inp = inp + "dump            dd all xyz 1 lammps.xyz\n"
       inp = inp + "variable        step equal step\n"
       inp = inp + "variable        pe equal pe\n"
       inp = inp + "fix             5 all print 1 \"$(v_step) $(v_pe)\" file lammps_energy\n"
       inp = inp + "dump            mq all custom 200 mq_lammps id mass q\n"
       inp = inp + "variable        ab equal cella\n"
       inp = inp + "variable        bb equal cellb\n"
       inp = inp + "variable        cb equal cellc\n"
       inp = inp + "fix             7 all print 1 \"$(v_ab) $(v_bb) $(v_cb)\" file lammps_cell\n"
       inp = inp + "velocity        all create " + str(mdparams["itemperature"]) + " ${seed} dist gaussian"
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

   def getMassCharge(self, rundir, col) :
       data = np.zeros( self.getNumberOfAtoms(rundir)[0] )
       f = open( rundir + "/mq_lammps", "r")
       inp, inmasses = f.read(), False
       f.close()
       for line in inp.splitlines() :
           if inmasses :
              words = line.split()
              data[int(words[0])-1] = float(words[col])
           elif "ATOMS id mass q" in line : inmasses = True
           
       return data

   def getCell( self, rundir ) :
       cell_size = np.loadtxt( rundir + "/lammps_cell")
       cell_final = np.zeros([cell_size.shape[0],9])
       for i in range(cell_size.shape[0]) : cell_final[i,0], cell_final[i,4], cell_final[i,8] = cell_size[i,0]/10, cell_size[i,1]/10, cell_size[i,2]/10
       return cell_final 

   def getMasses( self, rundir ) :
       return self.getMassCharge( rundir, 1 ) 

   def getCharges( self, rundir ) :
       return self.getMassCharge( rundir, 2 ) 

   def getEnergy( self, rundir ) :
       return  np.loadtxt(rundir + "/lammps_energy")[:,1]*4.184 
