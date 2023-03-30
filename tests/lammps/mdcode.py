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
  
   def runMD( self, mdparams ) :
       of = open("in.peptide-plumed","w+")
       of.write("units           real\n")
       of.write("atom_style      full\n")
       of.write("pair_style      lj/charmm/coul/long 8.0 10.0 10.0\n")
       of.write("bond_style      harmonic\n")
       of.write("angle_style     charmm\n")
       of.write("dihedral_style  charmm\n")
       of.write("improper_style  harmonic\n")
       of.write("kspace_style    pppm 0.0001\n")
       of.write("read_data       data.peptide\n")
       of.write("neighbor        2.0 bin\n")
       of.write("neigh_modify    delay 5\n")
       of.write("timestep        " + str(mdparams["tstep"]) + "\n")
       of.write("group           peptide type <= 12\n")
       of.write("group           one id 2 4 5 6\n")
       of.write("group           two id 80 82 83 84\n")
       of.write("group           ref id 37\n")
       of.write("group           colvar union one two ref\n")
       of.write("fix             2 all plumed plumedfile plumed.dat outfile p.log\n")
       of.write("fix             1 all nvt temp  " + str(mdparams["temperature"]) + " " + str(mdparams["temperature"]) + " " + mdparams["friction"] + " tchain 1\n")
       of.write("fix             2a ref setforce 0.0 0.0 0.0\n")
       of.write("fix             4 all shake 0.0001 10 100 b 4 6 8 10 12 14 18 a 31\n")
       of.write("thermo_style    custom step temp etotal pe ke epair ebond f_2\n")
       of.write("thermo          10\n")
       of.write("dump            dd all xyz 10 lammps.xyz\n")
       of.write("variable        step equal step\n")
       of.write("variable        pe equal pe\n")
       of.write("fix             5 all print 10 "$(v_step) $(v_pe)" file lammps_energy\n")
       of.write("dump            mq all custom 200 mq_lammps id mass q\n")
       of.write("run             " + str(mdparams["nsteps"]) + "\n")
       of.close()

