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
       inp = "<simulation verbosity='high'>\n"
       inp = inp + "  <output prefix='tut1'>\n"
       inp = inp + "    <properties filename='md' stride='1'> [step, time{picosecond}, conserved{kelvin}, temperature{kelvin}, potential{j/mol}, kinetic_cv{j/mol}] </properties>\n"
       inp = inp + "    <properties filename='force' stride='20'> [atom_f{piconewton}(atom=0;bead=0)] </properties>\n"
       inp = inp + "    <trajectory filename='pos' stride='1' format='xyz' cell_units='angstrom'> positions{angstrom} </trajectory>\n"
       inp = inp + "    <checkpoint filename='checkpoint' stride='1000' overwrite='True'/>\n"
       inp = inp + "  </output>\n"
       inp = inp + "  <total_steps> " + str(mdparams["nsteps"]) + " </total_steps>\n"
       inp = inp + "  <ffsocket mode='inet' name='driver'>\n"
       inp = inp + "    <address>localhost</address>\n"
       inp = inp + "    <port> 31415 </port>\n"
       inp = inp + "  </ffsocket>\n"
       inp = inp + "  <ffplumed name='plumed'>\n"
       inp = inp + "   <file mode='xyz'> structure.xyz </file>\n"
       inp = inp + "   <plumeddat> plumed.dat </plumeddat>\n"
       inp = inp + "  </ffplumed>\n"
       inp = inp + "  <system>\n"
       inp = inp + "    <initialize nbeads='1'>\n"
       inp = inp + "      <file mode='xyz'> structure.xyz </file>\n"
       inp = inp + "      <velocities mode='thermal' units='kelvin'> " + str(mdparams["temperature"]) + " </velocities>\n"
       inp = inp + "    </initialize>\n"
       inp = inp + "    <forces>\n"
       inp = inp + "      <force forcefield='driver'/>\n"
       inp = inp + "    </forces>\n"
       inp = inp + "    <ensemble>\n"
       inp = inp + "      <temperature units='kelvin'> " + str(mdparams["temperature"]) + " </temperature>\n"
       if mdparams["ensemble"]=="npt" : inp = inp + "      <pressure units='bar'>" + str(mdparams["pressure"]) + "</pressure>\n"
       inp = inp + "      <bias>\n"
       inp = inp + "        <force forcefield='plumed' nbeads='1'></force>\n"
       inp = inp + "      </bias>\n"
       inp = inp + "    </ensemble>\n"
       inp = inp + "    <motion mode='dynamics'>\n"
       if mdparams["ensemble"]=="npt" : 
          inp = inp + "      <dynamics mode='npt'>\n"
          inp = inp + "        <barostat mode='isotropic'>\n"
          inp = inp + "        <thermostat mode='langevin'>\n"
          inp = inp + "          <tau units='femtosecond'> " + str(mdparams["prelaxtime"]) + " </tau>\n"
          inp = inp + "        </thermostat>\n"
          inp = inp + "        <tau units='femtosecond'> " + str(mdparams["prelaxtime"]) + " </tau>\n"
          inp = inp + "        </barostat>\n"
       else : inp = inp + "      <dynamics mode='nvt'>\n"
       inp = inp + "        <thermostat mode='pile_g'>\n"
       inp = inp + "          <tau units='femtosecond'> " + str(mdparams["relaxtime"]) + " </tau>\n"
       inp = inp + "        </thermostat>\n"
       inp = inp + "        <timestep units='femtosecond'> " + str(mdparams["tstep"]) + " </timestep>\n"
       inp = inp + "      </dynamics>\n"
       inp = inp + "    </motion>\n"
       inp = inp + "  </system>\n"
       inp = inp + "  <smotion mode='metad'>\n"
       inp = inp + "      <metad> <metaff> [ plumed ] </metaff> <use_energy> True </use_energy> </metad>\n"
       inp = inp + "  </smotion>\n"
       inp = inp + "</simulation>\n"
       of = open("input.xml","w+")
       of.write(inp)
       of.close()
       # Work out the script that will run ipi for us
       executible = mdparams["executible"] 
       # Now run the calculation using subprocess
       with open("stdout","w") as stdout:
        with open("stderr","w") as stderr:
           out = subprocess.run([executible], text=True, input=inp, stdout=stdout, stderr=stderr )
       # Just pause for a second after run to make sure there is time to write out all files
       time.sleep(5)
       return out.returncode

   def getTimestep( self ) :
       return 0.001  # Convert timestep in fs to ps

   def getNumberOfAtoms( self, rundir ) :
       natoms, traj, fnum = [], mda.coordinates.XYZ.XYZReader( rundir + "/tut1.pos_0.xyz"), 0
       for frame in traj.trajectory : 
           if fnum>0 : natoms.append( frame.positions.shape[0] )
           fnum = fnum + 1
       return natoms
       
   def getPositions( self, rundir ) :
       fnum, traj = 0, mda.coordinates.XYZ.XYZReader( rundir + "/tut1.pos_0.xyz")
       for frame in traj.trajectory :
          if fnum==1 : pos, first = frame.positions / 10, False
          elif fnum>0 : pos = np.concatenate( (pos, frame.positions / 10), axis=0 )
          fnum = fnum + 1
       return pos

   def getCell( self, rundir ) :
       nframes = len( self.getNumberOfAtoms( rundir ) )
       cell = np.zeros([nframes,9])
       f = open( rundir + "/tut1.pos_0.xyz", "r" )
       lines = f.read().splitlines()
       f.close() 
       natoms = int( lines[0] )
       for i in range(1,nframes) : 
           cellstr = lines[i*(2+natoms)+1].split()
           cell[i,0], cell[i,4], cell[i,8] = float(cellstr[2]) / 10, float(cellstr[3]) / 10, float(cellstr[4]) / 10
       return cell

   def getMasses( self, rundir ) :
       raise Exception("No function to get masses yet")

   def getCharges( self, rundir ) :
       raise Exception("No function to get charges yet")
  
   def getEnergy( self, rundir ) :
       return np.loadtxt( rundir + "/tut1.md" )[1:,4] / 1000
