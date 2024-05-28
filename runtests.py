import os
import sys
import yaml
import getopt
import shutil
import subprocess
import numpy as np
import importlib
import MDAnalysis as mda
from datetime import date
from contextlib import contextmanager
from PlumedToHTML import test_plumed, get_html

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def processMarkdown( filename ) :
    if not os.path.exists(filename) :
       raise RuntimeError("Found no file called " + filename )
    f = open( filename, "r" )
    inp = f.read()
    f.close()

    ofile, inplumed, plumed_inp, ninputs = open( filename, "w+" ), False, "", 0
    for line in inp.splitlines() :
        # Detect and copy plumed input files 
        if "```plumed" in line :
           inplumed, plumed_inp, ninputs = True, "", ninputs + 1
        # Test plumed input files that have been found in tutorial 
        elif inplumed and "```" in line :
           inplumed = False
           solutionfile = "working" + str(ninputs) + ".dat"
           sf = open( solutionfile, "w+" )
           sf.write( plumed_inp )
           sf.close()
           # Test whether the input solution can be parsed
           success = success=test_plumed( "plumed", solutionfile )
           success_master=test_plumed( "plumed_master", solutionfile, printjson=True  )
           # Find the stable version 
           stable_version=subprocess.check_output('plumed info --version', shell=True).decode('utf-8').strip()
           # Use PlumedToHTML to create the input with all the bells and whistles
           html = get_html( plumed_inp, solutionfile, solutionfile, ("v"+ stable_version,"master"), (success,success_master), ("plumed","plumed_master"), usejson=(not success_master) )
           # Print the html for the solution
           ofile.write( "{% raw %}\n" + html + "\n {% endraw %} \n" )
        elif inplumed :
             if "__FILL__" in line : raise RuntimeError("Should not be incomplete files in this page")
             plumed_inp += line + "\n"
        # Just copy any line that isn't part of a plumed input
        elif not inplumed : ofile.write( line + "\n" )
    ofile.close()

def buildTestPages( directory ) :
   for page in os.listdir(directory) :
       if ".md" in page : processMarkdown( directory + "/" + page )

def runMDCalc( name, code, version, runner, params ) :
    # Get the name of the executible
    stram=open("tests/" + code + "/info.yml", "r") 
    params["executible"]=yaml.load(stram,Loader=yaml.BaseLoader)["executible"] + "_" + version
    stram.close()
    # Now test that the executable exists if it doesn't then the test is broken
    if shutil.which(params["executible"]) == None : return True
    # Copy all the input needed for the MD calculation
    shutil.copytree("tests/" + code + "/input", "tests/" + code + "/" + name + "_" + version ) 
    # Change to the directory to run the calculation
    with cd("tests/" + code + "/" + name + "_" + version ) :
       # Output the plumed file  
       of = open("plumed.dat","w+")
       of.write(params["plumed"])
       of.close() 
       # Now run the MD calculation
       broken = runner.runMD( params )
    # Make a zip archive that contains the input and output
    shutil.make_archive("tests/" + code + "/" + name + "_" + version, 'zip', "tests/" + code + "/" + name + "_" + version )
    return broken

def runTests(code,version,runner) :
   # Read in the information on the tests that should be run for this code
   stram=open("tests/" + code + "/info.yml", "r")
   info=yaml.load(stram,Loader=yaml.BaseLoader)["tests"]
   stram.close()

   fname, usestable = "testout.md", version=="stable"
   if version=="master" : fname = "testout_" + version + ".md"
   elif version!="stable" : ValueError("version should be master or stable")

   of = open("tests/" + code + "/" + fname, "w+")
   of.write("Testing " + code + "\n")
   of.write("------------------------\n \n")
   stable_version=subprocess.check_output('plumed info --version', shell=True).decode('utf-8').strip()
   of.write("The tests described in the following table were performed on __" + date.today().strftime("%B %d, %Y") + "__ to test whether the interface between " + code + " and ")
   if version=="stable" : version = "v" + stable_version
   of.write("the " + version + " version of PLUMED is working correctly.\n\n") 
   if info["virial"]=="no" : of.write("WARNING: " + code + " does not pass the virial to PLUMED and it is thus not possible to run NPT simulations with this code\n\n")

   params, basic_md_failed = runner.setParams(), True
   if info["positions"]=="yes" or info["timestep"]=="yes" or info["mass"]=="yes" or info["charge"]=="yes" : 
      params["plumed"] = "DUMPATOMS ATOMS=@mdatoms FILE=plumed.xyz\n"
      params["plumed"] = params["plumed"] + "c: CELL \n PRINT ARG=c.* FILE=cell_data\n"
      if info["mass"]=="yes" and info["charge"]=="yes" : params["plumed"] = params["plumed"] + "DUMPMASSCHARGE FILE=mq_plumed\n"
      elif info["mass"]=="yes" : params["plumed"] = params["plumed"] + "DUMPMASSCHARGE FILE=mq_plumed ONLY_MASSES\n"
      if info["timestep"]=="yes" : params["plumed"] = params["plumed"] + "t1: TIME\nPRINT ARG=t1 FILE=colvar\n"
      params["nsteps"], params["ensemble"] = 10, "npt"
      basic_md_failed = runMDCalc( "basic", code, version, runner, params )
  
   val1, val2 = 0.1, 0.1  
   of.write("| Description of test | Status | \n")
   of.write("|:--------------------|:------:| \n")
   if info["positions"]=="yes" :
      plumednatoms, codenatoms, codepos, plumedpos, codecell, plumedcell = [], [], np.ones(10), np.ones(10), np.ones(10), np.ones(10)
      if not basic_md_failed :
         # Get the trajectory that was output by PLUMED
         if os.path.exists("tests/" + code + "/basic_" + version + "/plumed.xyz") : 
            plumedtraj = mda.coordinates.XYZ.XYZReader("tests/" + code + "/basic_" + version + "/plumed.xyz")
            # Get the number of atoms in each frame from plumed trajectory
            plumednatoms, codenatoms = [], runner.getNumberOfAtoms( "tests/" + code + "/basic_" + version )
            for frame in plumedtraj.trajectory : plumednatoms.append( frame.positions.shape[0] )
            # Concatenate all the trajectory frames
            codepos, first = runner.getPositions( "tests/" + code + "/basic_" + version ), True
            for frame in plumedtraj.trajectory :
                if first : plumedpos, first = frame.positions, False
                else : plumedpos = np.concatenate( (plumedpos, frame.positions), axis=0 )
            codecell, plumedcell = runner.getCell( "tests/" + code + "/basic_" + version ), np.loadtxt("tests/" + code + "/basic_" + version + "/cell_data")[:,1:]
         else : basic_md_failed = True

      # Output results from tests on natoms
      writeReportPage( "natoms", code, version, basic_md_failed, ["basic"], codenatoms, plumednatoms, 0.01*np.ones(len(codenatoms)) ) 
      of.write("| MD code number of atoms passed correctly | " + getBadge( check(basic_md_failed, np.array(codenatoms), np.array(plumednatoms), 0.01*np.ones(len(codenatoms)) ), "natoms", code, version) + "| \n") 
      # Output results from tests on positions
      writeReportPage( "positions", code, version, basic_md_failed, ["basic"], codepos, plumedpos, 0.001*np.ones(plumedpos.shape) )
      of.write("| MD code positions passed correctly | " + getBadge( check(basic_md_failed, np.array(codepos), plumedpos, 0.001*np.ones(plumedpos.shape) ), "positions", code, version) + "| \n")
      writeReportPage( "cell", code, version, basic_md_failed, ["basic"], codecell, plumedcell, 0.001*np.ones(plumedcell.shape) )
      of.write("| MD code cell vectors passed correctly | " + getBadge( check(basic_md_failed, np.array(codecell), plumedcell, 0.001*np.ones(plumedcell.shape) ), "cell", code, version) + " | \n")
   if info["timestep"]=="yes" :
      md_tstep, plumed_tstep = 0.1, 0.1
      if not basic_md_failed :
         plumedtimes = np.loadtxt("tests/" + code + "/basic_" + version + "/colvar")[:,1]
         md_tstep, plumed_tstep = runner.getTimestep(), plumedtimes[1]-plumedtimes[0]
         for i in range(1,len(plumedtimes)) : 
             if plumedtimes[i]-plumedtimes[i-1]!=plumed_tstep : ValueError("Timestep should be the same for all MD steps")
      writeReportPage( "timestep", code, version, basic_md_failed, ["basic"], md_tstep, plumed_tstep, 0.0001 )
      of.write("| MD timestep passed correctly | " + getBadge( check(basic_md_failed, md_tstep, plumed_tstep, 0.0001), "timestep", code, version) + " | \n")
   if info["mass"]=="yes" : 
      md_masses, pl_masses = np.ones(10), np.ones(10)
      if not basic_md_failed : md_masses, pl_masses = runner.getMasses("tests/" + code + "/basic_" + version), np.loadtxt("tests/" + code + "/basic_" + version + "/mq_plumed")[:,1]
      writeReportPage( "mass", code, version, basic_md_failed, ["basic"], md_masses, pl_masses, 0.001*np.ones(pl_masses.shape) ) 
      of.write("| MD code masses passed correctly | " + getBadge( check( basic_md_failed, np.array(md_masses), pl_masses, 0.001*np.ones(pl_masses.shape) ), "mass", code, version) + " | \n")
   if info["charge"]=="yes" :
      md_charges, pl_charges = np.ones(10), np.ones(10)
      if not basic_md_failed : md_charges, pl_charges = runner.getCharges("tests/" + code + "/basic_" + version), np.loadtxt("tests/" + code + "/basic_" + version + "/mq_plumed")[:,2]
      writeReportPage( "charge", code, version, basic_md_failed, ["basic"], md_charges, pl_charges, 0.001*np.ones(pl_charges.shape) ) 
      of.write("| MD code charges passed correctly | " + getBadge( check( basic_md_failed, np.array(md_charges), pl_charges, 0.001*np.ones(pl_charges.shape) ), "charge", code, version) + " | \n")
   if info["forces"]=="yes" :
      # First run a calculation to find the reference distance between atom 1 and 2
      rparams = runner.setParams()
      rparams["nsteps"], rparams["ensemble"] = 2, "nvt"
      rparams["plumed"] = "dd: DISTANCE ATOMS=1,2 \nPRINT ARG=dd FILE=colvar"
      refrun, mdrun, plrun = runMDCalc("refres", code, version, runner, rparams ), True, True
      if not refrun : 
         # Get the reference distance betwene the atoms
         refdist = np.loadtxt("tests/" + code + "/refres_" + version + "/colvar")[0,1]
         # Run the calculation with the restraint applied by the MD code
         rparams["nsteps"], rparams["ensemble"] = 20, "nvt"
         rparams["restraint"] = refdist
         rparams["plumed"] = "dd: DISTANCE ATOMS=1,2 \nPRINT ARG=dd FILE=colvar"
         mdrun = runMDCalc("forces1", code, version, runner, rparams )
         # Run the calculation with the restraint applied by PLUMED
         rparams["restraint"] = -10
         rparams["plumed"] = "dd: DISTANCE ATOMS=1,2 \nRESTRAINT ARG=dd KAPPA=2000 AT=" + str(refdist) + "\nPRINT ARG=dd FILE=colvar"
         plrun = runMDCalc("forces2", code, version, runner, rparams )
      # And create our reports from the two runs
      md_failed, val1, val2 = mdrun or plrun, np.ones(1), np.ones(1) 
      if not md_failed : val1, val2 = np.loadtxt("tests/" + code + "/forces1_" + version + "/colvar")[:,1], np.loadtxt("tests/" + code + "/forces2_" + version + "/colvar")[:,1]
      writeReportPage( "forces", code, version, md_failed, ["forces1", "forces2"], val1, val2, 0.001*np.ones(val1.shape) )
      of.write("| PLUMED forces passed correctly | " + getBadge( check( md_failed, val1, val2, 0.001*np.ones(val1.shape) ), "forces", code, version) + " | \n")
   if info["virial"]=="yes" :
      params = runner.setParams()
      params["nsteps"], params["ensemble"] = 50, "npt"
      params["plumed"] = "vv: VOLUME \n PRINT ARG=vv FILE=volume"
      run1 = runMDCalc("virial1", code, version, runner, params )
      params["pressure"] = 1001*params["pressure"] 
      run3 = runMDCalc("virial3", code, version, runner, params )
      params["plumed"] = "vv: VOLUME \n RESTRAINT AT=0.0 ARG=vv SLOPE=-60.221429 \nPRINT ARG=vv FILE=volume"
      run2 = runMDCalc("virial2", code, version, runner, params )
      md_failed, val1, val2, val3 = run1 or run2 or run3, np.ones(1), np.ones(1), np.ones(1)
      if not md_failed : val1, val2, val3 = np.loadtxt("tests/" + code + "/virial1_" + version + "/volume")[:,1], np.loadtxt("tests/" + code + "/virial2_" + version + "/volume")[:,1], np.loadtxt("tests/" + code + "/virial3_" + version + "/volume")[:,1]
      writeReportPage( "virial", code, version, md_failed, ["virial1", "virial2"], val1, val2, np.abs(val3-val1) )
      of.write("| PLUMED virial passed correctly | " + getBadge( check( md_failed, val1, val2, np.abs(val3-val1) ), "virial", code, version) + " | \n")
   if info["energy"]=="yes" :
      params["nsteps"], params["plumed"] = 150, "e: ENERGY \nPRINT ARG=e FILE=energy"
      md_failed, md_energy, pl_energy = runMDCalc( "energy", code, version, runner, params ), np.ones(1), np.ones(1) 
      if not md_failed and os.path.exists("tests/" + code + "/energy_" + version + "/energy") : md_energy, pl_energy = runner.getEnergy("tests/" + code + "/energy_" + version), np.loadtxt("tests/" + code + "/energy_" + version + "/energy")[:,1]
      else : md_failed = True
      writeReportPage( "energy", code, version, md_failed, ["energy"], md_energy, pl_energy, 0.001*np.ones(len(md_energy)) )
      of.write("| MD code potential energy passed correctly | " + getBadge( check( md_failed, md_energy, pl_energy, 0.001*np.ones(len(md_energy)) ), "energy", code, version) + " | \n") 
      sqrtalpha = 1.1
      alpha = sqrtalpha*sqrtalpha
      if info["engforces"]=="yes" :
         params = runner.setParams()
         params["nsteps"], params["ensemble"] = 50, "nvt"
         params["plumed"] = "e: ENERGY\n PRINT ARG=e FILE=energy"
         run1 = runMDCalc("engforce1", code, version, runner, params )
         params["temperature"] = params["temperature"]*alpha
         params["relaxtime"] = params["relaxtime"] / sqrtalpha
         params["tstep"] = params["tstep"] / sqrtalpha
         run3 = runMDCalc("engforce3", code, version, runner, params )
         params["plumed"] = "e: ENERGY\n PRINT ARG=e FILE=energy \n RESTRAINT AT=0.0 ARG=e SLOPE=" + str(alpha - 1)
         run2 = runMDCalc("engforce2", code, version, runner, params )
         md_failed, val1, val2, val3 = run1 or run2 or run3, np.ones(1), np.ones(1), np.ones(1)
         if not md_failed : val1, val2, val3 = np.loadtxt("tests/" + code + "/engforce1_" + version + "/energy")[:,1], np.loadtxt("tests/" + code + "/engforce2_" + version + "/energy")[:,1], np.loadtxt("tests/" + code + "/engforce3_" + version + "/energy")[:,1]
         writeReportPage( "engforce", code, version, md_failed, ["engforce1", "engforce2"], val1, val2, np.abs(val1-val3) ) 
         of.write("| PLUMED forces on potential energy passed correctly | " + getBadge( check( md_failed, val1, val2, np.abs(val1-val3) ), "engforce", code, version) + " | \n") 
      if info["engforces"] and info["virial"]=="yes" :
         params = runner.setParams()
         params["nsteps"], params["ensemble"] = 150, "npt"
         params["plumed"] = "e: ENERGY\n v: VOLUME \n PRINT ARG=e,v FILE=energy"
         run1 = runMDCalc("engvir1", code, version, runner, params )
         params["temperature"] = params["temperature"]*alpha
         params["relaxtime"], params["prelaxtime"] = params["relaxtime"] / sqrtalpha, params["prelaxtime"] / sqrtalpha
         params["tstep"] = params["tstep"] / sqrtalpha
         run3 = runMDCalc("engvir3", code, version, runner, params )
         params["plumed"] = "e: ENERGY\n v: VOLUME \n PRINT ARG=e,v FILE=energy \n RESTRAINT AT=0.0 ARG=e SLOPE=" + str(alpha - 1)
         run2 = runMDCalc("engvir2", code, version, runner, params )
         md_failed, val1, val2, val3 = run1 or run2 or run3, np.ones(1), np.ones(1), np.ones(1)
         if not md_failed : val1, val2, val3 = np.loadtxt("tests/" + code + "/engvir1_" + version + "/energy")[:,1:], np.loadtxt("tests/" + code + "/engvir2_" + version + "/energy")[:,1:], np.loadtxt("tests/" + code + "/engvir3_" + version + "/energy")[:,1:]
         writeReportPage( "engvir", code, version, md_failed, ["engvir1", "engvir2"], val1, val2, np.abs(val1-val3) )
         of.write("| PLUMED contribution to virial due to force on potential energy passed correctly | " + getBadge( check( md_failed, val1, val2, np.abs(val1-val3) ), "engvir", code, version) + " | \n") 
   of.close()
   # Read output file to get status
   ifn, of = open("tests/" + code + "/" + fname, "r"), open("tests/" + code + "/info.yml", "a")
   inp = ifn.read()
   ifn.close()
   if "failed-red.svg" in inp : of.write("test_plumed" + version + ": broken \n")
   elif "%25-green.svg" in inp and ("%25-red.svg" in inp or "%25-yellow.svg" in inp) : of.write("test_plumed" + version + ": partial\n")
   elif "%25-yellow.svg" in inp : of.write("test_plumed" + version + ": partial\n")
   elif "%25-green.svg" in inp : of.write("test_plumed" + version + ": working \n")
   elif "%25-red.svg" in inp : of.write("test_plumed" + version + ": broken \n")
   else : raise Exception("Found no test badges in output for tests on " + code + " with " + version)
   of.close()

def getBadge( sucess, filen, code, version ) :
   badge = '[![tested on ' + version + '](https://img.shields.io/badge/' + version + '-'
   if sucess<0 : badge = badge + 'failed-red.svg'
   elif sucess<5 : badge = badge + 'fail ' + str(sucess) + '%25-green.svg'
   elif sucess<20 : badge = badge + 'fail ' + str(sucess) + '%25-yellow.svg'
   else : badge = badge + 'fail ' + str(sucess) + '%25-yellow.svg'
   return badge + ')](' + filen + '_' + version + '.html)'

def writeReportPage( filen, code, version, md_fail, zipfiles, ref, data, denom ) :
   # Read in the file
   f = open( "pages/" + filen + ".md", "r" )
   inp = f.read()
   f.close()
   # Now output the file 
   of = open("tests/" + code + "/" + filen + "_" + version + ".md", "w+" )
   for line in inp.splitlines() :
       if "Trajectory" in line and "#" in line : 
           if len(zipfiles)!=1 : ValueError("wrong number of trajectories")
           of.write(line + "\n")
           of.write("Input and output files for the test calculation are available in this [zip archive](" + zipfiles[0] + "_" + version + ".zip) \n\n")
       elif "Trajectories" in line and "#" in line :
           if len(zipfiles)!=2 : ValueError("wrong number of trajectories")
           of.write(line + "\n")
           of.write("1. Input and output files for the first calculation described above are available in this [zip archive](" + zipfiles[0] + "_" + version + ".zip) \n")
           of.write("2. Input and output files for the second calculation described above are available in this [zip archive](" + zipfiles[1] + "_" + version + ".zip) \n\n") 
       elif "Results" in line and "#" in line and md_fail :
           of.write(line + "\n")
           of.write("Calculations were not sucessful and no data was generated for comparison\n")  
       else : of.write(line + "\n")
   if not md_fail and hasattr(data, "__len__") : 
      if len(zipfiles)==1 : of.write("\n| MD code output | PLUMED output | % Difference | \n")
      else : of.write("\n| First result | Second result | % Difference | \n")
      of.write("|:-------------|:--------------|:--------------| \n")
      nlines = min( 20, len(ref) )
      percent_diff = np.divide( np.abs( ref - data ), denom, out=np.zeros_like(denom), where=denom!=0 )
      for i in range(nlines) : 
          if hasattr(ref[i], "__len__") :
             ref_strings = [ "%.4f" % x for x in ref[i] ] 
             data_strings = [ "%.4f" % x for x in data[i] ]
             pp_strings = [ "%.4f" % x for x in percent_diff[i] ]
             of.write("|" + " ".join(ref_strings) + " | " + " ".join(data_strings) + " | " + " ".join(pp_strings) + "| \n")
          else : of.write("|" + str(ref[i]) + " | " + str(data[i]) + " | " + str(pecent_diff[i]) + "| \n")
   elif not md_fail : 
      if len(zipfiles)==1 : of.write("\n| MD code output | PLUMED output | % Difference | \n")
      else : of.write("| First result | Second result | % Difference | \n")
      of.write("|:-------------|:--------------|:--------------| \n")
      of.write("| " + str(ref) + " | " + str(data) + " | " + str(percent_diff) + " | \n")
   of.close()

def check( md_failed, val1, val2, val3 ) :
   if md_failed : return -1
   if hasattr(val2, "__len__") and len(val1)!=len(val2) : return -1
   if hasattr(val2, "__len__") and len(val3)!=len(val2) : return -1
   percent_diff = np.divide( np.abs( val1 - val2 ), val3, out=np.zeros_like(val3), where=val3!=0 ) 
   return int(np.round( np.average( percent_diff ) )) 

if __name__ == "__main__" :
   code, version, argv = "", "", sys.argv[1:]
   try:
     opts, args = getopt.getopt(argv,"hc:v:",["code="])
   except:
       print('runtests.py -c <code> -v <version>')

   for opt, arg in opts:
       if opt in ['-h'] :
          print('runtests.py -c <code>')
          sys.exit()
       elif opt in ["-c", "--code"]:
          code = arg
       elif opt in ["-v", "--version"]:
          version = arg

   # Build all the pages that describe the tests for this code
   buildTestPages( "tests/" + code )
   # Copy the engforce file 
   shutil.copy("pages/engforce.md","pages/engvir.md")
   # Build the default test pages
   buildTestPages( "pages" )
   # Create an __init__.py module for the desired code
   ipf = open("tests/" + code + "/__init__.py", "w+")
   ipf.write("from .mdcode import mdcode\n")
   ipf.close()
   # Now import the module
   d = importlib.import_module( "tests." + code, "mdcode" )
   # And create the class that interfaces with the MD code output
   runner = d.mdcode()
   # Now run the tests 
   runTests( code, version, runner )
