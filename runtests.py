import os
import sys
import yaml
import getopt
import shutil
import subprocess
import numpy as np
import importlib
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

    ofile, inplumed, plumed_inp = open( filename, "w+" ), False, ""
    for line in inp.splitlines() :
        # Detect and copy plumed input files 
        if "```plumed" in line :
           inplumed, plumed_inp = True, ""
        # Test plumed input files that have been found in tutorial 
        elif inplumed and "```" in line :
           inplumed = False
           solutionfile = "this_input_should_work.dat"
           sf = open( solutionfile, "w+" )
           sf.write( plumed_inp )
           sf.close()
           # Test whether the input solution can be parsed
           success = success=test_plumed( "plumed", solutionfile )
           success_master=test_plumed( "plumed_master", solutionfile  )
           # Find the stable version 
           stable_version=subprocess.check_output('plumed info --version', shell=True).decode('utf-8').strip()
           # Use PlumedToHTML to create the input with all the bells and whistles
           html = get_html( plumed_inp, solutionfile, solutionfile, ("v"+ stable_version,"master"), (success,success_master), ("plumed","plumed_master") )
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

def runMDCalc( name, code, runner, params ) :
    # Copy all the input needed for the MD calculation
    shutil.copytree("tests/" + code + "/input", "tests/" + code + "/" + name + "_" + params["version"] ) 
    # Change to the directory to run the calculation
    with cd("tests/" + code + "/" + name + "_" + params["version"] ) :
       # Output the plumed file  
       of = open("plumed.dat","w+")
       of.write(params["plumed"])
       of.close() 
       # Now run the MD calculation
       broken = runner.runMD( params )
    # Make a zip archive that contains the input and output
    shutil.make_archive("tests/" + code + "/" + name + "_" + params["version"], 'zip', "tests/" + code + "/" + name + "_" + params["version"] )
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

   params = runner.setParams()
   params["version"] = version
   params["stable_version"], basic_md_failed = usestable, True
   if info["positions"]=="yes" or info["timestep"]=="yes" or info["mass"]=="yes" or info["charge"]=="yes" : 
      params["plumed"] = "DUMPATOMS ATOMS=@mdatoms FILE=plumed.xyz\n"
      if info["mass"]=="yes" and info["charge"]=="yes" : params["plumed"] = params["plumed"] + "DUMPMASSCHARGE FILE=mq_plumed\n"
      elif info["mass"]=="yes" : params["plumed"] = params["plumed"] + "DUMPMASSCHARGE FILE=mq_plumed ONLY_MASSES\n"
      if info["timestep"]=="yes" : params["plumed"] = params["plumed"] + "t1: TIME\nPRINT ARG=t1 FILE=colvar\n"
      params["nsteps"] = 10
      basic_md_failed = runMDCalc( "basic", code, runner, params )
  
   val1, val2 = 0.1, 0.1  
   of.write("| Description of test | Status | \n")
   of.write("|:--------------------|:------:| \n")
   if info["positions"]=="yes" :
      codepos = 0.1 # runner.getPositions()
      plumedpos = 0.1 # io.read("tests/" + code + "/plumed.xyz")
      writeReportPage( "natoms", code, version, basic_md_failed, ["basic"], "Number of atoms", codepos, plumedpos ) 
      of.write("| MD code number of atoms passed correctly | " + getBadge( check(basic_md_failed, codepos, plumedpos), "natoms", code, version) + "| \n") 
      writeReportPage( "positions", code, version, basic_md_failed, ["basic"], "Positions", codepos, plumedpos )
      of.write("| MD code positions passed correctly | " + getBadge( check(basic_md_failed, codepos, plumedpos), "pos", code, version) + "| \n")
      writeReportPage( "cell", code, version, basic_md_failed, ["basic"], "Cell vectors", codepos, plumedpos )
      of.write("| MD code cell vectors passed correctly | " + getBadge( check(basic_md_failed, codepos, plumedpos), "cell", code, version) + " | \n")
   if info["timestep"]=="yes" :
      md_tstep = runner.getTimestep()
      plumedtimes = np.loadtxt("tests/" + code + "/basic_" + version + "/colvar")[:,1]
      writeReportPage( "timestep", code, version, basic_md_failed, ["basic"], "Timestep", md_tstep, plumedtimes[1]-plumedtimes[0] )
      of.write("| MD timestep passed correctly | " + getBadge( check(basic_md_failed, md_tstep, plumedtimes[1]-plumedtimes[0]), "timestep", code, version) + " | \n")
   if info["mass"]=="yes" : 
      md_masses = 0.1 # runner.getMasses()
      pl_masses = 0.1 # np.loadtxt("tests/" + code + "/mq_plumed")[1]
      writeReportPage( "mass", code, version, basic_md_failed, ["basic"], "Masses", codepos, plumedpos ) 
      of.write("| MD code masses passed correctly | " + getBadge( check( basic_md_failed, md_masses, pl_masses ), "mass", code, version) + " | \n")
   if info["charge"]=="yes" :
      md_charges = 0.1 # runner.getCharges()
      pl_charges = 0.1 # np.loadtxt("tests/" + code + "/mq_plumed")[2]  
      writeReportPage( "charge", code, version, basic_md_failed, ["basic"], "Charges", codepos, plumedpos ) 
      of.write("| MD code charges passed correctly | " + getBadge( check( basic_md_failed, md_charges, pl_charges ), "charge", code, version) + " | \n")
   if info["forces"]=="yes" :
      md_failed = True
      of.write("| PLUMED forces passed correctly | " + getBadge( check( md_failed, val1, val2 ), "forces", code, version) + " | \n")
   if info["virial"]=="yes" :
      md_failed = True
      of.write("| PLUMED virial passed correctly | " + getBadge( check( md_failed, val1, val2 ), "virial", code, version) + " | \n")
   if info["energy"]=="yes" :
      md_failed = True 
      of.write("| MD code potential energy passed correctly | " + getBadge( check( md_failed, val1, val2 ), "energy", code, version) + " | \n") 
      if info["forces"]=="yes" :
         md_failed = True
         of.write("| PLUMED forces on potential energy passed correctly | " + getBadge( check( md_failed, val1, val2 ), "engforce", code, version) + " | \n") 
      if info["virial"]=="yes" :
         md_failed = True 
         of.write("| PLUMED contribution to virial due to force on potential energy passed correctly | " + getBadge( check( md_failed, val1, val2 ), "engvir", code, version) + " | \n") 
   of.close()
   # Read output file to get status
   ifn, of = open("tests/" + code + "/" + fname, "r"), open("tests/" + code + "/info.yml", "a")
   inp = ifn.read()
   ifn.close()
   if "passing-green.svg" in inp and "failed-red.svg" in inp : of.write("test_plumed" + version + ": partial\n")
   elif "passing-green.svg" in inp : of.write("test_plumed" + version + ": working \n")
   elif "failed-red.svg" in inp : of.write("test_plumed" + version + ": broken \n")
   of.close()

def getBadge( sucess, filen, code, version ) :
   badge = '[![tested on ' + version + '](https://img.shields.io/badge/' + version + '-'
   if sucess : badge = badge + 'passing-green.svg'
   else : badge = badge + 'failed-red.svg'
   return badge + ')](' + filen + '_' + version + '.html)'

def writeReportPage( filen, code, version, md_fail, zipfiles, description, ref, data ) :
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
           of.write("Input and output files for the first calculation described above are available in this [zip archive](" + zipfiles[0] + "_" + version + ".zip) \n")
           of.write("Input and output files for the second calculation described above are available in this [zip archive](" + zipfiles[1] + "_" + version + ".zip) \n\n") 
       elif "Results" in line and "#" in line and md_fail :
           of.write(line + "\n")
           of.write("Calculations were not sucessful and no data was generated for comparison\n")  
       else : of.write(line + "\n")
   if not md_fail and hasattr(data, "__len__") : 
      if len(zipfiles)==1 : of.write("\n| PLUMED output | MD code output | \n")
      else : of.write("| First result | Second result | \n")
      of.write("|:-------------|:--------------| \n")
      for d in range(len(ref)) : of.write("|" + str(ref[i]) + " | " + str(data[i]) + " | \n")
   elif not md_fail : 
      if len(zipfiles)==1 : of.write("\n| PLUMED output | MD code output | \n")
      else : of.write("| First result | Second result | \n")
      of.write("|:-------------|:--------------| \n")
      of.write("| " + str(ref) + " | " + str(data) + " | \n")
   of.close()

def check( md_failed, val1, val2 ) :
   if md_failed : return False
   return np.array_equal( val1, val2 )

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
   print("Running tests for ", runner.getName() )
   # Now run the tests 
   runTests( code, version, runner )
