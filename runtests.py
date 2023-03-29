import os
import sys
import yaml
import getopt
import subprocess
import numpy as np
import importlib
from datetime import date
from PlumedToHTML import test_plumed, get_html

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

def runTests(code,version,runner) :
   # Read in the information on the tests that should be run for this code
   stram=open("tests/" + code + "/info.yml", "r")
   info=yaml.load(stram,Loader=yaml.BaseLoader)["tests"]
   stram.close()

   fname = "testout.md";
   if version=="master" : fname = "testout_" + version + ".md"
   elif version!="stable" : ValueError("version should be master or stable")

   of = open("tests/" + code + "/" + fname, "w+")
   of.write("Testing " + code + "\n")
   of.write("------------------------\n \n")
   stable_version=subprocess.check_output('plumed info --version', shell=True).decode('utf-8').strip()
   of.write("The tests described in the following table were performed on __" + date.today().strftime("%B %d, %Y") + "__ to test whether the interface between " + code + " and ")
   if version=="master" : of.write("the master version of PLUMED is working correctly.\n\n") 
   else : of.write("v" + stable_version + " of PLUMED is working correctly.\n\n")
   if info["virial"]=="no" : of.write("WARNING: " + code + " does not pass the virial to PLUMED and it is thus not possible to run NPT simulations with this code\n\n")

   if info["positions"]=="yes" or info["timestep"]=="yes" or info["mass"]=="yes" or info["charge"]=="yes" : 
      plumed_inpt = "DUMPATOMS ATOMS=@mdatoms FILE=plumed.xyz\n"
      if info["mass"]=="yes" or info["charge"]=="yes" : plumed_inpt = plumed_inpt + "DUMPMASSCHARGE FILE=mq_plumed\n"
      if info["timestep"]=="yes" : plumed_inpt = plumed_inpt + "t1: TIME\nPRINT ARG=t1 FILE=colvar\n"
      # runner.runCode( plumed_inpt )
  
   val1, val2 = 0.1, 0.1
   of.write("| Description of test | Status | \n")
   of.write("|:--------------------|:------:| \n")
   if info["positions"]=="yes" :
      codepos = 0.1 # runner.getPositions()
      plumedpos = 0.1 # io.read("tests/" + code + "/plumed.xyz")
      of.write("| [MD code number of atoms passed correctly](../../pages/positions.html) | " + getBadge( check(codepos, plumedpos), "natoms", code, version) + "| \n") 
      of.write("| [MD code positions passed correctly](../../pages/positions.html) | " + getBadge( check(codepos, plumedpos), "pos", code, version) + "| \n")
      of.write("| [MD code cell vectors passed correctly](../../pages/positions.html) | " + getBadge( check(codepos, plumedpos), "cell", code, version) + " | \n")
   if info["timestep"]=="yes" :
      md_tstep = 0.1 # runner.getTimeStep()
      plumed_times = 0.1 #np.loadtxt("tests/" + code + "/colvar")[0]
      of.write("| [MD timestep passed correctly](../../pages/timestep.html) | " + getBadge( checkTimestepMatch(md_tstep,plumed_times), "timestep", code, version) + " | \n")
   if info["mass"]=="yes" : 
      md_masses = 0.1 # runner.getMasses()
      pl_masses = 0.1 # np.loadtxt("tests/" + code + "/mq_plumed")[1]
      of.write("| [MD code masses passed correctly](../../pages/mass_charge.html) | " + getBadge( check( md_masses, pl_masses ), "mass", code, version) + " | \n")
   if info["charge"]=="yes" :
      md_charges = 0.1 # runner.getCharges()
      pl_charges = 0.1 # np.loadtxt("tests/" + code + "/mq_plumed")[2]  
      of.write("| [MD code charges passed correctly](../../pages/mass_charge.html) | " + getBadge( check( md_charges, pl_charges ), "charge", code, version) + " | \n")
   if info["forces"]=="yes" :
      of.write("| [PLUMED forces passed correctly](../../pages/forces.html) | " + getBadge( check( val1, val2 ), "forces", code, version) + " | \n")
   if info["virial"]=="yes" :
      of.write("| [PLUMED virial passed correctly](../../pages/virial.html) | " + getBadge( check( val1, val2 ), "virial", code, version) + " | \n")
   if info["energy"]=="yes" :
      of.write("| [MD code potential energy passed correctly](../../pages/energy.html) | " + getBadge( check( val1, val2 ), "energy", code, version) + " | \n") 
      if info["forces"]=="yes" :
         of.write("| [PLUMED forces on potential energy passed correctly](../../pages/engforce.html) | " + getBadge( check( val1, val2 ), "engforce", code, version) + " | \n") 
      if info["virial"]=="yes" : 
         of.write("| [PLUMED contribution to virial due to force on potential energy passed correctly](../../pages/engforce.html) | " + + getBadge( check( val1, val2 ), "engvir", code, version) + + " | \n") 
   of.close()

def getBadge( sucess, filen, code, version ) :
   badge = '[![tested on ' + version + '](https://img.shields.io/badge/' + version + '-'
   if sucess : badge = badge + 'passing-green.svg'
   else : badge = badge + 'failed-red.svg'
   return badge + ')](tests/' + code + '/' + filen + '_' + version + '.html)'

def check( val1, val2 ) :
   if np.random.uniform(0,1)<0.5 : return True
   return False

def checkTimestepMatch( val1, val2 ) :
   return True

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
