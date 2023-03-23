import sys
import yaml
import getopt
import subprocess
from datetime import date
from buildTestPages import buildTestPages

def runTests(code) :
   # Read in the information on the tests that should be run for this code
   stram=open("tests/" + code + "/info.yml", "r")
   info=yaml.load(stram,Loader=yaml.BaseLoader)["tests"]
   stram.close()

   of = open("tests/" + code + "/testout.md", "w+")
   of.write("Testing " + code + "\n")
   of.write("------------------------\n \n")
   stable_version=subprocess.check_output('plumed info --version', shell=True).decode('utf-8').strip()
   of.write("The tests described in the following table were performed on " + date.today().strftime("%B %d, %Y") + " to test whether the interface between " + code + " and PLUMED is working correctly.\n\n") 
   if info["virial"]=="no" : of.write("WARNING: " + code + " does not pass the virial to PLUMED and it is thus not possible to run NPT simulations with this code\n\n")
   of.write("| Description of test | Status | \n")
   of.write("|:-------------------:|:------:| \n")
   if info["positions"]=="yes" : 
      of.write("| [MD code positions passed correctly](../../pages/positions.html) | BADGE | \n")
      of.write("| [MD code cell vectors passed correctly](../../pages/positions.html) | BADGE | \n")
   if info["timestep"]=="yes" :
      of.write("| [MD timestep passed correctly](../../pages/timestep.html) | BADGE | \n")
   if info["mass_charge"]=="yes" : 
      of.write("| [MD code masses passed correctly](../../pages/mass_charge.html) | BADGE | \n")
      of.write("| [MD code charges passed correctly](../../pages/mass_charge.html) | BADGE | \n")
   if info["forces"]=="yes" :
      of.write("| [PLUMED forces passed correctly](../../pages/forces.html) | BADGE | \n")
   if info["virial"]=="yes" :
      of.write("| [PLUMED virial passed correctly](../../pages/virial.html) | BADGE | \n")
   if info["energy"]=="yes" :
      of.write("| [MD code potential energy passed correctly](../../pages/energy.html) | BADGE | \n") 
      if info["forces"]=="yes" :
         of.write("| [PLUMED forces on potential energy passed correctly](../../pages/engforce.html) | BADGE | \n") 
      if info["virial"]=="yes" : 
         of.write("| [PLUMED contribution to virial due to force on potential energy passed correctly](../../pages/engforce.html) | BADGE | \n") 
   of.close()

if __name__ == "__main__" :
   code, argv = "", sys.argv[1:]
   try:
     opts, args = getopt.getopt(argv,"hc:",["code="])
   except:
       print('runtests.py -c <code>')

   for opt, arg in opts:
       if opt in ['-h'] :
          print('runtests.py -c <code>')
          sys.exit()
       elif opt in ["-c", "--code"]:
          code = arg

   # Build all the pages that describe the tests
   buildTestPages( "tests/" + code )
   # Now run the tests 
   runTests(code)
