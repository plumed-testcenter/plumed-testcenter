import os
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
