import os

def checkWorkflow() : 
   print("Checking Workflow")
   f = open(".github/workflows/main.yml","r")
   inp = f.read()
   f.close()

   for line in inp.splitlines() :
       if "replica:" in line :
          adirs = line.replace("replica:","").replace("[","").replace("]","").replace(",","").split()
          rdirs = os.listdir("tests")
          if rdirs != adirs :
             ValueError("Tests have not been updated. Run the command create_workflow.py in create_workflow directory and commit the changes to .github/actions/main.yml")

def buildBrowsePage( tested ) :
   print("Building browse page")
   f = open("browse.md","w+")
   f.write("Browse the tests\n")
   f.write("-----------------\n")
   f.write("The codes that are currently tested by PLUMED-TESTCENTER are listed below.  PLUMED-TESTCENTER monitors whether the current and development versions of the code can be used to complete the tests for each of these codes.\n")
   f.write("| Name of Program  | Compiles | Passes tests | \n")
   f.write("|:----------------:|:--------:|:------------:| \n")
   for code in os.listdir("tests") :
       of = open( "tests/" + code + "/install.md", "w+" )
       of.write("Compiling " + code + "\n")
       of.write("------------------------\n \n")
       of.write("To compile " + code + " and PLUMED the following bash script was used\n \n")
       of.write("```bash\n")
       sf = open("tests/" + code + "/install.sh" ,"r")
       inp = sf.read()
       for line in inp.splitlines() : of.write( line + "\n")
       of.write("```\n")
       of.close()
       compile_badge, test_badge = "", ""
       for i in range(len(tested)):
           compile_badge = compile_badge + ' [![tested on ' + tested[i] + '](https://img.shields.io/badge/' + tested[i] + '-'
           compile_badge = compile_badge + 'passing-green.svg'
           compile_badge = compile_badge + ')](tests/' + code + '/install.md)'
           test_badge = test_badge + ' [![tested on ' + tested[i] + '](https://img.shields.io/badge/' + tested[i] + '-'
           test_badge = test_badge + 'passing-green.svg'
           test_badge = test_badge + ')](www.youtube.com)'
       f.write("| [" + code + "](www.youtube.com) | " + compile_badge + " | " + test_badge + " | \n")  
   f.write(" \n")
   f.write("**Building PLUMED**\n")
   f.write(" \n")
   f.write("When the tests above are run PLUMED is built using the following bash script.\n")
   f.write(" \n")
   f.write("```bash\n")
   sf = open(".ci/install.plumed","r")
   inp = sf.read()
   for line in inp.splitlines() : f.write( line + "\n")
   f.write("```\n")
   f.close()

if __name__ == "__main__":
   # Check that the workflow matches with the directories
   checkWorkflow()
   # Build the page with all the MD codes
   stable_version = "2.7"
   buildBrowsePage( ("v"+ stable_version,"master") )
