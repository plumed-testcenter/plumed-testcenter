import os
import sys
import getopt
import zipfile

def zip(path):
    """ Zip a path removing the original file """
    with zipfile.ZipFile(path + ".zip", "w") as f_out:
        f_out.write(path)
    os.remove(path)

def buildInstallPage( code ) :
   # Zip all the logs
   zip("tests/" + code + "/stdout.txt.zip")
   zip("tests/" + code + "/stderr.txt.zip")
   zip("tests/" + code + "/stdout_master.txt.zip")
   zip("tests/" + code + "/stderr_master.txt.zip")
   of = open( "tests/" + code + "/install.md", "w+" )
   of.write("Compiling " + code + "\n")
   of.write("------------------------\n \n")
   of.write("To compile " + code + " and PLUMED the following bash script was used.  " + code + " was statically linked with the latest stable version of PLUMED. In a separate build, the master version of PLUMED was linked to " + code + " as a runtime library. \n \n")
   of.write("Build with stable version download: [zipped raw stdout](tests/" + code + "/stdout.txt.zip)  - [zipped raw stderr](tests/" + code + "/stderr.txt.zip) \n")
   of.write("Build with master version download: [zipped raw stdout](tests/" + code + "/stdout_master.txt.zip)  - [zipped raw stderr](tests/" + "/stderr_master.txt.zip) \n\n")
   of.write("```bash\n")
   sf = open("tests/" + code + "/install.sh" ,"r")
   inp = sf.read()
   for line in inp.splitlines() : of.write( line + "\n")
   of.write("```\n\n")
   of.close()

if __name__ == "__main__":
   code, argv = "", sys.argv[1:]
   try:
     opts, args = getopt.getopt(argv,"hc:",["code="])
   except:
       print('setupInstall.py -c <code>')

   for opt, arg in opts:
       if opt in ['-h'] :
          print('setupInstall.py -c <code>')
          sys.exit()
       elif opt in ["-c", "--code"]:
          code = arg
   # Setup compile page
   buildInstallPage(code)
