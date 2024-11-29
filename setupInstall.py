import os
import sys
import getopt
import subprocess
import zipfile

def zip(path):
    """ Zip a path removing the original file """
    with zipfile.ZipFile(path + ".zip", "w") as f_out:
        f_out.write(path)
    # I need the file in the ouput
    os.remove(path)

def buildInstallPage( code ) :
   # Zip all the logs
   zip("tests/" + code + "/stdout.txt")
   zip("tests/" + code + "/stderr.txt")
   zip("tests/" + code + "/stdout_master.txt")
   zip("tests/" + code + "/stderr_master.txt")
   of = open( "tests/" + code + "/install.md", "w+" )
   of.write("Compiling " + code + "\n")
   of.write("------------------------\n \n")
   stable_version=subprocess.check_output('plumed info --version', shell=True).decode('utf-8').strip()
   # Save the stable version of plumed to a file to pass to the update.py script (not classy Gareth)
   with open("stable_version.md", "w+") as vf:
    vf.write( stable_version )
   
   of.write("To compile " + code + " and PLUMED the following bash script was used.  " + code + " was statically linked with the v" + stable_version + " of PLUMED. In a separate build, the master version of PLUMED was linked to " + code + " as a runtime library. \n \n")
   of.write("Build with stable version download: [zipped raw stdout](stdout.txt.zip)  - [zipped raw stderr](stderr.txt.zip) \n \n")
   of.write("Build with master version download: [zipped raw stdout](stdout_master.txt.zip)  - [zipped raw stderr](stderr_master.txt.zip) \n\n")
   of.write("```bash\n")

   with open("tests/" + code + "/install.sh", "r") as sf:   
    inp = sf.read()
   for line in inp.splitlines() :
      of.write( line + "\n")
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
