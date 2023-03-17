import sys
import getopt

def buildInstallPage( code ) :
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

if __name__ == "__main__":
   code, argv = "", sys.argv[1:]
   try:
     opts, args = getopt.getopt(argv,"hc:",["code="])
   except:
       print('setupInstall.py -c <code>')
   # Setup compile page
   buildInstallPage(code)
