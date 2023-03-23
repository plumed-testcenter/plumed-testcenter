import getopt
from buildTestPages import buildTestPages

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
