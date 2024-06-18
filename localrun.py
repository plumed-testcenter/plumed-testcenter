import sys
import getopt
import shutil
import importlib
from runtests import buildTestPages, runTests

if __name__ == "__main__":
    code = "quantum_espresso"
    local = "local_"
    plumedToRun = [{"plumed": "plumed", "printJson": False}]
    print("Code: " + code)
    buildTestPages("tests/" + code, local, plumedToRun)
    buildTestPages("pages", local, plumedToRun)
    # Create an __init__.py module for the desired code
    with open("tests/" + code + "/__init__.py", "w") as ipf:
        ipf.write("from .mdcode import mdcode\n")
    d = importlib.import_module("tests." + code, "mdcode")
    # And create the class that interfaces with the MD code output
    runner = d.mdcode()
    # Now run the tests
    runTests(code, version, runner)

if __name__ == "__main_":
    code, version, argv = "", "", sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "hc:v:", ["code="])
    except:
        print("runtests.py -c <code> -v <version>")

    for opt, arg in opts:
        if opt in ["-h"]:
            print("runtests.py -c <code>")
            sys.exit()
        elif opt in ["-c", "--code"]:
            code = arg
        elif opt in ["-v", "--version"]:
            version = arg

    # Build all the pages that describe the tests for this code
    buildTestPages("tests/" + code)
    # Copy the engforce file
    shutil.copy("pages/engforce.md", "pages/engvir.md")
    # Build the default test pages
    buildTestPages("pages")
    # Create an __init__.py module for the desired code
    ipf = open("tests/" + code + "/__init__.py", "w+")
    ipf.write("from .mdcode import mdcode\n")
    ipf.close()
    # Now import the module
    d = importlib.import_module("tests." + code, "mdcode")
    # And create the class that interfaces with the MD code output
    runner = d.mdcode()
    # Now run the tests
    runTests(code, version, runner)
