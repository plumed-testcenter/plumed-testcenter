import sys
import getopt
import shutil
import importlib
from runtests import buildTestPages, runTests

if __name__ == "__main__":
    # to run PATH must contain the dir to plumed and the one to the executable of your code
    code = "quantum_espresso"
    prefix = "local_"
    plumedToRun = [{"plumed": "plumed", "printJson": False}]
    print("Code: " + code)
    print("Preparing pages")
    # Engforce and engvir share the same procedure
    shutil.copy("pages/engforce.md", "pages/engvir.md")
    # buildTestPages("tests/" + code, prefix, plumedToRun)
    buildTestPages("pages", prefix, plumedToRun)
    # Create an __init__.py module for the desired code
    with open("tests/" + code + "/__init__.py", "w") as ipf:
        ipf.write("from .mdcode import mdcode\n")
    d = importlib.import_module("tests." + code, "mdcode")
    # And create the class that interfaces with the MD code output
    runner = d.mdcode()
    # Now run the tests
    print("Running the tests on stable")
    # execNameChanged=False because in my case I have compiled qe without changing its suffix
    runTests(
        code,
        "stable",
        runner,
        prefix=prefix,
        settigsFor_runMDCalc=dict(execNameChanged=False, makeArchive=False),
    )
