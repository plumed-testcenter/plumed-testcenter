# formatted with ruff 0.6.4
import yaml
import os
from datetime import date


def isTest(path) -> bool:
    """
    Check if a directory can be considered a test directory.
    It is a test directory if it contains the three files:
    - info.yml
    - install.sh
    - mdcode.py
    """
    if os.path.isdir(path):
        for file in ["info.yml", "install.sh", "mdcode.py"]:
            if not os.path.isfile(path + "/" + file):
                return False
        return True
    return False


def checkWorkflow():
    print("Checking Workflow")
    # this works with the replica list written like this
    ###
    # replica:
    #  - "simplemd"
    #  - "lammps"
    #  - "quantum_espresso"
    #  - "i-pi"
    #  - "gromacs"
    #  - "gromacs-vanilla"
    ### or like this:
    # replica: ["simplemd", "lammps", "quantum_espresso", "i-pi", "gromacs", "gromacs-vanilla"]
    with open(".github/workflows/main.yml", "r") as f:
        workflow = yaml.load(f, Loader=yaml.BaseLoader)
        replicalist = workflow["jobs"]["build"]["strategy"]["matrix"]["replica"]

        testdirs = [d for d in os.listdir("tests") if isTest("tests/" + d)]

        replicalist = sorted(replicalist)
        testdirs = sorted(testdirs)
        if testdirs != replicalist:
            error = "Tests have not been updated.\n"
            if len(testdirs) > len(replicalist):
                error += "The following tests are missing from the workflow file:\n"
                for t in testdirs:
                    if t not in replicalist:
                        error += " - " + t + "\n"
            else:
                error += "The following test in the workflow do not have their own directory:\n"
                for t in replicalist:
                    if t not in testdirs:
                        error += " - " + t + "\n"
            raise ValueError(error)


def buildBrowsePage(stable_version, tested):
    print("Building browse page")

    browse = f"""## Browse the tests  
   
The codes listed below below were tested on __{date.today().strftime("%B %d, %Y")}__.
PLUMED-TESTCENTER tested whether the current and development versions of the code can be used to complete the tests for each of these codes.
| Name of Program  | Short description | Compiles | Passes tests |
|:-----------------|:------------------|:--------:|:------------:|
"""
    testdirs = [d for d in os.listdir("tests") if isTest("tests/" + d)]
    testdirs = sorted(testdirs)
    for code in testdirs:
        compile_badge = ""
        test_badge = ""

        with open("tmp/extract/tests/" + code + "/info.yml", "r") as stream:
            info = yaml.load(stream, Loader=yaml.BaseLoader)

        for version in tested:
            # building the compilation badge
            compile_badge += (
                f" [![tested on {version}](https://img.shields.io/badge/{version}-"
            )

            compile_status = info["install_plumed_" + version]
            if compile_status == "working":
                compile_badge += "passing-green.svg"
            elif compile_status == "broken":
                compile_badge += "failed-red.svg"
            else:
                raise ValueError(
                    f"found invalid compilation status for {code}['test_plumed{version}'] should be 'working' or 'broken', is '{compile_status}'"
                )
            compile_badge += ")](tests/" + code + "/install.html)"

            # building the tests badge
            test_badge += (
                f" [![tested on {version}](https://img.shields.io/badge/{version}-"
            )
            test_status = info["test_plumed_" + version]

            if test_status == "working":
                test_badge += "passing-green.svg"
            elif test_status == "partial":
                test_badge += "partial-yellow.svg"
            elif test_status == "broken":
                test_badge += "broken-red.svg"
            elif test_status == "failing":
                test_badge += "failed-red.svg"
            else:
                raise ValueError(
                    f"found invalid test status for {code}['test_plumed{version}'] should be 'working', 'partial', 'failing' or 'broken', is '{test_status}'"
                )
            test_badge += f")](tests/{code}/testout_{version}.html)"

        browse += f"| [{code}]({info['link']}) | {info['description']} | {compile_badge} | {test_badge} | \n"
    browse += " \n"
    browse += """#### Building PLUMED

When the tests above are run PLUMED is built using the install plumed action.
```yaml
- name: Install plumed
      uses: Iximiel/install-plumed@v1
      with:
         CC: "ccache mpicc"
         CXX: "ccache mpic++"
         suffix: "${{ steps.get-key.outputs.suffix }}"
         version: "${{ steps.get-key.outputs.branch }}"
         extra_options: --enable-boost_serialization --enable-fftw --enable-libtorch LDFLAGS=-Wl,-rpath,$LD_LIBRARY_PATH --disable-basic-warnings
```
"""
    # no more necessary
    # browse+=" \n"
    # browse+="```bash\n"
    # sf = open(".ci/install.plumed","r")
    # inp = sf.read()
    # for line in inp.splitlines() : f.write( line + "\n")
    # f.write("```\n")
    with open("browse.md", "w+") as f:
        f.write(browse)


if __name__ == "__main__":
    try:
        # Check that the workflow matches with the directories
        checkWorkflow()

        with open("tmp/extract/stable_version.md", "r") as vf:
            stable_version = vf.read()
        # Build the page with all the MD codes
        buildBrowsePage("v" + stable_version, ("v" + stable_version, "master"))
    except Exception as e:
        print(e)
        exit(1)
