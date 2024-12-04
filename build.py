# formatted with ruff 0.6.4
import yaml
import os
from datetime import date


def checkWorkflow():
    print("Checking Workflow")
    f = open(".github/workflows/main.yml", "r")
    inp = f.read()
    f.close()

    for line in inp.splitlines():
        if "replica:" in line:
            adirs = (
                line.replace("replica:", "")
                .replace("[", "")
                .replace("]", "")
                .replace(",", "")
                .split()
            )
            rdirs = os.listdir("tests")
            if rdirs != adirs:
                ValueError(
                    "Tests have not been updated. Run the command create_workflow.py in create_workflow directory and commit the changes to .github/actions/main.yml"
                )


def buildBrowsePage(stable_version, tested):
    print("Building browse page")

    browse = f"""## Browse the tests
   
   The codes listed below below were tested on __{date.today().strftime("%B %d, %Y")}__. """
    browse += """PLUMED-TESTCENTER tested whether the current and development versions of the code can be used to complete the tests for each of these codes.
   | Name of Program  | Short description | Compiles | Passes tests |
   |:-----------------|:------------------|:--------:|:------------:|"""

    for code in os.listdir("tests"):
        if os.path.isfile("tests/" + code):
            continue
        compile_badge, test_badge = "", ""

        with open("tmp/extract/tests/" + code + "/info.yml", "r") as stram:
            info = yaml.load(stram, Loader=yaml.BaseLoader)

        for i in range(len(tested)):
            compile_badge += (
                f" [![tested on {tested[i]}](https://img.shields.io/badge/{tested[i]}-"
            )

            compile_status = info["install_plumed_" + tested[i]]
            if compile_status == "working":
                compile_badge += "passing-green.svg"
            elif compile_status == "broken":
                compile_badge += "failed-red.svg"
            else:
                ValueError(
                    f"found invalid compilation status for {code} should be 'working' or 'broken', is '{compile_status}'"
                )
            compile_badge += ")](tests/" + code + "/install.html)"

            test_badge += (
                f" [![tested on {tested[i]}](https://img.shields.io/badge/{tested[i]}-"
            )
            test_status = info["test_plumed" + tested[i]]
            if test_status == "working":
                test_badge += "passing-green.svg"
            elif test_status == "partial":
                test_badge += "partial-yellow.svg"
            elif test_status == "broken":
                test_badge += "broken-red.svg"
            elif test_status == "failing":
                test_badge += "failed-red.svg"
            if tested[i] != stable_version:
                test_badge += f")](tests/{code}/testout_{tested[i]}.html)"

            else:
                test_badge = test_badge + ")](tests/" + code + "/testout.html)"
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
    # Check that the workflow matches with the directories
    checkWorkflow()
    # Build the page with all the MD codes
    vf = open("tmp/extract/stable_version.md", "r")
    stable_version = vf.read()
    vf.close()
    buildBrowsePage("v" + stable_version, ("v" + stable_version, "master"))
