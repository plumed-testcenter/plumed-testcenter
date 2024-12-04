# formatted with ruff 0.6.4
import os
import yaml
import shutil
import subprocess
import numpy as np
import importlib
from pathlib import Path
from MDAnalysis.coordinates.XYZ import XYZReader
from datetime import date
from contextlib import contextmanager
from PlumedToHTML import test_plumed, get_html
from runhelper import (
    writeReportForSimulations,
    dictToReport,
    dictToTestoutTableEntry,
)
from runhelper import SUCCESS as SUCCESS_THRESHOLD, PARTIAL as PARTIAL_THRESHOLD
from typing import TextIO, Literal

STANDARD_RUN_SETTINGS = [
    {"plumed": "plumed", "printJson": False},
    {"plumed": "plumed_master", "printJson": True, "version": "master"},
]
# tuple becasue I do't want to mutate it and keep the order
TEST_ORDER = (
    "natoms",
    "positions",
    "cell",
    "timestep",
    "mass",
    "charge",
    "forces",
    "virial",
    "energy",
    "engforces",
    "engvir",
)

# trying to be a little more declarative:
TEST_DESCRIPTIONS = {
    "natoms": "MD code number of atoms passed correctly",
    "positions": "MD code positions passed correctly",
    "cell": "MD code cell vectors passed correctly",
    "timestep": "MD timestep passed correctly",
    "mass": "MD code masses passed correctly",
    "charge": "MD code charges passed correctly",
    "forces": "PLUMED forces passed correctly",
    "virial": "PLUMED virial passed correctly",
    "energy": "MD code potential energy passed correctly",
    "engforces": "PLUMED forces on potential energy passed correctly",
    "engvir": "PLUMED contribution to virial due to force on potential energy passed correctly",
}


@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def yamlToDict(filename, **yamlOpts):
    """Simply opens a yaml file an returns an object with the parsed data"""
    with open(filename, "r") as stram:
        ymldata = yaml.load(stram, **yamlOpts)
        # the return closes the file :)
        return ymldata


def processMarkdown(filename, prefix="", runSettings=STANDARD_RUN_SETTINGS):
    if not os.path.exists(filename):
        raise RuntimeError("Found no file called " + filename)
    with open(filename, "r") as f:
        inp = f.read()

    if prefix != "":
        directory = prefix + os.path.dirname(filename)
        Path(f"./{directory}").mkdir(parents=True, exist_ok=True)
        shutil.copy(filename, prefix + filename)
        filename = prefix + filename
    with open(filename, "w+") as ofile:
        inplumed = False
        plumed_inp = ""
        ninputs = 0
        for line in inp.splitlines():
            # Detect and copy plumed input files
            if "```plumed" in line:
                inplumed = True
                plumed_inp = ""
                ninputs += 1
            # Test plumed input files that have been found in tutorial
            elif inplumed and "```" in line:
                inplumed = False
                solutionfile = f"working{ninputs}.dat"
                with open(solutionfile, "w+") as sf:
                    sf.write(plumed_inp)
                # preparing for get_html
                successes = []
                plumed_exec = []
                versions = []
                usejson = False
                for plmd in runSettings:
                    successes.append(
                        test_plumed(
                            plmd["plumed"], solutionfile, printjson=plmd["printJson"]
                        )
                    )
                    plumed_exec.append(plmd["plumed"])
                    # if we asked to print json and we suceed we can use json in the get_html
                    if not usejson and not successes[-1] and plmd["printJson"]:
                        usejson = True

                    # Get the version
                    if "version" in plmd:
                        versions.append(plmd["version"])
                    else:
                        versions.append(
                            subprocess.check_output(
                                f"{plmd['plumed']} info --version", shell=True
                            )
                            .decode("utf-8")
                            .strip()
                        )

                # Use PlumedToHTML to create the input with all the bells and whistles
                html = get_html(
                    plumed_inp,
                    solutionfile,
                    solutionfile,
                    versions,
                    successes,
                    plumed_exec,
                    usejson=usejson,
                )
                # Print the html for the solution
                ofile.write("{% raw %}\n" + html + "\n {% endraw %} \n")
            elif inplumed:
                if "__FILL__" in line:
                    raise RuntimeError("Should not be incomplete files in this page")
                plumed_inp += line + "\n"
            # Just copy any line that isn't part of a plumed input
            elif not inplumed:
                ofile.write(line + "\n")


def buildTestPages(directory, prefix="", runSettings=STANDARD_RUN_SETTINGS):
    for page in os.listdir(directory):
        if ".md" in page:
            processMarkdown(directory + "/" + page, prefix, runSettings)


def runMDCalc(
    name: str,
    code: str,
    version: str,
    runner,
    params: dict,
    *,
    prefix: str = "",
    execNameChanged: bool = True,
    makeArchive: bool = True,
):
    # Get the name of the executible
    basedir = f"tests/{code}"
    params["executible"] = yamlToDict(f"{basedir}/info.yml", Loader=yaml.BaseLoader)[
        "executible"
    ]
    if execNameChanged:
        params["executible"] += f"_{version}"

    print(f'Starting run "{name}"')
    # Now test that the executable exists if it doesn't then the test is broken
    if shutil.which(params["executible"]) is None:
        print(f"Executable {params['executible']} does not exist in current PATH.")
        return True
    # Copy all the input needed for the MD calculation
    wdir = f"{basedir}/{name}_{version}"
    if prefix != "":
        wdir = f"{prefix}{wdir}"
    shutil.copytree(f"{basedir}/input", f"{wdir}")
    # Change to the directory to run the calculation
    # print(f"{params=}")
    with cd(f"{wdir}"):
        # Output the plumed file
        with open("plumed.dat", "w+") as of:
            of.write(params["plumed"])
        # Now run the MD calculation
        mdExitCode = runner.runMD(params)
    # Make a zip archive that contains the input and output
    if makeArchive:
        shutil.make_archive(f"{wdir}", "zip", f"{wdir}")
    return mdExitCode


def runBasicTests(
    testout: TextIO, outdir: str, info: dict, runMDCalcSettings: dict, tolerance: float
) -> dict:
    """run the (eventual) MD test for position, timestep, mass, and charge"""
    params = runMDCalcSettings["runner"].setParams()
    results = {"mdruns": {}}
    basic_md_failed = True
    if (
        info["positions"] == "yes"
        or info["timestep"] == "yes"
        or info["mass"] == "yes"
        or info["charge"] == "yes"
    ):
        dumpMassesStr = ""
        if info["mass"] == "yes":
            dumpMassesStr = f"DUMPMASSCHARGE FILE=mq_plumed {'' if info['charge'] == 'yes' else 'ONLY_MASSES'}"
        timeStepStr = ""
        if info["timestep"] == "yes":
            timeStepStr = "t1: TIME\nPRINT ARG=t1 FILE=colvar"
        params["plumed"] = f"""DUMPATOMS ATOMS=@mdatoms FILE=plumed.xyz
c: CELL
PRINT ARG=c.* FILE=cell_data
{dumpMassesStr}
{timeStepStr}
"""
        params["nsteps"] = 10
        params["ensemble"] = "npt"
        basic_md_failed = runMDCalc("basic", params=params, **runMDCalcSettings)

    testout.write("| Description of test | Status | \n")
    testout.write("|:--------------------|:------:| \n")
    results["mdruns"]["basic"] = basic_md_failed
    basicSR = writeReportForSimulations(
        testout,
        runMDCalcSettings["code"],
        runMDCalcSettings["version"],
        basic_md_failed,
        ["basic"],
        prefix=runMDCalcSettings["prefix"],
    )
    basicDir = f"{outdir}/basic_{runMDCalcSettings['version']}"
    if info["positions"] == "yes":
        plumednatoms = np.empty(0)
        codenatoms = np.empty(0)
        codepos = np.ones(params["nsteps"])
        plumedpos = np.ones(params["nsteps"])
        codecell = np.ones(params["nsteps"])
        plumedcell = np.ones(params["nsteps"])
        if not basic_md_failed and os.path.exists(f"{basicDir}/plumed.xyz"):
            # Get the trajectory that was output by PLUMED
            plumedtraj = XYZReader(f"{basicDir}/plumed.xyz")
            # Get the number of atoms in each frame from plumed trajectory
            codenatoms = np.array(
                runMDCalcSettings["runner"].getNumberOfAtoms(f"{basicDir}")
            )
            plumednatoms = np.array(
                [frame.positions.shape[0] for frame in plumedtraj.trajectory]
            )
            # Concatenate all the trajectory frames
            codepos = np.array(runMDCalcSettings["runner"].getPositions(f"{basicDir}"))
            first = True

            for frame in plumedtraj.trajectory:
                if first:
                    first = False
                    plumedpos = frame.positions.copy()
                else:
                    plumedpos = np.concatenate((plumedpos, frame.positions), axis=0)
            codecell = np.array(runMDCalcSettings["runner"].getCell(f"{basicDir}"))
            plumedcell = np.loadtxt(f"{basicDir}/cell_data")[:, 1:]

        else:
            basicSR.md_failed = True
            basic_md_failed = True
        # Output results from tests on natoms
        results.update(
            basicSR.writeReportAndTable(
                "natoms",
                "MD code number of atoms passed correctly",
                codenatoms,
                plumednatoms,
                0.01 * np.ones(codenatoms.shape[0]),
            )
        )
        # Output results from tests on positions
        results.update(
            basicSR.writeReportAndTable(
                "positions",
                "MD code positions passed correctly",
                codepos,
                plumedpos,
                tolerance * np.ones(plumedpos.shape),
            )
        )
        # Output results from tests on cell
        results.update(
            basicSR.writeReportAndTable(
                "cell",
                "MD code cell vectors passed correctly",
                codecell,
                plumedcell,
                tolerance * np.ones(plumedcell.shape),
            )
        )

    if info["timestep"] == "yes":
        md_tstep = 0.1
        plumed_tstep = 0.1
        if not basic_md_failed:
            plumedtimes = np.loadtxt(f"{basicDir}/colvar")[:, 1]
            md_tstep = runMDCalcSettings["runner"].getTimestep()
            plumed_tstep = plumedtimes[1] - plumedtimes[0]

            for i in range(1, len(plumedtimes)):
                if plumedtimes[i] - plumedtimes[i - 1] != plumed_tstep:
                    ValueError("Timestep should be the same for all MD steps")

        # Output results from tests on timestep
        results.update(
            basicSR.writeReportAndTable(
                "timestep",
                "MD timestep passed correctly",
                md_tstep,
                plumed_tstep,
                0.0001,
            )
        )
    if info["mass"] == "yes":
        md_masses = np.ones(10)
        pl_masses = np.ones(10)
        if not basic_md_failed:
            md_masses = np.array(runMDCalcSettings["runner"].getMasses(f"{basicDir}"))
            pl_masses = np.loadtxt(f"{basicDir}/mq_plumed")[:, 1]

        # Output results from tests on mass
        results.update(
            basicSR.writeReportAndTable(
                "mass",
                "MD code masses passed correctly",
                md_masses,
                pl_masses,
                0.01 * np.ones(pl_masses.shape),
            )
        )

    if info["charge"] == "yes":
        md_charges = np.ones(10)
        pl_charges = np.ones(10)
        if not basic_md_failed:
            md_charges = np.array(runMDCalcSettings["runner"].getCharges(f"{basicDir}"))
            pl_charges = np.loadtxt(f"{basicDir}/mq_plumed")[:, 2]

        # Output results from tests on charge
        results.update(
            basicSR.writeReportAndTable(
                "charge",
                "MD code charges passed correctly",
                md_charges,
                pl_charges,
                tolerance * np.ones(pl_charges.shape),
            )
        )
    return results


def runForcesTest(
    testout: TextIO, outdir: str, runMDCalcSettings: dict, tolerance: float
) -> dict:
    # First run a calculation to find the reference distance between atom 1 and 2
    version = runMDCalcSettings["version"]
    rparams = runMDCalcSettings["runner"].setParams()
    rparams["nsteps"] = 2
    rparams["ensemble"] = "nvt"
    rparams["plumed"] = "dd: DISTANCE ATOMS=1,2 \nPRINT ARG=dd FILE=colvar"
    refrun_fail = runMDCalc("refres", params=rparams, **runMDCalcSettings)
    mdrun_fail = True
    plrun_fail = True
    results = {"mdruns": {}}
    results["mdruns"]["refres"] = refrun_fail
    if not refrun_fail:
        # Get the reference distance between the atoms
        refdist = np.loadtxt(f"{outdir}/refres_{version}/colvar")[0, 1]
        # Run the calculation with the restraint applied by the MD code
        rparams["nsteps"] = 20
        rparams["ensemble"] = "nvt"
        rparams["restraint"] = refdist
        rparams["plumed"] = "dd: DISTANCE ATOMS=1,2 \nPRINT ARG=dd FILE=colvar"
        mdrun_fail = runMDCalc("forces1", params=rparams, **runMDCalcSettings)
        results["mdruns"]["forces1"] = mdrun_fail
        # Run the calculation with the restraint applied by PLUMED
        rparams["restraint"] = -10
        rparams["plumed"] = (
            "dd: DISTANCE ATOMS=1,2\n"
            f"RESTRAINT ARG=dd KAPPA=2000 AT={refdist}\n"
            "PRINT ARG=dd FILE=colvar\n"
        )
        plrun_fail = runMDCalc("forces2", params=rparams, **runMDCalcSettings)
        results["mdruns"]["forces2"] = plrun_fail
    # And create our reports from the two runs
    md_failed = mdrun_fail or plrun_fail
    val1 = np.ones(1)
    val2 = np.ones(1)
    if not md_failed:
        val1 = np.loadtxt(f"{outdir}/forces1_{version}/colvar")[:, 1]
        val2 = np.loadtxt(f"{outdir}/forces2_{version}/colvar")[:, 1]
    results.update(
        writeReportForSimulations(
            testout,
            runMDCalcSettings["code"],
            version,
            md_failed,
            ["forces1", "forces2"],
            prefix=runMDCalcSettings["prefix"],
        ).writeReportAndTable(
            "forces",
            "PLUMED forces passed correctly",
            val1,
            val2,
            tolerance * np.ones(val1.shape),
        )
    )
    return results


def runVirialTest(
    testout: TextIO, outdir: str, runMDCalcSettings: dict, tolerance: float
) -> dict:
    version = runMDCalcSettings["version"]
    params = runMDCalcSettings["runner"].setParams()
    params["nsteps"] = 50
    params["ensemble"] = "npt"
    params["plumed"] = "vv: VOLUME \n PRINT ARG=vv FILE=volume"
    run1_fail = runMDCalc("virial1", params=params, **runMDCalcSettings)
    params["pressure"] = 1001 * params["pressure"]
    run3_fail = runMDCalc("virial3", params=params, **runMDCalcSettings)
    params["plumed"] = (
        "vv: VOLUME\n"
        "RESTRAINT AT=0.0 ARG=vv SLOPE=-60.221429\n"
        "PRINT ARG=vv FILE=volume\n"
    )
    run2_fail = runMDCalc("virial2", params=params, **runMDCalcSettings)
    results = {"mdruns": {}}
    results["mdruns"]["virial1"] = run1_fail
    results["mdruns"]["virial2"] = run2_fail
    results["mdruns"]["virial3"] = run3_fail
    md_failed = run1_fail or run2_fail or run3_fail
    val1 = np.ones(1)
    val2 = np.ones(1)
    val3 = np.ones(1)

    if not md_failed:
        val1 = np.loadtxt(f"{outdir}/virial1_{version}/volume")[:, 1]
        val2 = np.loadtxt(f"{outdir}/virial2_{version}/volume")[:, 1]
        val3 = np.loadtxt(f"{outdir}/virial3_{version}/volume")[:, 1]
    results.update(
        writeReportForSimulations(
            testout,
            runMDCalcSettings["code"],
            version,
            md_failed,
            ["virial1", "virial2", "virial3"],
            prefix=runMDCalcSettings["prefix"],
        ).writeReportAndTable(
            "virial",
            "PLUMED virial passed correctly",
            val1,
            val2,
            np.abs(val3 - val1),
            denominatorTolerance=tolerance,
        )
    )
    return results


def energyTest(
    testout: TextIO,
    outdir: str,
    title: str,
    nsteps: int,
    ensemble: str,
    sqrtalpha: float,
    docstring: str,
    runMDCalcSettings: dict,
    tolerance: float = 0.0,
    prerelaxtime: bool = False,
) -> dict:
    alpha = sqrtalpha * sqrtalpha
    prefix = runMDCalcSettings["prefix"]
    version = runMDCalcSettings["version"]
    params = runMDCalcSettings["runner"].setParams()
    params["nsteps"] = nsteps
    params["ensemble"] = ensemble
    params["plumed"] = "e: ENERGY\n" "v: VOLUME\n" "PRINT ARG=e,v FILE=energy\n"
    run1_fail = runMDCalc(f"{title}1", params=params, **runMDCalcSettings)
    params["temperature"] = params["temperature"] * alpha
    params["relaxtime"] = params["relaxtime"] / sqrtalpha
    if prerelaxtime:
        params["prelaxtime"] = params["prelaxtime"] / sqrtalpha
    params["tstep"] = params["tstep"] / sqrtalpha
    run3_fail = runMDCalc(f"{title}3", params=params, **runMDCalcSettings)
    params["plumed"] = (
        "e: ENERGY\n"
        "v: VOLUME\n"
        "PRINT ARG=e,v FILE=energy\n"
        f"RESTRAINT AT=0.0 ARG=e SLOPE={alpha - 1}\n"
    )
    run2_fail = runMDCalc(f"{title}2", params=params, **runMDCalcSettings)
    results = {"mdruns": {}}
    results["mdruns"][f"{title}1"] = run1_fail
    results["mdruns"][f"{title}2"] = run2_fail
    results["mdruns"][f"{title}3"] = run3_fail
    md_failed = run1_fail or run2_fail or run3_fail
    val1 = np.ones(1)
    val2 = np.ones(1)
    val3 = np.ones(1)

    if not md_failed:
        val1 = np.loadtxt(f"{outdir}/{title}1_{version}/energy")[:, 1:]
        val2 = np.loadtxt(f"{outdir}/{title}2_{version}/energy")[:, 1:]
        val3 = np.loadtxt(f"{outdir}/{title}3_{version}/energy")[:, 1:]
    results.update(
        writeReportForSimulations(
            testout,
            runMDCalcSettings["code"],
            version,
            md_failed,
            [f"{title}1", f"{title}2", f"{title}3"],
            prefix=prefix,
        ).writeReportAndTable(
            title,
            docstring,
            val1,
            val2,
            np.abs(val1 - val3),
            denominatorTolerance=tolerance,
        )
    )
    return results


def runEnergyTests(
    testout: TextIO, outdir: str, info: dict, runMDCalcSettings: dict, tolerance: float
) -> dict:
    code = runMDCalcSettings["code"]
    prefix = runMDCalcSettings["prefix"]
    version = runMDCalcSettings["version"]
    params = runMDCalcSettings["runner"].setParams()
    params["nsteps"] = 150
    params["ensemble"] = "npt"
    params["plumed"] = "e: ENERGY \nPRINT ARG=e FILE=energy"
    md_failed = runMDCalc("energy", params=params, **runMDCalcSettings)
    results = {"mdruns": {}}
    results["mdruns"]["energy"] = md_failed
    md_energy = np.ones(1)
    pl_energy = np.ones(1)

    if not md_failed and os.path.exists(f"{outdir}/energy_{version}/ene rgy"):
        md_energy = runMDCalcSettings["runner"].getEnergy(f"{outdir}/energy_{version}")
        pl_energy = np.loadtxt(f"{outdir}/energy_{version}/energy")[:, 1]

    else:
        md_failed = True
    results.update(
        writeReportForSimulations(
            testout,
            code,
            version,
            md_failed,
            ["energy"],
            prefix=prefix,
        ).writeReportAndTable(
            "energy",
            "MD code potential energy passed correctly",
            md_energy,
            pl_energy,
            tolerance * np.ones(len(md_energy)),
        )
    )
    # TODO:https://docs.python.org/3/library/string.html#template-strings
    # the .md files can be templated with this string built-in feature,
    # so in the engforces/engvir mds we can change sqrtalpha to an arbitray
    # value for each code (and postprocess the mds a second time here)

    sqrtalpha = 1.1
    if info["engforces"] == "yes":
        results.update(
            energyTest(
                testout,
                outdir,
                "engforce",
                50,
                "nvt",
                sqrtalpha,
                "PLUMED forces on potential energy passed correctly",
                runMDCalcSettings,
                tolerance,
            )
        )

    sqrtalpha = 1.1
    if info["engforces"] and info["virial"] == "yes":
        results.update(
            energyTest(
                testout,
                outdir,
                "engvir",
                150,
                "npt",
                sqrtalpha,
                "PLUMED contribution to virial due to force on potential energy passed correctly",
                runMDCalcSettings,
                tolerance,
                prerelaxtime=True,
            )
        )
    return results


def runTests(
    code: str,
    version: Literal["master", "stable"],
    runner,
    *,
    prefix: str = "",
    settingsFor_runMDCalc: dict = {},
):
    # Read in the information on the tests that should be run for this code
    basedir = f"tests/{code}"
    outdir = basedir
    if prefix != "":
        outdir = f"{prefix}{outdir}"
        Path(f"./{outdir}").mkdir(parents=True, exist_ok=True)
    ymldata = yamlToDict(f"{basedir}/info.yml", Loader=yaml.BaseLoader)
    info = ymldata["tests"]
    tolerance = float(ymldata["tolerance"])
    fname = "testout.md"

    if version == "master":
        fname = "testout_" + version + ".md"
    elif version != "stable":
        raise ValueError("version should be master or stable")

    with open(f"{outdir}/{fname}", "w+") as testout:
        testout.write(f"Testing {code}\n")
        testout.write("------------------------\n \n")
        stable_version = (
            subprocess.check_output("plumed info --version", shell=True)
            .decode("utf-8")
            .strip()
        )
        if version == "stable":
            version = "v" + stable_version
        # it looks strange, but strings do not need the + to be concatenated
        testout.write(
            f"The tests described in the following table were performed on "
            f"__{date.today().strftime('%B %d, %Y')}__ to test whether the "
            f"interface between {code} and "
            f"the {version} version of PLUMED is working correctly.\n\n"
        )
        if info["virial"] == "no":
            testout.write(
                f"WARNING: {code} does not pass the virial to PLUMED and it is thus "
                "not possible to run NPT simulations with this code\n\n"
            )
        if "warning" in ymldata.keys():
            for warn in ymldata["warning"]:
                testout.write(f"WARNING: {warn}\n\n")

        # sugar with the settings that are always the same for runMDCalc
        # note that if I modify directly the input `settingsFor_runMDCalc`,
        # I will change the default parameteres on subsequent calls!!!
        runMDCalcSettings = dict(
            code=code,
            version=version,
            runner=runner,
            prefix=prefix,
            **settingsFor_runMDCalc,
        )
        results = runBasicTests(testout, outdir, info, runMDCalcSettings, tolerance)
        mddict = results["mdruns"]
        # the next runs are not based on the basic run
        if info["forces"] == "yes":
            tmp = runForcesTest(testout, outdir, runMDCalcSettings, tolerance)
            mddict.update(tmp["mdruns"])
            results.update(tmp)

        if info["virial"] == "yes":
            tmp = runVirialTest(testout, outdir, runMDCalcSettings, tolerance)
            mddict.update(tmp["mdruns"])
            results.update(tmp)

        if info["energy"] == "yes":
            tmp = runEnergyTests(testout, outdir, info, runMDCalcSettings, tolerance)
            mddict.update(tmp["mdruns"])
            results.update(tmp)
        results["mdruns"] = mddict
        howbad = []
        for test in TEST_ORDER:
            if test in results.keys():
                dictToReport(results[test])
                howbad.append(results[test]["failure_rate"])
                testout.write(dictToTestoutTableEntry(results[test]))

        if any(val == -1 for val in howbad):
            # at least one test is broken
            test_result = "broken"
        else:
            if all(val < SUCCESS_THRESHOLD for val in howbad):
                # everithyng passes
                test_result = "working"
            elif any(val >= PARTIAL_THRESHOLD for val in howbad):
                # at least one test is failing
                test_result = "failing"
                # should I add a failing over total?
            else:
                # elif any(howbad<SUCCESS_THRESHOLD):
                # elif any(SUCCESS_THRESHOLD<val<PARTIAL_THRESHOLD for val in howbad):
                # no test are failing red (previous if) and
                # at least one test is partially failing
                # but there can be green tests
                # if we are here the if is redundant
                test_result = "partial"

    print(f"Test result for {code} with version {version}: {test_result}")
    with open(f"{outdir}/info.yml", "a") as infoOut:
        infoOut.write(f"test_plumed{version}: {test_result} \n")


if __name__ == "__main__":
    import sys
    import getopt

    code, version, argv = "", "", sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "hc:v:", ["code="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
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
    # Engforce and engvir share the same procedure
    shutil.copy("pages/engforce.md", "pages/engvir.md")
    # Build the default test pages
    buildTestPages("pages")
    # Create an __init__.py module for the desired code
    ipf = open(f"tests/{code}/__init__.py", "w+")
    ipf.write("from .mdcode import mdcode\n")
    ipf.close()
    # Now import the module
    myMDcode = importlib.import_module("tests." + code, "mdcode")
    # And create the class that interfaces with the MD code output
    runner = myMDcode.mdcode()
    # Now run the tests
    runTests(code, version, runner)
