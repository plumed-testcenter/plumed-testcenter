from typing import TextIO, Literal
import numpy as np
# formatted with ruff 0.6.4


# tuple becasue I don't want to mutate it and I want to keep the order
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


SUCCESS = 5
PARTIAL = 20


def successState(success: int) -> Literal["broken", "success", "partial", "failure"]:
    if success < 0:
        return "broken"
    if success < SUCCESS:
        return "success"
    elif success < PARTIAL:
        return "partial"
    return "failure"


def badgeColor(success: int) -> Literal["red", "yellow", "green"]:
    toret = {
        "broken": "red",
        "failure": "red",
        "partial": "yellow",
        "success": "green",
    }
    return toret[successState(success)]


def testOpinion(
    howbad: Literal["broken", "success", "partial", "failure"],
) -> Literal["broken", "failing", "partial", "working"]:
    if howbad is str:
        howbad = [howbad]
    if "broken" in howbad:
        # at least one test is broken
        return "broken"
    elif "failure" in howbad:
        # at least one test is failing
        return "failing"
        # should I add a failing over total?
    elif "partial" in howbad:
        # no test are failing red (previous if) but
        # at least one test is partially failing
        return "partial"
    elif all(val == "success" for val in howbad):
        # the if shoudl be redundant
        # everithyng passes
        return "working"
    raise RuntimeError(
        'testOpinion arguments should be a list of "broken", "failing", "partial", "success"\n'
        f"{howbad=}"
    )


def getBadge(success, filen, version: str):
    badge = f"[![tested on {version}](https://img.shields.io/badge/{version}-"
    if success < 0:
        badge += "failed-red.svg"
    else:
        color = badgeColor(success)
        badge += f"fail%20{success}%25-{color}.svg"
    return badge + f")]({filen}_{version}.html)"


def writeReportPage(
    filen, code, version, md_fail, zipfiles, ref, data, denom, *, prefix=""
):
    # Read in the file
    with open(f"{prefix}pages/{filen}.md", "r") as f:
        inp = f.read()

    # Now output the file
    with open(f"{prefix}tests/{code}/{filen}_{version}.md", "w+") as of:
        for line in inp.splitlines():
            if "Trajectory" in line and "#" in line:
                if len(zipfiles) != 1:
                    raise ValueError("wrong number of trajectories")
                of.write(
                    f"{line}\n"
                    "Input and output files for the test calculation are available in "
                    f"this [zip archive]({zipfiles[0]}_{version}.zip)\n\n"
                )
            elif "Trajectories" in line and "#" in line:
                if len(zipfiles) == 2:
                    of.write(
                        f"{line}\n"
                        "1. Input and output files for the calculation where the restraint "
                        "is applied by the MD code available in "
                        f"this [zip archive]({zipfiles[0]}_{version}.zip)\n"
                        "2. Input and output files for the calculation where the restraint "
                        "is applied by PLUMED are available in t"
                        f"his [zip archive]({zipfiles[1]}_{version}.zip)\n\n"
                    )
                elif len(zipfiles) == 3:
                    of.write(
                        f"{line}\n"
                        "1. Input and output files for the unpeturbed calculation "
                        "are available in "
                        f"this [zip archive]({zipfiles[0]}_{version}.zip)\n"
                        "2. Input and output files for the peturbed calculation "
                        "are available in "
                        f"this [zip archive]({zipfiles[2]}_{version}.zip)\n\n"
                        "3. Input and output files for the peturbed calculation in "
                        "which a PLUMED restraint is used to undo the effect of the "
                        "changed MD parameters are available in "
                        f"this [zip archive]({zipfiles[1]}_{version}.zip)\n\n"
                    )
                else:
                    raise ValueError("wrong number of trajectories")
            elif "Results" in line and "#" in line and md_fail:
                of.write(
                    f"{line}\n"
                    "Calculations were not successful and no data was generated for comparison\n"
                )
            else:
                of.write(line + "\n")
        if not md_fail and hasattr(data, "__len__"):
            if len(zipfiles) == 1:
                of.write(
                    "\n| MD code output | PLUMED output | Tolerance | % Difference | \n"
                )
            else:
                of.write(
                    "\n| Original result | Result with PLUMED | Effect of peturbation | % Difference | \n"
                )
            of.write(
                "|:-------------|:--------------|:--------------|:--------------| \n"
            )
            nlines = min(20, len(ref))
            percent_diff = 100 * np.divide(
                np.abs(ref - data), denom, out=np.zeros_like(denom), where=denom != 0
            )
            for i in range(nlines):
                if hasattr(ref[i], "__len__"):
                    ref_strings = " ".join([f"{x:.4f}" for x in ref[i]])
                    data_strings = " ".join([f"{x:.4f}" for x in data[i]])
                    denom_strings = " ".join([f"{x:.4f}" for x in denom[i]])
                    pp_strings = " ".join([f"{x:.4f}" for x in percent_diff[i]])
                    of.write(
                        f"| {ref_strings} | {data_strings} | {denom_strings} | {pp_strings} | \n"
                    )
                else:
                    # TODO:ask if also this needs formatting (just append ":.4f")
                    of.write(
                        f"| {ref[i]} | {data[i]} | {denom[i]} | {percent_diff[i]} |\n"
                    )
        elif not md_fail:
            if len(zipfiles) == 1:
                of.write(
                    "\n| MD code output | PLUMED output | Tolerance | % Difference | \n"
                )
            else:
                of.write(
                    "| Original result | Result with PLUMED | Effect of peturbation | % Difference | \n"
                )
            of.write(
                "|:-------------|:--------------|:--------------|:--------------| \n"
                f"| {ref} | {data} | {denom} | {100 * np.abs(ref - data) / denom} | \n"
            )


def check(
    md_failed: "int|bool",
    ref: "float|np.ndarray",
    data: "float|np.ndarray",
    denom: "float|np.ndarray",
    denominatorTolerance: float = 0.0,
) -> int:
    # this may be part of writeReportForSimulations
    if md_failed:
        return -1
    if hasattr(data, "__len__") and len(ref) != len(data):
        return -1
    if hasattr(data, "__len__") and len(denom) != len(data):
        return -1
    percent_diff = 100 * np.divide(
        np.abs(ref - data),
        denom,
        out=np.zeros_like(denom),
        where=denom > denominatorTolerance,
    )
    return int(np.round(np.average(percent_diff)))


class writeReportForSimulations:
    """helper class to write the report of the simulations"""

    testout: TextIO
    code: str
    version: Literal["master", "stable"]
    md_failed: "bool | int"
    simulations: list[str]
    prefix: str = ""

    def __init__(
        self,
        code: str,
        version: str,
        md_failed: "bool | int",
        simulations: list[str],
    ) -> None:
        self.code = code
        self.version = version
        self.md_failed = md_failed
        self.simulations = simulations

    def writeReportAndTable(
        self,
        kind: str,
        ref: "float|np.ndarray",
        data: "float|np.ndarray",
        denom: "float|np.ndarray",
        *,
        denominatorTolerance: float = 0.0,
    ) -> dict:
        report = {
            "filen": kind,
            "code": self.code,
            "version": self.version,
            "md_fail": self.md_failed,
            "zipfiles": self.simulations,
            "ref": ref,
            "data": data,
            "denom": denom,
        }
        failure_rate = check(
            self.md_failed,
            ref,
            data,
            denom,
            denominatorTolerance=denominatorTolerance,
        )
        report["failure_rate"] = failure_rate
        report["docstring"] = TEST_DESCRIPTIONS[kind]
        return report


def dictToReport(input: dict, *, prefix: str = ""):
    # isolates the needed data from the dictionary
    report = {
        "filen": input["filen"],
        "code": input["code"],
        "version": input["version"],
        "md_fail": input["md_fail"],
        "zipfiles": input["zipfiles"],
        "ref": input["ref"],
        "data": input["data"],
        "denom": input["denom"],
        "prefix": prefix,
    }
    writeReportPage(**report)


def dictToTestoutTableEntry(input: dict):
    docstring = input["docstring"]
    failure_rate = input["failure_rate"]
    kind = input["filen"]
    version = input["version"]
    return (
        f"| {docstring} | "
        + getBadge(
            failure_rate,
            kind,
            version,
        )
        + " |\n"
    )
