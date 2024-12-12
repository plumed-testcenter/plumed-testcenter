from typing import Literal
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
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
    if success < 0:
        badge = "failed-red.svg"
    else:
        color = badgeColor(success)
        badge = f"fail%20{success}%25-{color}.svg"
    return f"[![tested on {version}](https://img.shields.io/badge/{version}-{badge})]({filen}_{version}.html)"


def tabulate3x3(data, fmt="{x:.4f}"):
    mytable = r"$\begin{array}{ccc} "
    mytable += " & ".join([fmt.format(x=x) for x in data[0:3]])
    mytable += r" \\ "
    mytable += " & ".join([fmt.format(x=x) for x in data[3:6]])
    mytable += r" \\ "
    mytable += " & ".join([fmt.format(x=x) for x in data[6:9]])
    mytable += r" \end{array}$"
    return mytable


def figure_cell(ax, data, ref, tolerance):
    cmap = mpl.cm.viridis.with_extremes(over="k")
    diff = np.abs(data - ref)
    norm = mpl.colors.Normalize(vmin=0.0, vmax=tolerance)
    # tolerance is denom[0,0]
    im = ax.imshow(diff.T, cmap=cmap, norm=norm)
    ax.set_xlabel("Timestep")

    ax.set_yticks(
        range(9),
        labels=[f"[{x},{y}]" for x in [0, 1, 2] for y in [0, 1, 2]],
    )
    ax.set_ylabel("Cell diff (components)")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.figure.colorbar(
        im,  # cmap=mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
        ax=ax,
        # cax=ax.inset_axes([1.05, 0, 0.05, 1]),
        # orientation="horizontal",
        extend="max",
        label="Error up to tolerance",
    )
    ax.set_ylim(bottom=0)
    return ax


def figure_engforces(ax, xmd, xpl, xmd_xmdp):
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    x = np.arange(len(xmd), dtype=int)
    delta = np.abs(xmd - xpl)
    ax.plot(x, delta, label=r"$\vert x_{md} - x_{pl}\vert$")
    ax.plot(x, xmd_xmdp, label=r"$\vert x_{md}'-x_{md}\vert$ (perturbation)")
    ax2 = ax.twinx()
    divlabel = r"$100*\frac{\vert x_{md} - x_{pl} \vert}{\vert x_{md}'-x_{md}\vert}$"
    ax2.plot(
        x,
        100 * np.divide(delta, xmd_xmdp, out=np.zeros_like(delta), where=xmd_xmdp != 0),
        label=divlabel,
        color="tab:red",
    )
    ax2.set_ylabel(divlabel, color="tab:red")

    ax.legend(loc="upper left")
    return ax


def writeReportPage(
    filen, code, version, md_fail, zipfiles, ref, data, denom, *, prefix="", extra={}
):
    with_image = False
    output = {
        # this is a workaround for not modify plumed2html
        # the other brackets have been "doubled" {{}}
        "% raw %": "{% raw %}",
        "% endraw %": "{% endraw %}",
    }
    # Now prepare the file output

    # if the file template contains the wrong keyword the format function will throw
    if len(zipfiles) == 1:
        output["Trajectory"] = (
            f"Input and output files for the test calculation are available in"
            f"this [zip archive]({zipfiles[0]}_{version}.zip)"
        )

    elif len(zipfiles) == 2:
        output["Trajectories"] = (
            " 1. Input and output files for the calculation where the restraint"
            f" is applied by the MD code available in this [zip archive]({zipfiles[0]}_{version}.zip)\n\n"
            " 2. Input and output files for the calculation where the restraint"
            f" is applied by PLUMED are available in this [zip archive]({zipfiles[1]}_{version}.zip"
        )

    elif len(zipfiles) == 3:
        output["Trajectories"] = (
            " 1. Input and output files for the unpeturbed calculation are available in "
            f"this [zip archive]({zipfiles[0]}_{version}.zip)\n\n"
            " 2. Input and output files for the peturbed calculation are available in "
            f"this [zip archive]({zipfiles[2]}_{version}.zip)\n\n"
            " 3. Input and output files for the peturbed calculation in which a PLUMED restraint is used to undo the effect of the "
            "changed MD parameters are available in "
            f"this [zip archive]({zipfiles[1]}_{version}.zip)\n"
        )
    else:
        raise ValueError("wrong number of trajectories")

    if md_fail:
        output["Results"] = (
            "Calculations were not successful and no data was generated for comparison"
        )

    else:
        if len(zipfiles) == 1:
            col1 = "MD code output"
            col2 = "PLUMED output"
            col3 = "Tolerance"
        else:
            col1 = "Original"
            col2 = "With PLUMED"
            col3 = "Effect of peturbation"
        output["Results"] = f"| {col1} | {col2} | {col3} | % Difference | \n"
        output["Results"] += (
            "|:-------------|:--------------|:--------------|:--------------| \n"
        )
        if hasattr(data, "__len__"):
            # then I would like to set up an image
            nlines = min(20, len(ref))
            percent_diff = 100 * np.divide(
                np.abs(ref - data), denom, out=np.zeros_like(denom), where=denom != 0
            )
            if hasattr(ref[0], "__len__"):
                if filen == "cell":
                    with_image = True
                    fig, ax = plt.subplots()
                    figure_cell(ax, data, ref, denom[0, 0])
                    fig.tight_layout()
                    fig.savefig(f"{prefix}tests/{code}/{filen}_{version}.png")
                if filen == "engvir" or filen == "engforces":
                    with_image = True
                    fig, axes = plt.subplots(2, sharex=True)
                    figure_engforces(axes[0], data[:, 0], ref[:, 0], denom[:, 0])
                    axes[0].set_ylabel("Energy")
                    figure_engforces(axes[1], data[:, 1], ref[:, 1], denom[:, 1])
                    axes[1].set_ylabel("Volume")
                    axes[1].set_xlabel("Timesteps")
                    fig.tight_layout()
                    fig.savefig(f"{prefix}tests/{code}/{filen}_{version}.png")

                for i in range(nlines):
                    if ref.shape[1] == 9:
                        ref_strings = tabulate3x3(ref[i])
                        data_strings = tabulate3x3(data[i])
                        denom_strings = tabulate3x3(denom[i])
                        pp_strings = tabulate3x3(percent_diff[i])
                    else:
                        ref_strings = " ".join([f"{x:.4f}" for x in ref[i]])
                        data_strings = " ".join([f"{x:.4f}" for x in data[i]])
                        denom_strings = " ".join([f"{x:.4f}" for x in denom[i]])
                        pp_strings = " ".join([f"{x:.4f}" for x in percent_diff[i]])

                    output["Results"] += (
                        f"| {ref_strings} | {data_strings} | {denom_strings} | {pp_strings} | \n"
                    )

            else:
                with_image = True
                fig, ax = plt.subplots()
                ax.xaxis.set_major_locator(MaxNLocator(integer=True))
                x = np.arange(len(ref), dtype=int)
                diff = np.abs(np.array(ref) - np.array(data))
                ax.plot(x, diff, label=f"$|| (${col1}$) - (${col2}$)||$")
                ax.plot(x, denom, label=col3)
                ax2 = ax.twinx()
                ax2.plot(x, percent_diff, label="% Difference", color="tab:red")
                if ax2.get_ylim()[1] < 1:
                    ax2.set_ylim(0, 1)
                ax2.set_ylabel("% Difference", color="tab:red")

                fig.legend(loc="outside upper center", ncol=3)
                fig.savefig(f"{prefix}tests/{code}/{filen}_{version}.png")
                for i in range(nlines):
                    # TODO:ask if also this needs formatting (just append ":.4f")

                    output["Results"] += (
                        f"| {ref[i]} | {data[i]} | {denom[i]} | {percent_diff[i]} |\n"
                    )

        else:
            output["Results"] += (
                f"| {ref} | {data} | {denom} | {100 * np.abs(ref - data) / denom} | \n"
            )
    if "sqrtalpha" in extra:
        output["sqrtalpha"] = extra["sqrtalpha"]
    # Read in the file template
    with open(f"{prefix}pages/{filen}.md", "r") as f:
        inp = f.read()
    # And output it:
    with open(f"{prefix}tests/{code}/{filen}_{version}.md", "w+") as of:
        of.write(inp.format(**output))
        if with_image:
            of.write("\n### Graphical representation (_beta_)\n")
            of.write("A visualization of the table above:  \n")
            of.write(f"![{filen}_{version}](./{filen}_{version}.png)\n")


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
    writeReportPage(**report, extra=input)


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
