#!/usr/bin/env python3
"""Regression tests for PLUMED with GROMACS.

    regtest.py --gmx PREFIX/bin/gmx_mpi [--kernel PLUMED/lib/libplumedKernel.so] [--tests REGEX] [--out DIR]

There are no reference files: every test checks a property that must hold for any GROMACS version, any
PLUMED version and any way of patching GROMACS, so the same tests run on all their combinations:
- an independent reference: a PLUMED RESTRAINT on a distance must act as the GROMACS pull code with the
  same parameters (same forces, virial and trajectory);
- invariances: the same simulation with or without ENERGY in the PLUMED input, or on 1, 2 or 3 ranks, must
  give the same results;
- exact relations: ENERGY equals the GROMACS potential energy, a bias on ENERGY rescales the forces, a rerun
  gives the CVs of the run, PLUMED's lambda is the lambda state of GROMACS...;
- protocols: checkpoints, stop requests (COMMITTOR), replica exchange, multiple walkers, shells, errors.
Each test prints PASS, FAIL (with the numbers) or SKIP (when the combination does not support the feature).
"""
import argparse
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SYSTEMS = HERE / "systems"
TESTS = []


class Fail(Exception):
    pass


class Skip(Exception):
    pass


def test(name):
    """Register a test; its docstring is the description"""
    def deco(f):
        TESTS.append((name, f))
        return f
    return deco


def check(cond, msg):
    if not cond:
        raise Fail(msg)


def maxdiff(a, b):
    return max((abs(x - y) for x, y in zip(a, b)), default=0.0)


def force_diff(a, b):
    """Largest difference of two force arrays, relative to the largest force: gmx dump prints 6 digits"""
    return maxdiff(a, b) / max(max((abs(x) for x in a), default=0.0), 1.0)


# ------------------------------------------------------------------------------------------------- files

def write(path, text):
    Path(path).write_text(text)


def colvar(path):
    """Read a PLUMED output file into {field: [values]}"""
    fields, rows = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#! FIELDS"):
            fields = line.split()[2:]
        elif line.strip() and not line.startswith("#"):
            rows.append([float(x) for x in line.split()])
    if fields is None:
        return {}
    return {f: [r[i] for r in rows] for i, f in enumerate(fields)}


def mdp_text(system, overrides):
    """The mdp file of SYSTEM with OVERRIDES (None removes a key)"""
    keys, values = [], {}
    for line in (SYSTEMS / system / "md.mdp").read_text().splitlines():
        line = line.split(";")[0]
        if "=" in line:
            k, v = (s.strip() for s in line.split("=", 1))
            k = k.replace("_", "-")
            if k not in values:
                keys.append(k)
            values[k] = v
    for k, v in overrides.items():
        k = k.replace("_", "-")
        if v is None:
            values.pop(k, None)
            continue
        if k not in values:
            keys.append(k)
        values[k] = str(v)
    return "".join(f"{k} = {values[k]}\n" for k in keys if k in values)


# ------------------------------------------------------------------------------------------------- GROMACS

class Gromacs:
    """One combination of a patched GROMACS and a PLUMED kernel"""

    def __init__(self, gmx, kernel, mpirun, out, timeout):
        gmx = str(Path(shutil.which(gmx) or gmx).resolve())
        kernel = kernel and str(Path(kernel).resolve())
        self.gmx, self.kernel, self.mpirun, self.out, self.timeout = gmx, kernel, mpirun, Path(out).resolve(), timeout
        r = subprocess.run([gmx, "--version"], capture_output=True, text=True, env=self.env())
        ver = r.stdout
        m = re.search(r"GROMACS version:\s+(\d{4})\.?(\d*)", ver)
        if not m:
            raise RuntimeError(f"cannot run {gmx} --version: " + " ".join((ver + r.stderr).split()[-40:]))
        self.version = (int(m.group(1)), int(m.group(2) or 0))
        self.cuda = bool(re.search(r"GPU support:\s+CUDA", ver))
        h = subprocess.run([gmx, "mdrun", "-h"], capture_output=True, text=True, env=self.env())
        # the legacy patches add -hrex (shown as -[no]hrex), the native module does not
        self.legacy = bool(re.search(r"-(\[no\])?hrex", h.stdout + h.stderr))
        self.out.mkdir(parents=True, exist_ok=True)

    def env(self, extra=None):
        e = dict(os.environ, OMP_NUM_THREADS="1")
        if self.kernel:
            e["PLUMED_KERNEL"] = self.kernel
        e.pop("GMX_DISABLE_MODULAR_SIMULATOR", None)
        e.update(extra or {})
        return e

    def dir(self, name):
        d = self.out / name
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True)
        return d

    def grompp(self, d, system, mdp, conf=None, ndx=None, out="topol.tpr"):
        write(d / "grompp.mdp", mdp_text(system, mdp))
        conf = SYSTEMS / system / (conf or next(p.name for p in sorted((SYSTEMS / system).iterdir())
                                              if p.suffix in (".gro", ".g96")))
        top = next((SYSTEMS / system).glob("*.top"))
        cmd = [self.gmx, "grompp", "-f", "grompp.mdp", "-c", str(conf), "-p", str(top), "-o", out, "-maxwarn", "5"]
        if ndx:
            write(d / "index.ndx", ndx)
            cmd += ["-n", "index.ndx"]
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True, env=self.env())
        if r.returncode != 0:
            text = r.stdout + r.stderr
            err = re.findall(r"(?:ERROR \d+ \[[^\]]*\]:|Fatal error:)\n((?:.+\n){1,3})", text)
            raise Fail("grompp failed: " + " | ".join(" ".join(e.split()) for e in err)[:400])

    def mdrun(self, d, np=1, plumed="plumed.dat", args=(), env=None, multidir=None, timeout=None):
        """Run mdrun in D (or in MULTIDIR subdirectories of D); return the exit code (124 on timeout).
        With a GPU build, everything runs on the CPU unless ARGS choose otherwise: GPU sums are not
        reproducible, and most tests compare two runs."""
        cmd = [*self.mpirun, "-np", str(np), self.gmx, "mdrun", "-ntomp", "1", "-deffnm", "md", "-s", "topol.tpr"]
        if self.cuda and "-nb" not in args:
            cmd += ["-nb", "cpu", "-pme", "cpu", "-bonded", "cpu", "-update", "cpu"]
        if multidir:
            cmd += ["-multidir", *multidir]
        if plumed:
            cmd += ["-plumed", plumed]
        cmd += list(args)
        with open(d / "mdrun.out", "w") as out:
            # in its own session, so that killing a hanging run does not reach this script (the ranks of Open MPI
            # are in their own process groups, and stop when mpirun is killed)
            p = subprocess.Popen(cmd, cwd=d, stdout=out, stderr=subprocess.STDOUT, env=self.env(env),
                                 start_new_session=True)
            try:
                return p.wait(timeout=timeout or self.timeout)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
                return 124

    def output(self, d):
        return (d / "mdrun.out").read_text(errors="replace")

    def frames(self, trr):
        """Frames of a trr file: [{'step', 'box', 'x', 'f'}], coordinates as flat lists"""
        r = subprocess.run([self.gmx, "dump", "-f", str(trr)], capture_output=True, text=True, env=self.env())
        frames, cur = [], None
        num = r"([-+\d.eE]+)"
        for line in r.stdout.splitlines():
            if re.match(r"\S.* frame \d+:", line):
                cur = {"x": [], "f": [], "box": []}
                frames.append(cur)
                continue
            if cur is None:
                continue
            m = re.search(r"step=\s*(\d+)", line)
            if m and "natoms" in line:
                cur["step"] = int(m.group(1))
                continue
            m = re.match(r"\s+(box|x|v|f)\[\s*\d+\]=\{\s*" + num + r",\s*" + num + r",\s*" + num + r"\}", line)
            if m and m.group(1) in cur:
                cur[m.group(1)].extend(float(m.group(i)) for i in (2, 3, 4))
        return frames

    def energy(self, d, terms):
        """{term: [values]} from md.edr in D"""
        xvg = d / "energy.xvg"
        subprocess.run([self.gmx, "energy", "-f", "md.edr", "-o", xvg.name], cwd=d, input="\n".join(terms) + "\n\n",
                       capture_output=True, text=True, env=self.env())
        names, rows = [], []
        for line in xvg.read_text().splitlines() if xvg.exists() else []:
            m = re.match(r'@ s\d+ legend "(.*)"', line)
            if m:
                names.append(m.group(1))
            elif line.strip() and line[0] not in "#@":
                rows.append([float(x) for x in line.split()])
        return {n: [r[i + 1] for r in rows] for i, n in enumerate(names)}

    def checkpoint_steps(self, d):
        return [int(s) for s in re.findall(r"Writing checkpoint, step (\d+)", (d / "md.log").read_text())]


def dd(np):
    """mdrun arguments for a reproducible domain decomposition over NP ranks"""
    return [] if np == 1 else ["-dd", str(np), "1", "1", "-npme", "0", "-dlb", "no"]


def water_index(extra):
    natoms = 3636
    groups = {"System": list(range(1, natoms + 1)), **extra}
    # grompp reads index files line by line with a limited buffer: write 15 atoms per line
    lines = lambda a: "".join(" ".join(map(str, a[i:i + 15])) + "\n" for i in range(0, len(a), 15))
    return "".join(f"[ {g} ]\n" + lines(a) for g, a in groups.items())


PULL = {"pull": "yes", "pull-ngroups": 2, "pull-ncoords": 1, "pull-group1-name": "A1", "pull-group2-name": "A2",
        "pull-coord1-type": "umbrella", "pull-coord1-geometry": "distance", "pull-coord1-groups": "1 2",
        "pull-coord1-dim": "Y Y Y", "pull-coord1-init": 0.3, "pull-coord1-rate": 0, "pull-coord1-k": 5000,
        "pull-coord1-start": "no", "pull-nstxout": 0, "pull-nstfout": 0}
RESTRAINT = "d: DISTANCE ATOMS={a},{b}\nRESTRAINT ARG=d AT=0.3 KAPPA=5000\nPRINT ARG=d STRIDE=1 FILE=COLVAR FMT=%.8f\n"
SPECTATOR = "e: ENERGY\nPRINT ARG=e STRIDE={stride} FILE=ENERGY FMT=%.8e\n"
MTS = {"mts": "yes", "mts-levels": 2, "mts-level2-factor": 2, "mts-level2-forces": "longrange-nonbonded"}


def compare_runs(env, d1, d2, what, xtol=1e-4, ftol=2e-5, vtol0=0.1, vtol=1e-4, force_steps=(0,)):
    """Compare positions (all frames), forces (FORCE_STEPS) and the virial (energy frames) of two runs. The virial is
    compared within VTOL0 kJ/mol at step 0 (that of a restraint between two water oxygens is about 1.7e3 kJ/mol), and
    then within VTOL of the total (8e4 kJ/mol for the water box), as round-off makes the two runs drift apart."""
    f1, f2 = env.frames(d1 / "md.trr"), env.frames(d2 / "md.trr")
    check(len(f1) == len(f2) and f1, f"{what}: {len(f1)} and {len(f2)} frames")
    dx = max(maxdiff(a["x"], b["x"]) for a, b in zip(f1, f2))
    df = max((force_diff(a["f"], b["f"]) for a, b in zip(f1, f2) if a.get("step") in force_steps and a["f"]), default=0)
    terms = ["Vir-XX", "Vir-YY", "Vir-ZZ"]
    e1, e2 = env.energy(d1, terms), env.energy(d2, terms)
    dv0 = max(abs(e1[t][0] - e2[t][0]) for t in terms)
    dv = max(maxdiff(e1[t], e2[t]) for t in terms)
    vmax = max(abs(x) for t in terms for x in e1[t])
    msg = (f"{what}: max |dx| {dx:.2g} nm, |df|/max|f| at step 0 {df:.2g}, |dvirial| at step 0 {dv0:.2g} kJ/mol, "
           f"max {dv:.3g} kJ/mol")
    check(dx < xtol and df < ftol and dv0 < vtol0 and dv < vtol * vmax, msg + f" (virial ~ {vmax:.3g})")
    return msg


def gromacs_error(text):
    """The message of the GROMACS error box in TEXT, if any"""
    m = re.search(r"\n((?:Fatal error|Internal error \(bug\)|Feature not implemented|Error in user input|"
                  r"Inconsistency in user input)[^\n]*:\n(?:.+\n){1,4})", text)
    return " ".join(m.group(1).split()) if m else " ".join(text.split()[-40:])


def run_ok(env, d, rc, what):
    if rc == 124:
        raise Fail(f"{what}: mdrun did not finish and was killed after the timeout")
    if rc != 0:
        raise Fail(f"{what}: mdrun exit {rc}: {gromacs_error(env.output(d))[:400]}")


# ------------------------------------------------------------------------------------------------- tests

@test("energy")
def t_energy(env):
    """ENERGY equals the potential energy of GROMACS at every step, on 1 and 2 ranks"""
    msgs = []
    for np in (1, 2):
        d = env.dir(f"energy-np{np}")
        env.grompp(d, "water", {"nsteps": 20, "nstcalcenergy": 1, "nstenergy": 1, "nstfout": 0})
        write(d / "plumed.dat", "e: ENERGY\nPRINT ARG=e STRIDE=1 FILE=COLVAR FMT=%.10e\n")
        run_ok(env, d, env.mdrun(d, np, args=dd(np)), f"{np} ranks")
        pl, md = colvar(d / "COLVAR")["e"], env.energy(d, ["Potential"])["Potential"]
        check(len(pl) == len(md) == 21, f"{np} ranks: {len(pl)} ENERGY values, {len(md)} energies")
        rel = max(abs(a - b) / abs(b) for a, b in zip(pl, md))
        check(rel < 1e-6, f"{np} ranks: ENERGY differs from the potential energy by {rel:.2g} (relative)")
        msgs.append(f"{np} ranks: {rel:.1g}")
    return "relative difference " + ", ".join(msgs)


def restraint_vs_pull(env, tag, np=1, plumed_extra="", mdp_extra=None, atoms=(1, 5)):
    mdp = {"nsteps": 20, "nstxout": 1, "nstfout": 1, "nstcalcenergy": 5, "nstenergy": 5, "nstlog": 0, **(mdp_extra or {})}
    ndx = water_index({"A1": [atoms[0]], "A2": [atoms[1]]})
    d1, d2 = env.dir(f"{tag}-plumed"), env.dir(f"{tag}-pull")
    env.grompp(d1, "water", mdp, ndx=ndx)
    env.grompp(d2, "water", {**mdp, **PULL}, ndx=ndx)
    write(d1 / "plumed.dat", RESTRAINT.format(a=atoms[0], b=atoms[1]) + plumed_extra)
    run_ok(env, d1, env.mdrun(d1, np, args=dd(np)), "PLUMED")
    run_ok(env, d2, env.mdrun(d2, np, plumed=None, args=dd(np)), "pull code")
    return compare_runs(env, d1, d2, f"PLUMED against pull code, {np} rank{'s' * (np > 1)}")


@test("pull")
def t_pull(env):
    """a RESTRAINT on the distance of two water oxygens acts as the GROMACS pull code: same positions over 20
    steps, forces at step 0, and virial every 5 steps (the virial is computed on those steps only)"""
    return restraint_vs_pull(env, "pull-np1")


@test("pull-dd")
def t_pull_dd(env):
    """as pull, on 2 ranks"""
    return restraint_vs_pull(env, "pull-np2", np=2)


@test("pull-energy")
def t_pull_energy(env):
    """as pull, with ENERGY in the PLUMED input (PLUMED then runs when the energy is complete)"""
    return restraint_vs_pull(env, "pull-energy", plumed_extra=SPECTATOR.format(stride=1))


@test("pull-mts")
def t_pull_mts(env):
    """as pull, with multiple time stepping and ENERGY needed on the slow steps"""
    if env.version < (2021, 0):
        raise Skip("no multiple time stepping before GROMACS 2021")
    return restraint_vs_pull(env, "pull-mts", plumed_extra=SPECTATOR.format(stride=2),
                             mdp_extra={**MTS, "nstcalcenergy": 4, "nstenergy": 4, "nstfout": 2})


@test("mts-fast-step")
def t_mts_fast(env):
    """with multiple time stepping, mdrun stops with an error if PLUMED needs the energy on a fast step"""
    if env.version < (2021, 0):
        raise Skip("no multiple time stepping before GROMACS 2021")
    d = env.dir("mts-fast-step")
    env.grompp(d, "water", {**MTS, "nsteps": 4, "nstcalcenergy": 2, "nstenergy": 2, "nstfout": 0})
    write(d / "plumed.dat", SPECTATOR.format(stride=1))
    rc = env.mdrun(d)
    check(rc != 0 and "multiple time stepping" in env.output(d), f"mdrun exit {rc}, with no error about MTS")
    return "stops with an error"


def path_invariance(env, tag, atoms, np):
    """the same run with and without ENERGY in the PLUMED input"""
    mdp = {"nsteps": 20, "nstxout": 1, "nstfout": 1, "nstcalcenergy": 5, "nstenergy": 5, "nstlog": 0}
    d1, d2 = env.dir(f"{tag}-plain"), env.dir(f"{tag}-energy")
    for d, extra in ((d1, ""), (d2, SPECTATOR.format(stride=1))):
        env.grompp(d, "water", mdp)
        write(d / "plumed.dat", RESTRAINT.format(a=atoms[0], b=atoms[1]) + extra)
        run_ok(env, d, env.mdrun(d, np, args=dd(np)), tag)
    return compare_runs(env, d1, d2, f"with ENERGY against without, {np} rank{'s' * (np > 1)}")


@test("vsite")
def t_vsite(env):
    """a RESTRAINT between two virtual sites (TIP4P M sites) acts the same with and without ENERGY in the PLUMED
    input (with ENERGY, PLUMED runs after GROMACS spread the forces on the virtual sites), on 1 and 2 ranks"""
    return "; ".join(path_invariance(env, f"vsite-np{np}", (4, 8), np) for np in (1, 2))


@test("energy-bias")
def t_energy_bias(env):
    """a bias s*ENERGY multiplies all the forces by 1+s, also with virtual sites"""
    d0, d1 = env.dir("energy-bias-ref"), env.dir("energy-bias")
    for d in (d0, d1):
        env.grompp(d, "water", {"nsteps": 0, "nstfout": 1})
    write(d0 / "plumed.dat", "e: ENERGY\nPRINT ARG=e FILE=COLVAR\n")
    write(d1 / "plumed.dat", "e: ENERGY\nRESTRAINT ARG=e AT=0 KAPPA=0 SLOPE=-0.5\nPRINT ARG=e FILE=COLVAR\n")
    run_ok(env, d0, env.mdrun(d0), "reference")
    run_ok(env, d1, env.mdrun(d1), "bias")
    f0, f1 = env.frames(d0 / "md.trr")[0]["f"], env.frames(d1 / "md.trr")[0]["f"]
    err = force_diff([0.5 * a for a in f0], f1)
    check(err < 2e-5, f"forces differ from 0.5 times the unbiased ones by up to {err:.3g} of the largest force")
    return f"forces are 0.5 times the unbiased ones within {err:.2g} of the largest force"


@test("energy-bias-dd")
def t_energy_bias_dd(env):
    """with a bias on ENERGY, the virial is the same on 1, 2 and 3 ranks"""
    vir = {}
    for np in (1, 2, 3):
        d = env.dir(f"energy-bias-dd-np{np}")
        env.grompp(d, "alanine", {"nsteps": 0}, conf="conf.0.gro")
        write(d / "plumed.dat", "e: ENERGY\nRESTRAINT ARG=e AT=0 KAPPA=0 SLOPE=-0.5\n")
        run_ok(env, d, env.mdrun(d, np, args=dd(np)), f"{np} ranks")
        e = env.energy(d, ["Vir-XX", "Vir-YY", "Vir-ZZ"])
        vir[np] = [e[t][0] for t in ("Vir-XX", "Vir-YY", "Vir-ZZ")]
    dv = max(maxdiff(vir[1], vir[np]) for np in (2, 3))
    check(dv < 0.05, "virial at step 0 on 1/2/3 ranks: " + " / ".join(f"{vir[n][0]:.4g}" for n in (1, 2, 3)) + " (XX)")
    return f"same virial within {dv:.2g} kJ/mol"


def em_distances(env, tag, atoms, np, extra=""):
    d = env.dir(tag)
    env.grompp(d, "water", {"integrator": "steep", "nsteps": 20, "nstxout": 0, "nstfout": 0, "tcoupl": "no",
                            "gen-vel": "no", "nstcalcenergy": 1, "nstenergy": 1})
    write(d / "plumed.dat", RESTRAINT.format(a=atoms[0], b=atoms[1]) + extra)
    run_ok(env, d, env.mdrun(d, np, args=dd(np)), tag)
    return colvar(d / "COLVAR")["d"]


@test("em")
def t_em(env):
    """energy minimization with a restraint: the same steps on 1 and 2 ranks, and with and without ENERGY
    (restraints between two oxygens and between two virtual sites)"""
    msgs = []
    for name, atoms in (("oxygens", (1, 5)), ("vsites", (4, 8))):
        ref = em_distances(env, f"em-{name}-np1", atoms, 1)
        for tag, np, extra in ((f"em-{name}-np2", 2, ""), (f"em-{name}-energy", 1, SPECTATOR.format(stride=1)),
                               (f"em-{name}-energy-np2", 2, SPECTATOR.format(stride=1))):
            dist = em_distances(env, tag, atoms, np, extra)
            diff = maxdiff(ref, dist)
            # the minimizer takes a step depending on energies that 1 and 2 ranks sum in a different order: the
            # distances can differ by 1e-4 nm, the bugs this checks give other atoms or NaN
            check(len(dist) == len(ref) and diff < 1e-3,
                  f"{name}, {tag[len('em-' + name) + 1:]}: {len(dist)} steps, distances differ by {diff:.3g} nm from 1 rank")
        msgs.append(f"{name}: same {len(ref)} distances")
    return "; ".join(msgs)


@test("data")
def t_data(env):
    """PLUMED gets the positions, box, masses, charges and time of GROMACS"""
    d = env.dir("data")
    env.grompp(d, "alanine", {"nsteps": 10, "nstxout": 10}, conf="conf.0.gro")
    write(d / "plumed.dat", "DUMPATOMS ATOMS=1-22 FILE=pos.xyz STRIDE=10 PRECISION=6\nDUMPMASSCHARGE FILE=mq\n"
                            "t: TIME\nc: CELL\nPRINT ARG=t,c.ax,c.by,c.cz STRIDE=1 FILE=COLVAR FMT=%.6f\n")
    run_ok(env, d, env.mdrun(d), "run")
    fr = env.frames(d / "md.trr")
    xyz = [l.split() for l in (d / "pos.xyz").read_text().splitlines()]
    pos = [float(v) for l in xyz if len(l) == 4 for v in l[1:]]
    n = len(fr[0]["x"])
    dx = max(maxdiff(pos[:n], fr[0]["x"]), maxdiff(pos[n:2 * n], fr[-1]["x"]))
    check(dx < 1e-4, f"positions differ from the trajectory by {dx:.3g} nm")
    cv = colvar(d / "COLVAR")
    dt = maxdiff(cv["t"], [0.002 * i for i in range(11)])
    db = maxdiff([cv["c.ax"][0], cv["c.by"][0], cv["c.cz"][0]], [fr[0]["box"][0], fr[0]["box"][4], fr[0]["box"][8]])
    check(dt < 1e-6 and db < 1e-5, f"time differs by {dt:.3g} ps, box by {db:.3g} nm")
    dump = subprocess.run([env.gmx, "dump", "-s", "topol.tpr"], cwd=d, capture_output=True, text=True,
                          env=env.env()).stdout
    tpr = [(float(m.group(1)), float(m.group(2))) for m in re.finditer(r"atom\[\s*\d+\]=\{type=.*?m=\s*([-\d.e+]+),\s*q=\s*([-\d.e+]+)", dump)]
    mq = [l.split() for l in (d / "mq").read_text().splitlines() if not l.startswith("#")]
    check(len(mq) == 22 and len(tpr) >= 22, f"{len(mq)} masses/charges from PLUMED, {len(tpr)} atoms in the tpr")
    dm = max(abs(float(a[1]) - b[0]) for a, b in zip(mq, tpr))
    dq = max(abs(float(a[2]) - b[1]) for a, b in zip(mq, tpr))
    check(dm < 1e-4 and dq < 1e-5, f"masses differ by {dm:.3g}, charges by {dq:.3g}")
    return "positions, box, time, masses and charges match"


@test("checkpoint")
def t_checkpoint(env):
    """PLUMED writes its state (OPES STATE_WFILE with STORE_STATES) at the steps where GROMACS writes a
    checkpoint, on 1 and 2 ranks"""
    msgs = []
    for np in (1, 2):
        d = env.dir(f"checkpoint-np{np}")
        env.grompp(d, "alanine", {"nsteps": 20000}, conf="conf.0.gro")
        write(d / "plumed.dat", "phi: TORSION ATOMS=5,7,9,15\nopes: OPES_METAD ARG=phi PACE=1 BARRIER=20 "
                                "SIGMA=0.3 STATE_WFILE=STATE STORE_STATES\n")
        run_ok(env, d, env.mdrun(d, np, args=dd(np) + ["-cpt", "0.001"]), f"{np} ranks")
        cpt = env.checkpoint_steps(d)
        counters = [int(float(m.group(1))) for m in re.finditer(r"#! SET counter\s+(\S+)", (d / "STATE").read_text())]
        check(len(cpt) >= 2 and counters, f"{np} ranks: {len(cpt)} checkpoints, {len(counters)} states")
        # with PACE=1 the counter grows by one per step: the last state is written at the last step
        last = 20000
        steps = [last - (counters[-1] - c) for c in counters]
        missing = [s for s in steps if s not in cpt + [last]]
        check(not missing, f"{np} ranks: states at steps {steps}, checkpoints at {cpt}")
        msgs.append(f"{np} ranks: {len(steps)} states at the {len(cpt)} checkpoint steps")
    return "; ".join(msgs)


# COMMITTOR on the time commits at a known step, whatever the trajectory
COMMITTOR = ("t: TIME\nphi: TORSION ATOMS=5,7,9,15\n"
             "COMMITTOR ARG=t STRIDE=1 BASIN_LL1={t} BASIN_UL1=1e6 FILE=COMMIT\n"
             "PRINT ARG=phi STRIDE=1 FILE=COLVAR\n")


def first_commit_step(d):
    """Step of the first commit written by COMMITTOR in D, or None"""
    for f in sorted(d.glob("COMMIT*")):
        for line in f.read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                return round(float(line.split()[0]) / 0.002)
    return None


def last_step(path):
    t = colvar(path).get("time", [])
    return round(t[-1] / 0.002) if t else None


@test("committor")
def t_committor(env):
    """COMMITTOR stops mdrun soon after committing, with a final checkpoint, on 1 and 2 ranks"""
    msgs = []
    for np in (1, 2):
        d = env.dir(f"committor-np{np}")
        env.grompp(d, "alanine", {"nsteps": 20000}, conf="conf.0.gro")
        write(d / "plumed.dat", COMMITTOR.format(t=0.4999))  # step 250 (0.5 ps is not exact)
        run_ok(env, d, env.mdrun(d, np, args=dd(np), timeout=300), f"{np} ranks")
        c, last = first_commit_step(d), last_step(d / "COLVAR")
        check(c == 250, f"{np} ranks: committed at step {c}, instead of 250")
        check(last < 20000 and last - c <= 200, f"{np} ranks: committed at step {c}, ran to step {last}")
        cpt = env.checkpoint_steps(d)
        check(cpt and cpt[-1] == last, f"{np} ranks: last checkpoint at {cpt[-1:]} , last step {last}")
        msgs.append(f"{np} ranks: committed at {c}, stopped at {last} with a checkpoint")
    return "; ".join(msgs)


@test("committor-walkers")
def t_committor_walkers(env):
    """with two simulations (-multidir), all of them stop when one commits"""
    d = env.dir("committor-walkers")
    for i in (0, 1):
        (d / f"w{i}").mkdir()
        env.grompp(d / f"w{i}", "alanine", {"nsteps": 20000}, conf=f"conf.{i}.gro")
        # only the second simulation commits
        write(d / f"w{i}" / "plumed.dat", COMMITTOR.format(t=0.4999 if i else 1e5))
    rc = env.mdrun(d, 2, multidir=["w0", "w1"], timeout=120)
    check(rc != 124, "mdrun did not stop (killed after 120 s)")
    run_ok(env, d, rc, "walkers")
    lasts = [last_step(next((d / f"w{i}").glob("COLVAR*"))) for i in (0, 1)]
    check(lasts[0] == lasts[1] and lasts[0] < 20000, f"last steps {lasts}")
    return f"both stopped at step {lasts[0]}"


@test("rerun")
def t_rerun(env):
    """a rerun of a trajectory gives the CVs that PLUMED computed in the run, on 1 and 2 ranks; a stop request
    (COMMITTOR) ends the rerun"""
    d0 = env.dir("rerun-md")
    env.grompp(d0, "alanine", {"nsteps": 200, "nstxout": 20, "nstvout": 0, "nstfout": 0}, conf="conf.0.gro")
    cvs = "phi: TORSION ATOMS=5,7,9,15\npsi: TORSION ATOMS=7,9,15,17\nPRINT ARG=phi,psi STRIDE={s} FILE=COLVAR FMT=%.8f\n"
    write(d0 / "plumed.dat", cvs.format(s=20))
    run_ok(env, d0, env.mdrun(d0), "run")
    ref = colvar(d0 / "COLVAR")
    msgs = []
    for np in (1, 2):
        d = env.dir(f"rerun-np{np}")
        shutil.copy(d0 / "topol.tpr", d)
        write(d / "plumed.dat", cvs.format(s=1))
        run_ok(env, d, env.mdrun(d, np, args=dd(np) + ["-rerun", str(d0 / "md.trr")]), f"rerun on {np} ranks")
        cv = colvar(d / "COLVAR")
        diff = max(maxdiff(ref["phi"], cv.get("phi", [])), maxdiff(ref["psi"], cv.get("psi", [])))
        check(len(cv.get("phi", [])) == len(ref["phi"]) and diff < 1e-5,
              f"rerun on {np} ranks: {len(cv.get('phi', []))} of {len(ref['phi'])} frames, CVs differ by {diff:.3g}")
        msgs.append(f"{np} ranks: same CVs")
    d = env.dir("rerun-stop")
    shutil.copy(d0 / "topol.tpr", d)
    write(d / "plumed.dat", COMMITTOR.format(t=0.0))
    run_ok(env, d, env.mdrun(d, 1, args=["-rerun", str(d0 / "md.trr")]), "rerun with COMMITTOR")
    t = colvar(d / "COLVAR").get("time", [])
    commit = [float(line.split()[0]) for line in (d / "COMMIT").read_text().splitlines()
              if line.strip() and not line.startswith("#")] if (d / "COMMIT").exists() else []
    check(commit and t and t[-1] == commit[0],
          f"rerun with COMMITTOR: committed at time {commit[0] if commit else None}, stopped after {len(t)} of 11 "
          f"frames (last time {t[-1] if t else None})")
    return "; ".join(msgs) + "; COMMITTOR stops it"


def hrex_run(env, tag, system, nsteps, plumed, np_per_replica=1, replex=10):
    d = env.dir(tag)
    for i in (0, 1):
        (d / f"w{i}").mkdir()
        conf = "conf.0.gro" if system == "alanine" else None  # identical replicas
        env.grompp(d / f"w{i}", system, {"nsteps": nsteps, "nstfout": 0, "nstxout": 0}, conf=conf)
        write(d / f"w{i}" / "plumed.dat", plumed)
    rc = env.mdrun(d, 2 * np_per_replica, multidir=["w0", "w1"], args=["-replex", str(replex), "-hrex"] + (
        ["-dd", str(np_per_replica), "1", "1", "-npme", "0", "-dlb", "no"] if np_per_replica > 1 else []))
    run_ok(env, d, rc, tag)
    log = (d / "w0" / "md.log").read_text()
    dpl = [float(x) for x in re.findall(r"dplumed =\s*(\S+)", log)]
    return d, dpl


@test("hrex")
def t_hrex(env):
    """-hrex with two identical replicas: the PLUMED energy difference is 0 (1 and 2 ranks per replica), and
    PLUMED does not calculate again on the exchange steps (one COLVAR line per step, the expected hills)"""
    if not env.legacy:
        raise Skip("-hrex is only in the legacy patches")
    msgs = []
    for tag, system, nsteps, npr in (("hrex-np1", "alanine", 100, 1), ("hrex-np2", "water", 120, 2)):
        d, dpl = hrex_run(env, tag, system, nsteps, "d: DISTANCE ATOMS=1,5\nPRINT ARG=d STRIDE=10 FILE=COLVAR\n", npr,
                          replex=10 if npr == 1 else 15)
        check(dpl and max(abs(x) for x in dpl) < 1e-3, f"{npr} rank(s) per replica: dplumed {dpl}")
        msgs.append(f"{npr} rank(s)/replica: max |dplumed| {max(abs(x) for x in dpl):.1g}")
    plumed = ("phi: TORSION ATOMS=5,7,9,15\nmetad: METAD ARG=phi PACE=3 HEIGHT=1.0 SIGMA=0.3 FILE=HILLS\n"
              "PRINT ARG=phi,metad.bias STRIDE=1 FILE=COLVAR\n")
    d, dpl = hrex_run(env, "hrex-metad", "alanine", 40, plumed)
    times = colvar(next((d / "w0").glob("COLVAR*")))["time"]
    hills = colvar(next((d / "w0").glob("HILLS*")))["time"]
    check(len(times) == len(set(times)) == 41 and len(hills) == 13,
          f"METAD, exchanges every 10 steps: {len(times)} COLVAR lines ({len(set(times))} distinct), {len(hills)} hills")
    return "; ".join(msgs) + "; no extra PLUMED calculation"


@test("shells")
def t_shells(env):
    """shell (polarizable) particles: PLUMED is applied once per step, also in a rerun, and mdrun warns when PLUMED
    needs the energy, which is the one before the shells are relaxed (or mdrun refuses shells)"""
    d = env.dir("shells")
    env.grompp(d, "shell", {})
    write(d / "plumed.dat", "d: DISTANCE ATOMS=1,6\nmetad: METAD ARG=d PACE=3 HEIGHT=1.0 SIGMA=0.05 FILE=HILLS\n"
                            "PRINT ARG=d,metad.bias STRIDE=1 FILE=COLVAR\n")
    rc = env.mdrun(d)
    if env.legacy:
        check(rc != 0 and "shell" in env.output(d), f"mdrun exit {rc}: PLUMED is not applied with shells "
              f"({len(colvar(d / 'COLVAR').get('time', [])) if (d / 'COLVAR').exists() else 0} COLVAR lines)")
        return "refused with an error"
    run_ok(env, d, rc, "shells")
    t = colvar(d / "COLVAR")["time"]
    h = colvar(d / "HILLS")["time"]
    check(len(t) == len(set(t)) == 21 and len(h) == 6, f"{len(t)} COLVAR lines ({len(set(t))} distinct), {len(h)} hills")
    # the starting structure, whose shells are not relaxed: the rerun relaxes them in a few iterations
    d2 = env.dir("shells-rerun")
    shutil.copy(d / "topol.tpr", d2)
    write(d2 / "plumed.dat", "d: DISTANCE ATOMS=1,6\nPRINT ARG=d STRIDE=1 FILE=COLVAR\n")
    run_ok(env, d2, env.mdrun(d2, args=["-rerun", str(SYSTEMS / "shell" / "sw-dimer.g96")]), "rerun")
    n = len(colvar(d2 / "COLVAR")["time"])
    check(n == 1, f"rerun of one frame: {n} COLVAR lines")
    d3 = env.dir("shells-energy")
    env.grompp(d3, "shell", {"nsteps": 2})
    write(d3 / "plumed.dat", "e: ENERGY\nPRINT ARG=e STRIDE=1 FILE=COLVAR\n")
    run_ok(env, d3, env.mdrun(d3), "with ENERGY")
    check("before the shells are relaxed" in " ".join((d3 / "md.log").read_text().split()),
          "with ENERGY: no warning that PLUMED gets the energy before the shells are relaxed")
    return "one PLUMED update per step, also in a rerun; a warning with ENERGY"


EXPANDED = {"integrator": "md-vv", "nsteps": 200, "nstcalcenergy": 10, "free-energy": "expanded",
            "couple-moltype": "Protein", "couple-lambda0": "vdw-q", "couple-lambda1": "vdw", "couple-intramol": "yes",
            "init-lambda-state": 0, "coul-lambdas": "0.0 0.005 0.01 0.015", "nstdhdl": 10, "nstexpanded": 10,
            "lmc-stats": "no", "lmc-move": "metropolis", "init-lambda-weights": "0 0 0 0", "nstlog": 10,
            "lmc-seed": 1234}


@test("expanded")
def t_expanded(env):
    """expanded ensemble: PLUMED's lambda is the lambda state of GROMACS; with lmc-stats other than 'no',
    GROMACS still updates the weights"""
    if not env.legacy:
        raise Skip("lambda dynamics with PLUMED is only in the legacy patches")
    d = env.dir("expanded")
    env.grompp(d, "alanine", EXPANDED, conf="conf.0.gro")
    write(d / "plumed.dat", "l: EXTRACV NAME=lambda\nPRINT ARG=l STRIDE=10 FILE=COLVAR\n")
    run_ok(env, d, env.mdrun(d, env={"GMX_DISABLE_MODULAR_SIMULATOR": "1"}), "lmc-stats = no")
    state = [round(r) for r in env_dhdl_states(d)]
    lam = [round(x) for x in colvar(d / "COLVAR")["l"]]
    check(state and len(state) == len(lam) and len(set(state)) > 1, f"{len(state)} states {set(state)}, {len(lam)} lambdas")
    bad = sum(a != b for a, b in zip(state, lam))
    check(bad == 0, f"PLUMED's lambda differs from GROMACS's state in {bad} of {len(state)} samples")
    d = env.dir("expanded-weights")
    env.grompp(d, "alanine", {**EXPANDED, "lmc-stats": "metropolis-transition"}, conf="conf.0.gro")
    write(d / "plumed.dat", "l: EXTRACV NAME=lambda\nPRINT ARG=l STRIDE=10 FILE=COLVAR\n")
    run_ok(env, d, env.mdrun(d, env={"GMX_DISABLE_MODULAR_SIMULATOR": "1"}), "metropolis-transition")
    blocks = (d / "md.log").read_text().split("MC-lambda information")
    weights = [float(l.split()[-2]) for l in blocks[-1].splitlines()[2:6] if len(l.split()) > 3]
    check(any(abs(w) > 0 for w in weights), f"weights with metropolis-transition stay {weights}")
    return f"same lambda in {len(state)} samples; GROMACS updates the weights"


def env_dhdl_states(d):
    states = []
    for line in (d / "md.xvg").read_text().splitlines():
        if line.strip() and line[0] not in "#@":
            states.append(float(line.split()[1]))
    return states


@test("md-vv")
def t_mdvv(env):
    """md-vv runs with PLUMED (the legacy patches use the legacy simulator, which supports it)"""
    d = env.dir("md-vv")
    env.grompp(d, "alanine", {"integrator": "md-vv", "nsteps": 50}, conf="conf.0.gro")
    write(d / "plumed.dat", "phi: TORSION ATOMS=5,7,9,15\nPRINT ARG=phi STRIDE=10 FILE=COLVAR\n")
    run_ok(env, d, env.mdrun(d), "md-vv")
    return "runs"


@test("replex-nex")
def t_nex(env):
    """replica exchange with -nex, which ignores the PLUMED bias, is refused"""
    if not env.legacy:
        raise Skip("replica exchange with a PLUMED bias is only in the legacy patches")
    d = env.dir("replex-nex")
    for i in (0, 1):
        (d / f"w{i}").mkdir()
        env.grompp(d / f"w{i}", "alanine", {"nsteps": 50}, conf=f"conf.{i}.gro")
        write(d / f"w{i}" / "plumed.dat", "phi: TORSION ATOMS=5,7,9,15\nPRINT ARG=phi STRIDE=10 FILE=COLVAR\n")
    rc = env.mdrun(d, 2, multidir=["w0", "w1"], args=["-replex", "10", "-nex", "4"])
    check(rc != 0, "mdrun ran with -nex and PLUMED")
    return "refused"


@test("replex-bias")
def t_replex_bias(env):
    """temperature replica exchange with a PLUMED bias: the legacy patches add the bias to the exchange
    criterion; the native module refuses it with a clear error"""
    d = env.dir("replex-bias")
    for i, temp in ((0, 300), (1, 320)):
        (d / f"w{i}").mkdir()
        env.grompp(d / f"w{i}", "alanine", {"nsteps": 200, "ref-t": temp, "gen-temp": temp}, conf=f"conf.{i}.gro")
        write(d / f"w{i}" / "plumed.dat", "phi: TORSION ATOMS=5,7,9,15\nmetad: METAD ARG=phi PACE=5 HEIGHT=1.0 "
                                          "SIGMA=0.3 FILE=HILLS\n")
    rc = env.mdrun(d, 2, multidir=["w0", "w1"], args=["-replex", "20"])
    out = env.output(d)
    if env.legacy:
        run_ok(env, d, rc, "replica exchange")
        n = len(re.findall(r"dplumed", (d / "w0" / "md.log").read_text()))
        check(n > 0, "no PLUMED bias in the exchange criterion")
        return f"{n} exchange attempts with the PLUMED bias"
    check(rc != 0, "mdrun ran replica exchange without the PLUMED bias in the exchange criterion, and without an error")
    check("not compatible" in out and "Internal error" not in out, f"mdrun exit {rc}: {gromacs_error(out)[:300]}")
    return "refused with a clear error"


@test("gpu-buffer-ops")
def t_gpu(env):
    """with nonbonded and PME on the GPU, the restraint acts the same with GPU force buffer operations"""
    if not env.cuda:
        raise Skip("not a CUDA build")
    res = {}
    for tag, extra in (("cpu-buffer", {}), ("gpu-buffer", {"GMX_USE_GPU_BUFFER_OPS": "1"})):
        d = env.dir(f"gpu-{tag}")
        env.grompp(d, "water", {"nsteps": 50, "nstcalcenergy": 100, "nstenergy": 100, "nstfout": 0, "nstlog": 0})
        write(d / "plumed.dat", RESTRAINT.format(a=1, b=5))
        run_ok(env, d, env.mdrun(d, 1, env=extra, args=["-nb", "gpu", "-pme", "gpu", "-update", "cpu"]), tag)
        res[tag] = colvar(d / "COLVAR")["d"]
    diff = maxdiff(res["cpu-buffer"], res["gpu-buffer"])
    check(diff < 1e-3, f"restrained distance after 50 steps: {res['cpu-buffer'][-1]:.4f} nm, "
                       f"{res['gpu-buffer'][-1]:.4f} nm with GPU buffer ops")
    return f"same distances within {diff:.2g} nm"


# ------------------------------------------------------------------------------------------------- main

MPIRUN = "mpirun --oversubscribe --bind-to none"


def run(gmx, kernel=None, out="regtests", tests=".", mpirun=MPIRUN, timeout=120, log=print):
    """Run the tests (those whose name matches the regular expression TESTS) in the directory OUT, with GMX (the
    gmx_mpi of a patched GROMACS, or a command that runs it with its arguments) and the PLUMED KERNEL (None: the
    one that GMX loads by itself). Return {test: {"status", "message", "description", "seconds"}}."""
    try:
        env = Gromacs(gmx, kernel, mpirun.split(), out, timeout)
    except Exception as e:
        return {"gromacs": {"status": "FAIL", "message": str(e), "description": "GROMACS runs", "seconds": 0}}
    log(f"# GROMACS {env.version[0]}.{env.version[1]}, {'legacy patch' if env.legacy else 'native module'}"
        f"{', CUDA' if env.cuda else ''}; {gmx}; kernel {kernel or '(loaded by GROMACS)'}")
    results = {}
    for name, f in TESTS:
        if not re.search(tests, name):
            continue
        t0 = time.time()
        try:
            status, msg = "PASS", f(env)
        except Skip as e:
            status, msg = "SKIP", str(e)
        except Fail as e:
            status, msg = "FAIL", str(e)
        except Exception as e:  # a bug of the harness or an unexpected output: report it as a failure
            status, msg = "FAIL", f"{type(e).__name__}: {e}"
        results[name] = {"status": status, "message": msg, "description": " ".join(f.__doc__.split()),
                         "seconds": round(time.time() - t0, 1)}
        log(f"{status:4s}  {name:18s} {msg}  [{results[name]['seconds']} s]")
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gmx", help="the gmx_mpi of a patched GROMACS (a path, or a command in the PATH)")
    ap.add_argument("--kernel", help="libplumedKernel.so of the PLUMED to test (default: PLUMED_KERNEL, or the "
                                     "PLUMED that GROMACS is linked to)")
    ap.add_argument("--mpirun", default=MPIRUN, help="MPI launcher")
    ap.add_argument("--tests", default=".", help="regular expression selecting the tests")
    ap.add_argument("--out", default="regtests", help="directory for the runs and results.json")
    ap.add_argument("--timeout", type=int, default=120, help="seconds before a run is killed as hanging")
    ap.add_argument("--list", action="store_true", help="list the tests and exit")
    a = ap.parse_args()
    if a.list:
        for name, f in TESTS:
            print(f"{name:18s} {' '.join(f.__doc__.split())}")
        return 0
    if not a.gmx:
        ap.error("--gmx is required")
    if os.environ.get("MPI_HOME"):  # an MPI that is not in the PATH
        mpi = os.environ["MPI_HOME"]
        os.environ["PATH"] = f"{mpi}/bin:{os.environ['PATH']}"
        os.environ["LD_LIBRARY_PATH"] = f"{mpi}/lib:{os.environ.get('LD_LIBRARY_PATH', '')}"
    results = run(a.gmx, a.kernel, a.out, a.tests, a.mpirun, a.timeout, log=lambda m: print(m, flush=True))
    Path(a.out).mkdir(parents=True, exist_ok=True)
    (Path(a.out) / "results.json").write_text(json.dumps(results, indent=1))
    return 1 if any(r["status"] == "FAIL" for r in results.values()) else 0


if __name__ == "__main__":
    sys.exit(main())
