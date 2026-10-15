# Regression tests for PLUMED with GROMACS

The general tests of the testcenter check that GROMACS passes positions, masses, charges, the timestep and the
energy to PLUMED, and that a bias on the volume gives the right pressure. The GROMACS patches are large and differ
for each GROMACS version, and these tests check the rest of the interface: forces on atoms and their virial,
virtual sites, multiple time stepping, domain decomposition, checkpoints, stop requests, reruns, energy
minimization, replica exchange, multiple walkers, expanded ensemble and shells.

The testcenter runs them (`runRegtests` in `mdcode.py`) on `gromacs`, patched with the patch of each PLUMED
version, and on `gromacs-vanilla`, the PLUMED module of GROMACS without patches. `gromacs-vanilla` fails the
tests of the features that only the patches of PLUMED add, e.g. `ENERGY`, `COMMITTOR` and multiple walkers.

## Approach

No test stores reference results: every test checks a property that must hold for any GROMACS version, any
PLUMED version and any way of patching GROMACS. The same tests then run on all their combinations.

- **Independent reference.** A PLUMED `RESTRAINT` on the distance between two atoms must act as the GROMACS
  pull code (an umbrella on the same distance, with the same constant): the same trajectory, the same forces and
  the same virial (`pull*` tests). The pull code is part of GROMACS and does not depend on PLUMED or on the patch.
- **Invariances.** The same simulation must give the same results with and without `ENERGY` in the PLUMED input
  (with it, PLUMED runs after the energy is complete, a different code path: `pull-energy`, `vsite`, `em`), and
  on 1, 2 or 3 ranks (`pull-dd`, `energy-bias-dd`, `em`, `rerun`, `hrex`).
- **Exact relations.** `ENERGY` equals the potential energy of GROMACS; a bias s·`ENERGY` multiplies all the
  forces by 1+s; a rerun gives the CVs that PLUMED computed during the run; PLUMED's `lambda` is the lambda state
  of GROMACS; two identical replicas have a zero exchange energy.
- **Protocols.** PLUMED writes its state at the checkpoints of GROMACS; `COMMITTOR` stops mdrun, with a final
  checkpoint, also with several simulations; shells, `-nex` and energies on MTS fast steps are refused with an
  error where they are not supported.

Comparisons use tolerances that cover single precision and the 6 digits printed by `gmx dump`; the bugs these
tests were written for give differences that are orders of magnitude larger (e.g. a double-counted virial of
1.7e3 kJ/mol against a tolerance of 0.1 kJ/mol at the first step). A test is skipped when the combination does not
have the feature, e.g. `-hrex` with the native PLUMED module of GROMACS 2025 and later.

`./regtest.py --list` lists the tests and what each one checks.

## Running the tests on your computer

The tests need a GROMACS built with MPI (`gmx_mpi`) and `mpirun` with `--oversubscribe` and `--bind-to none` (as in
Open MPI; `--mpirun` sets another launcher). They work with any GROMACS version that PLUMED can patch (2020 and
later, with the patches in runtime mode or linked to PLUMED) and with the PLUMED module of GROMACS 2025 and later.

```bash
./regtest.py --list                       # the tests and what they check
./regtest.py --gmx /path/to/gromacs/bin/gmx_mpi --kernel /path/to/plumed/lib/libplumedKernel.so
./regtest.py --gmx ... --kernel ... --tests 'pull|vsite' --out mytest
```

`--kernel` sets `PLUMED_KERNEL` (not needed when GROMACS is linked to PLUMED), `--out` is the directory of the runs
(the results are also written to `results.json` there) and `MPI_HOME` selects an MPI that is not in the `PATH`.
`python localrun.py tests/gromacs` (from the root of the repository) runs them with all the other tests, on the
`gromacs` script that `install.sh` writes, which must be in the `PATH` (without arguments it runs the general tests,
with arguments it runs `gmx_mpi` with them). Each test takes from a fraction of a second to a few seconds (small systems: alanine
dipeptide in vacuum, a TIP4P water box of 3636 atoms, a polarizable water dimer); all of them take one or two
minutes, unless a run hangs (it is killed after 2 minutes, or 5 for the longer run of the `committor` test).

To test a branch of PLUMED, build it, patch GROMACS with it in runtime mode, and run the tests with its kernel and
with the kernel of the branch it starts from: in runtime mode the same GROMACS works with both.

### GPU

The `gpu-buffer-ops` test needs GROMACS built with CUDA (`-DGMX_GPU=CUDA`) and a GPU, so it does not run in the
testcenter (GitHub's free runners have no GPU); without a CUDA build it is skipped. With CUDA 13 and GROMACS
2024 or older, add `-DCUDA_NVCC_FLAGS=-static-global-template-stub=false` to the cmake command. With a CUDA build,
all the other tests run on the CPU: sums on the GPU are not reproducible, and most tests compare two runs.

## Adding a test

A test is a function in `regtest.py` decorated with `@test("name")`, whose docstring says what it checks. It
prepares runs with `env.dir()`, `env.grompp()` (an input of `systems/` with mdp overrides) and `env.mdrun()`,
reads the results with `env.frames()` (trr), `env.energy()` (edr) and `colvar()`, and raises `Fail` (with the
numbers) or `Skip` (feature not available in this combination), or returns a message. Prefer properties that
hold everywhere over reference values.

## Not covered

- thread-MPI builds (the patches support only one rank with thread-MPI, which they check);
- GPU runs other than the force buffer test;
- performance, and simulations long enough to compare averages.
