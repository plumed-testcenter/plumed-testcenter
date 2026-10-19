Regression tests of {code} with PLUMED {version}
------------------------

These tests are specific to {code}: they check parts of its interface with PLUMED that the general tests do not cover.
Each test checks a property that must hold for any version of PLUMED and of {code}, such as a PLUMED bias that must act as
the same bias applied by {code}, or the same simulation that must give the same results in two ways of running it, so that
no reference results are stored.
A test is skipped when this version of {code} does not have the feature it tests.
The tests, and how to run them on your own computer, are in the
[directory of the code](https://github.com/plumed-testcenter/plumed-testcenter/tree/main/tests/{code}) in the repository of
the testcenter.

{results}

The text files of the runs (inputs, logs and output of PLUMED) are available in this [zip archive](regtests_{version}.zip).
