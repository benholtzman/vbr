## 0_tests

This directory contains test scripts that should be run periodically when making changes.

`run_tests.m` is the main driver. The default is to run all `.m` files in this directory.

When creating new code functionality, add a test function here.

### running tests

To run a test:

```
cd vbr/0_tests
test_results=run_tests;
```

The code will execute all the `.m` files. If any fail to run, a warning message will print to screen. After the tests run, if there are any failed tests you should run those tests individually to see the full matlab error. The output variable `test_results` is a structure with each test function as a field with each set to 0 or 1 to indicate failure or success. The `run_test` function does include the option for specifying a `test_type`, e.g., `run_test(test_type)` but the only current type is `full_test`. You could add a new type that is a sub-set of the test functions if you want (e.g., all `.m` files containing some keyword).

Example output containing a failed test:

```
Running test_type: full_test

initializing with vbr_version: VBR_v0p95
VBR calculator initialized
Starting full_test

    **** Running test_000_vbrcore ****
    test_000_vbrcore failed :(
    please run test_000_vbrcore() and debug.

Testing complete.

Displaying failed test functions. Please run each one and debug:
    test_000_vbrcore
```

You would then need to run `test_000_vbrcore()` to debug the problem.
