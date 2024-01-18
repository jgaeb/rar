## Resubmission

This is a resubmission. In this version, I have:

1. Fixed the spelling of `arXiv` in the `DESCRIPTION` file.
2. Ensured that no more than two cores are used in examples, vignettes, etc.
  - The only exception is one test which checks that the proper warning message
    is issued when a user attempts to use more threads than they have available
    cores (`tests/testthat/test-sens.R:279`). However, I have modified this
    test so that it will not run on CRAN.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
