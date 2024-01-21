## R CMD check results

0 errors | 0 warnings | 1 note

## LTO Optimization

* This version differs from 0.0.2 only in that the signature of
  `run_testthat_test` declared in `src/cpp11.cpp` has been fixed to match the
  signature declared in the `testthat` header. This should fix the LTO bug that
  occurred in that version.
* The note comes from the fact that the most recent version of this package was
  submitted only several days ago. However, in light of the deadline for fixing
  the issue (Feb. 02), I am trying to resolve the issue as quickly as I am
  able.
