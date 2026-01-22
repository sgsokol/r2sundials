News for r2sundials package
===========================

Version 7.2.1-4
---------------

* 2026-01-21
 - removed xref to archived ast2ast (signaled by CRAN team)

Version 7.2.1-3
---------------

* 2025-04-10
 - upgraded sundials/cvodes to 7.2.1
 - fixed the presence of sprintf, vprintf, vsprintf, stderr, stdout and abort signaled by CRAN team

Version 6.5.0-5
---------------

* 2023-12-08
 - fixed gcc-14 issue signaled by CRAN team
 - fixed Robertson test failure on some platforms

Version 6.5.0-4
---------------

* 2023-12-08
 - fixed warnings "-Wdeprecated-declarations" from RcppArmadillo (reported by CRAN team)
 - added .Rbuildignore
 - fixed html links

Version 6.5.0-3
---------------

* 2023-03-21
 - passed to v6.5.0 of SUNDIALS library
 - fixed warnings from "-Wstrict-prototypes" (reported by R-core)
 - updated CIATION format

Version 5.0.0-10
----------------

* 2021-05-17
 - fixed errors from "-fno-common" compile option relative to global variables (reported by R-core)

Version 5.0.0-9
---------------

* 2021-03-13
 - fixed erroneous memory freeing reported by valgrind check (run by R-core)
 - fixed https links

Version 5.0.0-7
---------------

* 2020-03-12
 - fixed too long line in examples

Version 5.0.0-6
---------------

* 2020-03-09
 - fixed clan-UBSAN error about function signatures mismatching
 - fixed clan-UBSAN error about memory leak
 - manual page for r2cvodes() is updated for C++ function signatures passed as callback arguments

Version 5.0.0-5
---------------

* 2020-03-06
 - fixed gcc-UBSAN error about misaligned memory access
 - fixed clan-ASAN error about stack-use-after-scope
 - fixed unit test failure on solaris

Version 5.0.0-4
---------------

* 2020-01-06
 - fixed Copyright holders in DESCRIPTION

Version 5.0.0-3
---------------

* 2019-12-18
 - fixed doc and DESCRIPTION

Version 5.0.0-2
---------------

* 2019-12-12
 - integrated CVODES-5.0.0 library into the source distribution

Version 5.0.0-1
---------------

* 2019-11-18
 - moved to CVODES v5.0.0 (but should be still compatible with v4.1.0)
 - renamed exported function from cvodes() to r2cvodes()
 - changed interaction with rmumps to pass by CCallable

Version 4.1.0-1
---------------

* 2019-09-03
 - first release. It includes cvodes() for ODE and optionally sensitivity problem solving
