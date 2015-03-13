#SIMPSON

SIMPSON is a widely used simulation program for solid-state NMR spectroscopy.

The original homepage with releases and examples is at http://nmr.au.dk/software/simpson

This site is intended for development and pre-release work.

Required libraries:

* CBLAS and LAPACK – can be replaced by any equivalent library (ATLAS, GotoBLAS, Intel MKL when flag –DINTEL_MKL is used,...)
* Tcl up to version 8.5 (with threads enabled like for the ActiveStateTcl?, on linux distros threads are typically disabled)
* pthread
* FFTW3
* NFFT3 (can be avoided when flag –DNO_NFFT is used but then no FWT interpolation will be possible)
* MPI when flag –DMPI is used (parallel version)
