# EnKF_sampling
Sampling of 1D and 2D pseudo random  curves and fields.

This code is used to test the sampling of pseudo random fields
on different architectures.  It currently supports different IBM,
SGI, LINUX DEC, and CRAY fft libraries.

The algorithm is the one derived in Evensen (1994, 2003, 2004) for
the generation of the pseudo random fields and the improved sampling.

Note that for IBM the ESSL library is used.

For LINUX we have decided to use the FFTW3 library which is 
available for many unix systems from http://www.fftw.org/

Before compiling define your achitecture in the file MODEL.CPP

Installation:
    install the scripts in EnKF_sampling/bin in your path

    cd EnKF_sampling/lib; make

To test (depends on EnKF_analysis routines): 

    git clone EnKF_analysis

    cd EnKF_analysis/lib

    change build in EnKF_analysis/lib/makefile to point to EnKF_sampling/build

    make
   
    cd EnKF_sampling/test; make

    cd build; testsampling
