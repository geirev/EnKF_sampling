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

You are probably using the mkdepend.pl script found elsewhere on 
this site. 

This code is developed by Geir.Evensen@hydro.com and/or Geir.Evensen@nersc.no
The code can be freely used but improvements, upgrades and bugfixes
should be reported,  and as everything free, it comes with no guarantee.  

