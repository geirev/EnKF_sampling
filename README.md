# EnKF_sampling
Sampling of 1D and 2D pseudo random  curves and fields.

The algorithm is the one derived in Evensen (1994, 2003, 2004) for
the generation of the pseudo random fields and the improved sampling.

However, the latest verison avoids the use of the newton solver by giving analytic 
expressions for r1=sqrt(3)/rx and ry=sqrt(3)/ry

We have decided to use the FFTW3 library which is 
available for many unix systems from http://www.fftw.org/

Installation:
install the scripts in EnKF_sampling/bin in your path
```
git clone git@github.com:geirev/EnKF_sampling.git
cd EnKF_sampling/lib
make
```

To test:
```
cd EnKF_sampling/test
make
../build/testsampling
```
