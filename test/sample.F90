program sample
! Code for testing sampling of pseudo random fields on different architectures.
   use m_pseudo1D
   use m_pseudo2D
   use m_fixsample1D
   use m_fixsample2D
   use m_tecfld
   use m_set_random_seed2
   implicit none


   real :: xlength=10000.0                ! domain size
   real :: ylength=10000.0                ! domain size
   integer, parameter :: nrens=250        ! ensemble size 
   integer, parameter :: nx=100           ! gridsize in x direction
   integer, parameter :: ny=100           ! gridsize in y direction
   real :: cor1=1500.0                    ! decorrelation length in princepal direction
   real :: cor2=600.0                     ! decorrelation length normal to princepal direction
   real :: dir=30.0                       ! princepal direction

   integer n1                             ! x-dimension of grid
   integer n2                             ! y-dimension of grid

   integer i

   real dx,dy

   real A(nx,ny,nrens),B(nx,nrens)        ! 2d and 1D samples


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set a variable random seed
   call set_random_seed2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Correction for Gaussian variogram
   cor1=cor1/sqrt(3.0)
   cor2=cor2/sqrt(3.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   dx=xlength/float(nx-1)
   dy=ylength/float(ny-1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sample 2D pseudo random fields
   n1=nint(real(nx)*1.2)
   n2=nint(real(ny)*1.2)

   call pseudo2D(A,nx,ny,nrens,cor1,cor2,dx,dy,n1,n2,dir,.true.)
   call fixsample2D(A,nx,ny,nrens)

   call tecfld('ens2D',nx,ny,min(10,nrens),A)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sample all random fields and interpolate to model grid
   call pseudo1D(B,nx,nrens,cor1,dx,n1)
   call fixsample1D(B,nx,nrens)


   open(10,file='tec_ens1D.dat')
   do i=1,nx
      write(10,'(i3,20f12.2)')i,(i-1)*dx,B(i,1:min(nrens,10))
   enddo
   close(10)

end program
