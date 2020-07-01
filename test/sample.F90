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
   integer, parameter :: nrens=1000       ! ensemble size 
   integer, parameter :: nx=100           ! gridsize in x direction
   integer, parameter :: ny=100           ! gridsize in y direction
   real :: cor1=2000.0                     ! decorrelation length in princepal direction
   real :: cor2=600.0                     ! decorrelation length normal to princepal direction
   real :: dir=60.0                       ! princepal direction

   integer n1                             ! x-dimension of grid
   integer n2                             ! y-dimension of grid

   integer i


   real dx,dy

   real A(nx,ny,nrens)                    ! Samples in 2D case
   real B(nx,nrens)                       ! Samples in 1D case    
   real R(nx,nx)                          ! Error covariance matrix in 1D case


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
   n1=2*nint(real(nx)*cor1/(2.0*dx)) ; print '(2(a,i5))','nx=',nx,', n1=',n1
   n2=2*nint(real(ny)*cor2/(2.0*dy)) ; print '(2(a,i5))','ny=',ny,', n2=',n2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sample 2D pseudo random fields
   call pseudo2D(A,nx,ny,nrens,cor1,cor2,dx,dy,n1,n2,dir,.true.)
   call fixsample2D(A,nx,ny,nrens)
!A(1,1,1)=0.0
!dir=0.0

   call tecfld('ens2D',nx,ny,min(10,nrens),A)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sample all random fields and interpolate to model grid
   print '(a)','Calling pseudo1D'
   call pseudo1D(B,nx,nrens,cor1,dx,n1,.true.)

   print '(a)','Calling fixsample1D'
   call fixsample1D(B,nx,nrens)

   print '(a)','Computing covariance matrix'
   R=matmul(B,transpose(B))/real(nrens-1)  

   print '(a)','Dumping outputs'
   open(10,file='tec_cov1D.dat')
   do i=1,nx
      write(10,'(2000f13.5)')real(i-1)*dx,sum(B(i,:))/real(nrens),exp(-(real(i-1)*dx/cor1)**2),R(i,1:nx)
   enddo
   close(10)

   open(10,file='tec_ens1D.dat')
   do i=1,nx
      write(10,'(i3,2000f13.5)')i,(i-1)*dx,B(i,1:min(nrens,1000))
   enddo
   close(10)

end program
