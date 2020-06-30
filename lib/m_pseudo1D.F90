module m_pseudo1D
contains
subroutine pseudo1D(A,nx,nrfields,rx,dx,n1,lmean)
   use m_random
   use mod_fftw3
   implicit none
   integer, intent(in) :: nx               ! Dimension of random vector
   integer, intent(in) :: nrfields         ! Number of random vectors
   real, intent(out)   :: A(nx,nrfields)   ! The random vectors
   real, intent(in)    :: rx               ! characteristic lengthscale
   real, intent(in)    :: dx               ! delta x  ( >1 )
   integer, intent(in) :: n1               ! Dimension of random vector
   logical, intent(inout), optional :: lmean  ! random mean 
! n1 should visely be selected as n1=nx*rx/dx unless periodic fields
! For periodic fields n1=nx. If rx=n*dx and n > 0.1*nx there is a risk for negative
! long range correlations.
! Setting lmean=.false. generates random fields with mean approx equal zero for each realization
! by setting the first fourier amplitude equal to zero.

   integer l,i
   real c
   real pi2,deltak,summ
   logical lm

   real fampl(2,0:n1/2)
   real phi(0:n1/2)

   real tt
   real, parameter :: pi=3.141592653589

   integer(kind=8) plan
   real y(n1)

   if (present(lmean)) then
      lm=lmean    ! optional override or confirm default
   else
      lm=.true.   ! default is true
   endif

   pi2=2.0*pi
   deltak=pi2/(real(n1)*dx)

   call dfftw_plan_dft_c2r_1d(plan,n1,fampl,y,FFTW_ESTIMATE)

! Ensuring variance of samples = 1.0
   summ=0.0
   do l=-n1/2+1,n1/2
      summ=summ+exp( -( (pi*rx*real(l)) / (real(n1)*dx) )**2 )
   enddo
   if (.not.lm) summ=summ-1.0   ! Subtract one if mean is not sampled.
   c=sqrt(1.0/(deltak*summ))

   do i=1,nrfields
      call random_number(phi)
      phi=pi2*phi

      do l=0,n1/2 
         tt=0.5*((pi*rx*real(l))/(n1*dx))**2
         fampl(1,l)=exp(-tt)*cos(phi(l))*sqrt(deltak)*c
         fampl(2,l)=exp(-tt)*sin(phi(l))*sqrt(deltak)*c
      enddo
      if (.not.lm) fampl(1,0)=0.0 
      fampl(2,0)=0.0 ! This is always zero

      call dfftw_execute_dft_c2r(plan, fampl, y)
      
      A(1:nx,i)=y(1:nx)

   enddo

   call dfftw_destroy_plan(plan)

end subroutine pseudo1D
end module m_pseudo1D
