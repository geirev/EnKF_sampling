
program main
   use ISO_C_BINDING
   use mod_fftw3
   implicit none

   real, parameter :: pi=3.14159265359
   real, allocatable :: test_in(:)
   integer(8) :: plan_c2r, plan_r2c
   integer :: n, i
   real :: rd=5.0
   real :: dx=0.5


   n = 100
   allocate(test_in(0:n+1))

   test_in = 0.

   call dfftw_plan_dft_r2c_1d(plan_r2c, n, test_in, test_in, FFTW_MEASURE)
   call dfftw_plan_dft_c2r_1d(plan_c2r, n, test_in, test_in, FFTW_MEASURE)

   do i=0,n/2-1
!      test_in(i) = (1.0/sqrt(2.0*pi*rd**2)) * exp(-0.5*(real(i))**2/rd**2)
      test_in(i) = (dx/sqrt(pi*rd**2)) * exp(-(real(i)*dx)**2/rd**2)
      test_in(n-i) = test_in(i)
   enddo

   open(10,file='in.dat')
   do i=0,n-1
      write(10,'(3f13.5)')real(i),test_in(i)
   enddo
   close(10)

   call dfftw_execute_dft_r2c(plan_r2c, test_in, test_in)

   open(10,file='out.dat')
   do i=0,n-1,2
      write(10,'(3f13.5)')real(i),test_in(i)
   enddo
   close(10) 

   open(10,file='anaFFT.dat')
   print *,'2pi/n=',2.0*pi/real(n),0.07979
   do i=0,n-1,2
!      write(10,'(3f13.5)')real(i),exp(-2.0*(0.5*pi*rd*real(i)/real(n))**2) 
      write(10,'(3f13.5)')real(i),exp(     -( (0.5*pi*rd*real(i))/(real(n)*dx) )**2) 
   enddo
   close(10)

   call dfftw_execute_dft_c2r(plan_c2r, test_in, test_in)

   test_in=test_in/real(n)

   open(10,file='back.dat')
   do i=0,n-1
      write(10,'(2f13.5)')real(i),test_in(i)
   enddo
   close(10)


   call dfftw_destroy_plan(plan_r2c)
   call dfftw_destroy_plan(plan_c2r)

   deallocate(test_in)


end program
