module m_set_random_seed2
contains
subroutine set_random_seed2
! Sets a random seed based on the system and wall clock time
   implicit none 

   logical ex
   real x
   INTEGER :: i, n, clock
   INTEGER, DIMENSION(:), ALLOCATABLE :: seed

   CALL RANDOM_SEED(size = n)
   ALLOCATE(seed(n))

   inquire(file='seed.dat',exist=ex)
   if (ex) then
      open(10,file='seed.dat')
         read(10,*)seed
      close(10)
      CALL RANDOM_SEED(PUT = seed)
      call random_number(x)
      print '(a,f20.18)','Random seed set using stored seed from seed.dat:',x

   else
      CALL SYSTEM_CLOCK(COUNT=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)
      call random_number(x)
      print '(a,f20.18)','New random seed set using system clock and stored to seed.dat:',x
      open(10,file='seed.dat')
         write(10,*)seed
      close(10)
   endif

   DEALLOCATE(seed)

end subroutine set_random_seed2
end module m_set_random_seed2
