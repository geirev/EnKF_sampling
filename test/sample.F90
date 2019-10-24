program sample
! Samples the measurement perturbation matrix E for the ert-Eclipse applications.
! The program reads the ert/observations/obshist.txt file and the schedulefile defined below.
! Only configured for opr, wpr, gpr observations.
! Predefined well names.

   use m_pseudo1D
   use m_fixsample1D
   use m_tecfld
   use m_set_random_seed2
   implicit none

   integer, parameter :: nrens=100        ! ensemble size 
! read observations info assumed stored as follows:
   integer, parameter :: nx=36            ! Number of dates keywords
   integer, parameter :: nrwells=5        ! number of wells
   integer, parameter :: nrdata=3         ! number of variables per well
! We will then sample nrwells*nrdata time series of length nx and store these as one realization
! A total of nrens such realizations are computed.

   real :: xlength=35.0                   ! length of time series is 36 months
   real :: cor1=15.0                      ! decorrelation length in months
   character(len=100) :: obshistfile='obshistnew.txt' 
   character(len=100) :: schedulefile='history.sch' 
   integer n1                             ! x-dimension of grid
   integer i,j,k,ic
   integer idat
   integer iwell
   integer ipos(4)
   real dx

   real B(nx,nrwells*nrdata*nrens)                       !  1D samples
   real A(nx*nrwells*nrdata,nrens)
   real C(nx,nrwells*nrdata,nrens)
   real obs(nx,nrwells,nrdata)
   real obs1(nx,nrwells*nrdata)
   real stddev(nx,nrwells*nrdata)

   type wellrates
      character(len=12) date
      character(len=4)  wellname(nrwells)
      real opr(nrwells)
      real gpr(nrwells)
      real wpr(nrwells)
   end type
   type (wellrates) rates(nx)

   character(len=100) cline
   character(len=4) var
   character(len=4) well
   real rel_err
   real min_err

   integer status
   character(len=200) string1  
   character(len=200) string2  
   string1="grep -A1 -e 'DATES' -e 'OPEN RESV' "//trim(schedulefile)// " > rates.txt"
   string2="sed -i -e '/^\/$/d' -e '/--/d' -e 's?/??' -e 's/ OPEN RESV //' -e 's/ /#/g'  rates.txt"
   
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set a variable random seed
   call set_random_seed2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Correction for Gaussian variogram
   cor1=cor1/sqrt(3.0)

   dx=xlength/float(nx-1)
   n1=nint(real(nx)*1.5)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sample all random time series
   call pseudo1D(B,nx,nrwells*nrdata*nrens,cor1,dx,n1)
   call fixsample1D(B,nx,nrwells*nrdata*nrens)

! testplot
   A=reshape(B,(/ nx*nrwells*nrdata, nrens /) )
   open(10,file='tec_ens1D.dat')
   do i=1,nx*nrwells*nrdata
      write(10,'(i3,f12.2,100f12.2)')i,(i-1)*dx,A(i,1:min(nrens,100))
   enddo
   close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read rates from schedule file
   print '(a,a)','string1=',trim(string1)
   status = system( string1 )
   print '(a,a)','string2=',trim(string2)
   status = system( string2 )
   status = system( 'head rates.txt' )

   do i=1,nx
      do j=1,nrwells
         write(rates(i)%wellname(j),'(a,i1)')'OP_',j
         rates(i)%opr(j)=0.0
         rates(i)%wpr(j)=0.0
         rates(i)%gpr(j)=0.0
         rates(i)%date='xxxxxxxxxxxx'
      enddo
   enddo

   open (10,file='rates.txt')
   idat=1
   do i=1,10000
      read(10,'(a)',end=200)cline
      if (cline(1:1) == 'D') then
         read(10,'(a)',end=200)cline
         rates(idat)%date=cline(3:15)
         idat=idat+1
      else
         read(cline(6:6),*)iwell
         ic=0
         do j=3,len_trim(cline)
            if(cline(j:j)=='#') then
               ic=ic+1
               ipos(ic)=j
            endif
         enddo
         read(cline(ipos(1)+1:ipos(2)-1),*)rates(idat)%opr(iwell)
         read(cline(ipos(2)+1:ipos(3)-1),*)rates(idat)%wpr(iwell)
         read(cline(ipos(3)+1:ipos(4)-1),*)rates(idat)%gpr(iwell)
      endif
   enddo
   200 close(10)
   open(10,file='output.dat')
      do i=1,nx
         do j=1,nrwells
            write(10,'(a12,tr1,a4,tr2,2f12.3,e15.3)')&
              rates(i)%date,rates(i)%wellname(j),rates(i)%opr(j),rates(i)%wpr(j),rates(i)%gpr(j)
         enddo
      enddo
   close(10)

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! scale with standard deviations



   C=reshape(B,(/ nx,nrwells*nrdata, nrens /) )
   do i=1,nx
      do j=1,nrwells
         obs(i,j,1) = rates(i)%opr(j)
         obs(i,j,2) = rates(i)%wpr(j)
         obs(i,j,3) = rates(i)%gpr(j)
      enddo
   enddo
   obs1=reshape(obs,(/nx,nrwells*nrdata /) )


   open(10,file=trim(obshistfile),err=100)
   do i=1,nrwells*nrdata
      read(10,'(tr21,a3,tr1,a4,tr11,f3.1,tr36,f6.1)',end=100) var,well,rel_err,min_err
      write(*,'(a,a3,a,a4,a,f3.1,a,f6.1)') ' var=',var,' well=',well,' rel_err=',rel_err,' min_err=',min_err
      do k=1,nx
         if (obs1(k,i)/=0.0) then
            stddev(k,i)=max(abs(rel_err*obs1(k,i)),min_err)
         else
            stddev(k,i)=0.0
         endif
      enddo
   enddo
   100 close(10)

   do j=1,nrens
      do i=1,nrwells*nrdata
         do k=1,nx
            C(k,i,j)=stddev(k,i)*C(k,i,j)
         enddo
      enddo
   enddo

! testplot
   A=reshape(C,(/ nx*nrwells*nrdata, nrens /) )
   open(10,file='tec_ens1D_B.dat')
   do i=1,nx*nrwells*nrdata
      write(10,'(i3,f12.2,100f12.2)')i,(i-1)*dx,A(i,1:min(nrens,100))
   enddo
   close(10)


   open(10,file='E.dat')
   do j=1,nrens
   do i=1,nx*nrwells*nrdata
      write(10,'(i6,i4,e15.6)')i,j,A(i,j)
   enddo
   enddo
   close(10)




end program
