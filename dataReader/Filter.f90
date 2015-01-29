program filter
   implicit none
   integer(kind=8) :: planf, planb
   real(kind=8),allocatable :: fftw_rin(:)   
   complex(kind=8),allocatable :: fftw_cout(:)
   integer :: n_days,istat
   real(kind=8),allocatable :: cdoms(:),kds(:), pars(:),mxds(:)
   real(kind=8) :: dmxds
   integer :: cutoff,i,n
   integer :: i_count,data_unitnum
   character(len=72) :: Aline
! read in data

   data_unitnum = 30
   open( unit = data_unitnum, file = 'data.csv', action="read")
   i_count = 0
   do
      read( data_unitnum, '(A)', iostat = istat ) Aline
      if( istat /= 0 ) exit
         i_count = i_count + 1
   enddo 
   n_days = i_count - 3
   
print*,"n_days", n_days
   
   allocate(cdoms(n_days))
   allocate(pars(n_days))
   allocate(kds(n_days))
   allocate(mxds(n_days))
   allocate(fftw_rin(n_days))
   allocate(fftw_cout(n_days/2+1))

   close(data_unitnum)

   open( unit = data_unitnum, file = 'data.csv', action="read")
   do i = 1,3
      read( data_unitnum, '(A)', iostat = istat ) Aline
   enddo
   do i = 1, n_days
      read( data_unitnum, '(A)', iostat = istat ) Aline
      do n = 1,72
        if (Aline(n:n) == ',') then
           Aline(n:n)= ' '
        endif
      enddo
      read( Aline,*)cdoms(i),kds(i),pars(i),mxds(i),dmxds
      write(Aline,*)' '
   enddo

   close(data_unitnum)

! prepare fftw

   call dfftw_plan_dft_r2c_1d(planf,n_days,fftw_rin, &
     &     fftw_cout,0)

   call dfftw_plan_dft_c2r_1d(planb,n_days,fftw_cout, &
     &     fftw_rin,0)

! transform to frequency space
! cuttoff numMonth = 12*(n_days/365) 
!  any frequency >= numMonth will be filtered out
   cutoff = 12*n_days/365
   print*,"cutoff",cutoff

! cdoms
   fftw_cout = 0.0d0
   call dfftw_execute_dft_r2c(planf,cdoms,fftw_cout)
   fftw_cout(cutoff+1:) = 0.0d0
   cdoms = 0.0d0
   call dfftw_execute_dft_c2r(planb,fftw_cout,cdoms)
! kds
   fftw_cout = 0.0d0
   call dfftw_execute_dft_r2c(planf,kds,fftw_cout)
   fftw_cout(cutoff+1:) = 0.0d0
   kds = 0.0d0
   call dfftw_execute_dft_c2r(planb,fftw_cout,kds)

! pars
   fftw_cout = 0.0d0
   call dfftw_execute_dft_r2c(planf,pars,fftw_cout)
   fftw_cout(cutoff+1:) = 0.0d0
   pars = 0.0d0
   call dfftw_execute_dft_c2r(planb,fftw_cout,pars)

! mxds
   fftw_cout = 0.0d0
   call dfftw_execute_dft_r2c(planf,mxds,fftw_cout)
   fftw_cout(cutoff+1:) = 0.0d0
   mxds = 0.0d0
   call dfftw_execute_dft_c2r(planb,fftw_cout,mxds)

! normalized

  cdoms = cdoms/n_days
  kds = kds/n_days
  pars= pars/n_days
  mxds = mxds/n_days

! write to files

  open( 1, file = 'low_passed_data.csv', &
            action='write', &
            status = 'unknown' )

   write(1,*) 'frequecies higher than month is filtered out'
   write(1,*) 'day_start:',1, ', day_end:',n_days
   write(1,*) 'cdoms,Kds,pars,mxls, dmxldt_day'
   do i = 1,n_days

      if (i/=n_days) then
         dmxds=mxds(i+1)-mxds(i)
      endif

      write(1,106)cdoms(i),Kds(i),pars(i),mxds(i),dmxds

   enddo
   close(1)
106    format(5(1x,e11.4))

end program
