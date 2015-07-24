program main
!     program to calculate the autocorrelation functions for the SeaWiFS data
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
   use MXDReader_module
   use COLORReader_module

   implicit none
!
   character(len=150):: mxlfolder_in
   character(len=150):: cdomfolder_in
   character(len=150):: kdfolder_in
   character(len=150):: parfolder_in

   character(len=70) :: tmpfile_in
   character(len=70) :: parfile_in
   character(len=70) :: Kdfile_in
   character(len=70) :: mxlfile_in

   character(len=50) :: par_suffix
   character(len=50) :: Kd_suffix
   character(len=50) :: mxl_prefix

   integer  :: i,imonth,iday,iyear,julian_day,jday,julday
   integer  :: year_start,year_end,kth
   integer  :: day_start, day_end
   character (len=8) :: ch8
   character (len=4) :: ch4
   character (len=2) :: ch2_m,ch2_d      

   character(len=300)::  CMD
   character(len=250)::  filename
   integer :: ioend,NDIMS,days 
   real :: BermudaLat, BermudaLon
   type(MXDReader) :: mxlIO
   type(COLORReader) :: cdomIO,KdIO,parIO
   real :: cdom, par,mxl,Kd,filledValue,dmxldt_day
   real,allocatable :: cdoms(:),pars(:),mxls(:),Kds(:)

   interface 
      subroutine fillgap(f,startday,endday,filledValue)
!TODO : to make it simple, make sure the first and the end are not filled
      implicit none
      integer,intent(in) :: startday,endday
      real,intent(inout) :: f(:)
      real,intent(in) :: filledValue
      end subroutine
   end interface
!create a new directory listing

   year_start = 1998
  ! year_end   = 2010
   year_end   = 2009

   days = 365*(year_end-year_start+1)
   allocate(cdoms(1:days))
   allocate( pars(1:days))
   allocate( mxls(1:days))
   allocate(  Kds(1:days))
 
   filledValue = -32766.0

   cdoms = filledValue
   pars  = filledValue
   mxls  = filledValue
   Kds   = filledValue
  
   BermudaLat = 32+20/60.0
   BermudaLon = 360-(64+45/60.0)
   NDIMS = 3
   mxlIO = newMXDReader(NDIMS)
   NDIMS = 2
   cdomIO = newCOLORReader(NDIMS)
   KdIO = newCOLORReader(NDIMS)
   parIO = newCOLORReader(NDIMS)

   mxl_prefix = 'MXLDEPTH.1440x720.'
   Kd_suffix  = '.L3m_DAY_KD490_Kd_490_9km.bz2'
   par_suffix ='.L3m_DAY_PAR_par_9km.bz2' 

   cdomfolder_in = '/Volumes/HD-LBU3/Gpdata/ocean_color_data/'// &
      &   'SeaWiFS/Mapped/Daily/9km/cdom/'

   Kdfolder_in = '/Volumes/HD-LBU3/Gpdata/ocean_color_data/'// &
      &   'SeaWiFS/Mapped/Daily/9km/Kd/'
   
   parfolder_in = '/Volumes/HD-LBU3/Gpdata/ocean_color_data/'// &
      &   'SeaWiFS/Mapped/Daily/9km/par/'

   mxlfolder_in='/Volumes/HD-LBU3/Gpdata/ECCO2/ecco2.jpl.nasa.gov/'// &
      &   'data1/cube/cube92/lat_lon/quart_90S_90N/MXLDEPTH.nc/'

   kth = 0
   do iyear =  year_start, year_end

      write(ch4(1:4),'(i4)') iyear

      CMD='ls ' // trim(cdomfolder_in)//ch4//'/' //' > files.list'
      call system(CMD)
      open(20,file='files.list',status='old',form='formatted')
!
      do
         print*,"iyear,days", iyear,kth
         kth = kth+1
         read  (20,'(A)',IOSTAT=ioend) tmpfile_in
         if (ioend .NE. 0 ) exit
         read(tmpfile_in(1:8),'(A)')ch8

! get cdom file and read

         filename=trim(cdomfolder_in)//ch4//'/'//trim(tmpfile_in)
         CMD='bunzip2 -k -c '// trim(filename) // ' > tmp_file.nc'
         call system(CMD)
         cdom = cdomIO%readDATA('tmp_file.nc',BermudaLon,BermudaLat)

! get kd file and read

         kdfile_in = ch8//trim(kd_suffix)
         filename=trim(Kdfolder_in)//ch4//'/'//trim(Kdfile_in)
         CMD='bunzip2 -k -c '// trim(filename) // ' > tmp_file.nc'
         call system(CMD)
         Kd = KdIO%readDATA('tmp_file.nc',BermudaLon,BermudaLat)

! get par file and read

         parfile_in = ch8//trim(par_suffix)
         filename=trim(parfolder_in)//ch4//'/'//trim(parfile_in)
         CMD='bunzip2 -k -c '// trim(filename) // ' > tmp_file.nc'
         call system(CMD)
         par = parIO%readDATA('tmp_file.nc',BermudaLon,BermudaLat)

! get mxl file and read

         read(ch8(6:8),'(i3)') jday
         julian_day=julday(12,31,iyear-1)+jday
         kth = (iyear - year_start)*365+jday

         imonth=0
         iday=0
         DO i=1,11
            IF (julday(i,1,iyear) .le. julian_day .and.  julday(i+1,1,iyear) .gt. julian_day) THEN
               imonth=i
               iday=julian_day-julday(i,1,iyear)+1
            END IF
         END DO
         IF (julday(12,1,iyear) .le. julian_day) THEN
            imonth=12
            iday=julian_day-julday(12,1,iyear)+1
         END IF

! convert the iyear,imonth,iday values to character strings

         !WRITE(ch4(1:4),'(i4)') iyear

         IF (imonth .lt. 10) THEN
            ch2_m(1:1)='0'
            WRITE(ch2_m(2:2),fmt='(i1)') imonth
         ELSE IF (imonth .ge. 10) THEN
            WRITE(ch2_m(1:2),'(i2)') imonth
         ENDIF

         IF (iday .lt. 10) THEN
            ch2_d(1:1)='0'
            WRITE(ch2_d(2:2),fmt='(i1)') iday
         ELSE IF (iday .ge. 10) THEN
            WRITE(ch2_d(1:2),'(i2)') iday
         ENDIF
!
        filename =trim(mxlfolder_in)//trim(mxl_prefix)//ch4//ch2_m//ch2_d//'.nc'
        mxl = mxlIO%readDATA( trim(filename),BermudaLon,BermudaLat)

        cdoms(kth)=cdom
        Kds(kth)  =Kd
        Pars(kth) =par
        mxls(kth) =mxl 
      enddo ! go  to file
      close(20)
   enddo ! go to year

   do i = 1,days
      if ( cdoms(i) > 0 .and. &
         &    Kds(i) > 0 .and. &
         &   pars(i) > 0 .and. &
         &  mxls(i) > 0 ) then
         day_start = i
         exit
      endif
   enddo 

   do i = days,1,-1
      if ( cdoms(i) > 0 .and. &
         &    Kds(i) > 0 .and. &
         &   pars(i) > 0 .and. &
         &  mxls(i) > 0 ) then
         day_end = i
         exit
      endif
   enddo 

   filledValue = 0.0

  print*,day_start,day_end

   call fillgap(cdoms,day_start,day_end,filledValue)
   call fillgap(Kds,  day_start,day_end,filledValue)
   call fillgap(pars, day_start,day_end,filledValue)
   call fillgap(mxls, day_start,day_end,filledValue)

   open( 1, file = 'data.csv', &
            action='write', & 
            status = 'unknown' )

   write(1,*) 'year_start:',year_start, ',year_end:',year_end
   write(1,*) 'day_start:',day_start, ', day_end:',day_end
   write(1,*) 'cdoms,Kds,pars,mxls, dmxldt_day'
   do i = 1,days

      if (i/=days) then
         dmxldt_day=mxls(i+1)-mxls(i)
         !if(dmxldt_day <=0) dmxldt_day =0.0d0
      endif

      write(1,106)cdoms(i),Kds(i),pars(i),mxls(i),dmxldt_day

   enddo
   close(1)
106    format(5(1x,e11.4))    
end program

subroutine fillgap(f,startday,endday,filledValue)
!TODO : to make it simple, make sure the first and the end are not filled
   implicit none
   integer,intent(in) :: startday,endday
   real,intent(inout) :: f(:)
   real,intent(in) :: filledValue
   integer :: i,j,k,n,k1,k2,k0 

   k1 = 0
   k2 = 0
   k0 = 0

   do i = startday,endday
      if(f(i) <=filledValue ) then
         k0= i
         do j = k0+1,endday
            if(f(j) > filledValue ) then
               k2 = j
               exit
            endif
         enddo
         f(k0) = f(k1)+(f(k2)-f(k1))/(k2-k1)*(k0-k1)
      endif
      k1=i
   enddo
end subroutine fillgap

function julday(mm,id,iyyy) result(jul)
   implicit none
   integer ,intent(in) :: mm,id,iyyy
   integer  :: jul
   integer  :: igreg
   integer  :: iyr,jy,jm,ja

   igreg=(15+31*(10+12*1582))
   iyr=iyyy
   IF (iyr .eq. 0) THEN
      STOP 'There is no year zero.'
   END IF
   IF (iyr .lt. 0) iyr=iyr+1

   IF (mm .gt. 2) THEN
      jy=iyr
      jm=mm+1
   ELSE
      jy=iyr-1
      jm=mm+13
   END IF

   jul=int(365.25*jy)+int(30.6001*jm)+id+1720995

   IF (id+31*(mm+12*iyr) .ge. igreg) THEN
      ja=int(0.01*jy)
      jul=jul+2-ja+int(0.25*ja)
   END IF
end function julday
