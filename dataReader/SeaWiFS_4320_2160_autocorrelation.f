      program SeaWiFS_4320_2160_autocorrelation
!     program to calculate the autocorrelation functions for the SeaWiFS data
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      implicit none
!
!     set the parameters
      integer nx,ny,nloop,nlag
      real spval
      parameter(nx=4320,ny=2160,nloop=10000,nlag=60,
     &          spval=-9999.)
!
      integer dfacc_read,dfnt_int32
      parameter(dfacc_read=1,dfnt_int32=24)
!
      real chla_A(nx,ny),chla_B(nx,ny),bar(nx,ny),
     &     SSE(nx,ny),cnt(nx,ny),rho(nx,ny)
      integer(kind=2) schla(nx,ny)
!
      real amin,amax
      real value,tmin,tmax,aslope,atercept,offset
      real alon,alat,date,cff
      real*4 base,slope,intercept,data_minimum,data_maximum
      integer it,i,j,k,mo,imo,F_year,F_jday,iclim,iloop,ilag
      integer irv,igv,ibv,icnt,icff,iday,ndays,iyear,lday,lyear
      integer imin,imax
      integer bmax,bmin
!
      real station_latitude,station_longitude
      integer*2 period_start_year,period_start_day,period_end_year,
     &        period_end_day
      character*17 start_time
      character*17 end_time
      integer*2 start_year,start_day,end_year,end_day
      integer*4 start_millisec,end_millisec,orbit,start_orbit,end_orbit
      character*24 map_projection
      character*14 latitude_units
      character*13 longitude_units
      real*4 northernmost_latitude,southernmost_latitude,
     &       westernmost_longitude,easternmost_longitude,
     &       latitude_step,longitude_step,
     &       SW_point_latitude,SW_point_longitude
      integer*4 data_bins,number_of_lines,number_of_columns
!
      integer sfstart,sfselect,sfrdata,sfendacc,sfend,sffinfo,
     &        sffattr,sfrattr,sffrattr,sfgainfo
      integer attr_index,data_type,count
      integer dpgpal
      integer julday
      character*20 attr_name
!
      integer sd_id,sds_id,sds_index,ndatasets,ngl_atts,status
!
      integer start(2),edges(2),stride(2)
!
      logical echo,go_A,go_B
!
      character*27 file_nm
      character*23 out_file_nm
      character*60 tmp_file_nm
      character*27 cfilnam,ifilnm_A,ifilnm_B
      character*150 filnam
      character*26 product_name
      character*26 title
      character*7 sensor_name
      character*2 product_type
      character*9 replacement_flag
      character*7 software_name
      character*4 software_version
      character*17 processing_time
      character*43 processing_control
      character*371 input_parameters
      character*22 input_files
      character*152 L2_flag_names
      character*1 data_center,station_name
      character*1 mission
      character*1 sensor,sensor_char
      character*1 ch1
      character*2 ch2
      character*3 ch3
      character*4 ch4
      character*200 CMD
      character*5 chmin5,chmax5
      character*28 parameter_name
      character*5 measure
      character*8 units
      character*12 scaling
      character*55 scaling_equation
!
      echo=.true.   ! .true. ==> write out to screen the programs progress
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!     create a new directory listing
      CMD='rm -f files.list'
      write(*,*) CMD
      call system(CMD)
      CMD='ls /Volumes/san0/data/mirror/seawifs/l3m_day/version_05/'//
     &    ' > files.list'
      write(*,*) CMD
      call system(CMD)
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     calculate the bar
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      do i=1,nx
        do j=1,ny
          bar(i,j)=0.
          cnt(i,j)=0.
        enddo
      enddo
!
      filnam='files.list'
      open(20,file=filnam,status='old',form='formatted')
!
      do iloop=1,nloop
        if (iloop .eq. nloop) then
          write(*,*) 'ERROR: nloop value too low'
          stop
        endif
!
        read(20,1,end=100) ch1
 1      format(a1)
        if (ch1 .eq. 'S') then
          backspace(20)
          read(20,2,end=100) file_nm
 2        format(a27)
!
          write(*,*) iloop,iclim,file_nm
          filnam='/Volumes/san0/data/mirror/seawifs/'//
     &           'l3m_day/version_05/'//
     &           file_nm
!
          CMD='rm -f tmp_file.bz2'
          call system(CMD)
          CMD='rm -f tmp_file'
          call system(CMD)
!
          CMD='cp '//filnam//' tmp_file.bz2'
          call system(CMD)
!
          CMD='bunzip2 tmp_file.bz2'
          call system(CMD)
!
          tmp_file_nm='tmp_file'
!
          write(*,*) 'going to sfstart ',tmp_file_nm
!         open the file and initialize the sd interface.
          sd_id=sfstart(tmp_file_nm,dfacc_read)
          write(*,*) 'Given "sd_id" number = ',sd_id
          write(*,*) '"filnam" = ',tmp_file_nm
          if (sd_id .eq. -1) then
            write(*,*) 'sfstart failed',sd_id
            pause
          endif
!
          attr_index=sffattr(sd_id,'Base')
          status = sfgainfo(sd_id,attr_index,attr_name,
     &                      data_type,count)
          status = sfrattr(sd_id,attr_index,base)
!
          attr_index=sffattr(sd_id,'Slope')
          status = sfgainfo(sd_id,attr_index,attr_name,
     &                      data_type,count)
          status = sfrattr(sd_id,attr_index,slope)
!
          attr_index=sffattr(sd_id,'Intercept')
          status = sfgainfo(sd_id,attr_index,attr_name,
     &                      data_type,count)
          status = sfrattr(sd_id,attr_index,intercept)
          write(*,*) 'Base/Slope/Intercept: ',base,
     &                slope,intercept
!
          sds_id=sffinfo(sd_id,ndatasets,ngl_atts)
          write(*,*) sds_id,sd_id,ndatasets,ngl_atts
!         select the first data set.
          sds_index=0
          write(*,*) 'going to sfselect'
          sds_id=sfselect(sd_id,sds_index)
          write(*,*) 'sfselect "sds_id" =',sds_id
!
!         set elements of the array start to 0, 
!         elements of the array edges to sds dimensions,
!         and elements of the array stride to 1 
!         to read the entire data.
          start(1)=0
          start(2)=0
          edges(1)=nx
          edges(2)=ny
          stride(1)=1
          stride(2)=1
!
!         read entire data into bchla array. note that sfrdata is used
!         to read the numeric data. 
          write(*,*) 'going to sfrdata',sds_id
          status=sfrdata(sds_id,start,stride,edges,schla)
          write(*,*) 'sfrdata "status" = ',status
!
!         terminate access to the data set.
          status=sfendacc(sds_id)
!
!         terminate access to the sd interface and close the file.
          status=sfend(sd_id)
!
!         convert the data over from integer to real values
!
!         i)   convert from signed byte to integer [(-32768 <==> 32767) to (0<==>65536)]
!         ii)  flip the array to obtain correct image 
!         iii) swap data array halves to lay on the correct grid
!         iv)  move data to ichla array
          do i=1,nx
            do j=1,ny
              if (i .le. 2160) then   ! issue of swapping halves
                if (schla(i,j) .lt. 0) then  ! issue of signed integer to integer
                  chla_A(i+2160,ny-j+1)=float(schla(i,j))+65536.
                else
                  chla_A(i+2160,ny-j+1)=float(schla(i,j))
                endif
              elseif (i .gt. 2160) then
                if (schla(i,j) .lt. 0) then
                  chla_A(i-2160,ny-j+1)=float(schla(i,j))+65536.
                else
                  chla_A(i-2160,ny-j+1)=float(schla(i,j))
                endif
              endif
            enddo
          enddo
!
!         move the data into the data arrays
          do i=1,nx
            do j=1,ny
!             REM: 65535 is the land/cloud/ice mask.
              if (chla_A(i,j) .lt. 65535.) then ! set chla value [mg Chla m-3]
                chla_A(i,j)=base**((chla_A(i,j)*slope)+intercept)
              elseif (chla_A(i,j) .eq. 65535.) then ! set spval
                chla_A(i,j)=spval
              endif
            enddo
          enddo
!
!         move the data into the data arrays
          do i=1,nx
            do j=1,ny
!             REM: 255 is the land/cloud/ice mask.
              if (chla_A(i,j) .ne. spval) then ! set chla value [mg Chla m-3]
                bar(i,j)=bar(i,j)+chla_A(i,j)
                cnt(i,j)=cnt(i,j)+1.
              endif
            enddo
          enddo
!
        endif ! data file available
!
      enddo ! go to next image file
 100  close(20)
!
!     complete the calculation
      do i=1,nx
        do j=1,ny
          if (cnt(i,j) .gt. 0.) then
            bar(i,j)=bar(i,j)/cnt(i,j)
          endif
        enddo
      enddo
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     calculate the SSE
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      do iloop=1,nloop
        if (iloop .eq. nloop) then
          write(*,*) 'ERROR: nloop value too low'
          stop
        endif
!
        filnam='files.list'
        open(20,file=filnam,status='old',form='formatted')
!
        read(20,1,end=110) ch1
        if (ch1 .eq. 'S') then
          backspace(20)
          read(20,2,end=110) file_nm
!
          write(*,*) iloop,iclim,file_nm
          filnam='/Volumes/san0/data/mirror/seawifs/'//
     &           'l3m_day/version_05/'//
     &           file_nm
!
          CMD='rm -f tmp_file.bz2'
          call system(CMD)
          CMD='rm -f tmp_file'
          call system(CMD)
!
          CMD='cp '//filnam//' tmp_file.bz2'
          call system(CMD)
!
          CMD='bunzip2 tmp_file.bz2'
          call system(CMD)
!
          tmp_file_nm='tmp_file'
!
          write(*,*) 'going to sfstart ',tmp_file_nm
!         open the file and initialize the sd interface.
          sd_id=sfstart(tmp_file_nm,dfacc_read)
          write(*,*) 'Given "sd_id" number = ',sd_id
          write(*,*) '"filnam" = ',tmp_file_nm
          if (sd_id .eq. -1) then
            write(*,*) 'sfstart failed',sd_id
            pause
          endif
!
          attr_index=sffattr(sd_id,'Base')
          status = sfgainfo(sd_id,attr_index,attr_name,
     &                      data_type,count)
          status = sfrattr(sd_id,attr_index,base)
!
          attr_index=sffattr(sd_id,'Slope')
          status = sfgainfo(sd_id,attr_index,attr_name,
     &                      data_type,count)
          status = sfrattr(sd_id,attr_index,slope)
!
          attr_index=sffattr(sd_id,'Intercept')
          status = sfgainfo(sd_id,attr_index,attr_name,
     &                      data_type,count)
          status = sfrattr(sd_id,attr_index,intercept)
          write(*,*) 'Base/Slope/Intercept: ',base,
     &                slope,intercept
!
          sds_id=sffinfo(sd_id,ndatasets,ngl_atts)
          write(*,*) sds_id,sd_id,ndatasets,ngl_atts
!         select the first data set.
          sds_index=0
          write(*,*) 'going to sfselect'
          sds_id=sfselect(sd_id,sds_index)
          write(*,*) 'sfselect "sds_id" =',sds_id
!
!         set elements of the array start to 0, 
!         elements of the array edges to sds dimensions,
!         and elements of the array stride to 1 
!         to read the entire data. 
          start(1)=0
          start(2)=0
          edges(1)=nx
          edges(2)=ny
          stride(1)=1
          stride(2)=1
!
!         read entire data into bchla array. note that sfrdata is used
!         to read the numeric data. 
          write(*,*) 'going to sfrdata',sds_id
          status=sfrdata(sds_id,start,stride,edges,schla)
          write(*,*) 'sfrdata "status" = ',status
!
!         terminate access to the data set.
          status=sfendacc(sds_id)
!
!         terminate access to the sd interface and close the file.
          status=sfend(sd_id)
!
!         convert the data over from integer to real values
!
!         i)   convert from signed byte to integer [(-32768 <==> 32767) to (0<==>65536)]
!         ii)  flip the array to obtain correct image 
!         iii) swap data array halves to lay on the correct grid
!         iv)  move data to ichla array
          do i=1,nx
            do j=1,ny
              if (i .le. 2160) then   ! issue of swapping halves
                if (schla(i,j) .lt. 0) then  ! issue of signed integer to integer
                  chla_A(i+2160,ny-j+1)=float(schla(i,j))+65536.
                else
                  chla_A(i+2160,ny-j+1)=float(schla(i,j))
                endif
              elseif (i .gt. 2160) then
                if (schla(i,j) .lt. 0) then
                  chla_A(i-2160,ny-j+1)=float(schla(i,j))+65536.
                else
                  chla_A(i-2160,ny-j+1)=float(schla(i,j))
                endif
              endif
            enddo
          enddo
!
!         move the data into the data arrays
          do i=1,nx
            do j=1,ny
!             REM: 65535 is the land/cloud/ice mask.
              if (chla_A(i,j) .lt. 65535.) then ! set chla value [mg Chla m-3]
                chla_A(i,j)=base**((chla_A(i,j)*slope)+intercept)
              elseif (chla_A(i,j) .eq. 65535.) then ! set spval
                chla_A(i,j)=spval
              endif
            enddo
          enddo
!
!         move the data into the data arrays
          do i=1,nx
            do j=1,ny
!             REM: 255 is the land/cloud/ice mask.
              if (chla_A(i,j) .ne. spval) then ! set chla value [mg Chla m-3]
                SSE(i,j)=SSE(i,j)+((chla_A(i,j)-bar(i,j))**2)
              endif
            enddo
          enddo
!
        endif ! data file available
!
      enddo ! go to next image file
 110  close(20)
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     calclate the autocorrelation for a set number of lags
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      do ilag=0,nlag
!
        write(unit=ch2,fmt='(i2)') ilag
        filnam='/tmp/SeaWiFS_L3m_DAY_CHLO_9_autocorr_lag_'//ch2//'.bin'
        open(unit=30,file=filnam,form='unformatted',status='new')
!
!       reset the rho array to 0.
        do i=1,nx
          do j=1,ny
            rho(i,j)=0.
            cnt(i,j)=0.
          enddo
        enddo
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       loop through the entire data set
!       file name format: 'S2002265.L3m_DAY_CHLO_9.bz2'
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
        do iyear=2002,2008
!
!         how many days in this year?
          ndays=julday(1,1,iyear+1)-julday(1,1,iyear)
!
          do iday=1,ndays
!
            if (iday+ilag .le. ndays) then
              lday=iday+ilag
              lyear=iyear
            else
              lday=iday+ilag-ndays
              lyear=iyear+1
            endif
!
!           find the name of the first file
            write(unit=ch4,fmt='(i4)') iyear
            if (iday .le. 9) then
              write(unit=ch1,fmt='(i1)') iday
              ifilnm_A='S'//ch4//'00'//ch1//'.L3m_DAY_CHLO_9.bz2'
            elseif (iday .ge. 10 .and. iday .le. 99) then
              write(unit=ch2,fmt='(i2)') iday
              ifilnm_A='S'//ch4//'0'//ch2//'.L3m_DAY_CHLO_9.bz2'
            elseif (iday .ge. 100) then
              write(unit=ch3,fmt='(i3)') iday
              ifilnm_A='S'//ch4//ch3//'.L3m_DAY_CHLO_9.bz2'
            endif
!
!           find the name of the second file at some lag in time
            write(unit=ch4,fmt='(i4)') lyear
            if (lday .le. 9) then
              write(unit=ch1,fmt='(i1)') lday
              ifilnm_B='S'//ch4//'00'//ch1//'.L3m_DAY_CHLO_9.bz2'
            elseif (lday .ge. 10 .and. lday .le. 99) then
              write(unit=ch2,fmt='(i2)') lday
              ifilnm_B='S'//ch4//'0'//ch2//'.L3m_DAY_CHLO_9.bz2'
            elseif (lday .ge. 100) then
              write(unit=ch3,fmt='(i3)') lday
              ifilnm_B='S'//ch4//ch3//'.L3m_DAY_CHLO_9.bz2'
            endif
            write(*,*) ifilnm_A
            write(*,*) ifilnm_B
!
!           search the 'files.list' file for the presence of both files
!
            filnam='files.list'
            open(20,file=filnam,status='old',form='formatted')
            go_A=.false.
            do iloop=1,nloop
              read(20,2,end=120) file_nm
              if (file_nm .eq. ifilnm_A) go_A=.true.
            enddo
 120        rewind(20)
            close(20)
!
            filnam='files.list'
            open(20,file=filnam,status='old',form='formatted')
            go_B=.false.
            do iloop=1,nloop
              read(20,2,end=130) file_nm
              if (file_nm .eq. ifilnm_B) go_B=.true.
            enddo
 130        rewind(20)
            close(20)
!
            if (go_A .and. go_B) then
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!             get and read in the first data set
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
              write(*,*) iyear,iday,ifilnm_A
              filnam='/Volumes/san0/data/mirror/seawifs/'//
     &               'l3m_day/version_05/'//
     &               ifilnm_A
!
              CMD='rm -f tmp_file.bz2'
              call system(CMD)
              CMD='rm -f tmp_file'
              call system(CMD)
!
              CMD='cp '//filnam//' tmp_file.bz2'
              call system(CMD)
!
              CMD='bunzip2 tmp_file.bz2'
              call system(CMD)
!
              tmp_file_nm='tmp_file'
!
              write(*,*) 'going to sfstart ',tmp_file_nm
!             open the file and initialize the sd interface.
              sd_id=sfstart(tmp_file_nm,dfacc_read)
              write(*,*) 'Given "sd_id" number = ',sd_id
              write(*,*) '"filnam" = ',tmp_file_nm
              if (sd_id .eq. -1) then
                write(*,*) 'sfstart failed',sd_id
                pause
              endif
!
              attr_index=sffattr(sd_id,'Base')
              status = sfgainfo(sd_id,attr_index,attr_name,
     &                          data_type,count)
              write(*,*) count,attr_index
              status = sfrattr(sd_id,attr_index,base)
!
              attr_index=sffattr(sd_id,'Slope')
              status = sfgainfo(sd_id,attr_index,attr_name,
     &                          data_type,count)
              write(*,*) count,attr_index
              status = sfrattr(sd_id,attr_index,slope)
!
              attr_index=sffattr(sd_id,'Intercept')
              status = sfgainfo(sd_id,attr_index,attr_name,
     &                          data_type,count)
              write(*,*) count,attr_index
              status = sfrattr(sd_id,attr_index,intercept)
              write(*,*) 'Base/Slope/Intercept: ',base,
     &                    slope,intercept
!
              sds_id=sffinfo(sd_id,ndatasets,ngl_atts)
              write(*,*) sds_id,sd_id,ndatasets,ngl_atts
!             select the first data set.
              sds_index=0
              write(*,*) 'going to sfselect'
              sds_id=sfselect(sd_id,sds_index)
              write(*,*) 'sfselect "sds_id" =',sds_id
!
!             set elements of the array start to 0, 
!             elements of the array edges to sds dimensions,
!             and elements of the array stride to 1 
!             to read the entire data. 
              start(1)=0
              start(2)=0
              edges(1)=nx
              edges(2)=ny
              stride(1)=1
              stride(2)=1
!
!             read entire data into bchla array. note that sfrdata is used
!             to read the numeric data. 
              write(*,*) 'going to sfrdata',sds_id
              status=sfrdata(sds_id,start,stride,edges,schla)
              write(*,*) 'sfrdata "status" = ',status
!
!             terminate access to the data set.
              status=sfendacc(sds_id)
!
!             terminate access to the sd interface and close the file.
              status=sfend(sd_id)
!
!             convert the data over from integer to real values
!             i)   convert from signed byte to integer [(-32768 <==> 32767) to (0<==>65536)]
!             ii)  flip the array to obtain correct image 
!             iii) swap data array halves to lay on the correct grid
!             iv)  move data to ichla array
              do i=1,nx
                do j=1,ny
                  if (i .le. 2160) then   ! issue of swapping halves
                    if (schla(i,j) .lt. 0) then  ! issue of signed integer to integer
                      chla_A(i+2160,ny-j+1)=float(schla(i,j))+65536.
                    else
                      chla_A(i+2160,ny-j+1)=float(schla(i,j))
                    endif
                  elseif (i .gt. 2160) then
                    if (schla(i,j) .lt. 0) then
                      chla_A(i-2160,ny-j+1)=float(schla(i,j))+65536.
                    else
                      chla_A(i-2160,ny-j+1)=float(schla(i,j))
                    endif
                  endif
                enddo
              enddo
!
!             move the data into the data arrays
              do i=1,nx
                do j=1,ny
!                 REM: 65535 is the land/cloud/ice mask.
                  if (chla_A(i,j) .lt. 65535.) then ! set chla value [mg Chla m-3]
                    chla_A(i,j)=base**((chla_A(i,j)*slope)+intercept)
                  elseif (chla_A(i,j) .eq. 65535.) then ! set spval
                    chla_A(i,j)=spval
                  endif
                enddo
              enddo
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!             get and read in the second data set
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
              write(*,*) iyear,iday,ifilnm_B
              filnam='/Volumes/san0/data/mirror/seawifs/'//
     &               'l3m_day/version_05/'//
     &               ifilnm_B
!
              CMD='rm -f tmp_file.bz2'
              call system(CMD)
              CMD='rm -f tmp_file'
              call system(CMD)
!
              CMD='cp '//filnam//' tmp_file.bz2'
              call system(CMD)
!
              CMD='bunzip2 tmp_file.bz2'
              call system(CMD)
!
              tmp_file_nm='tmp_file'
!
              write(*,*) 'going to sfstart ',tmp_file_nm
!             open the file and initialize the sd interface.
              sd_id=sfstart(tmp_file_nm,dfacc_read)
              write(*,*) 'Given "sd_id" number = ',sd_id
              write(*,*) '"filnam" = ',tmp_file_nm
              if (sd_id .eq. -1) then
                write(*,*) 'sfstart failed',sd_id
                pause
              endif
!
              attr_index=sffattr(sd_id,'Base')
              status = sfgainfo(sd_id,attr_index,attr_name,
     &                          data_type,count)
              write(*,*) count,attr_index
              status = sfrattr(sd_id,attr_index,base)
!
              attr_index=sffattr(sd_id,'Slope')
              status = sfgainfo(sd_id,attr_index,attr_name,
     &                          data_type,count)
              write(*,*) count,attr_index
              status = sfrattr(sd_id,attr_index,slope)
!
              attr_index=sffattr(sd_id,'Intercept')
              status = sfgainfo(sd_id,attr_index,attr_name,
     &                          data_type,count)
              write(*,*) count,attr_index
              status = sfrattr(sd_id,attr_index,intercept)
              write(*,*) 'Base/Slope/Intercept: ',base,
     &                    slope,intercept
!
              sds_id=sffinfo(sd_id,ndatasets,ngl_atts)
              write(*,*) sds_id,sd_id,ndatasets,ngl_atts
!             select the first data set.
              sds_index=0
              write(*,*) 'going to sfselect'
              sds_id=sfselect(sd_id,sds_index)
              write(*,*) 'sfselect "sds_id" =',sds_id
!
!             set elements of the array start to 0, 
!             elements of the array edges to sds dimensions,
!             and elements of the array stride to 1 
!             to read the entire data. 
              start(1)=0
              start(2)=0
              edges(1)=nx
              edges(2)=ny
              stride(1)=1
              stride(2)=1
!
!             read entire data into bchla array. note that sfrdata is used
!             to read the numeric data. 
               write(*,*) 'going to sfrdata',sds_id
              status=sfrdata(sds_id,start,stride,edges,schla)
              write(*,*) 'sfrdata "status" = ',status
!
!             find the min/max values
              imin=10000.
              imax=-10000.
              do i=1,nx
                do j=1,ny
                  if (imin .gt. schla(i,j)) imin=schla(i,j)
                  if (imax .lt. schla(i,j)) imax=schla(i,j)
                enddo
              enddo
              write(*,*) 'Region Chla min/max = ',imin,imax
!
!             terminate access to the data set.
              status=sfendacc(sds_id)
!
!             terminate access to the sd interface and close the file.
              status=sfend(sd_id)
!
!             convert the data over from integer to real values
!             i)   convert from signed byte to integer [(-32768 <==> 32767) to (0<==>65536)]
!             ii)  flip the array to obtain correct image 
!             iii) swap data array halves to lay on the correct grid
!             iv)  move data to ichla array
              do i=1,nx
                do j=1,ny
                  if (i .le. 2160) then   ! issue of swapping halves
                    if (schla(i,j) .lt. 0) then  ! issue of signed integer to integer
                      chla_B(i+2160,ny-j+1)=float(schla(i,j))+65536.
                    else
                      chla_B(i+2160,ny-j+1)=float(schla(i,j))
                    endif
                  elseif (i .gt. 2160) then
                    if (schla(i,j) .lt. 0) then
                      chla_B(i-2160,ny-j+1)=float(schla(i,j))+65536.
                    else
                      chla_B(i-2160,ny-j+1)=float(schla(i,j))
                    endif
                  endif
                enddo
              enddo
!
!             move the data into the data arrays
              do i=1,nx
                do j=1,ny
!                 REM: 65535 is the land/cloud/ice mask.
                  if (chla_B(i,j) .lt. 65535.) then ! set chla value [mg Chla m-3]
                    chla_B(i,j)=base**((chla_B(i,j)*slope)+intercept)
                  elseif (chla_B(i,j) .eq. 65535.) then ! set spval
                    chla_B(i,j)=spval
                  endif
                enddo
              enddo
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
            endif ! set of data files available
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!           carry out the correlation for this time link
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
              do i=1,nx
                do j=1,ny
                  if (chla_A(i,j) .ne. spval .and.
     &                chla_B(i,j) .ne. spval) then
                    rho(i,j)=rho(i,j)+
     &                ((chla_A(i,j)-bar(i,j))*(chla_B(i,j)-bar(i,j)))
                    cnt(i,j)=cnt(i,j)+1.
                  endif
                enddo
              enddo
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
          enddo ! go to the next yearday
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
        enddo ! go to the next year
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       complete the calculation for the correlation at the specific ilag
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
        do i=1,nx
          do j=1,ny
            if (cnt(i,j) .gt. 50) then
               rho(i,j)=rho(i,j)/SSE(i,j)
            else
              rho(i,j)=spval
            endif
          enddo
        enddo
!
        write(30) rho,cnt
        close(30)
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      enddo  ! go to the next ilag
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     stop end
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      stop
      end
!234567890123456789012345678901234567890123456789012345678901234567890
      function julday(mm,id,iyyy)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      parameter (igreg=15+31*(10+12*1582))
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      iyr=iyyy
      if (iyr .eq. 0) pause 'there is no year zero.'
      if (iyr .lt. 0) iyr=iyr+1
!
      if (mm .gt. 2) then
        jy=iyr
        jm=mm+1
      else
        jy=iyr-1
        jm=mm+13
      endif
!
      julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
!
      if (id+31*(mm+12*iyr).ge.igreg) then
        ja=int(0.01*jy)
        julday=julday+2-ja+int(0.25*ja)
      endif
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      return
      end
