!
!latitude =(-90+1/8):(1/4):90
!longitude =(1/8):(1/4):360
!
module MXDReader_module
   use netcdf
   use NetcdfReader_module

   implicit none
!   include 'netcdf.inc'

   public :: MXDReader
   public :: newMXDReader

   type,extends(NetcdfReader) :: MXDReader
      integer :: NDIMS !number of dimension
      character (len = 20) :: VNAME
   contains
      procedure :: readData
   end type
  
   interface newMXDReader
      module procedure newMXDReader
   end interface

   contains
      function newMXDReader(ndims) result(io)
         integer,intent(in) :: ndims !number of dimension
         type(MXDReader) :: io

         io%NDIMS     =ndims
         io%VNAME = "MXLDEPTH"

      end function newMXDReader

      function readDATA(this,file_name,xLon,yLat) result (d) 
         class(MXDReader),intent(inout) :: this
         character (len= *),intent(in) :: file_name
         real,intent(in) :: xLon,yLat
         real :: d 
         integer :: lonStart,latStart
         integer :: ncid
         ! The start and count arrays will tell the netCDF library where to
         ! Read our data.
         integer,allocatable :: start(:), count(:),stride(:) ! NDIMS

         real,allocatable :: lats(:), lons(:)
         integer :: lon_varid, lat_varid,time_varid,dimids(3)
         integer :: mxl_varid
         ! Program variables to hold the data we will read in. We will only
         ! need enough space to hold one timestep of data; one record.
         real,allocatable :: mxl_in(:,:)! (NLONS,NLATS)
         character(len=:),allocatable :: lon_name, lat_name, time_name,time_origin
         integer :: nlon,nlat,ntime

         lon_name = '             '
         lat_name = '             '
         time_name = '            '
         time_origin = '                                '

         ! Open the file. 
         call check( nf90_open(trim(file_name), nf90_nowrite, ncid) )
         call check( nf90_inq_varid(ncid, this%VNAME, mxl_varid) )
         call check(nf90_inquire_variable(ncid, mxl_varid, dimids=dimids(1:3)))
         call check(nf90_inquire_dimension(ncid, dimids(1), &
                &      name=lon_name,len=nlon))
         call check(nf90_inquire_dimension(ncid, dimids(2), &
               &       name=lat_name,len=nlat))
         call check(nf90_inquire_dimension(ncid, dimids(3), &
               &       name=time_name,len=ntime))

         allocate(start(this%NDIMS),count(this%NDIMS))
         allocate(lons(nlon),lats(nlat))
         allocate(mxl_in(2,2))

         call check( nf90_inq_varid(ncid, lat_name, lat_varid) )
         call check( nf90_inq_varid(ncid, lon_name, lon_varid) )
         call check( nf90_inq_varid(ncid, time_name, time_varid) )
         call check( nf90_get_att(ncid,time_varid,"time_origin",time_origin))

         call check( nf90_get_var(ncid, lat_varid, lats) )
         call check( nf90_get_var(ncid, lon_varid, lons) )

         ! Get the point needed

         lonStart = binarysearch(nlon,lons, xLon)        
         latStart = binarysearch(nlat,lats, yLat)        

         count = (/2,2,1/)
         start = (/lonStart,latStart,1 /)
         stride = (/1,1,1/)

         ! Read the mixed layer length data from the file
         call check( nf90_get_var(ncid, mxl_varid, mxl_in, start = start, &
                              count = count) )
         ! Close the file. This frees up any internal netCDF resources
         ! associated with the file.
         call check( nf90_close(ncid) )
         d = interpolate(lons(lonStart),lons(lonStart+1),lats(latStart),&
             lats(latStart+1),mxl_in, xLon, yLat)

      contains
         subroutine check(status)
            integer, intent ( in) :: status
    
            if(status /= nf90_noerr) then 
                print *, trim(nf90_strerror(status))
                stop "Stopped"
            end if
         end subroutine check
   end function
 
end module MXDReader_module
