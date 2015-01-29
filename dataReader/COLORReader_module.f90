module COLORReader_module
   use netcdf
   use NetcdfReader_module

   implicit none
!   include 'netcdf.inc'
   public :: COLORReader
   public :: newCOLORReader

   type,extends(NetcdfReader) ::  COLORReader
      integer :: NDIMS !number of dimension
      character (len = 20) :: VNAME 
      contains
         procedure :: readDATA
   end type
  
   interface newCOLORReader
      module procedure newCOLORReader
   end interface

   contains
      function newCOLORReader(ndims) &
              result(io)
         integer,intent(in) :: ndims !number of dimension
         type(COLORReader) :: io

         io%NDIMS     =ndims
         io%VNAME = "l3m_data"

      end function newCOLORReader

      function readDATA(this,file_name,xLon,yLat) result(d)
         class(COLORReader),intent(inout) :: this
         character(len=*),intent(in) :: file_name
         real,intent(in) :: xLon,yLat
         real :: d
         integer :: lonStart,latStart
         integer :: ncid
         integer,allocatable :: start(:), count(:) ,stride(:)! NDIMS
         real,allocatable :: lats(:), lons(:)
         integer :: varid
         ! Program variables to hold the data we will read in. We will only
         ! need enough space to hold one timestep of data; one record.
         real,allocatable :: d_in(:,:)! (NLONS,NLATS)
         ! Loop indices
         integer :: lat, lon,i,j,dimids(2)
         real :: longitude_step,latitude_step
         real :: slope,intercept
         
         character(Len=:),allocatable :: lon_name
         character(Len=:),allocatable :: lat_name

         lon_name = '        '
         lat_name = '        '
         ! Open the file. 
         call check( nf90_open(trim(file_name), nf90_nowrite, ncid) )
         call check( nf90_inq_varid(ncid, this%VNAME, varid) )
         call check(nf90_inquire_variable(ncid, varid, dimids=dimids(1:2)))

         call check(nf90_inquire_dimension(ncid, dimids(1), &
                &      name=lon_name,len=lon))
         call check(nf90_inquire_dimension(ncid, dimids(2), &
               &       name=lat_name,len=lat))
         allocate(lons(lon),lats(lat))
         allocate(d_in(2,2))
         call check( nf90_get_att(ncid, varid,"Slope", slope ))
         call check( nf90_get_att(ncid, varid,"Intercept", intercept))
         call check( nf90_get_att(ncid, nf90_global,"Latitude Step", Latitude_step))
         call check( nf90_get_att(ncid, nf90_global,"Longitude Step", Longitude_step))

         do i = 1, lon
            lons(i)= (i-1)*longitude_step
         enddo
         do i = 1, lat
            lats(i)= -90.0+(i-1)*latitude_step
         enddo

         lonStart = binarysearch(lon,lons, xLon)        
         latStart = binarysearch(lat,lats, yLat)        

         count = (/2,2/)
         start = (/lonStart,latStart /)
         stride =(/1,1/)

         call check( nf90_get_var(ncid, varid, d_in, start = start, &
                              count = count,stride=stride) )
         call check( nf90_close(ncid) )

         d = interpolate(lons(lonStart),lons(lonStart+1),lats(latStart),&
                 lats(latStart+1),d_in, xLon, yLat)

         d = slope*d+intercept

      contains
         subroutine check(status)
            integer, intent ( in) :: status
    
            if(status /= nf90_noerr) then 
                print *, trim(nf90_strerror(status))
                stop "Stopped"
            end if
         end subroutine check
   end function  
end module COLORReader_module
