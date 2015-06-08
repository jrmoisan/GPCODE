module NetcdfReader_module
   use netcdf
   implicit none
!   include 'netcdf.inc'
   public :: NetcdfReader

   character (len = 5), parameter :: UNITS = "units"
   character (len = 1), parameter :: MXL_UNITS = "m"
   character (len = 13), parameter :: LAT_UNITS = "degrees_north"
   character (len = 12), parameter :: LON_UNITS = "degrees_east"

   type,abstract ::  NetcdfReader
      contains
         procedure(readD),deferred :: readDATA
   end type
  
   abstract interface
      function readD(this,file_name,xlon,ylat) result (d) 
         import NetcdfReader
         class(NetcdfReader), intent(inout) :: this
         character(len=*),intent(in) :: file_name
         real,intent(in) :: xlon
         real,intent(in) :: ylat
         real :: d
      end function
   end interface

contains

   function binarysearch(length, array, value, delta)
   ! Given an array and a value, returns the index of the element that
   ! is closest to, but less than, the given value.
   ! Uses a binary search algorithm.
   ! "delta" is the tolerance used to determine if two values are equal
   ! if ( abs(x1 - x2) <= delta) then
   !    assume x1 = x2
   ! endif
      implicit none
      integer, intent(in) :: length
      real, dimension(length), intent(in) :: array
      real, intent(in) :: value
      real, intent(in), optional :: delta

      integer :: binarysearch

      integer :: left, middle, right
      real :: d

      if (present(delta)) then
         d = delta
      else
         d = 1e-9
      endif
      left = 1
      right = length
      do
         if (left > right) then
            exit
         endif
         middle = nint((left+right) / 2.0)
         if ( abs(array(middle) - value) <= d) then
            binarySearch = middle
            return
         else if (array(middle) > value) then
                 right = middle - 1
              else
                 left = middle + 1
         end if
      end do
      binarysearch = right

   end function binarysearch

   real function interpolate(x1, x2, y1, y2, f, x, y)
   ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation
      implicit none
      real, dimension(2, 2), intent(in) :: f
      real, intent(in) :: x1,x2,y1,y2,x,y

      real :: denom
      integer :: i,j

      denom = (x2 - x1)*(y2 - y1)

      i=1
      j=1
      interpolate = (f(i,j)*(x2-x)*(y2-y) + f(i+1,j)*(x-x1)*(y2-y) + &
      f(i,j+1)*(x2-x)*(y-y1) + f(i+1, j+1)*(x-x1)*(y-y1))/denom

   end function interpolate

end module NetcdfReader_module
