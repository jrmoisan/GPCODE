subroutine spline ( x, y, n, yp1, ypn, y2 )



    use kinds_mod 

    implicit none

    integer(kind=i4b), intent(in) :: n
    real(kind=r8b), intent(in) :: x(n), y(n), yp1, ypn
    real(kind=r8b), intent(out) :: y2(n)
    integer(kind=i4b), parameter :: NMAX = 500
    integer(kind=i4b) :: i
    real(kind=r8b) :: p, qn, sig, un, u(NMAX)

    if ( yp1 .gt. 0.99D+30 ) then
        y2(1) = 0.D+0
        u(1) = 0.D+0
    else
        y2(1) = -0.5D+0
        u(1) = (3.D+0 / (x(2)-x(1))) * ((y(2)-y(1)) / (x(2)-x(1)) - yp1)
    endif

    do i = 2,n-1
        sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
        p = sig*y2(i-1) + 2.D+0
        y2(i) = (sig - 1.D+0) / p
        u(i) = ( 6.D+0 * ((y(i+1)-y(i))/(x(i+1)-x(i)) - &
                 (y(i)-y(i-1))/(x(i)-x(i-1))) / (x(i+1)-x(i-1)) - sig*u(i-1)) / p
    enddo

    if ( ypn .gt. 0.99D+30 ) then
        qn = 0.D+0
        un = 0.D+0
    else
        qn = 0.5D+0
        un = (3.D+0 / (x(n)-x(n-1))) * (ypn - (y(n)-y(n-1))/(x(n)-x(n-1)))
    endif

    y2(n) = (un - qn*u(n-1)) / (qn*y2(n-1) + 1.D+0)

    do i = n-1,1,-1
        y2(i) = y2(i) * y2(i+1) + u(i)
    enddo

    return
end subroutine spline

subroutine splint2 ( xa, ya, y2a, n, x, y, dydx )
    implicit none

    integer (kind = 4) :: n
    real (kind = 8), intent(in) :: x, xa(n), y2a(n), ya(n)
    real (kind = 8), intent(out) :: y, dydx
    integer (kind = 4) :: k, khi, klo
    real (kind = 8) :: a, b, h

    klo=1
    khi=n

    do while ( khi - klo .gt. 1 )
        k = (khi + klo) / 2
        if ( xa(k) .gt. x ) then
            khi = k
        else
            klo = k
        endif
    enddo

    h = xa(khi) - xa(klo)
    if ( h .eq. 0 ) write(*,*) 'bad xa input in splint'
    a = (xa(khi) - x) / h
    b = (x - xa(klo)) / h
    y = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi))*(h**2)/6
    dydx = (ya(khi) - ya(klo))/(xa(khi) - xa(klo)) - ((3.D+0*a**2 - 1.D+0) * &
        (xa(khi) - xa(klo))*y2a(klo) - (3.D+0*b**2-1.D+0)*(xa(khi) - xa(klo)) * &
        y2a(khi)) / 6.D+0

    return
end subroutine splint2
