real function decay(t,b1,b2) result(y)
    real, intent(in) :: t,b1,b2
    y = b1 * exp(-b2 * t)
end function decay

subroutine random_stduniform(u)
    implicit none
    real,intent(out) :: u
    real :: r
    call random_number(r)
    u = 1 - r
 end subroutine random_stduniform

 subroutine random_stdnormal(x)
    implicit none
    real,intent(out) :: x
    real,parameter :: pi= 4 * atan (1.0_8)
    real :: u1,u2
    call random_stduniform(u1)
    call random_stduniform(u2)
    x = sqrt(-2*log(u1))*cos(2*pi*u2)
 end subroutine random_stdnormal

program generate
    implicit none
    integer, parameter :: n = 2000
    real, dimension(n) :: x
    real, dimension(n) :: y
    real, dimension(n) :: noise
    real, external :: decay

    integer :: i, ios, u=20
    real :: b1 = 5. , b2 = 0.6

    do i = 1, n
        x(i) = 0.01 * i
        y(i) = decay(x(i),b1,b2)
    end do

    do i = 1, n
        call random_stdnormal(noise(i))
    end do
    y = y + 0.1*noise


    open(unit=u, iostat=ios, file='example_data.txt', status='replace', action='write')
    if (ios /= 0) then
      print*, 'error opening the file.'
      stop
    else

            
        write(u,'(a1,a35)') '#', 'example data for exponential decay fit'
        write(u,'(a1,a15,a16)') '#', 'x', 'y' 
        do i = 1, n
            write(u,'(2f16.4)') x(i), y(i)
        end do
        close(u)
    end if

end program generate