module ode

    implicit none

    CONTAINS

    ! euler solver for a system of first order ODEs
    subroutine euler(f, par, t0, x0, dt, tf, xout, tout, saveat)
        
        implicit NONE

        interface

        ! we allow x to be a vector of length "dim" and the function f also returns a vector of derivatives
            function f(t,x,par)
                real(8), intent(in) :: t
                real(8), dimension(:), intent(in) :: x 
                real(8), dimension(:), intent(in) :: par ! additional parameters

                real(8), dimension(size(x)) :: f
            end function f


        end interface

        real(8), dimension(:), intent(in) :: par ! additional parameters passed to f function


        real(8), intent(in) :: t0, dt, tf
        real(8), dimension(:), intent(in) :: x0
        real(8), dimension(:), allocatable, intent(out) :: tout
        real(8), dimension(:,:), allocatable, intent(out) :: xout
        real(8), optional, intent(in) :: saveat

        real(8) :: saveat_ ! a copy of the optional argument
        integer :: save_every 

        integer :: i, j, dim, N
        real(8) :: t
        real(8), dimension(:), allocatable :: x

        if (present(saveat)) then
            saveat_ = saveat
        else
            saveat_ = dt
        end if

        save_every = nint(saveat_ / dt)
        
        N = ceiling((tf-t0)/saveat_) ! number of saves

        dim = size(x0)

        allocate(x(dim))

        allocate(xout(N+1, dim))
        allocate(tout(N+1))


        t = t0
        x = x0

        xout(1,:) = x0
        tout(1) = t0

        do i = 1, N

            do j = 1, save_every
                ! fortran does the update element-wise so this line is the same as for single-variable case.
                x = x + f(t,x,par) * dt
                t = t + dt
            end do
            ! print*, t, x
            xout(i+1,:) = x
            tout(i+1) = t
        end do

    end subroutine euler 

    ! midpoint solver for a system of first order ODEs
    subroutine midpoint(f, par, t0, x0, dt, tf, xout, tout, saveat)
        
        implicit NONE

        interface

            function f(t,x,par)
                real(8), intent(in) :: t
                real(8), dimension(:), intent(in) :: x 
                real(8), dimension(:), intent(in) :: par ! additional parameters
                real(8), dimension(size(x)) :: f
            end function f

        end interface

        real(8), dimension(:), intent(in) :: par ! additional parameters passed to f function

        
        real(8), intent(in) :: t0, dt, tf
        real(8), dimension(:), intent(in) :: x0
        real(8), dimension(:), allocatable, intent(out) :: tout
        real(8), dimension(:,:), allocatable, intent(out) :: xout
        real(8), optional, intent(in) :: saveat

        real(8) :: saveat_
        integer :: save_every 
        integer :: i, j, dim, N
        real(8) :: t, tm
        real(8), dimension(:), allocatable :: x, xm

        if (present(saveat)) then
            saveat_ = saveat
        else
            saveat_ = dt
        end if

        save_every = nint(saveat_ / dt)
        
        N = ceiling((tf-t0)/saveat_)

        dim = size(x0)

        allocate(x(dim))
        allocate(xm(dim))

        allocate(xout(N+1, dim))
        allocate(tout(N+1))


        t = t0
        x = x0

        xout(1,:) = x0
        tout(1) = t0

        do i = 1, N

            do j = 1, save_every
                ! first find t and x at the mid-point: tm and xm
                ! these are not used for updating t and x directly!!!
                tm = t + 0.5*dt
                xm = x+f(t,x,par)*0.5*dt
                ! then find the derivative at the mid point f(tm,xm) for the update
                x = x + f(tm,xm,par) * dt
                t = t + dt
            end do

            ! print*, t, x
            xout(i+1,:) = x
            tout(i+1) = t
        end do

    end subroutine midpoint 



    ! verlet solver for a conservative system where acc(x,par) is the acceleration function
    subroutine verlet(acc, par, t0, x0, v0, dt, tf, xout, vout, tout, saveat)
        
        implicit NONE

        interface

            function acc(x, par)
                real(8), dimension(:), intent(in) :: x 
                real(8), dimension(:), intent(in) :: par ! additional parameters
                real(8), dimension(size(x)) :: acc
            end function acc

        end interface

        real(8), dimension(:), intent(in) :: par ! additional parameters passed to acc function
        
        real(8), intent(in) :: t0, dt, tf
        real(8), dimension(:), intent(in) :: x0, v0
        real(8), dimension(:), allocatable, intent(out) :: tout
        real(8), dimension(:,:), allocatable, intent(out) :: xout, vout
        real(8), optional, intent(in) :: saveat


        real(8) :: saveat_
        integer :: save_every 
        integer :: N ! number of saves

        integer :: i, j, dim
        real(8) :: t
        real(8), dimension(:), allocatable :: x, v


        if (present(saveat)) then
            saveat_ = saveat
        else
            saveat_ = dt
        end if

        save_every = nint(saveat_ / dt)
        
        N = ceiling((tf-t0)/saveat_)

        dim = size(x0)

        allocate(x(dim))
        allocate(v(dim))

        allocate(xout(N+1, dim))
        allocate(vout(N+1, dim))
        allocate(tout(N+1))


        t = t0
        x = x0
        v = v0

        xout(1,:) = x0
        vout(1,:) = v0
        tout(1) = t0

        do i = 1, N

            do j = 1, save_every

                v = v + 0.5 * acc(x,par) * dt
                x = x + v * dt
                v = v + 0.5 * acc(x,par) * dt

                t = t + dt

            end do


            tout(i+1) = t
            vout(i+1,:) = v
            xout(i+1,:) = x

        end do

    end subroutine verlet 


    ! subroutine for writing the results in a file
    subroutine write_file(filename, tout, xout)
        character(len=*), intent(in) :: filename
        real(8), dimension(:), intent(in) :: tout
        real(8), dimension(:,:), intent(in) :: xout

        integer :: myunit, i, ios
        open(newunit = myunit, iostat=ios, file=filename,&
            status='replace', action='write')
            if(ios==0) then
                do i = 1, size(tout)

                    write(myunit,*) tout(i), xout(i,:)

                end do
                print*, 'written to ', filename


            else

               print*, 'error creating file.'
            end if
        close(myunit)
    end subroutine write_file

end module ode