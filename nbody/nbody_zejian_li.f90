program nbody
    use ode
    implicit NONE

    real(8) :: t0, dt, tf, saveat

    ! for verlet method
    integer :: n_body ! number of bodies
    real(8), dimension(:), allocatable :: t_out
    real(8), dimension(:), allocatable :: m_body
    real(8), dimension(:), allocatable :: x0, v0 ! vectors in 3d
    real(8), dimension(:,:), allocatable :: x_out, v_out

    
    integer :: u,i,ios
    open(unit=u, iostat=ios, file='nbody_in.txt',&
    status='old', action='read')

        if (ios==0) then
        
            read(u,*)
            read(u,*) n_body
            allocate(m_body(n_body), x0(n_body*3), v0(n_body*3))

            read(u,*)
            do i = 1, n_body
                read(u,*) m_body(i)
            end do

            read(u,*)
            do i = 1, n_body * 3
                read(u,*) x0(i)
            end do

            read(u,*)
            do i = 1, n_body * 3
                read(u,*) v0(i)
            end do

            read(u,*)
            read(u,*) t0
            
            read(u,*)
            read(u,*) tf

            read(u,*)
            read(u,*) dt

            read(u,*)
            read(u,*) saveat
            
        else
            print*, 'Input file not found. Using default parameters for 3body problem.'

                !!!!!!!!!!!!!!!!!!!!!!! input parameters
                t0 = 0.0
                dt = 0.01
                tf = 100.0
                saveat = 0.1


                n_body = 3 ! 3 body problem
                allocate(m_body(n_body), x0(n_body*3), v0(n_body*3))
                m_body(1) = 1. ! equal masses for all bodies
                m_body(2) = 1.
                m_body(3) = 1.

                ! initial conditions: http://www.jstor.org/stable/2661357
                ! A remarkable periodic solution of the three-body problem in the case of equal masses
                ! Pages 881-901 from Volume 152 (2000), Issue 3 by Alain Chenciner, Richard Montgomery

                x0(1:3) = (/ 0.97000436,  -0.24308753, 0.0 /)
                x0(4:6) = -x0(1:3)
                x0(7:9) = 0.0

                v0(7:9) = (/ -0.93240737, -0.86473146, 0.0 /) 
                v0(1:3) = -0.5 * v0(7:9)
                v0(4:6) = v0(1:3)

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end if

    close(u)


    call verlet(acc, m_body, t0, x0, v0, dt, tf, x_out,v_out, t_out, saveat)
    call write_file('nbody_verlet.txt', t_out, x_out)

    CONTAINS 

    ! norm function for length of a vector
    real(8) function norm(r)

        real(8), dimension(:), intent(in) :: r
        norm = sqrt(sum(r**2))

    end function norm

    ! acceleration vector for newtonian gravity
    function acc_pair(r)
        implicit none
        real(8), dimension(:), intent(in) :: r
        real(8), dimension(size(r)) :: acc_pair

        acc_pair = - r / ( norm(r)**3 )

    end function acc_pair

    function acc(x,par)
        implicit none
        real(8), dimension(:), intent(in) :: x
        real(8), dimension(:), intent(in) :: par ! masses
        real(8), dimension(size(x)) :: acc
        integer :: i, j
        integer :: n_body
        real(8), dimension(3) :: acc_ij, r_ij

        n_body = size(x) / 3

        acc = 0.0

        do i = 1, n_body

            do j = i+1, n_body

                r_ij = x(3*i-2:3*i) - x(3*j-2:3*j)
                acc_ij = acc_pair(r_ij) * par(i) * par(j)
                acc(3*i-2:3*i) = acc(3*i-2:3*i) + acc_ij
                acc(3*j-2:3*j) = acc(3*j-2:3*j) - acc_ij

            end do

        end do

    end function acc

end program nbody