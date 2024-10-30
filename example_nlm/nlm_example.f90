! example nonlinear least-square fitting program
! for fitting an exponential decay
! with f(x,beta) = beta(1) * exp( -beta(2) * x )

module util ! helper functions

    implicit none
    contains
    ! determinant of 2x2 matrix
    real function det22(m) result(res)
  
        real, dimension(2,2), intent(in) :: m
  
        res =  m(1,1)*m(2,2) - m(1,2)*m(2,1) 
  
    end function det22

    ! inverse of 2x2 matrix
    function inv22(m) 

        real, dimension(2,2) :: inv22
        real, dimension(2,2), intent(in) :: m
        inv22(1,1) = m(2,2)
        inv22(1,2) = -m(1,2)
        inv22(2,1) = -m(2,1)
        inv22(2,2) = m(1,1)
        inv22 = inv22 / det22(m)

    end function inv22

    ! the model function
    real function f(x,beta)
        real, dimension(2), intent(in) :: beta
        real, intent(in) :: x
        real :: b1, b2
        b1 = beta(1)
        b2 = beta(2)

        f = b1 * exp(-b2 * x)
    end function f
    
    ! the gradient function (partial derivatives)
    function df(x,beta)

        real, dimension(2) :: df 
        real, intent(in) :: x
        real, dimension(2), intent(in) :: beta

        real :: b1, b2, df1, df2
        b1 = beta(1)
        b2 = beta(2)

        df1 = exp(-b2*x)
        df2 = -x * b1 * exp(-b2*x)
        df = (/ df1, df2/)


    end function df

end module util

program nlm_example

    use util
    implicit none
    integer, parameter :: nr = 2000
    integer ::  u, ios, i, maxiter, iter
    real, dimension(2) ::  dbeta, beta, beta0
    real, dimension(nr) :: x, y, yfit, dy
    real, dimension(nr,2) :: J  ! Jacobian
    real, dimension(2,2) :: var ! covariance matrix

    
    real :: epsilon, alpha ! tolerance and step length
    real :: sres, sigsq ! sum of squared residuals and sigma squared
    u = 10
    open(unit=u, iostat=ios, file='example_data.txt',&
            status='old', action='read')

        if (ios==0) then
            ! skip the two lines of headers
            read (u, *)
            read (u, *)
            do i = 1, nr

                read (u, *) x(i), y(i)
                
            end do

            beta0 = (/1.0,1.0/)

            maxiter = 1000

            epsilon = 1.0e-6

            alpha = 0.1

            beta = beta0

            do iter = 1, maxiter

                do i = 1, nr

                    yfit(i) = f(x(i),beta)

                    J(i,:) = df(x(i),beta)
    
                end do

                dy = y - yfit


                dbeta = matmul( inv22( matmul( transpose(J),J ) ) ,matmul(transpose(J),dy) )

                beta = beta + alpha * dbeta

                if ( sum(dbeta**2) <= epsilon ) exit

            end do

            if (iter < maxiter) then
                print*, 'converged after', iter, 'iterations.'
            else
                print*, 'warning: maxiter reached without convergence!'
            end if


            print*, 'fit for beta : ', beta

            sres = 0.
            do i= 1,nr
                sres = sres + ( y(i) - f(x(i),beta) )**2
            end do

            sigsq = sres / (nr - size(beta))

            var = sigsq * inv22( matmul(transpose(J),J) )

            print*, 'standard error of beta_1 : ', sqrt(var(1,1))
            print*, 'standard error of beta_2 : ', sqrt(var(2,2))

        else
            print '(a25)', 'error: file not opened.'
        end if

    close(u)

end program nlm_example