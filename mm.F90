program main
    
    implicit none

    integer :: n,i,j,k
    real(8),dimension(3,3) :: Q,C,Z

    Q(1,1) = 0.316
    Q(1,2) =-0.95
    Q(1,3) = 0
    Q(2,1) = -0.95
    Q(2,2) = -0.32
    Q(2,3) = 0
    Q(3,1) = 0
    Q(3,2) = 0
    Q(3,3) = 1

    C(1,1) = -0.36
    C(1,2) = 0.79
    C(1,3) = -0.5
    C(2,1) = 0.65
    C(2,2) = -0.18
    C(2,3) = -0.74
    C(3,1) = 0.67
    C(3,2) = 0.59
    C(3,3) = 0.45

    do i = 1,3
        do j = 1,3
            Z(i,j) = 0
            do k = 1,3
                Z(i,j) = Z(i,j) + Q(i,k)*C(k,j)
            end do
        end do
    end do

    print*, ""

    do i = 1,size(Z,1)

      do j = 1,size(Z,2)
        write(*, fmt='(1X, A, F0.2)', advance="no") " ", z(i,j)
      end do

      print*, ""
    end do
end program 