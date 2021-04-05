program main
    use cudafor
    use utils
    use dsytrd_gpu
    implicit none

    integer                                         :: N, lwork, nb
    real(8), dimension(:,:), allocatable, device    :: A
    real(8), dimension(:), allocatable, device      :: work, d, e, tau
    real(8), dimension(:,:), allocatable            :: A_h
    real(8), dimension(:), allocatable              :: work_h, d_h, e_h, tau_h

    N = 5
    lwork = 2*N
    nb = N

    allocate(A_h(N,N))
    ! A_h(1, 1) = 1
    ! A_h(1, 2) = 2
    ! A_h(1, 3) = 3
    ! A_h(2, 1) = 2
    ! A_h(2, 2) = 4
    ! A_h(2, 3) = 1
    ! A_h(3, 1) = 3
    ! A_h(3, 2) = 1
    ! A_h(3, 3) = 1

    ! A_h(1, 1) = 1
    ! A_h(1, 2) = 2
    ! A_h(1, 3) = 3
    ! A_h(1, 4) = 4
    ! A_h(2, 1) = 2
    ! A_h(2, 2) = 4
    ! A_h(2, 3) = 1
    ! A_h(2, 4) = 1
    ! A_h(3, 1) = 3
    ! A_h(3, 2) = 1
    ! A_h(3, 3) = 1
    ! A_h(3, 4) = 2
    ! A_h(4, 1) = 4
    ! A_h(4, 2) = 1
    ! A_h(4, 3) = 2
    ! A_h(4, 4) = 2

    A_h(1, 1) = 1
    A_h(1, 2) = 2
    A_h(1, 3) = 3
    A_h(1, 4) = 4
    A_h(1, 5) = 5
    A_h(2, 1) = 2
    A_h(2, 2) = 4
    A_h(2, 3) = 1
    A_h(2, 4) = 1
    A_h(2, 5) = 5
    A_h(3, 1) = 3
    A_h(3, 2) = 1
    A_h(3, 3) = 1
    A_h(3, 4) = 2
    A_h(3, 5) = 3
    A_h(4, 1) = 4
    A_h(4, 2) = 1
    A_h(4, 3) = 2
    A_h(4, 4) = 2
    A_h(4, 5) = 3
    A_h(5, 1) = 5
    A_h(5, 2) = 5
    A_h(5, 3) = 3
    A_h(5, 4) = 3
    A_h(5, 1) = 1

    allocate(A, source=A_h)
    allocate(work_h(lwork))
    allocate(work, source=work_h)
    allocate(d_h(N))
    allocate(d, source=d_h)
    allocate(e_h(N-1))
    allocate(e, source=e_h)
    allocate(tau_h(N-1))
    allocate(tau, source=tau_h)

    call dsytrd_gpu('U', N, A, N, d, e, tau, work, lwork, nb)

end program