program main
    use cudafor
    use utils
    use dsyevd_gpu
    implicit none

    integer                                         :: N, lwork, liwork_h, info, il, iu
    real(8), dimension(:,:), allocatable, device    :: A,Z
    real(8), dimension(:), allocatable, device      :: work, w
    real(8), dimension(:,:), allocatable            :: A_h,Z_h
    real(8), dimension(:), allocatable              :: work_h, w_h
    integer, dimension(:), allocatable              :: iwork_h

    N = 3
    lwork = 1 + 6*N + 2*N*N
    liwork_h = 10*N
    iu = n
    il = 1

    allocate(A_h(N,N))
    allocate(Z_h(N,N))
    A_h(1, 1) = 1
    A_h(1, 2) = 2
    A_h(1, 3) = 3
    A_h(2, 1) = 2
    A_h(2, 2) = 4
    A_h(2, 3) = 1
    A_h(3, 1) = 3
    A_h(3, 2) = 1
    A_h(3, 3) = 1

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

    ! A_h(1, 1) = 1
    ! A_h(1, 2) = 2
    ! A_h(1, 3) = 3
    ! A_h(1, 4) = 4
    ! A_h(1, 5) = 5
    ! A_h(2, 1) = 2
    ! A_h(2, 2) = 4
    ! A_h(2, 3) = 1
    ! A_h(2, 4) = 1
    ! A_h(2, 5) = 5
    ! A_h(3, 1) = 3
    ! A_h(3, 2) = 1
    ! A_h(3, 3) = 1
    ! A_h(3, 4) = 2
    ! A_h(3, 5) = 3
    ! A_h(4, 1) = 4
    ! A_h(4, 2) = 1
    ! A_h(4, 3) = 2
    ! A_h(4, 4) = 2
    ! A_h(4, 5) = 3
    ! A_h(5, 1) = 5
    ! A_h(5, 2) = 5
    ! A_h(5, 3) = 3
    ! A_h(5, 4) = 3
    ! A_h(5, 1) = 1

    Z_h = A_h
    allocate(A, source=A_h)
    allocate(Z, source=Z_h)
    allocate(work_h(lwork))
    allocate(work, source=work_h)
    allocate(w_h(N))
    allocate(w, source=w_h)
    allocate(iwork_h(liwork_h))

    call dsyevd_gpu('V', 'U', il, iu, N, A, N, Z, N, w, work, lwork, &
        work_h, lwork, iwork_h, liwork_h, Z_h, N, w_h, info)
    !call dsyevd('V', 'U', N, A_h, N, w_h, work_h, lwork, &
    !iwork_h, liwork_h, info)
    Z_h = Z
    !call print_matrix(A_h)
    call print_matrix(Z_h)
end program