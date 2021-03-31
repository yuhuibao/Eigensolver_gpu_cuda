module dsytd2_gpu
  contains
    attributes(global) subroutine dsytd2_gpu(n,a,lda,d,e,tau)
      use cudafor
      implicit none
      integer, value    :: lda
      real(8),device    :: a(lda,32),tau(32)
      real(8),device    :: d(32),e(32)
      real(8),shared    :: a_s(32,32)
      real(8),shared    :: a_s_1(32,32)
      real(8),shared    :: alpha
      real(8),shared    :: taui
      real(8)           :: beta
      real(8)           :: alphar
      real(8)           :: xnorm,x,y,z,w
      real(8)           :: wc
      integer, value    :: n
      integer           :: tx,ty,tl,i,j,ii,stride

      tx=threadIdx%x
      ty=threadIdx%y
      ! Linear id of the thread (tx,ty)
      tl=tx+ blockDim%x*(ty-1)

      ! Load a_d in shared memory
      if (tx <= N .and. ty <= N) then
         a_s(tx           ,ty           )=a(tx           ,ty)
      endif
      call syncthreads()

      ! Symmetric matrix from upper triangular
      if (tx >ty) then
         a_s(tx,ty)=a_s(ty,tx)
      end if
      call syncthreads()

      ! For each column working backward
      ! do i=n-1,1,-1
      do i = n-1, 1, -1
        ! Generate elementary reflector
        ! Sum the vectors above the diagonal, only one warp active
        ! Reduce in a warp
        if (tl <=32) then
          if (tl <i) then
            a_s_1(tl,i+1)=a_s(tl,i+1)*a_s(tl,i+1)
          else 
            a_s_1(tl,i+1)=0._8
          endif
        end if
        call syncthreads()
        !   !  xnorm=__shfl_down(w,1)
        !   !  w=w+xnorm
        !   !  xnorm=__shfl_down(w,2)
        !   !  w=w+xnorm
        !   !  xnorm=__shfl_down(w,4)
        !   !  w=w+xnorm
        !   !  xnorm=__shfl_down(w,8)
        !   !  w=w+xnorm
        !   !  xnorm=__shfl_down(w,16)
        !   !  w=w+xnorm

        stride=2
        do while (stride <= 32)
          if (tl <= 32 .and. mod(tl-1, stride)==0) then 
            a_s_1(tl,i+1) = a_s_1(tl,i+1) + a_s_1(tl+stride/2, i+1)
          endif 
          stride = 2*stride
          call syncthreads()
        end do

        ! if(tl <= 32) then
        !   stride=2
        !   do while (stride<=32)
        !     if (mod(tl-1, stride) == 0) then
        !       a_s_1(tl,i+1) = a_s_1(tl,i+1) + a_s_1(tl+stride/2, i+1)
        !     end if
        !     call syncthreads()
        !     stride = 2*stride
        !   end do
        ! end if
        ! call syncthreads()

        print*, "75: taui when tx=",tx,",ty=",ty, taui

        if(tl==1) then
          w = a_s_1(tl,i+1)
          print*, "w= ",w
          alpha=a_s(i,i+1)
          alphar=dble(alpha)
          xnorm=dsqrt(w)
          
          if (xnorm .eq. 0_8) then
          ! H=1
            taui= 0._8
            alpha = 1.d0 ! To prevent scaling by dscal in this case
          else
            !Compute sqrt(alphar^2+xnorm^2) with  dlapy2(alphar,xnorm)
            x=abs(alphar)
            y=abs(xnorm)
            w=max(x,y)
            z=min(x,y)

            if (z .eq. 0.d0) then
              beta=-sign(w, alphar)
            else
              beta=-sign(w*sqrt(1.d0 + (z/w)**2), alphar)
            endif

            taui= (beta-alphar)/beta
            alpha = 1.d0/(alphar - beta) ! scale factor for dscal

          end if
        end if

        call syncthreads()

        print*, "109: taui when tx=",tx,",ty=",ty, taui

        ! dscal
        if (tl<i) then
          a_s(tl,i+1)=a_s(tl,i+1)*alpha
        end if

        if (tl==1) then 
          if (xnorm .ne. 0_8) then
            alpha=beta
          else
            alpha=a_s(i,i+1) ! reset alpha to original value
          endif

          e(i)=alpha
        end if

        print*, "121: taui when tx=",tx,",ty=",ty, taui

        if(taui.ne.(0.d0,0.d0)) then
          a_s(i,i+1)=1.d0
          call syncthreads()
          if(tl<=i) then
            tau(tl)=0.d0
            do j=1,i
              tau(tl)=tau(tl)+taui*a_s(tl,j)*a_s(j,i+1)
            end do
            print*, "tau(tl)=", tau(tl)
          end if

          call syncthreads()
       
          if (tl <=32) then
            if (tl <=i) then
              a_s_1(tl,i+1)=-.5d0*taui*tau(tl)*a_s(tl,i+1)
              print*,a_s_1(tl,i+1)
            else
              a_s_1(tl,i+1)=0._8
            endif
          end if
          call syncthreads() 

          if (tl <=32) then
            stride=2
            do while(stride<=32)
              if (mod(tl-1,stride) == 0) then
                a_s_1(tl,i+1) =a_s_1(tl,i+1) + a_s_1(tl+stride/2,i+1) 
              endif
              call syncthreads()
              stride = stride*2
            end do
            x = a_s_1(1,i+1)
            !print*,x
          end if
            ! z=__shfl_xor(x,1)
            ! x=x+z
            ! z=__shfl_xor(x,2)
            ! x=x+z
            ! z=__shfl_xor(x,4)
            ! x=x+z
            ! z=__shfl_xor(x,8)
            ! x=x+z
            ! z=__shfl_xor(x,16)
            ! x=x+z

          ! end if

          ! call syncthreads()
          ! if (tl == 1) then
          !   print*,x
          ! endif

          if (tl <=i) then
            tau(tl)=tau(tl)+x*a_s(tl,i+1)
            !print*,tau(tl)
          end if

          if(tl==1) alpha=x

          call syncthreads()

          print*, "Before Check i=", i
          if(tx<=i .and. ty<=i) then
            print*, "check:", tx, ty, a_s(tx,ty)
            ! a_s(tx,i+1),tau(ty),a_s(ty,i+1),tau(tx)
            ! a_s(tx,ty)=a_s(tx,ty)-a_s(tx,i+1)*tau(ty)-a_s(ty,i+1)*tau(tx)
          end if
          call syncthreads()

        endif

        if (tl==1) then
          a_s(i,i+1)=e(i)
          d(i+1)=a_s(i+1,i+1)
          tau(i)=taui
        end if

        call syncthreads()

      end do

      if (tl==1) then
        d(1) = a_s(1,1)
      endif

      call syncthreads()

      ! Back to device memory
      if (tx <= N .and. ty <= N) then
        print*, "Back to device memory"
        print*, "a_s(", tx, ",", ty, ")=", a_s(tx,ty) 
        a(tx,ty)=a_s(tx,ty)
      endif


    end subroutine dsytd2_gpu
end module dsytd2_gpu
  
