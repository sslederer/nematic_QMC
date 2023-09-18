subroutine OperateT_right(jtau, dir, jmat, is_inv, dtau, mat, mc)
    use global
    implicit none
    integer(4) dir, is_inv, jmat, ind1, ind2, jtau
    $G_type$ mat
    type(McConfig) mc
    real(8) dtau, alpha
    $thop_type$ thop
    $thop_type$, allocatable :: t(:,:)
    integer(4) j, k
    real clock_start, clock_finish
    integer(8) sclock_start,sclock_finish, rate

    call cpu_time(clock_start)
    call system_clock(sclock_start,rate)

	alpha =mc%par%alpha
	if((jmat.eq.1).or.(jmat.eq.2)) alpha=-alpha

    allocate(t(2*mc%kin_mat(jmat)%num_pairs,2))
    do j = 1, mc%kin_mat(jmat)%num_pairs
        thop = is_inv * dtau *mc%kin_mat(jmat)%thop(j)*(1.d0+alpha*mc%o_p(jmat)%eta(j,jtau))
        t(2*j-1,1) =  cosh(abs(thop))
        if (abs(thop).gt.1.d-10) then
			t(2*j-1,2) =  sinh(abs(thop)) * thop/abs(thop)
        else
			t(2*j-1,2) = 0.d0
		endif
        t(2*j,2) = t(2*j-1,1)
        t(2*j,1) = myconj(t(2*j-1,2))
    end do
    if (dir.eq.RIGHT) then
        do j = 1, mc%kin_mat(jmat)%num_pairs
            ! indices of sites to be multiplied
            ind1 = mc%kin_mat(jmat)%ind(j,1)
            ind2 = mc%kin_mat(jmat)%ind(j,2)
            mat%m(:,(/ind1,ind2/)) = matmul(mat%m(:,(/ind1,ind2/)),t((/2*j-1,2*j/),:))
        end do
    endif
    deallocate(t)
    call cpu_time(clock_finish)
    call system_clock(sclock_finish)
    mc%prof%t(2) = mc%prof%t(2) + clock_finish - clock_start
    mc%prof_s%t(2)=mc%prof_s%t(2)+dble(sclock_finish-sclock_start)/dble(rate)
end subroutine OperateT_right

subroutine OperateT_left(jtau, dir, jmat, is_inv, dtau, mat, mc)
    use global
    implicit none
    integer(4) dir, is_inv, jmat, ind1, ind2,jtau
    $G_type$ mat
    type(McConfig) mc
    real(8) dtau, alpha
	$thop_type$ thop
    $thop_type$, allocatable :: t(:,:)
    integer(4) j, k
    real clock_start, clock_finish
    integer(8) sclock_start,sclock_finish, rate

    call cpu_time(clock_start)
    call system_clock(sclock_start,rate)

	alpha =mc%par%alpha
	if((jmat.eq.1).or.(jmat.eq.2)) alpha=-alpha


    allocate(t(2*mc%kin_mat(jmat)%num_pairs,2))
    do j = 1, mc%kin_mat(jmat)%num_pairs
		thop = is_inv * dtau * mc%kin_mat(jmat)%thop(j)*(1.d0+alpha*mc%o_p(jmat)%eta(j,jtau))
        t(2*j-1,1) =  cosh(abs(thop))
        if (abs(thop).gt.1.d-10) then
			t(2*j-1,2) =  sinh(abs(thop)) * thop/abs(thop)
        else
			t(2*j-1,2) = 0.d0
		endif
        t(2*j,2) = t(2*j-1,1)
        t(2*j,1) = myconj(t(2*j-1,2))

    end do
    if(dir.eq.LEFT) then
        do k = 1,mc%lat_dim,2  !! This extra sum compared to operateT_right doesn't seem to matter at all. Maybe it's faster this way. Yoni.
            do j = 1, mc%kin_mat(jmat)%num_pairs
                ! indices of sites to be multiplied
                ind1 = mc%kin_mat(jmat)%ind(j,1)
                ind2 = mc%kin_mat(jmat)%ind(j,2)
                mat%m((/ind1,ind2/),(/k,k+1/)) = matmul(t((/2*j-1,2*j/),:),mat%m((/ind1,ind2/),(/k,k+1/)))
            end do
        end do
    endif
    deallocate(t)

    call cpu_time(clock_finish)
    call system_clock(sclock_finish)
    mc%prof%t(1) = mc%prof%t(1) + clock_finish - clock_start
    mc%prof_s%t(1)=mc%prof_s%t(1)+dble(sclock_finish-sclock_start)/dble(rate)
end subroutine OperateT_left

!!!! Operate with exp(-is_inv*dtau * K) which can be expressed as a product of factors of the form
!!!! exp(-is_inv*dtau*thop*sigmax) = cosh(nudtau*thop) - sinh(nudtau*thop)*sigmax
!!!! Operation from direction dir.
subroutine OperateT(jtau, dir, jmat, is_inv, dtau, mat, mc)
    use global
    implicit none
    integer(4) dir, is_inv, jmat, jtau
    type(CMatrix) mat
    type(McConfig) mc
    real(8) dtau

	if (mc%is_verbose.eq.1) then
		write(6,*) 'OperateT, jtau, dir, jmat, is_inv'
		write(6,*) jtau, dir, jmat, is_inv
	endif
    if (dir.eq.RIGHT) then
        call OperateT_right(jtau, dir, jmat, is_inv, dtau, mat, mc)
    elseif(dir.eq.LEFT) then
        call OperateT_left(jtau, dir, jmat, is_inv, dtau, mat, mc)
    endif
end subroutine OperateT


!OperateB(dir = 'R' or 'L', is_inv, mat, nt1, nt2, mc)
!    Operate with B matrix from left or from right
!    mat: matrix of dimension lat_dim
!    nt1<nt2: time indices of B matrices to apply
!    mc: structure with Monte Carlo configuration
!    (kinetic energy matrices, current \eta configuration)
subroutine OperateB(dir, is_inv, mat, nt1, nt2, mc)
    use global
    implicit none
    integer(4) dir, is_inv, nt1, nt2, j, jt, nt,jmat, time_dir
    $G_type$ mat
    type(McConfig) mc
    real(8) muExp


	if (mc%is_verbose.eq.1) then
		write(6,*) 'OperateB, dir, is_inv, nt1, nt2'
		write(6,*) dir, is_inv, nt1, nt2
	endif

    if(nt2.ge.nt1) then
    	time_dir=1
    elseif(nt2.lt.nt1) then
    	time_dir=-1
    end if
    muExp=exp(is_inv*mc%dtau*mc%par%mu)

	do jt = 0, nt2-nt1,time_dir
		if (dir.eq.LEFT) then
		    nt = nt1 + jt
		elseif(dir.eq.RIGHT) then
		    nt = nt2 - jt
		else
		    stop 'From OperateB. dir should be either LEFT or RIGHT'
		endif


		if (-dir * is_inv * (2 * mod(nt,2) - 1).eq.1) then
			do jmat=1,4
				call OperateT(nt, dir, jmat, is_inv, mc%dtau, mat, mc)
			end do
		else
			do jmat=4,1,-1
				call OperateT(nt, dir, jmat, is_inv, mc%dtau, mat, mc)
			end do
		end if

	    mat%m=muExp*mat%m;
    end do
end subroutine OperateB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Initialize Green's function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form G = (1 + B(Ntau) * ... *B(1))^-1
!!!!!!! where B(n) = exp(-dtau*V(n)) T_1 T_2 (Right orientation) !!!!!!!!!!!
!!!!!!! Every Ntau_mul matrix multiplies, perform an SVD stabilization!!!!!
!!!!!!! Store intermediate U,S,Vd matrices !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitG(mc)
    use global
    implicit none
    type(McConfig) mc
    $G_type$ A
    type(SVD) USV_dummy
    integer(4) dir, is_inv, j, j1, j2
    character(*), parameter :: message = 'Initialize Green function'

    if (mc%is_verbose.eq.1) then
		write(6,*) 'initG'
	endif

	dir=LEFT
	is_inv=NO_INV

    call allocate_matrix(A,mc%lat_dim)

    do j = 1, mc%Nsvd
        call eye(A)     ! Initialize the matrix A as a unity matrix

    	call OperateB(dir, is_inv, A, (j-1)*mc%Ntau_mul+1, min(j*mc%Ntau_mul,mc%Ntau), mc)

        if (j.eq.1) then
            ! Perform SVD on A. store result in mc%USV(1)
            call get_svd(A, mc%USV(1)%U, mc%USV(1)%S, mc%USV(1)%Vd)
            mc%USV(1)%n1 = 1; mc%USV(1)%n2 = min(mc%Ntau_mul,mc%Ntau)
        else
            ! Advance previous SVD: mc%USV(j-1). Store the result in mc%USV(j)
            call AdvanceSVD(LEFT, A, mc%USV(j-1), mc%USV(j))
            mc%USV(j)%n1 = 1; mc%USV(j)%n2 = min(j*mc%Ntau_mul,mc%Ntau)
        endif
	end do
    !!!!!! Form initial G by stabilized matrix inversion !!!!!!

    call allocate_usv(USV_dummy,mc%lat_dim)   ! Choice of USV_dummy found to be numerically stable: U = U(Ns), Vd = U(Ns)^{+}, S =1
    USV_dummy%S=1.d0
    USV_dummy%U%m  =  mc%USV(mc%Nsvd)%U%m
    USV_dummy%Vd%m =  mc%USV(mc%Nsvd)%U%m
    USV_dummy%Vd%m=myconj(transpose(USV_dummy%Vd%m))

    ! G = (1 + (U0*S0*V0d)*(Un*Sn*Vnd))^-1, where U0=Vd0=S0=1
    call construct_G(mc, USV_dummy,mc%USV(mc%Nsvd))

    mc%n  = mc%Ntau
    mc%ns = mc%Nsvd
    mc%dir = RIGHT

    call deallocate_usv(USV_dummy)
    call deallocate_matrix(A)

    call profiling(message)
end subroutine InitG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Compute G = (1 + (U1*S1*Vd1)*(U2*S2*Vd2))^-1 , and reset mc%det_change to 1.
subroutine construct_G(mc, USV1, USV2)
    use global
    implicit none
    type(McConfig) mc
    type(svd) USV1,USV2, USV
    integer(4) d

	d=mc%lat_dim
	call construct_G_helper(mc%USV_G,mc%S,d,USV1,USV2)

	call USV_to_mat(mc%G,mc%USV_G,NO_INV)

    mc%det_change = 1.d0

end subroutine construct_G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Calculate USV=(1+(U1*S1*Vd1)*(U2*S2*Vd2))^-1 , and S_inv=1/S.
subroutine construct_G_helper(USV,S_inv,d,USV1,USV2)
	use global
	implicit none
	integer(4) :: d, j
	type(SVD) USV, USV1, USV2
	real(8) S_inv(d)
    $G_type$ mat


    call allocate_matrix(mat,d)

        ! mat = V1d*U2
    call get_matmul(mat,USV1%Vd,USV2%U, 'N','N')
        ! mat = S1*mat
    do j=1,d; mat%m(j,:) = USV1%S(j)*mat%m(j,:); end do
        ! mat = mat*S2
    do j=1,d; mat%m(:,j) = mat%m(:,j)*USV2%S(j); end do
        ! mat = U1^{+}*Vd2^{+} + mat
    if(is_complex.eq.1) then
        call zgemm('C','C',d,d,d,(1.d0,0.d0), USV1%U%m,d, USV2%Vd%m,d, (1.d0,0.d0), mat%m,d)
    else
        call dgemm('T','T',d,d,d,1.d0, USV1%U%m,d, USV2%Vd%m,d, 1.d0, mat%m,d)
    endif
        ! mat = U*S*Vd
    call get_svd(mat,USV%U,USV%S,USV%Vd)
        ! G = Vd2^{+}*Vd^{+}*S^-1*U^{+}*U1^{+}
    call get_matmul(mat,USV2%Vd,USV%Vd, 'C','C')
    call get_matmul(USV%Vd,USV%U,USV1%U, 'C','C')
    USV%U%m=mat%m

	S_inv=USV%S !!! Save inverse singular values of G
    USV%S=1.d0/USV%S

    call deallocate_matrix(mat)

end subroutine construct_G_helper

!!calculates the product mat=U*S*Vd or mat=(USV)^(-1)
subroutine USV_to_mat(mat, USV,is_inv)
	use global
	implicit none
	$G_type$ mat1,mat
	type(SVD) USV
	integer(4) j,is_inv

	call allocate_matrix(mat1,USV%U%d)
	if(is_inv.eq.NO_INV) then
		do j=1,USV%U%d; mat1%m(:,j) = USV%U%m(:,j)*USV%S(j); end do
		call get_matmul(mat,mat1,USV%Vd, 'N','N')
	elseif(is_inv.eq.INV) then !!V S^(-1) Ud = (S^-1 *Vd)^(dagger) *Ud
		do j=1,USV%U%d; mat1%m(j,:) = USV%Vd%m(j,:)/USV%S(j); end do
		call get_matmul(mat,mat1,USV%U, 'C','C')
	end if

    call deallocate_matrix(mat1)

end subroutine


!!Safely multiplies two SVDs USV=(U1*S1*Vd1)*(U2*S2*Vd2) or their inverses.
subroutine multiply_USV(USV,USV1,USV2,is_inv_1,is_inv_2)
	use global
	implicit none
	type(SVD) USV1,USV2,USV
	$G_type$ mat
	integer(4) d,is_inv_1,is_inv_2,j,i

	d=USV1%U%d;
    call allocate_matrix(mat,d)

	if(is_inv_2.eq.NO_INV) then
		if(is_inv_1.eq.NO_INV) then
			!mat=Vd1*U2
			call get_matmul(mat,USV1%Vd,USV2%U,'N','N')
			!mat=S1*mat*S2
			do j=1,d; do i=1,d
				mat%m(i,j)=USV1%S(i)*mat%m(i,j)*USV2%S(j)
			end do; end do
			!U*S*Vd=mat
			call get_svd(mat,USV%U,USV%S,USV%Vd)
			!U=U1*U
			call get_matmul(mat,USV1%U,USV%U,'N','N')
			USV%U%m=mat%m
			!Vd=Vd*Vd2
			call get_matmul(mat,USV%Vd,USV2%Vd,'N','N')
			USV%Vd%m=mat%m
		elseif(is_inv_1.eq.INV) then
			!mat=Ud1*U2
			call get_matmul(mat,USV1%U,USV2%U,'C','N')
			!mat=S1^(-1)*mat1*S2
			do j=1,d; do i=1,d
				mat%m(i,j)=mat%m(i,j)*USV2%S(j)/USV1%S(i)
			end do; end do
			!mat=U*S*Vd
			call get_svd(mat,USV%U,USV%S,USV%Vd)
			!U=V1*U
			call get_matmul(mat,USV1%Vd,USV%U,'C','N')
			USV%U%m=mat%m
			!Vd=Vd*Vd2
			call get_matmul(mat,USV%Vd,USV2%Vd,'N','N')
			USV%Vd%m=mat%m
		end if
	elseif(is_inv_2.eq.INV) then
		if(is_inv_1.eq.NO_INV) then
			!mat=Vd1*V_2
			call get_matmul(mat,USV1%Vd,USV2%Vd,'N','C')
			!mat=S1*mat1/S2
			do j=1,d; do i=1,d
				mat%m(i,j)=USV1%S(i)*mat%m(i,j)/USV2%S(j)
			end do; end do
			!mat=U*S*Vd
			call get_svd(mat,USV%U,USV%S,USV%Vd)
			!U=U1*U
			call get_matmul(mat,USV1%U,USV%U,'N','N')
			USV%U%m=mat%m
			!Vd=Vd*U2^(-1)
			call get_matmul(mat,USV%Vd,USV2%U,'N','C')
			USV%Vd%m=mat%m
		elseif(is_inv_1.eq.INV) then
			!mat=Ud1*V_2
			call get_matmul(mat,USV1%U,USV2%Vd,'C','C')
			!mat=S1^(-1)*mat1/S2
			do j=1,d; do i=1,d
				mat%m(i,j)=mat%m(i,j)/USV2%S(j)/USV1%S(i)
			end do; end do
			!mat=U*S*Vd
			call get_svd(mat,USV%U,USV%S,USV%Vd)
			!U=V1*U
			call get_matmul(mat,USV1%Vd,USV%U,'C','N')
			USV%U%m=mat%m
			!Vd=Vd*U2^(-1)
			call get_matmul(mat,USV%Vd,USV2%U,'N','C')
			USV%Vd%m=mat%m
		end if
	end if

    call deallocate_matrix(mat)
end subroutine multiply_USV

subroutine allocate_matrix(M,d)
    use global
    implicit none
    $G_type$ M
    integer(4) d

    if(M%is_complex.eq.1) then
        call allocate_cmatrix(M,d)
    else
        call allocate_rmatrix(M,d)
    endif
end subroutine allocate_matrix

subroutine deallocate_matrix(M)
    use global
    implicit none
    $G_type$ M

    M%d = 0; deallocate(M%m)
end subroutine deallocate_matrix

subroutine allocate_cmatrix(M, d)
    use global
    implicit none
    type(CMatrix) M
    integer(4) d
    M%d = d;    allocate(M%m(d,d))
end subroutine allocate_cmatrix

subroutine allocate_rmatrix(M, d)
    use global
    implicit none
    type(rMatrix) M
    integer(4) d
    M%d = d;    allocate(M%m(d,d))
end subroutine allocate_rmatrix

subroutine deallocate_cmatrix(M)
    use global
    implicit none
    type(CMatrix) M
    M%d = 0; deallocate(M%m)
end subroutine deallocate_cmatrix

subroutine deallocate_rmatrix(M)
    use global
    implicit none
    type(RMatrix) M
    M%d = 0; deallocate(M%m)
end subroutine deallocate_rmatrix


subroutine allocate_cmatrix_nonsq(M, d1, d2)
    use global
    implicit none
    type(CMatrix_nonsq) M
    integer(4) d1, d2
    M%d1 = d1; M%d2 = d2;   allocate(M%m(d1,d2))
end subroutine allocate_cmatrix_nonsq

subroutine deallocate_cmatrix_nonsq(M)
    use global
    implicit none
    type(CMatrix_nonsq) M
    M%d1 = 0; M%d2 = 0; deallocate(M%m)
end subroutine deallocate_cmatrix_nonsq

subroutine allocate_usv(USV,d)
    use global
    implicit none
    type(SVD) USV
    integer(4) d
    if (USV%is_alloc.eq.0) then
        call allocate_matrix(USV%U, d)
        call allocate_matrix(USV%Vd,d)
        allocate(USV%S(d))
        USV%is_alloc = 1
    endif
end subroutine allocate_usv

subroutine deallocate_usv(USV)
    use global
    implicit none
    type(SVD) USV
    integer(4) d
    if (USV%is_alloc.eq.1) then
        call deallocate_matrix(USV%U)
        call deallocate_matrix(USV%Vd)
        deallocate(USV%S)
        USV%is_alloc = 0
    endif
end subroutine deallocate_usv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Generate identity matrix (assumes that A is allocated and A%d is defined)
subroutine eye(A)
    use global
    implicit none
    $G_type$ A
    integer(4) j

    if(A%is_complex.eq.1) then
        A%m = (0.d0,0.d0)
        do j = 1,A%d; A%m(j,j) = (1.d0,0.d0); end do
    else
        A%m = 0.d0
        do j = 1,A%d; A%m(j,j) = 1.d0; end do
    endif
end subroutine eye



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! stabilized matrix multiply by A from left/right using SVD !!!!!!
!!! Result in SVD form   !!!!!!!!!!!!!!!!!!!!!!
!!! dir=Left:  Unew Snew Vdnew = A* U S Vd  !!!!!!!!!!!!!
!!! dir=RIGHT: Unew Snew Vdnew = U S Vd *A !!!!!!!!!!!!!
subroutine AdvanceSVD(dir, A, USV, USVnew)
    use global
    implicit none

    type(McConfig) mc
    integer(4) dir, d, j, k
    $G_type$ A, Atmp, Mtmp !, U, Vd, Unew, Vdnew
    type(SVD) USV, USVnew

    !real(8) S(A%d), Snew(A%d)

    d = A%d
	call allocate_matrix(Atmp,d)
	call allocate_matrix(Mtmp,d)

    if (dir.eq.LEFT) then
        call get_matmul(Atmp,A,USV%U,'N','N')           ! Atmp = A*U
        do j=1,d                                                     ! Atmp = Atmp*S
            Atmp%m(:,j) = Atmp%m(:,j)*USV%S(j)
        end do

        call get_svd(Atmp, USVnew%U, USVnew%S, USVnew%Vd)        ! SVD Atmp: Atmp = Unew Snew Vdnew

        call get_matmul(Mtmp,USVnew%Vd,USV%Vd,'N','N')    ! Mtmp = Vdnew*Vd
        USVnew%Vd%m = Mtmp%m

    elseif(dir.eq.RIGHT) then
        call get_matmul(Atmp,USV%Vd,A,'N','N')           ! Atmp = Vd*A
        do k=1,d
            do j=1,d                                         ! Atmp = S*Atmp
                Atmp%m(j,k) = USV%S(j)*Atmp%m(j,k)
            end do
        end do
        call get_svd(Atmp, USVnew%U, USVnew%S, USVnew%Vd)        ! SVD Atmp: Atmp = Unew Snew Vdnew

        call get_matmul(Mtmp,USV%U,USVnew%U,'N','N')           ! Mtmp = U*Unew
        USVnew%U%m = Mtmp%m
    endif

    deallocate(Atmp%m, Mtmp%m)

end subroutine AdvanceSVD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Advance G from n to n+1 or n-1
!!!!!!!! By calculating G(n+1) = B(n+1)  G(n) B(n+1)^-1  (LEFT)
!!!!!!!!                G(n-1) = B(n)^-1 G(n)   B(n)     (RIGHT)
subroutine rotateG(dir, mc)
    use global
    implicit none
    integer(4) dir, n_new
    type(McConfig) mc

	if (mc%is_verbose.eq.1) then
		write(6,*) 'RotateG'
	endif

    if (dir.eq.LEFT) then
        n_new = mod(mc%n,mc%Ntau)+1
        call OperateB(LEFT, NO_INV,  mc%G, n_new, n_new, mc)
        call OperateB(RIGHT,   INV,  mc%G, n_new, n_new, mc)
        mc%n = n_new
    elseif(dir.eq.RIGHT) then
        call OperateB(LEFT,   INV,  mc%G, mc%n, mc%n, mc)
        call OperateB(RIGHT, NO_INV, mc%G, mc%n, mc%n, mc)
        n_new = mod(mc%n-2+mc%Ntau, mc%Ntau)+1
        mc%n = n_new
    else
        stop 'From rotateG. dir must be either LEFT or RIGHT'
    endif
end subroutine rotateG


!!!! Stabilized reconstruction of G(ns*Ntau_mul) = [1 + B(ns*Nm)...B(1)*B(Ns*Nm)...B(ns*Nm+1)]^-1
! If dir=RIGHT, it is assumed that B((ns+1)*Nm)... B(ns*Nm+1) where changed recently, and B(ns*Nm)   will be changed next.
! If dir=LEFT,                     B(ns*Nm)... B((ns-1)*Nm+1) where changed recently, and B(ns*Nm+1) will be changed next.
!
!!!! dir=RIGHT:
!           Assumes knowledge of:     USV(ns)    =  B(ns*Nm)...B(1)
!                                     USV(ns+2)  =  B(Ns*Nm)...B((ns+1)*Nm+1)
!           Algorithm:
!               Calculate:                   A =  B((ns+1)*Nm)... B(ns*Nm+1)
!               AdvanceSVD:            USV(ns+1) = USV(ns+2)*A = B(Ns*Nm)... B(ns*Nm+1)
!               construct_G:      G(ns*Nm) = [1 + USV(ns)*USV(ns+1)]^-1
!               Works for ns=0
!
!!!! dir=LEFT:
!           Assumes knowledge of:     USV(ns-1)  =  B((ns-1)*Nm)...B(1)
!                                     USV(ns+1)  =  B(Ns*Nm)...B(ns*Nm+1)
!           Algorithm:
!               Calculate:                   A =  B(ns*Nm)... B((ns-1)*Nm+1)
!               AdvanceSVD:            USV(ns) = A*USV(ns-1) = B(ns*Nm)...B(1)
!               construct_G:          G(ns*Nm) = [1 + USV(ns)*USV(ns+1)]^-1
!               Works for ns=Ns
subroutine stab_reconstructG(dir,mc)
    use global
    implicit none
    integer(4) ns,dir,n1,n2
    type(McConfig) mc
    $G_type$  A
    real clock_start, clock_finish
    integer(8) sclock_start,sclock_finish, rate

    if (mc%is_verbose.eq.1) then
		write(6,*) 'stab_reconstruct_G'
	endif

    call cpu_time(clock_start)
    call system_clock(sclock_start,rate)
    call allocate_matrix(A,mc%lat_dim)

    if (dir.eq.LEFT) then    ! increase ns by 1
        ns = mod(mc%ns,mc%Nsvd)+1
        n1 = (ns-1)*mc%Ntau_mul+1
        n2 = min(ns*mc%Ntau_mul,mc%Ntau)
        call eye(A)     ! Initialize the matrix A as a unity matrix
        call OperateB(LEFT, NO_INV, A, n1, n2, mc)
        if ((ns.eq.1).or.(ns.eq.mc%Nsvd)) then
            call get_svd(A, mc%USV(ns)%U, mc%USV(ns)%S, mc%USV(ns)%Vd)  ! SVD A
            mc%USV(ns)%n1 = n1; mc%USV(ns)%n2 = n2
        else
            call AdvanceSVD(LEFT, A, mc%USV(ns-1), mc%USV(ns))
            mc%USV(ns)%n1 = mc%USV(ns-1)%n1; mc%USV(ns)%n2 = n2
        endif
        if (ns.lt.mc%Nsvd) then
            call construct_G(mc, mc%USV(ns), mc%USV(ns+1))
        else
            call construct_G(mc, mc%USV(ns), mc%USV(ns-1))
        endif
    elseif(dir.eq.RIGHT) then      ! decrease ns by 1
        ns = mc%ns-1
        if (ns.lt.0) stop 'From stab_reconstructG. cannot be smaller than 0'
        n1 = ns*mc%Ntau_mul+1
        n2 = min((ns+1)*mc%Ntau_mul, mc%Ntau)
        call eye(A)     ! Initialize the matrix A as a unity matrix
        call OperateB(LEFT, NO_INV, A, n1, n2, mc)
        if ((ns.eq.(mc%Nsvd-1)).or.(ns.eq.0)) then
            call get_svd(A, mc%USV(ns+1)%U, mc%USV(ns+1)%S, mc%USV(ns+1)%Vd)  ! SVD A
            mc%USV(ns+1)%n1 = n1; mc%USV(ns+1)%n2 = n2
        else
            call AdvanceSVD(RIGHT, A, mc%USV(ns+2), mc%USV(ns+1))
            mc%USV(ns+1)%n1 = n1; mc%USV(ns+1)%n2 = mc%USV(ns+2)%n2
        endif
        if (ns.gt.0) then
            call construct_G(mc, mc%USV(ns), mc%USV(ns+1))
        else
            call construct_G(mc, mc%USV(2), mc%USV(1))
        endif
    else
        stop 'stab_From rotateG. dir must be either LEFT or RIGHT'
    endif
    if (ns.gt.0) then
        mc%ns = ns;
    else
        mc%ns = mc%Nsvd
    endif
    mc%n  = min(mc%ns*mc%Ntau_mul,mc%Ntau)
    call deallocate_matrix(A)

    call cpu_time(clock_finish)
    call system_clock(sclock_finish)
    if (dir.eq.LEFT) then
        mc%prof%t(3) = mc%prof%t(3) + clock_finish - clock_start
        mc%prof_s%t(3)=mc%prof_s%t(3)+dble(sclock_finish-sclock_start)/dble(rate)
    else
        mc%prof%t(4) = mc%prof%t(4) + clock_finish - clock_start
        mc%prof_s%t(4)=mc%prof_s%t(4)+dble(sclock_finish-sclock_start)/dble(rate)
    endif
end subroutine stab_reconstructG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Update determinant due to change in eta !!!!!!!!!
subroutine det_ratio(rdet, mc, jbond, jmat, eta_old, eta_new)
    use global
    implicit none
    type(McConfig) mc
    integer(4) jbond, ind1, ind2, jmat, j
    $op_type$ eta_old, eta_new
    real(8) rdet
    $thop_type$ c_rdet, factor
    $thop_type$, allocatable ::bmat(:,:),gmat(:,:)
    real(8) rdet_old

	if((is_local_int.ne.1).and.(eta_old.ne.1).and.(eta_old.ne.-1)) then
		write(6,*) 'From det ratio, bad value for eta old', eta_old ;
		stop
	end if
	allocate(bmat(2,2),gmat(2,2))

	!! Form Bnew * B^-1 - 1 =  (B^-1) * Bnew - 1 . In the nematic B,B' commute.
	factor = -(eta_old - eta_new) * mc%dtau * mc%kin_mat(jmat)%thop(jbond) * mc%par%alpha;
	if ((jmat.eq.1).or.(jmat.eq.2)) factor=-factor


	if (abs(factor).gt.1.d-10) then
		bmat(1,1)=cosh(abs(factor))-1.d0
		bmat(2,2)=bmat(1,1)
		bmat(1,2)= sinh(abs(factor)) * factor/abs(factor)
        bmat(2,1)=myconj(bmat(1,2))
		!! multiply by 1-G
		ind1 = mc%kin_mat(jmat)%ind(jbond,1)
		ind2 = mc%kin_mat(jmat)%ind(jbond,2)

		gmat(1,1)=1.d0-mc%G%m(ind1,ind1)
		gmat(1,2)=-mc%G%m(ind1,ind2)
		gmat(2,1)=-mc%G%m(ind2,ind1)
		gmat(2,2)=1.d0-mc%G%m(ind2,ind2)

		bmat = matmul(gmat,bmat)

		!! Add unity matrix
		do j=1,2; bmat(j,j) = bmat(j,j) + 1.d0; end do

		!! calculate determinant:
        c_rdet = bmat(1,1)*bmat(2,2)-bmat(1,2)*bmat(2,1)

		rdet=(abs(c_rdet))**2
	else
		rdet=1.d0
	endif

		deallocate(bmat,gmat)

        ! !!!testing:
        ! call det_ratio_old(rdet_old, mc, jbond, jmat, eta_old)
        ! if(abs(rdet-rdet_old).ge.0.00000001) then
        !     write(6,*) 'bad det ratio?', abs(rdet-rdet_old)/rdet
        ! else
        !     write(6,*) 'good det ratio', abs(rdet-rdet_old)/rdet
        ! endif

end subroutine det_ratio

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Update determinant due to change in eta !!!!!!!!!
subroutine det_ratio_old(rdet, mc, jbond, jmat, eta_old)
    use global
    implicit none
    type(McConfig) mc
    integer(4) jbond, ind1, ind2, jmat, eta_old, j
    real(8) rdet
    complex(8) c_rdet, factor
    complex(8), allocatable ::bmat(:,:),gmat(:,:)

	if((eta_old.ne.1).and.(eta_old.ne.-1)) then
		write(6,*) 'From det ratio, bad value for eta old', eta_old ;
		stop
	end if
	allocate(bmat(2,2),gmat(2,2))

	!! Form Bnew * B^-1 - 1 =  (B^-1) * Bnew - 1 . In the nematic B,B' commute.
	factor = -2 * mc%dtau * mc%kin_mat(jmat)%thop(jbond) * mc%par%alpha;
	if ((jmat.eq.1).or.(jmat.eq.2)) factor=-factor


	if (abs(factor).gt.1.d-10) then
		bmat(1,1)=cosh(abs(factor))-1.d0
		bmat(2,2)=bmat(1,1)
		bmat(1,2)=eta_old * sinh(abs(factor)) * factor/abs(factor)
		bmat(2,1)=bmat(1,2)

		!! multiply by 1-G
		ind1 = mc%kin_mat(jmat)%ind(jbond,1)
		ind2 = mc%kin_mat(jmat)%ind(jbond,2)

		gmat(1,1)=1.d0-mc%G%m(ind1,ind1)
		gmat(1,2)=-mc%G%m(ind1,ind2)
		gmat(2,1)=-mc%G%m(ind2,ind1)
		gmat(2,2)=1.d0-mc%G%m(ind2,ind2)

		bmat = matmul(gmat,bmat)

		!! Add unity matrix
		do j=1,2; bmat(j,j) = bmat(j,j) + 1.d0; end do

		!! calculate determinant:
		call cmatrix_det2X2(bmat, c_rdet)
		rdet=(abs(c_rdet))**2
	else
		rdet=1.d0
	endif

		deallocate(bmat,gmat)

end subroutine det_ratio_old

subroutine update_G(dir, mc, jbond, jmat, eta_old, eta_new)
    use global
    implicit none
    type(McConfig) mc
    integer(4) dir, jmat, jbond
    $op_type$ eta_old, eta_new

    if(is_complex.eq.1) then
        call update_G_complex(dir, mc, jbond, jmat, eta_old, eta_new)
    else
        call update_G_real(dir, mc, jbond, jmat, eta_old, eta_new)
    endif
end subroutine update_G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Update Green's function due to change in eta !!!!!!!!!
!!! dir=LEFT, G(n) is updated by the change in eta(n-1)
!!! dir=RIGHT, G(n) is updated by the change in eta(n)
subroutine update_G_real(dir, mc, jbond, jmat, eta_old,eta_new)
    use global
    implicit none
    type(McConfig) mc
    $op_type$ eta_old, eta_new
    integer(4) dir,j, d, jmat, jbond, ind1,ind2
    real(8), allocatable:: u(:,:), vd(:,:), m(:,:), mp(:,:), dG(:,:)
    real(8) factor
    real clock_start, clock_finish
    integer(8) sclock_start,sclock_finish, rate

    if((is_local_int.ne.1).and.(eta_old.ne.1).and.(eta_old.ne.-1)) then
		write(6,*) 'From update_G, bad value for eta old', eta_old
		stop
	end if
    call cpu_time(clock_start)
    call system_clock(sclock_start,rate)

    allocate(u(2,2), m(2,2))

  	!! Form u= Bnew * B^-1 - 1 =  (B^-1) * Bnew - 1 . In the nematic B,B' commute.

    factor = -(eta_old - eta_new) * mc%dtau * mc%kin_mat(jmat)%thop(jbond) * mc%par%alpha;
    if ((jmat.eq.1).or.(jmat.eq.2)) factor=-factor

    if (abs(factor).gt.1.d-10) then
        u(1,1)=cosh(abs(factor))-1.d0
        u(2,2)=u(1,1)
        u(1,2)=sinh(abs(factor)) * factor/abs(factor)
        u(2,1)=u(1,2)

        d = mc%lat_dim

    	ind1 = mc%kin_mat(jmat)%ind(jbond,1)
    	ind2 = mc%kin_mat(jmat)%ind(jbond,2)

        if (dir.eq.LEFT) then !G'=G*(1+(Bn'*Bn^(-1)-1)*(1-G))^(-1)=G-G*u*(1+vd*u)^(-1)*vd
    	    allocate(vd(2,d))
    	    !!! Form vd=(1-G)((/ind1,ind2/),:)
    	    vd=-mc%G%m((/ind1,ind2/),:)
    	    vd(1,ind1)=vd(1,ind1)+1.d0
    	    vd(2,ind2)=vd(2,ind2)+1.d0

    	    !! Form m=u*inv(1+vd(:,(/ind1,ind2/))*u)
    	    m=matmul(vd(:,(/ind1,ind2/)),u)
    	    do j=1,2; m(j,j)=m(j,j)+ 1.d0; end do
    	    call Rmatrix_inv(m,2)
    	    m=matmul(u,m)

    	    !! Form mp=m*Vd
    	    allocate(mp(2,d))
    	    call dgemm('n', 'n', 2, d, 2, 1.d0, m, 2, vd, 2, 0.d0, mp, 2)

    	    !! G=G-G(:,(/ind1,ind2/))*mp
    	    allocate(dG(d,d))
    	    call dgemm('n', 'n' ,d, d, 2, -1.d0, mc%G%m(:,(/ind1,ind2/)), d, mp, 2, 0.d0, dG, d)
    	    mc%G%m=mc%G%m+dG

        elseif(dir.eq.RIGHT) then !G'=G*(1+(1-G)*(Bn^(-1) Bn'-1))^(-1)=G-vd(1+u*vd)^(-1)*u*G
            allocate(vd(d,2))
            !!! Form vd=(1-G)(:,(/ind1,ind2/))
    	    vd=-mc%G%m(:,(/ind1,ind2/))
    	    vd(ind1,1)=vd(ind1,1)+1.d0
    	    vd(ind2,2)=vd(ind2,2)+1.d0

    	    !! Form m=inv(1+vd(:,(/ind1,ind2/))*u)*u
    	    m=matmul(u,vd((/ind1,ind2/),:))
    	    do j=1,2; m(j,j)=m(j,j)+ 1.d0; end do
    	    call Rmatrix_inv(m,2)
    	    m=matmul(m,u)

    	    !! Form mp=Vd*m
    	    allocate(mp(d,2))
    	    call dgemm('n', 'n', d, 2, 2, 1.d0, vd, d, m, 2, 0.d0, mp, d)

    	    !! G=G-mp*G((/ind1,ind2/),:)
    	    allocate(dG(d,d))
    	    call dgemm('n', 'n', d, d, 2, -1.d0, mp, d, mc%G%m((/ind1,ind2/),:), 2, 0.d0, dG, d)
    	    mc%G%m=mc%G%m+dG

        else
            stop '///!!! from update_G. dir must be LEFT or RIGHT !!!///'
        endif

    deallocate(dG,m,mp,vd,u)
    endif

    call cpu_time(clock_finish)
    call system_clock(sclock_finish)
    if (dir.eq.LEFT) then
        mc%prof%t(5) = mc%prof%t(5) + clock_finish - clock_start
        mc%prof_s%t(5)=mc%prof_s%t(5)+dble(sclock_finish-sclock_start)/dble(rate)
    else
        mc%prof%t(6) = mc%prof%t(6) + clock_finish - clock_start
        mc%prof_s%t(6)=mc%prof_s%t(6)+dble(sclock_finish-sclock_start)/dble(rate)
    endif
end subroutine update_G_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Update Green's function due to change in eta !!!!!!!!!
!!! dir=LEFT, G(n) is updated by the change in eta(n-1)
!!! dir=RIGHT, G(n) is updated by the change in eta(n)
subroutine update_G_complex(dir, mc, jbond, jmat, eta_old, eta_new)
    use global
    implicit none
    type(McConfig) mc
    integer(4) dir,j, d, jmat, jbond, ind1,ind2
    $op_type$ eta_old, eta_new
    complex(8), allocatable:: u(:,:), vd(:,:), m(:,:), mp(:,:), dG(:,:), G_old(:,:), G_old2(:,:)
    complex(8) factor
    real clock_start, clock_finish
    integer(8) sclock_start,sclock_finish, rate

	if((is_local_int.ne.1).and.(eta_old.ne.1).and.(eta_old.ne.-1)) then
		write(6,*) 'From update_G, bad value for eta old', eta_old
		stop
	end if
    call cpu_time(clock_start)
    call system_clock(sclock_start,rate)

    allocate(u(2,2), m(2,2))
! !!!! testing
!     allocate(G_old(mc%lat_dim, mc%lat_dim), G_old2(mc%lat_dim, mc%lat_dim))
!     G_old = mc%G%m
!     call update_G_old(dir, mc, jbond, jmat, eta_old)
!     G_old2 = mc%G%m
!     mc%G%m = G_old
! !!!!end this part of testing

  	!! Form u= Bnew * B^-1 - 1 =  (B^-1) * Bnew - 1 . In the nematic B,B' commute.
    factor = -(eta_old - eta_new) * mc%dtau * mc%kin_mat(jmat)%thop(jbond) * mc%par%alpha;
	if ((jmat.eq.1).or.(jmat.eq.2)) factor=-factor

	if (abs(factor).gt.1.d-10) then
		u(1,1)=cosh(abs(factor))-1.d0
		u(2,2)=u(1,1)
		u(1,2)=sinh(abs(factor)) * factor/abs(factor)
        u(2,1)=myconj(u(1,2))

		d = mc%lat_dim

		ind1 = mc%kin_mat(jmat)%ind(jbond,1)
		ind2 = mc%kin_mat(jmat)%ind(jbond,2)

		if (dir.eq.LEFT) then !G'=G*(1+(Bn'*Bn^(-1)-1)*(1-G))^(-1)=G-G*u*(1+vd*u)^(-1)*vd
			allocate(vd(2,d))
			!!! Form vd=(1-G)((/ind1,ind2/),:)
			vd=-mc%G%m((/ind1,ind2/),:)
			vd(1,ind1)=vd(1,ind1)+1.d0
			vd(2,ind2)=vd(2,ind2)+1.d0

			!! Form m=u*inv(1+vd(:,(/ind1,ind2/))*u)
			m=matmul(vd(:,(/ind1,ind2/)),u)
			do j=1,2; m(j,j)=m(j,j)+ 1.d0; end do
            if(is_complex.eq.1) then
                call matrix_inv(m,2)
            else
                call rmatrix_inv(m,2)
            endif
			m=matmul(u,m)

			!! Form mp=m*Vd
			allocate(mp(2,d))
            call zgemm('n', 'n', 2, d, 2, (1.d0,0.d0), m, 2, vd, 2, (0.d0,0.d0), mp, 2)

			!! G=G-G(:,(/ind1,ind2/))*mp
			allocate(dG(d,d))
			call zgemm('n', 'n' ,d, d, 2, (-1.d0,0.d0), mc%G%m(:,(/ind1,ind2/)), d, mp, 2, (0.d0,0.d0), dG, d)

            mc%G%m=mc%G%m+dG

		elseif(dir.eq.RIGHT) then !G'=G*(1+(1-G)*(Bn^(-1) Bn'-1))^(-1)=G-vd(1+u*vd)^(-1)*u*G
		    allocate(vd(d,2))
		    !!! Form vd=(1-G)(:,(/ind1,ind2/))
			vd=-mc%G%m(:,(/ind1,ind2/))
			vd(ind1,1)=vd(ind1,1)+1.d0
			vd(ind2,2)=vd(ind2,2)+1.d0

			!! Form m=inv(1+vd(:,(/ind1,ind2/))*u)*u
			m=matmul(u,vd((/ind1,ind2/),:))
			do j=1,2; m(j,j)=m(j,j)+ 1.d0; end do
			call matrix_inv(m,2)
			m=matmul(m,u)

			!! Form mp=Vd*m
			allocate(mp(d,2))
			call zgemm('n', 'n', d, 2, 2, (1.d0,0.d0), vd, d, m, 2, (0.d0,0.d0), mp, d)

			!! G=G-mp*G((/ind1,ind2/),:)
			allocate(dG(d,d))
			call zgemm('n', 'n', d, d, 2, (-1.d0,0.d0), mp, d, mc%G%m((/ind1,ind2/),:), 2, (0.d0,0.d0), dG, d)
			mc%G%m=mc%G%m+dG

		else
		    stop '///!!! from update_G. dir must be LEFT or RIGHT !!!///'
		endif

		deallocate(dG,m,mp,vd,u)
	endif

    call cpu_time(clock_finish)
    call system_clock(sclock_finish)
    if (dir.eq.LEFT) then
        mc%prof%t(5) = mc%prof%t(5) + clock_finish - clock_start
        mc%prof_s%t(5)=mc%prof_s%t(5)+dble(sclock_finish-sclock_start)/dble(rate)
    else
        mc%prof%t(6) = mc%prof%t(6) + clock_finish - clock_start
        mc%prof_s%t(6)=mc%prof_s%t(6)+dble(sclock_finish-sclock_start)/dble(rate)
    endif

    !!!!! testing
    ! if (sum(abs(mc%G%m-G_old2)).ge.0.0000001) then
    !     write(6,*) 'large difference in update G', sum(abs(mc%G%m-G_old2))
    ! else
    !     write(6,*) 'small difference in update G', sum(abs(mc%G%m-G_old2))
    ! endif
    ! deallocate(G_old, G_old2)
!!!!!!!!!
end subroutine update_G_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Update Green's function due to change in eta !!!!!!!!!
!!! dir=LEFT, G(n) is updated by the change in eta(n-1)
!!! dir=RIGHT, G(n) is updated by the change in eta(n)
subroutine update_G_old(dir, mc, jbond, jmat, eta_old)
    use global
    implicit none
    type(McConfig) mc
    integer(4) dir,j, d, jmat, jbond, ind1,ind2, eta_old
    complex(8), allocatable:: u(:,:), vd(:,:), m(:,:), mp(:,:), dG(:,:)
    complex(8) factor
    real clock_start, clock_finish
    integer(8) sclock_start,sclock_finish, rate

	if((eta_old.ne.1).and.(eta_old.ne.-1)) then
		write(6,*) 'From update_G, bad value for eta old', eta_old
		stop
	end if
    call cpu_time(clock_start)
    call system_clock(sclock_start,rate)

    allocate(u(2,2), m(2,2))

  	!! Form u= Bnew * B^-1 - 1 =  (B^-1) * Bnew - 1 . In the nematic B,B' commute.
	factor = -2 * mc%dtau * mc%kin_mat(jmat)%thop(jbond) * mc%par%alpha;
	if ((jmat.eq.1).or.(jmat.eq.2)) factor=-factor

	if (abs(factor).gt.1.d-10) then
		u(1,1)=cosh(abs(factor))-1.d0
		u(2,2)=u(1,1)
		u(1,2)=eta_old*sinh(abs(factor)) * factor/abs(factor)
		u(2,1)=myconj(u(1,2))

		d = mc%lat_dim

		ind1 = mc%kin_mat(jmat)%ind(jbond,1)
		ind2 = mc%kin_mat(jmat)%ind(jbond,2)

		if (dir.eq.LEFT) then !G'=G*(1+(Bn'*Bn^(-1)-1)*(1-G))^(-1)=G-G*u*(1+vd*u)^(-1)*vd
			allocate(vd(2,d))
			!!! Form vd=(1-G)((/ind1,ind2/),:)
			vd=-mc%G%m((/ind1,ind2/),:)
			vd(1,ind1)=vd(1,ind1)+1.d0
			vd(2,ind2)=vd(2,ind2)+1.d0

			!! Form m=u*inv(1+vd(:,(/ind1,ind2/))*u)
			m=matmul(vd(:,(/ind1,ind2/)),u)
			do j=1,2; m(j,j)=m(j,j)+ 1.d0; end do
			call matrix_inv(m,2)
			m=matmul(u,m)

			!! Form mp=m*Vd
			allocate(mp(2,d))
			call zgemm('n', 'n', 2, d, 2, (1.d0,0.d0), m, 2, vd, 2, (0.d0,0.d0), mp, 2)

			!! G=G-G(:,(/ind1,ind2/))*mp
			allocate(dG(d,d))
			call zgemm('n', 'n' ,d, d, 2, (-1.d0,0.d0), mc%G%m(:,(/ind1,ind2/)), d, mp, 2, (0.d0,0.d0), dG, d)
			mc%G%m=mc%G%m+dG

		elseif(dir.eq.RIGHT) then !G'=G*(1+(1-G)*(Bn^(-1) Bn'-1))^(-1)=G-vd(1+u*vd)^(-1)*u*G
		    allocate(vd(d,2))
		    !!! Form vd=(1-G)(:,(/ind1,ind2/))
			vd=-mc%G%m(:,(/ind1,ind2/))
			vd(ind1,1)=vd(ind1,1)+1.d0
			vd(ind2,2)=vd(ind2,2)+1.d0

			!! Form m=inv(1+vd(:,(/ind1,ind2/))*u)*u
			m=matmul(u,vd((/ind1,ind2/),:))
			do j=1,2; m(j,j)=m(j,j)+ 1.d0; end do
			call matrix_inv(m,2)
			m=matmul(m,u)

			!! Form mp=Vd*m
			allocate(mp(d,2))
			call zgemm('n', 'n', d, 2, 2, (1.d0,0.d0), vd, d, m, 2, (0.d0,0.d0), mp, d)

			!! G=G-mp*G((/ind1,ind2/),:)
			allocate(dG(d,d))
			call zgemm('n', 'n', d, d, 2, (-1.d0,0.d0), mp, d, mc%G%m((/ind1,ind2/),:), 2, (0.d0,0.d0), dG, d)
			mc%G%m=mc%G%m+dG

		else
		    stop '///!!! from update_G. dir must be LEFT or RIGHT !!!///'
		endif

		deallocate(dG,m,mp,vd,u)
	endif

    call cpu_time(clock_finish)
    call system_clock(sclock_finish)
    if (dir.eq.LEFT) then
        mc%prof%t(5) = mc%prof%t(5) + clock_finish - clock_start
        mc%prof_s%t(5)=mc%prof_s%t(5)+dble(sclock_finish-sclock_start)/dble(rate)
    else
        mc%prof%t(6) = mc%prof%t(6) + clock_finish - clock_start
        mc%prof_s%t(6)=mc%prof_s%t(6)+dble(sclock_finish-sclock_start)/dble(rate)
    endif


end subroutine update_G_old



!! calculates the SVD of B(ntau)... B1 from direction dir, storing intermediate SVD's.
!! If dir=LEFT, is_inv=NO_INV: BNtau...(B3*(B2*B1)). USV(n) is the SVD of  B(n*Ntau_mul)*...*B(2)*B(1) = B(n*Ntau_mul,0) for n=1...Nsvd
!! If dir=RIGHT, is_inv=NO_INV: ((BNtau*BNtau-1)*BNtau-2)...B1. USV(n) is the SVD of B(Ntau)*B(Ntau-1)*...B((n-1)*Ntau_mul+1) = B(Ntau,(n-1)*Ntau_mul) for n=1..Nsvd
!! If dir=LEFT, is_inv=INV: B1^(-1)...(BNtau-1^(-1) BNtau^(-1))). USV(n) is the SVD of B((n-1)*Ntau_mul+1)^(-1)....B(Ntau)^(-1) = B(Ntau,(n-1)*Ntau_mul)^(-1)
!! If dir=RIGHT, is_inv=INV: (B1^(-1)* B2^(-1))...... BNtau^(-1))). USV(n) is the SVD of B(1)^(-1)....B(n*Ntau_mul)^(-1) = B(n*Ntau_mul,0)^(-1)
subroutine calc_B_USV(mc,USV,dir,is_inv)
   use global
    implicit none
    type(McConfig) mc
    $G_type$ A
    type(SVD) USV(mc%Nsvd)
    integer(4) dir, j, k, n1,n2, is_inv

    call allocate_matrix(A,mc%lat_dim)

    if (is_inv.eq.NO_INV) then
		if (dir.eq.LEFT) then
			n1=1
			n2=mc%Nsvd
		elseif(dir.eq.RIGHT) then
			n1=mc%Nsvd
			n2=1
		endif
	elseif (is_inv.eq.INV) then
		if (dir.eq.RIGHT) then
			n1=1
			n2=mc%Nsvd
		elseif(dir.eq.LEFT) then
			n1=mc%Nsvd
			n2=1
		endif
	endif

	if(is_inv.eq.NO_INV) then
		do j = n1, n2, -dir
			call eye(A)     ! Initialize the matrix A as a unity matrix
			call OperateB(LEFT, NO_INV, A, (j-1)*mc%Ntau_mul+1, min(j*mc%Ntau_mul,mc%Ntau), mc)
			if (j.eq.n1) then
				! Perform SVD on A. store result in USV(j)
				call get_svd(A, USV(j)%U, USV(j)%S, USV(j)%Vd)
			else
				! Advance previous SVD: mc%USV(j+dir). Store the result in mc%USV(j)
				call AdvanceSVD(dir, A, USV(j+dir), USV(j))
			endif

			if (dir.eq.LEFT) then
				USV(j)%n1 = 1; USV(j)%n2 = min(j*mc%Ntau_mul,mc%Ntau)
			elseif(dir.eq.RIGHT) then
				USV(j)%n1 = (j-1)*mc%Ntau_mul+1; USV(j)%n2 = mc%Ntau
			endif
		end do

	elseif(is_inv.eq.INV) then
		do j = n1 , n2, dir
			call eye(A)     ! Initialize the matrix A as a unity matrix
			call OperateB(RIGHT, INV, A, min(j*mc%Ntau_mul,mc%Ntau), (j-1)*mc%Ntau_mul+1, mc)

			if (j.eq.n1) then
				! Perform SVD on A. store result in USV(j)
				call get_svd(A, USV(j)%U, USV(j)%S, USV(j)%Vd)
			else
				! Advance previous SVD: mc%USV(j+dir). Store the result in mc%USV(j)
				call AdvanceSVD(dir, A, USV(j-dir), USV(j))
			endif

			if (dir.eq.RIGHT) then
				USV(j)%n1 = 1; USV(j)%n2 = min(j*mc%Ntau_mul,mc%Ntau)
			elseif(dir.eq.LEFT) then
				USV(j)%n1 = (j-1)*mc%Ntau_mul+1; USV(j)%n2 = mc%Ntau
			endif
		end do
	endif

    call deallocate_matrix(A)

end subroutine calc_B_USV


!! Returns the SVD of C= (A^(is_inverse_A)+B^(is_inverse_B))^(-1), where A  and B are given in USV form. The algorithm follows Scallettar et al., shown here for is_inverse=NO_INV:
!! define SAp = max (A%S, 1) ,SAm = min (A%S, 1) , so that S%A = SAm*SAp, and similarly for B.
!! Then C= (B%Vd)^dagger * ( 1/SBp * (SAm * A%Vd * (B%Vd)^dagger /SBp + 1/(SAp) * (A%U)^dagger * B%U * SBm) ^(-1) ) / SAp ) * (A%U)^dagger
!! Where the intermediate inversion is calculated using an additional SVD, and another SVD is taken after the final multiplication by 1/SAp, 1/SBp.
subroutine stable_inverse_sum(C, A, B, is_inverse_A, is_inverse_B)
	use global
	implicit none
	type(SVD) A, B, C, X, Y
	integer(4) is_inverse_A, is_inverse_B, d, j,k
	real(8), allocatable:: SXp(:), SXm(:), SYp(:), SYm(:)
	$G_type$ mat1, mat2

	d=A%U%d

!! X=A^(is_inverse_A), Y=B^(is_inverse_B) . This can be optimized, but I don't think it matters.
	call allocate_usv(X,d)
	call allocate_usv(Y,d)

	call invert_svd(X,A, is_inverse_A)
	call invert_svd(Y,B, is_inverse_B)

	allocate(SXp(d), SXm(d), SYp(d), SYm(d))
	call allocate_matrix(mat1,d)
	call allocate_matrix(mat2,d)

!SXp = max (X%S, 1) ,SXm = min (X%S, 1) and similarly for Y.
	do j=1,d
		SXp(j)=max(X%S(j),1.d0)
		SXm(j)=min(X%S(j),1.d0)
		SYp(j)=max(Y%S(j),1.d0)
		SYm(j)=min(Y%S(j),1.d0)
	end do

!mat1= SXm * X%Vd * (Y%Vd)^dagger /SYp
	call get_matmul(mat1,X%Vd,Y%Vd,'N','C')
	do j=1,d; do k=1,d
		mat1%m(j,k)=mat1%m(j,k)*SXm(j)/SYp(k)
	end do; end do

!mat2 = 1/(SXp) * (X%U)^dagger * Y%U * SYm
	call get_matmul(mat2,X%U,Y%U,'C','N')
	do j=1,d; do k=1,d
		mat2%m(j,k)=mat2%m(j,k)*SYm(k)/SXp(j)
	end do; end do

!mat1 = mat1+mat2
	mat1%m=mat1%m+mat2%m

!invert mat1: mat1=mat1^(-1)
	call get_svd(mat1, C%U, C%S, C%Vd)
	call usv_to_mat(mat1,C,INV)

!mat1 = 1/SYp * mat1 /SXp
	do j=1,d; do k=1,d
		mat1%m(j,k)=mat1%m(j,k) / SYp(j) / SXp(k)
	end do; end do

!mat1 = U S Vd
	call get_svd(mat1, C%U,C%S,C%Vd)

! U = (Y%Vd)^dagger * U , Vd = Vd * (X%U)^dagger
	call get_matmul(mat1,Y%Vd,C%U,'C','N')
	call get_matmul(mat2,C%Vd,X%U,'N','C')
	C%U%m=mat1%m
	C%Vd%m=mat2%m


	deallocate(SXp,SXm,SYp,SYm)
	call deallocate_USV(X)
	call deallocate_USV(Y)
	call deallocate_matrix(mat1)
	call deallocate_matrix(mat2)

end subroutine stable_inverse_sum


!! Returns C= (1+B)^(-1), where  B is given in SVD form. The algorithm follows Scallettar et al.:
!! define Sp = max (B%S, 1) ,Sm = min (B%S, 1) , so that S%B = Sm*Sp.
!! Then C= (B%Vd)^dagger * ( 1/Sp * ((B%Vd)^dagger /Sp + B%U * Sm) ^(-1) ) )
!! Where the intermediate inversion is calculated using an additional SVD, and another SVD is taken after the final multiplication by 1/SAp.
subroutine stable_invert_B(C, B, is_inverse_B)
	use global
	implicit none
	type(SVD) B, C, X
	integer(4) is_inverse_B, d, j,k
	real(8), allocatable:: Sp(:), Sm(:)
	$G_type$ mat1, mat2
	real(8) avgS

	d=B%U%d
	allocate(Sp(d), Sm(d))
	call allocate_matrix(mat1,d)
	call allocate_matrix(mat2,d)

	call allocate_usv(X,d)
	! X = B^(is_inverse_B)
	call invert_svd(X,B, is_inverse_B)

!Sp = max (X%S, 1) ,Sm = min (X%S, 1)
	do j=1,d
		Sp(j)=max(X%S(j),1.d0)
		Sm(j)=min(X%S(j),1.d0)
	end do

!mat1= (X%Vd)^dagger /Sp
	do j=1,d; do k=1,d
		mat1%m(j,k)=myconj(X%Vd%m(k,j))/Sp(k)
	end do; end do


!mat2 = X%U * Sm
	do j=1,d; do k=1,d
		mat2%m(j,k)=X%U%m(j,k)*Sm(k)
	end do; end do

!mat1 = mat1+mat2
	mat1%m=mat1%m+mat2%m

!invert mat1: mat1=mat1^(-1)
	call get_svd(mat1, C%U, C%S, C%Vd)
	call usv_to_mat(mat1,C,INV)

!mat1 = 1/Sp * mat1
	do j=1,d; do k=1,d
		mat1%m(j,k)=mat1%m(j,k) / Sp(j)
	end do; end do

!mat1 = U S Vd
	call get_svd(mat1, C%U,C%S,C%Vd)

! U = (X%Vd)^dagger * U
	call get_matmul(mat1,X%Vd,C%U,'C','N')
	C%U%m=mat1%m

	deallocate(Sp,Sm)
	call deallocate_matrix(mat1)
	call deallocate_matrix(mat2)

end subroutine stable_invert_B


! if is_invert=INV, return USV_inv= USV^(-1). If is_invert=NO_INV return a copy USV_inv=USV.
subroutine invert_svd(USV_inv, USV, is_invert)
	use global
	implicit none
	type(svd) USV, USV_inv
	integer(4) is_invert

	if (is_invert.eq.INV) then
		USV_inv%U%m= myconj(transpose(USV%Vd%m))
		USV_inv%Vd%m= myconj(transpose(USV%U%m))
		USV_inv%S= 1.d0/USV%S
	elseif (is_invert.eq.NO_INV) then
		USV_inv%U%m=USV%U%m
		USV_inv%Vd%m=USV%Vd%m
		USV_inv%S = USV%S
	end if


end subroutine invert_svd


subroutine test_Ntau_mul(mc)
    use global
    implicit none
    type(McConfig) mc , mc2
    real(8) det_ratio
    integer(4) j

    call copy_param(mc2, mc)
    mc2%ntau_mul = 2

    write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(6,*) 'Comparing with ntau_mul=',mc2%ntau_mul

    call allocate_mc_config(mc2)
    call initG(mc2)
    det_ratio=1.d0/mc%det_change
    do j = 1,mc%lat_dim
        det_ratio = det_ratio * (mc2%S(j)/mc%S(j))**2
    end do
    write(6,*) 'det ratio=', det_ratio
    write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!'
end subroutine test_Ntau_mul
