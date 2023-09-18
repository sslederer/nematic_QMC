!!!! Interface for lapack zgesvd
!!!! Compute A = U S V_dag. Note that A is destroyed in the process
subroutine zgesvd_(A, U, S, Vdag)
    use global
    implicit none

    type(CMatrix) A, U, Vdag
    real(8) S(A%d)
    complex(8), allocatable :: work(:)
    real(8) rwork(5*A%d)
    integer(4)    iwork(8*A%d)

    integer    lwork, info

!    allocate(rwork(5*A%d), iwork(8*A%d))

    lwork = -1
    allocate(work(4*A%d))
    call zgesvd( 'A', 'A', A%d, A%d, A%m, A%d, S, U%m, U%d, Vdag%m, Vdag%d, &
                    &     work, lwork, rwork, info )

!!!!!!!!!! zgesdd uses too much memory. use zgesvd
!    call zgesdd( 'A', A%d, A%d, A%m, A%d, S, U%m, U%d, Vdag%m, Vdag%d, &
!                    &     work, lwork, rwork, iwork, info )

    lwork = int(work(1))
    deallocate(work); allocate(work(lwork))

    call zgesvd( 'A', 'A', A%d, A%d, A%m, A%d, S, U%m, U%d, Vdag%m, Vdag%d, &
                    &     work, lwork, rwork, info )
!    call zgesdd( 'A', A%d, A%d, A%m, A%d, S, U%m, U%d, Vdag%m, Vdag%d, &
!                    &     work, lwork, rwork, iwork, info )

    IF( INFO.GT.0 ) THEN
        WRITE(*,*)'The algorithm computing SVD failed to converge.'
    END IF

    deallocate(work)

end subroutine zgesvd_

!!!! Interface for lapack zgesvd
!!!! Compute A = U S V_dag. Note that A is destroyed in the process
subroutine dgesvd_(A, U, S, Vdag)
    use global
    implicit none
    type(RMatrix) A, U, Vdag
    real(8) S(A%d)
    real(8), allocatable :: work(:)
    real(8) rwork(5*A%d)
    integer(4)    iwork(8*A%d)
    integer    lwork, info

!    allocate(rwork(5*A%d), iwork(8*A%d))

    lwork = -1
    allocate(work(4*A%d))
    call dgesvd( 'A', 'A', A%d, A%d, A%m, A%d, S, U%m, U%d, Vdag%m, Vdag%d, &
                    &     work, lwork, info )

!!!!!!!!!! zgesdd uses too much memory. use zgesvd
!    call zgesdd( 'A', A%d, A%d, A%m, A%d, S, U%m, U%d, Vdag%m, Vdag%d, &
!                    &     work, lwork, rwork, iwork, info )

    lwork = int(work(1))
    deallocate(work); allocate(work(lwork))

    call dgesvd( 'A', 'A', A%d, A%d, A%m, A%d, S, U%m, U%d, Vdag%m, Vdag%d, &
                    &     work, lwork, info )
!    call zgesdd( 'A', A%d, A%d, A%m, A%d, S, U%m, U%d, Vdag%m, Vdag%d, &
!                    &     work, lwork, rwork, iwork, info )

    IF( INFO.GT.0 ) THEN
        WRITE(*,*)'The algorithm computing SVD failed to converge.'
    END IF

    deallocate(work)

end subroutine dgesvd_

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Call zgemm to get A = B*C, all square matrices !!!!!!!
subroutine zsq_matmul(A,B,C,is_HC_B, is_HC_C)
    use global
    implicit none
    integer(4) d
    character*1 is_HC_B, is_HC_C
    type(CMatrix) A,B,C
    d = A%d
    call zgemm(is_HC_B,is_HC_C,d,d,d,(1.d0,0.d0), B%m,d, C%m,d, (0.d0,0.d0), A%m,d)
end subroutine zsq_matmul

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Call dgemm to get A = B*C, all square, real matrices !!!!!!!
subroutine dsq_matmul(A,B,C,is_HC_B, is_HC_C)
    use global
    implicit none
    integer(4) d
    character*1 is_HC_B, is_HC_C
    type(RMatrix) A,B,C
    d = A%d
    call dgemm(is_HC_B,is_HC_C,d,d,d,1.d0, B%m,d, C%m,d, 0.d0, A%m,d)
end subroutine dsq_matmul

subroutine get_matmul(A,B,C,is_HC_B, is_HC_C)
    use global
    implicit none
    integer(4) d
    character*1 is_HC_B, is_HC_C, is_Hc_B2, is_Hc_C2
    $G_type$ A,B,C
    d = A%d
    if(A%is_complex.eq.1) then
        call zgemm(is_HC_B,is_HC_C,d,d,d,(1.d0,0.d0), B%m,d, C%m,d, (0.d0,0.d0), A%m,d)
    else
        is_Hc_B2 = is_Hc_B
        is_Hc_C2 = is_Hc_C
        if (is_HC_B.eq.'C') is_Hc_B2 = 'T'
        if (is_HC_C.eq.'C') is_Hc_C2 = 'T'
        call dgemm(is_HC_B2,is_HC_C2,d,d,d,1.d0, B%m,d, C%m,d, 0.d0, A%m,d)
    endif
end subroutine get_matmul


! invert matrix using LAPACK functions. note that mat is overwritten by its inverse
subroutine matrix_inv(mat, d)
    implicit none
    integer(4) d, lwork, info
    integer ipiv(d)
    complex(8) mat(d,d)
    complex(8), allocatable :: work(:)

    !!!! LDU decomposition
    call zgetrf( d, d, mat, d, ipiv, info )

    !!!! Inquire about optimal lwork
    allocate(work(4*d))
    call zgetri( d, mat, d, ipiv, work, -1, info )
    lwork = int(work(1))

    deallocate(work); allocate(work(lwork))
    call zgetri( d, mat, d, ipiv, work, lwork, info )

    IF( INFO.GT.0 ) THEN
        WRITE(*,*)'The algorithm computing matrix inversion failed to converge.'
    END IF

    deallocate(work)
end subroutine matrix_inv

! invert real matrix using LAPACK functions. note that mat is overwritten by its inverse
subroutine Rmatrix_inv(mat, d)
    implicit none
    integer(4) d, lwork, info
    integer ipiv(d)
    real(8) mat(d,d)
    real(8), allocatable :: work(:)

    !!!! LDU decomposition
    call dgetrf( d, d, mat, d, ipiv, info )

    !!!! Inquire about optimal lwork
    allocate(work(4*d))
    call dgetri( d, mat, d, ipiv, work, -1, info )
    lwork = int(work(1))

    deallocate(work); allocate(work(lwork))
    call dgetri( d, mat, d, ipiv, work, lwork, info )

    IF( INFO.GT.0 ) THEN
        WRITE(*,*)'The algorithm computing matrix inversion failed to converge.'
    END IF

    deallocate(work)
end subroutine Rmatrix_inv


!!!! determinant of a matrix (non-destructive)
subroutine matrix_det(mat, d, det)
    implicit none
    integer(4) d, info, j
    integer ipiv(d)
    complex(8) mat(d,d), mat1(d,d), det

    !!!! LDU decomposition
    mat1 = mat
    call zgetrf( d, d, mat1, d, ipiv, info )
    det = (1.d0,0.d0)
    do j = 1,d
        det = det*mat1(j,j)
    end do
end subroutine matrix_det

!!!! determinant of a real matrix (non-destructive)
subroutine Rmatrix_det(mat, d, det)
    implicit none
    integer(4) d, info, j
    integer ipiv(d)
    real(8) mat(d,d), mat1(d,d), det

    !!!! LDU decomposition
    mat1 = mat
    call dgetrf( d, d, mat1, d, ipiv, info )
    det = 1.d0
    do j = 1,d
        det = det*mat1(j,j)
    end do
end subroutine Rmatrix_det

subroutine Rmatrix_det2X2(mat,det)
	implicit none
	real(8) mat(2,2)
	real(8) det

	det= mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)

end subroutine Rmatrix_det2X2

subroutine Cmatrix_det2X2(mat,det)
	implicit none
	complex(8) mat(2,2)
	complex(8) det

	det= mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)

end subroutine Cmatrix_det2X2


! test invert mat
!subroutine test_inv
!    use global
!    implicit none
!    complex(8) mat(2,2), det

!    mat = sigmay + 2*sigmaz
!    call matrix_inv(mat,2)
!!    call PRINT_MATRIX( 'inv=', 2, 2, mat, 2 )
!    call matrix_det(mat,2,det)
!end subroutine test_inv

! Determinant of a complex 4x4 matrix
subroutine det4x4(A, det)
    implicit none
    complex(8) A(4,4), det
    det = A(1,1)*A(2,2)*A(3,3)*A(4,4) - A(1,1)*A(2,2)*A(3,4)*A(4,3) - A(1,1)*A(3,2)*A(2,3)*A(4,4) + A(1,1)*A(3,2)*A(2,4)*A(4,3) &
      & + A(1,1)*A(4,2)*A(2,3)*A(3,4) - A(1,1)*A(4,2)*A(2,4)*A(3,3) - A(2,1)*A(1,2)*A(3,3)*A(4,4) + A(2,1)*A(1,2)*A(3,4)*A(4,3) &
      & + A(2,1)*A(3,2)*A(1,3)*A(4,4) - A(2,1)*A(3,2)*A(1,4)*A(4,3) - A(2,1)*A(4,2)*A(1,3)*A(3,4) + A(2,1)*A(4,2)*A(1,4)*A(3,3) &
      & + A(3,1)*A(1,2)*A(2,3)*A(4,4) - A(3,1)*A(1,2)*A(2,4)*A(4,3) - A(3,1)*A(2,2)*A(1,3)*A(4,4) + A(3,1)*A(2,2)*A(1,4)*A(4,3) &
      & + A(3,1)*A(4,2)*A(1,3)*A(2,4) - A(3,1)*A(4,2)*A(1,4)*A(2,3) - A(4,1)*A(1,2)*A(2,3)*A(3,4) + A(4,1)*A(1,2)*A(2,4)*A(3,3) &
      & + A(4,1)*A(2,2)*A(1,3)*A(3,4) - A(4,1)*A(2,2)*A(1,4)*A(3,3) - A(4,1)*A(3,2)*A(1,3)*A(2,4) + A(4,1)*A(3,2)*A(1,4)*A(2,3)
end subroutine det4x4



subroutine dgesvj_(A, U, S, Vd)
	use global
	implicit none
	type(rmatrix) A, U, Vd
	real(8) S(A%d)
	integer(4) d, lwork, info
	real(8), allocatable :: work(:)

	lwork = max(4, 2*A%d)
	allocate(work(lwork))


!	call dgesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, info)

	call dgesvj('G', 'U', 'V', A%d, A%d, A%m, A%d, S, 0, Vd%m, A%d, work, lwork, info)
	if (info.ne.0) then
		write(6,*) 'From dgesvj_: SVD did not converge, with info code', info
	endif
	U%m=A%m
	if(abs((work(1)-1.d0)).ge.(0.00001)) then
		write(6,*) 'From dgesvj_: overflow or underflow of some singular values. Scaling factor', work(1)
	endif

	deallocate(work)

end subroutine dgesvj_

subroutine dgejsv_(A, U, S, Vd)
	use global
	implicit none
	type(rmatrix) A, U, Vd
	real(8) S(A%d)
	real(8), allocatable :: work(:)
	integer(4), allocatable :: iwork(:)
	integer(4) lwork, info

	lwork=2*A%d*(6+2*A%d) ! I added a factor of 2 over the minimal lwork.
	allocate(iwork(4*A%d))
	allocate(work(lwork))


!	call dgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, work, lwork, iwork, info)

	call dgejsv('F', 'F', 'V', 'R', 'T', 'P', A%d, A%d, A%m, A%d, S, U%m, A%d, Vd%m, A%d, work, lwork, iwork, info)

	if(info.ne.0) then
		write(6,*) 'From dgejsv_: SVD did not converge, with error code', info
	endif

	if(abs(work(1)/work(2)-1.d0).ge.(0.00001)) then
		write(6,*) 'From dgesvj_: overflow or underflow of some singular values. Scaling factor', work(1)/work(2)
	end if


	deallocate(work)
	deallocate(iwork)

end subroutine dgejsv_


subroutine test_svd(A, U, S, Vdag)
    use global
    implicit none
    type(RMatrix) A, U, Vdag, B
    real(8) S(A%d)
	real(8), allocatable :: S2(:), S3(:)

	allocate(S2(A%d))
	allocate(S3(A%d))

	call allocate_rmatrix(B, A%d)

	B%m=A%m
	call dgejsv_(B, U, S3, Vdag)
	B%m=A%m
	call dgesvj_(B, U, S2, Vdag)
	B%m=A%m
	call dgesvd_(B, U, S, Vdag)

	write(6,*) 'dgesvd/dgesvj', sum(abs(S/S2-1.d0))
	write(6,*) 'dgesvd/dgejsv',sum(abs(S/S3-1.d0))
	write(6,*) 'dgesvj/dgejsv',sum(abs(S2/S3-1.d0))

	deallocate(S2)
	deallocate(S3)
	call deallocate_rmatrix(B)
end subroutine


subroutine get_svd(A, U, S, Vd)
	use global
	implicit none
	$G_type$ A, U, Vd
	real(8) S(A%d)

    if(A%is_complex.eq.1) then
        ! write(6,*) 'svd, complex'
        call zgesvd_(A, U, S, Vd)
    else
        ! write(6,*) 'svd, real'
        call dgesvd_(A, U, S, Vd)
    endif

end subroutine get_svd
