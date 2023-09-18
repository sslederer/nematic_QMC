subroutine test_loop_rotation(mc)
	use global
	implicit none
	type(mcconfig) mc
	type(Rmatrix) G_code
	integer(4) dir,nt,jmat,j,unitnum,is_even, is_inv
	integer(4) jmat_i,jmat_f, jmat_dir

	open(newunit=unitnum, file=trim(adjustl(mc%run_name))//'_tests.out' ,status='replace')
	call allocate_Rmatrix(G_code,mc%lat_dim)

	dir=LEFT
	call initG(mc)
	call rotateG(-dir, mc)

	
	if (dir.eq.LEFT) then
	    nt = mc%n;
	elseif (dir.eq.RIGHT) then 
	    nt = mod(mc%n,mc%Ntau)+1
	endif
    
    is_even=1-2*mod(nt,2);

	jmat_dir=-dir*is_even    
    	
	if (jmat_dir.eq.1) then
		jmat_i=1
		jmat_f=4
	else
		jmat_i=4
		jmat_f=1
	endif

	if (dir.eq.LEFT) then 
		is_inv=INV
	elseif (dir.eq.RIGHT) then
		is_inv=NO_INV
	endif

	do jmat=jmat_i,jmat_f,jmat_dir
		! rotate the current kinetic energy matrix to the other side
		call OperateT(nt, LEFT, jmat, is_inv, mc%dtau, mc%G, mc)
		call OperateT(nt, RIGHT, jmat, -is_inv, mc%dtau, mc%G, mc)
	end do
	!now change mc%n:
	mc%n=mc%n+dir
	
	G_code%m=mc%g%m
	
	call initG(mc)
	call rotateG(-dir,mc)
	call rotateG(-dir,mc)
	
	write(unitnum,*) sum(abs(mc%G%m-G_code%m))
	
end subroutine test_loop_rotation


!!!! Test subroutine det_ratio. Make sure you uncomment the calculation of mc%det in constructG.
subroutine test_rdet(mc)
    use global
    implicit none
    type(McConfig) mc
    real(8)  rdet, det1, det2, rdet_exact, error
    integer(4) jmat, jbond, jtau, eta_old, stat,unitnum, ntau, dir
    type(cmatrix) G0

	call allocate_rmatrix(G0, mc%lat_dim)
	
	open(newunit=unitnum, file=trim(adjustl(mc%run_name))//'_tests.out' ,status='replace')

	call initG(mc)
	det1=mc%det
	
	dir=LEFT

	jtau=mc%Ntau
	jtau=1
	jmat=1
	jbond=1
	eta_old=mc%o_p(jmat)%eta(jbond,jtau)
	call det_ratio(rdet, mc, jbond, jmat, eta_old)

	mc%o_p(jmat)%eta(jbond,jtau)=-eta_old
	call initG(mc)
	rdet_exact=mc%det/det1
	
	write(unitnum,*) rdet_exact
	write(unitnum,*) rdet
	
	close(unitnum)
	
	call deallocate_cmatrix(G0)
end subroutine test_rdet


subroutine test_update_G(mc)
	use global
	implicit none
	type(mcconfig) mc
	integer(4) unitnum, jmat, jtau, jbond, eta_old, dir,j
	type(rmatrix) G_new, G_old, G_new2
	real(8) error,error2
	
	call allocate_rmatrix(G_new, mc%lat_dim)
	call allocate_rmatrix(G_old, mc%lat_dim)
	call allocate_rmatrix(G_new2, mc%lat_dim)
	open(newunit=unitnum, file=trim(adjustl(mc%run_name))//'_tests.out' ,status='replace')
	
	dir = RIGHT
	call initG(mc)
!	call rotateG(-dir, mc)
	G_old%m=mc%G%m

!	jtau=mc%ntau
	jtau=1
	jmat=1
	jbond=1
	eta_old=mc%o_p(jmat)%eta(jbond,jtau)
	call update_G(dir, mc, jbond, jmat, eta_old)
	G_new%m=mc%G%m
	
	mc%o_p(jmat)%eta(jbond,jtau)=-eta_old
	call InitG(mc)
!	call rotateG(-dir,mc)
	error=sum(abs(mc%G%m-G_new%m))

	write(unitnum, *) (error)
	
	mc%o_p(jmat)%eta(jbond,jtau)=eta_old
	mc%G%m=G_old%m

	call deallocate_cmatrix(G_new)
	call deallocate_cmatrix(G_old)
	close(unitnum)
	
end subroutine

!subroutine test_operateT(mc)
!	use global
!	implicit none
!	type(mcconfig) mc
!	type(cmatrix) mat1,mat2
!	integer(4) jtau,dir,jmat,is_inv
!	
!	call allocate_cmatrix(mat1,mc%lat_dim)
!	call allocate_cmatrix(mat2,mc%lat_dim)
!	call initG(mc)
!	mat1%m=mc%G%m
!	mat2%m=mc%G%m
!	
!	is_inv=-1; jmat=4; dir=LEFT; jtau=1;
!	call OperateT_left(jtau, dir, jmat, is_inv, mc%dtau/2.d0, mat1, mc)
!	call OperateT_left_new(jtau, dir, 4, is_inv, mc%dtau/2.d0, mat2, mc)

!	write(6,*) '------------------------------'
!	write(6,*) sum(abs(mat1%m-mat2%m))
!	write(6,*) sum(abs(mat1%m))
!	write(6,*) sum(abs(mat2%m))
!	write(6,*) '------------------------------'
!	call deallocate_cmatrix(mat1)
!	call deallocate_cmatrix(mat2)

!end subroutine


!subroutine test_update_G(mc)
!    use global
!    implicit none 
!    type(McConfig) mc
!!    complex(8) mat(mc%lat_dim,mc%lat_dim)!, vmat1(4,4), vmat2(4,4), vmat3(4,4)
!    real(8) phi_new(3), rdet, det_updated
!    integer(4) jx, jy, j

!    phi_new = (/2.d0,1.d0,0.d0/)
!    phi_new = phi_new/sqrt(5.d0)

!    jx=1;jy=1
!    call det_ratio(rdet, LEFT, mc, jx, jy, phi_new)
!    det_updated = mc%det*rdet
!    
!    !! form_vmat test
!!    call form_vmat(vmat1, jx, jy,    phi_new(:),   NO_INV, mc)
!!    call form_vmat(vmat2, jx, jy,    phi_new(:),      INV, mc)
!!    vmat3 = matmul(vmat1,vmat2)

!    do j = 1, mc%num_kin_mat, 1
!        call OperateT(LEFT,  j, NO_INV, mc%dtau/2.d0, mc%G, mc)
!    end do                
!    do j = 1, mc%num_kin_mat, 1
!        call OperateT(RIGHT, j,    INV, mc%dtau/2.d0, mc%G, mc)
!    end do                

!    call savemat(mc%lat_dim, mc%G%m, 'G1.dat')

!    do jx=1,1
!        do jy=1,1
!            call update_G(LEFT, mc, jx, jy, phi_new)
!        end do    
!    end do
!    call savemat(mc%lat_dim, mc%G%m, 'Gu.dat')

!    do jx=1,1; do jy=1,1
!        mc%phi(jx, jy, mc%Ntau,:) = phi_new
!    end do; end do

!    call InitG(mc); !call rotateG_TTdag(mc)

!    call savemat(mc%lat_dim, mc%G%m, 'G2.dat')
!end subroutine

!subroutine test_reconstructG(mc)
!    use global
!    implicit none 
!    type(McConfig) mc   
!    integer(4) j
!    
!    call savemat(mc%lat_dim, mc%G%m, 'G.dat')

!    do j=1,3*mc%Ntau_mul; call rotateG(RIGHT, mc); end do
!    call savemat(mc%lat_dim, mc%G%m, 'G1.dat')

!    do j=1,3; call stab_reconstructG(RIGHT,mc); end do             ! ns: Nsvd   -> Nsvd-1
!    call savemat(mc%lat_dim, mc%G%m, 'G2.dat')
!    
!    do j=1,3; call stab_reconstructG(LEFT, mc); end do             ! ns  Nsvd-1 -> Nsvd 
!    call savemat(mc%lat_dim, mc%G%m, 'G3.dat')
!end subroutine

!subroutine test_rotateG(mc)
!    use global
!    implicit none 
!    type(McConfig) mc   
!    integer(4) j
!    
!    call savemat(mc%lat_dim, mc%G%m, 'G.dat')

!    do j=1,mc%Ntau_mul; call rotateG(RIGHT, mc); end do
!    call savemat(mc%lat_dim, mc%G%m, 'G1.dat')

!    do j=1,mc%Ntau_mul; call rotateG(LEFT, mc); end do             ! ns: Nsvd   -> Nsvd-1
!    call savemat(mc%lat_dim, mc%G%m, 'G2.dat')
!    
!end subroutine

!subroutine test_S_phi(mc)
!    use global
!    implicit none 
!    type(McConfig) mc   
!    real(8) S_phi, S_phi_new, dS1, dS2
!    real(8) phi_new(3)
!    integer(4) jx,jy,jt

!    phi_new = (/2.d0,1.d0,0.d0/)
!    phi_new = phi_new/sqrt(5.d0)

!    call calc_S_phi(S_phi,mc)
!    
!    jx=mc%Lx; jy=mc%Ly; jt=mc%Ntau
!    call delta_S_phi(dS1, jx, jy, jt, phi_new, mc)
!    
!    mc%phi(jx,jy,jt,:) = phi_new(:);
!    call calc_S_phi(S_phi_new,mc)
!    dS2 = S_phi_new - S_phi    
!end subroutine

!subroutine test_Gk(mc)
!    use global
!    implicit none 
!    type(McConfig) mc   
!    call initG(mc)
!    call center_G(mc)
!    call write_mean_Gk(mc%G, mc)
!end subroutine 

!subroutine test_global_update(mc)
!    use global
!    implicit none 
!    type(McConfig) mc   
!    call rescale_phi(mc)
!end subroutine


subroutine test_Gk(mc)
	use global
	implicit none
	type(mcConfig) mc
	type(Rmatrix) G, K(4)
	integer(4) jmat
	

end subroutine

