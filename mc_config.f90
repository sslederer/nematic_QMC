!!!! Initialize lattice sizes, lattice indices, allocate kin_mat structures (not their matrices data), and storage SVD
subroutine allocate_mc_config(mc)
    use global
    implicit none
    type(McConfig) mc
    integer(4) lat_dim, j,num_sites, ind,j_x,j_y

    num_sites = mc%Lx*mc%Ly
    lat_dim = num_sites
    mc%lat_dim = lat_dim
    mc%opt%num_trial = 0
    mc%opt%num_trial_global=0
    mc%opt%num_accept = 0
    mc%opt%size_global_tot=0
    mc%opt%num_add_rho=0
    mc%opt%num_adjust_mu=50
    mc%prof%t = 0.d0
    mc%prof_s%t = 0.d0

    allocate(mc%kin_mat(mc%num_kin_mat))
    allocate(mc%o_p(mc%num_kin_mat))
    allocate(mc%i(mc%Lx,mc%Ly), mc%ix(lat_dim), mc%iy(lat_dim) )

    !!! Initialize index matrices
    ind = 0
    do j_y = 1,mc%Ly
        do j_x = 1,mc%Lx
            ind = ind+1
            mc%i(j_x, j_y) = ind
            mc%ix(ind)   = j_x
            mc%iy(ind)   = j_y
        end do
    end do

    !!! Initialize SVD storage matrices
    mc%Nsvd = int(dble(mc%Ntau)/dble(mc%Ntau_mul))
    if (mod(mc%Ntau,mc%Ntau_mul).gt.0) mc%Nsvd = mc%Nsvd+1

    !!!! Allocate SVD storage !!!!!!!
    allocate(mc%USV(mc%Nsvd))
    do j = 1,(mc%Nsvd)
        call allocate_usv(mc%USV(j), lat_dim)
    end do

    !!! allocate Green's function
    if(is_complex.eq.1) then
        call allocate_cmatrix(mc%G, lat_dim)
    else
        call allocate_rmatrix(mc%G, lat_dim)
    endif

    !!! allocate Green's function USV
    call allocate_usv(mc%USV_G, lat_dim)


    allocate(mc%S(lat_dim))
end subroutine allocate_mc_config



subroutine deallocate_mc_config(mc)
    use global
    implicit none
    type(McConfig) mc
    integer(4) j

    do j = 1,mc%num_kin_mat

        deallocate(mc%kin_mat(j)%thop) !why was this commented out?
        deallocate(mc%kin_mat(j)%ind) !why was this commented out?
        deallocate(mc%kin_mat(j)%pair_ind)

        deallocate(mc%o_p(j)%eta) !NEWSL
        deallocate(mc%o_p(j)%ind) !NEWSL

        deallocate(mc%o_p(j)%J) !NEWSL
        deallocate(mc%o_p(j)%Q) !NEWSL
    end do

    deallocate(mc%kin_mat)
    deallocate(mc%o_p) !NEWSL
    deallocate(mc%i, mc%ix, mc%iy )

    do j = 1,mc%Nsvd
        deallocate(mc%USV(j)%U%m, mc%USV(j)%Vd%m, mc%USV(j)%S)
        mc%USV(j)%is_alloc = 0
    end do
	deallocate(mc%USV_G%U%m, mc%USV_G%Vd%m, mc%USV_G%S)
    mc%USV_G%is_alloc = 0
    deallocate(mc%USV)
    deallocate(mc%G%m)
    deallocate(mc%S)
end subroutine deallocate_mc_config

!!! Assign values to kinetic energy matrix number jmat, in direction dir=HORIZONTAL or VERTICAL, with shift=0,1
subroutine assign_kin_mat(jmat, dir, shift, thop, mc)
    use global
    implicit none
    integer(4) jmat, jx,jy,jx1,jx2,jy1,jy2,step_x,step_y,dx,dy, jp, ind1, ind2, j, jt, dir, shift_x, shift_y,shift
    type(McConfig) mc
    real(8) thop
    complex(8) thop2
    real (8) a ! 0<=a<=1 determines the gauge choice

	mc%kin_mat(jmat)%dir = dir
	mc%kin_mat(jmat)%shift = shift

    jp = 0
    if (dir.eq.VERTICAL) then
        step_x = 1; step_y = 2; dx = 0; dy = 1; shift_x=0;     shift_y=shift
    elseif(dir.eq.HORIZONTAL) then
        step_x = 2; step_y = 1; dx = 1; dy = 0; shift_x=shift; shift_y=0
    endif

	do jx = (1+shift_x),mc%Lx, step_x
        do jy = (1+shift_y),mc%Ly, step_y
            jp = jp+1
            jx1 = jx; jx2 = mod(jx+dx-1, mc%Lx)+1
            jy1 = jy; jy2 = mod(jy+dy-1, mc%Ly)+1
            ind1 = mc%i(jx1,jy1)
            ind2 = mc%i(jx2,jy2)
            mc%kin_mat(jmat)%ind(jp,1) = ind1
            mc%kin_mat(jmat)%ind(jp,2) = ind2
            mc%kin_mat(jmat)%pair_ind(ind1,1) = ind2
            mc%kin_mat(jmat)%pair_ind(ind1,2) = jp
            mc%kin_mat(jmat)%pair_ind(ind1,3) = 1
            mc%kin_mat(jmat)%pair_ind(ind2,1) = ind1
            mc%kin_mat(jmat)%pair_ind(ind2,2) = jp
            mc%kin_mat(jmat)%pair_ind(ind2,3) = 2


!!! Landau gauge
!~          if(dir.eq.VERTICAL) mc%kin_mat(jmat)%thop(jp) =  thop
!~ 			if(dir.eq.HORIZONTAL) mc%kin_mat(jmat)%thop(jp) =  thop * exp(-2 * Pi * (0.d0,1.d0) * dble(mc%Fz)/dble(mc%Lx * mc%Ly)  * (jy1-1))
!~ 			if ((dir.eq.VERTICAL).and.(jy1.eq.mc%Ly)) then
!~ 						mc%kin_mat(jmat)%thop(jp) = mc%kin_mat(jmat)%thop(jp) * exp((0.d0,1.d0) * 2 * Pi * (mc%Fy + mc%Fz * dble(jx1-1)/dble(mc%Lx))) !Verified the sign of Fz.
!~ 			endif
!~ 			if ((dir.eq.HORIZONTAL).and.(jx1.eq.mc%Lx)) then
!~ 				mc%kin_mat(jmat)%thop(jp) = mc%kin_mat(jmat)%thop(jp) * exp((0.d0,1.d0) * 2 * Pi * mc%Fx)
!~ 			endif


!!! gauge determined by a:
			a = mc%gauge_factor
			if(dir.eq.VERTICAL) thop2 =  thop * exp( 2 * ( 1.d0- a) * Pi * (0.d0,1.d0) * dble(mc%Fz)/dble(mc%Lx * mc%Ly)  * (jx1-1))
			if(dir.eq.HORIZONTAL) thop2 =  thop * exp(- 2 * a * Pi * (0.d0,1.d0) * dble(mc%Fz)/dble(mc%Lx * mc%Ly)  * (jy1-1))

			if ((dir.eq.VERTICAL).and.(jy1.eq.mc%Ly)) then
				thop2 = thop2 * exp((0.d0,1.d0) * 2 * Pi * (mc%Fy + a * mc%Fz * dble(jx1-1)/dble(mc%Lx)))

			endif
			if ((dir.eq.HORIZONTAL).and.(jx1.eq.mc%Lx)) then
				thop2 = thop2 * exp((0.d0,1.d0) * 2 * Pi * (mc%Fx - (1.d0 -a) * mc%Fz * dble(jy1-1)/dble(mc%Ly)))
			endif

            if (is_complex.eq.1) then
                mc%kin_mat(jmat)%thop(jp) = thop2
            else
                mc%kin_mat(jmat)%thop(jp) = dble(thop2)
                if (abs(imag(thop2)).ge.1.d-8) write(6,*) 'Imaginary thop while working with real numbers. Im(thop)=', imag(thop2)
            endif

!~             !!!!  Periodic/Antiperiodic BC in x/y direction?  !!!!!
!~             if ( ((dir.eq.VERTICAL)  .and.(jy1.eq.mc%Ly).and.(mc%bcy.eq.'antiperiodic')).or. &
!~             &    ((dir.eq.HORIZONTAL).and.(jx1.eq.mc%Lx).and.(mc%bcx.eq.'antiperiodic')) ) then
!~                 mc%kin_mat(jmat)%thop(jp) = -cmplx(thop)
!~             else
!~                 mc%kin_mat(jmat)%thop(jp) =  cmplx(thop)
!~             endif
        end do
    end do


	!set indices for order parameter field NEWSL
	jp=0
	do jx = (1+shift_x),mc%Lx, step_x
		do jy = (1+shift_y),mc%Ly, step_y
			jp = jp+1
			mc%o_p(jmat)%J(jx,jy)=jp
			jx1 = jx; jx2 = mod(jx+dx-1, mc%Lx)+1
			jy1 = jy; jy2 = mod(jy+dy-1, mc%Ly)+1
		    ind1 = mc%i(jx1,jy1)
		    ind2 = mc%i(jx2,jy2)
		    mc%o_p(jmat)%ind(jp,1) = ind1
		    mc%o_p(jmat)%ind(jp,2) = ind2
		end do
	end do


	!Lookup table for nearest neighbor bonds of a given bond NEWSL

	do j=1, mc%o_p(jmat)%num_pairs
		jx=mc%ix(mc%o_p(jmat)%ind(j,2))
		jy=mc%iy(mc%o_p(jmat)%ind(j,2))

	!the (j,n,1) component of Q specifies which OP matrix corresponds to the nearest neighbor bond n of bond j (where n runs from 1 to 4, and the initial site
	!of bond 1 is the second site of bond j, proceeding clockwise
		if (dir.eq.HORIZONTAL) then
			if (mod(jy,2).eq.1) then
				mc%o_p(jmat)%Q(j,1,1)=3; mc%o_p(jmat)%Q(j,2,1)=4 ; mc%o_p(jmat)%Q(j,3,1)=4; mc%o_p(jmat)%Q(j,4,1)=3
			else
				mc%o_p(jmat)%Q(j,1,1)=4; mc%o_p(jmat)%Q(j,2,1)=3 ; mc%o_p(jmat)%Q(j,3,1)=3; mc%o_p(jmat)%Q(j,4,1)=4
			end if
		else
			if (mod(jx,2).eq.1) then
				mc%o_p(jmat)%Q(j,1,1)=1; mc%o_p(jmat)%Q(j,2,1)=1 ; mc%o_p(jmat)%Q(j,3,1)=2; mc%o_p(jmat)%Q(j,4,1)=2
			else
				mc%o_p(jmat)%Q(j,1,1)=2; mc%o_p(jmat)%Q(j,2,1)=2 ; mc%o_p(jmat)%Q(j,3,1)=1; mc%o_p(jmat)%Q(j,4,1)=1
			end if
		end if
	end do

end subroutine assign_kin_mat

!NEWSL
subroutine assign_o_p(jmat,mc)
    use global
    implicit none
    integer(4) jmat, jx,jy, j,jxm,jym
    type(McConfig) mc


	do j=1, mc%o_p(jmat)%num_pairs
		jx=mc%ix(mc%o_p(jmat)%ind(j,2))
		jy=mc%iy(mc%o_p(jmat)%ind(j,2))

		jxm=mod(mc%Lx+jx-2,mc%Lx)+1
		jym=mod(mc%Ly+jy-2,mc%Ly)+1

	! the (j,n,2) component of Q specifies the bond index
		mc%o_p(jmat)%Q(j,1,2)=mc%o_p(mc%o_p(jmat)%Q(j,1,1))%J(jx,jy)
		mc%o_p(jmat)%Q(j,2,2)=mc%o_p(mc%o_p(jmat)%Q(j,2,1))%J(jx,jym)
		mc%o_p(jmat)%Q(j,3,2)=mc%o_p(mc%o_p(jmat)%Q(j,3,1))%J(jxm,jym)
		mc%o_p(jmat)%Q(j,4,2)=mc%o_p(mc%o_p(jmat)%Q(j,4,1))%J(jxm,jy)
	end do

end subroutine assign_o_p

subroutine assign_bosonic_action(mc)
    use global
    implicit none
    type(mcconfig) mc
    real(8), allocatable:: G0(:,:,:,:), A(:,:), temp(:,:,:,:),S(:)
    type(rmatrix) B, U, Vd
    integer(4) jmat1, jmat2, jp1, jp2, jmat3, jp3, jnmat, num_pairs
    real(8) tr

    num_pairs = mc%Lx * mc%Ly /2

    allocate(mc%S0(num_pairs, 4, num_pairs, 4 ), G0(num_pairs, 4, num_pairs, 4 ), A(4*num_pairs, 4*num_pairs))
    G0 = 0.d0

    do jmat1=1,4
        do jp1=1,num_pairs
            G0(jp1, jmat1, jp1, jmat1) = G0(jp1, jmat1, jp1, jmat1) + mc%par%u2
            do jnmat =1,4
                jmat2 = mc%o_p(jmat1)%Q(jp1,jnmat,1)
                jp2 = mc%o_p(jmat1)%Q(jp1,jnmat,2)
                G0(jp1, jmat1, jp2, jmat2) = G0(jp1, jmat1, jp2, jmat2) - mc%par%u1  ! the negative sign here means the nematic phase is a 'ferromagnetic' configuration of the etas.
            enddo
        enddo
    enddo

    G0 = G0 / (mc%par%u2 + 4.d0 * mc%par%u1 )
    A = reshape(G0, (/4 * num_pairs, 4* num_pairs/))
    ! write(6,*) '!!!!!!!!!!!!'
    ! do jmat1=1,4
    !     do jmat2 = 1,jmat1
    !         write(6,*) '!!!!!!!!!!!!'
    !         write(6,*) 'jmat1,jmat2=', jmat1, jmat2
    !         write(6,*) G0(:,jmat1, :,jmat2)
    !     enddo
    ! enddo
    call Rmatrix_inv(A, 4 * num_pairs)


    mc%S0 = reshape(A, (/num_pairs, 4, num_pairs, 4/))

    ! call allocate_rmatrix(B, 4*num_pairs)
    ! call allocate_rmatrix(U, 4*num_pairs)
    ! call allocate_rmatrix(Vd, 4*num_pairs)
    ! allocate(S(4 * num_pairs))
    ! B%m = A
    ! call dgesvd_(B, U, S, Vd)
    ! write(6,*) S

    ! !!! testing to make sure I got the inversion right:
    ! allocate(temp(num_pairs, 4, num_pairs, 4))
    ! temp = 0.d0
    ! do jmat1 = 1,4; do jp1 = 1,num_pairs
    !     do jmat2 = 1,4; do jp2 = 1,num_pairs
    !         do jmat3 = 1,4 ; do jp3 = 1, num_pairs
    !             temp(jp1, jmat1, jp2, jmat2) = temp(jp1, jmat1, jp2, jmat2) + G0(jp1, jmat1, jp3, jmat3) * mc%S0(jp3, jmat3, jp2, jmat2)
    !         enddo; enddo
    !     enddo; enddo
    !     temp(jp1, jmat1, jp1, jmat1) = temp(jp1, jmat1, jp1, jmat1) - 1.0
    ! enddo; enddo
    ! write(6,*) 'delta(S0^-1 * S0 -1)=', sum(abs(temp))
    ! temp = 0.d0
    ! do jmat1 = 1,4; do jp1 = 1,num_pairs
    !     do jmat2 = 1,4; do jp2 = 1,num_pairs
    !         temp(jp2, jmat2, jp1, jmat1) = mc%S0(jp1, jmat1, jp2, jmat2)
    !     enddo; enddo
    ! enddo; enddo
    ! write(6,*) 'delta (S0 - transpose(S0))' , sum(abs(temp-mc%S0))
    ! deallocate(temp)

    ! mc%S0 = mc%S0 / 4.0 ! just a matter of definition...
    deallocate(G0, A)

end subroutine

!!!!! Allocate kin_mat, initialize num_pairs
subroutine allocate_kin_mat(jmat, num_pairs,mc)
    use global
    implicit none
    integer(4) jmat, num_pairs
    type(McConfig) mc

    mc%kin_mat(jmat)%num_pairs = num_pairs
    allocate(mc%kin_mat(jmat)%ind(num_pairs,2), mc%kin_mat(jmat)%thop(num_pairs), mc%kin_mat(jmat)%pair_ind(mc%lat_dim,3))
end subroutine

!!!! Allocate order parameter, initialize num_pairs (here differs between kin_mat and o_p) NEWSL
subroutine allocate_o_p(jmat, num_pairs, Ntau,mc)
    use global
    implicit none
    integer(4) jmat, num_pairs, Ntau
    type(McConfig) mc

	mc%o_p(jmat)%num_pairs = num_pairs
	allocate(mc%o_p(jmat)%ind(num_pairs,2), mc%o_p(jmat)%eta(num_pairs,Ntau), mc%o_p(jmat)%J(mc%Lx,mc%Ly), mc%o_p(jmat)%Q(num_pairs,4,2))

end subroutine

subroutine assert_L(Lx,Ly,Ntau,Ntau_mul)
	implicit none
    integer(4) Lx,Ly,Ntau,Ntau_mul
    if ((mod(Lx,2).ne.0).or.(mod(Ly,2).ne.0)) then
        stop '/////    Lx and Ly must currently be even. ///////'
    endif
    if ((Lx.eq.0).or.(Ly.eq.0)) then
        stop '///// Lx and Ly must not be 0. /////'
    endif
    if (mod(Ntau,2).ne.0) then
		stop '////// Ntau must be even. /////'
	end if
	if (Ntau.eq.Ntau_mul) then
		stop '////// Ntau_mul must be different from Ntau. /////'
	end if
	if (mod(Ntau,Ntau_mul).ne.0) then
		stop '////// Ntau_mul should be a divisor of Ntau //////'
	end if

end subroutine

!!!!!! Initialize parameters
subroutine init_par(mc)
    use global
    use random_numbers
    implicit none
    integer(4) Lx, Ly, num_kin_mat, jp, ind1, ind2, k, j, jt, jx, jy, jmat
    integer(4) Ntau, Ntau_mul, seed, max_mc_it, n_write_config, n_write_out
    integer(4) is_overwrite, is_cont_run, num_global_update, num_check_global, min_cluster_size,n_eq
    integer(4) is_rand_initial, is_run_fermions
    type(McConfig) mc
    real(8) thop, beta, mu, alpha, V, h, hz, U1, U2, Vt,ht,dVt,dht, Fx, Fy, gauge_factor
    integer(4) Fz
    character(100) run_name, config_file_name, run_name_prev, run_name_new, fnum, bcx, bcy, purpose

	integer(4) is_measure_G, is_measure_Gk, is_measure_ETCF, is_measure_ETPC
    integer(4) is_measure_tdgf, is_measure_tdgf_K, is_measure_tdcf
	integer(4) is_measure_cccf, is_measure_tdpc, is_measure_A, is_measure_E
	integer(4) is_measure_chi, is_measure_vertex, n_write_meas

    logical :: file_exists
	integer(4) unitnum

    namelist/mc_param/   Lx, Ly, Ntau, Ntau_mul, seed, max_mc_it,&
                       & run_name, n_write_config, n_write_out, is_overwrite,&
                       & is_cont_run, config_file_name, is_rand_initial, num_global_update, num_check_global, min_cluster_size,&
                       & purpose, is_run_fermions, n_eq
    namelist/model_param/beta, thop, mu, alpha, Fx, Fy, Fz, gauge_factor, V, h, U1, U2, hz,Vt,ht,dVt,dht
    namelist/meas_param/	is_measure_G, is_measure_Gk, is_measure_ETCF, is_measure_ETPC, is_measure_tdgf,&
							& is_measure_tdgf_K, is_measure_tdcf, is_measure_cccf, is_measure_tdpc, is_measure_A,&
							& is_measure_E, is_measure_chi, is_measure_vertex, n_write_meas

    !!!!!!! Define PI !!!!!!!!
    PI = 4.d0*atan(1.d0);

    !!!!!!!! defaults !!!!!!!!!!
    purpose = 'dqmc'
    n_write_config = 100
    n_write_out    = 500
    is_cont_run = 0
    num_global_update = 0
    config_file_name = 'none'
    bcx = 'periodic'
    bcy = 'periodic'
    is_rand_initial = 1
    is_run_fermions= 1
    num_check_global=10
	min_cluster_size=1

	is_measure_G = 0
	is_measure_Gk = 1
	is_measure_E = 1
	is_measure_ETCF = 0
	is_measure_ETPC = 0
	is_measure_tdgf = 0
	is_measure_tdgf_K = 0
	is_measure_tdcf = 0
	is_measure_cccf = 0
	is_measure_tdpc = 0
	is_measure_A = 0
    is_measure_chi = 0
    is_measure_vertex = 0
	n_write_meas = 10

	Fx = 0.d0
	Fy = 0.d0
	Fz = 0
    gauge_factor = 0.5d0

!!!!!
    U1 = 0.0
    U2 = 0.0
!!!!!

    open(12,file='dqmc_input.txt',status='old') !status='old' means don't overwrite
    read(12,nml=mc_param)
    read(12,nml=model_param)
    read(12,nml=meas_param)
    close(12)

!    if cont_run.eq.1, then use <run_name>_config.out. If it doesn't exist, use externally given config file.
!    then, if is_overwrite=0, results will be appended to old .out, .in files.
    if (is_cont_run.eq.1) then
        if (config_file_name.eq.'none') then
            inquire(file=trim(adjustl(run_name))//'_config.out', exist=file_exists)
            if (file_exists) then
                config_file_name = trim(adjustl(run_name))//'_config.out'
            endif
        endif
    endif

    ! Write input
    write(6,nml=mc_param)
    write(6,nml=model_param)
    write(6,nml=meas_param)
    if (is_overwrite.eq.1) then
        open(15,file=trim(adjustl(run_name))//'.in' ,status='replace')
    else
        open(15,file=trim(adjustl(run_name))//'.in' ,status='unknown', access='append')
    endif
    write(15,nml=mc_param)
    write(15,nml=model_param)
    write(15,nml=meas_param)
    close(15)

    call assert_L(Lx,Ly,Ntau,Ntau_mul)

	mc%Lx= Lx
	mc%Ly= Ly
    mc%Ntau = Ntau
    mc%Ntau_mul = Ntau_mul
    mc%seed = seed
    mc%max_mc_it = max_mc_it
    mc%run_name = run_name
    mc%n_write_config = n_write_config
    mc%n_write_out = n_write_out
    mc%num_global_update = num_global_update
    mc%purpose = purpose
    mc%is_run_fermions = is_run_fermions
    mc%n_eq=n_eq
    mc%opt%num_check_global= num_check_global
	mc%opt%min_cluster_size= min_cluster_size
    mc%dtau = beta/dble(mc%Ntau)

    mc%num_kin_mat = 4

    if(is_local_int.eq.1) then
        alpha = 2.0 * sqrt(mc%par%u2 + 4.d0 * mc%par%u1 )
    endif


    call allocate_mc_config(mc)

    mc%par%thop = thop
    mc%par%alpha = alpha
    mc%par%V = V
    mc%par%h = h
    mc%par%U1 = U1
    mc%par%U2 = U2
    mc%par%Vt = Vt
    mc%par%ht = ht
    mc%par%dVt = dVt
    mc%par%dht = dht
    mc%par%hz = hz
    mc%par%mu = mu
!~     mc%bcx = bcx
!~     mc%bcy = bcy
	mc%Fx = Fx
    mc%Fy = Fy
    mc%Fz = Fz
    mc%gauge_factor = gauge_factor



    mc%K_cl=0.d0
    if (mc%par%h.ne.0) then
        !mc%K_cl=h*mc%dtau
        mc%K_cl=-0.5d0*log(tanh(mc%par%h*mc%dtau))
    end if
    mc%K_clt=0.d0
    if (mc%par%ht.ne.0) then
        !mc%K_clt=h*mc%dtau
        mc%K_clt=-0.5d0*log(tanh(abs(mc%par%ht*mc%dtau)))
    end if

	mc%ps=1.0d0-exp(-2.0d0*mc%par%Vt*mc%dtau) !Set the probabilities of adding spins to the cluster according to the fictitious couplings Vt and ht
	mc%pt=1.0d0-exp(-2.0d0*mc%K_clt)
    !the array 'weights' calculates up front the appropriate Exp (delta S) value for flipping a spin
    !depending on the number of like nearest neighbor spins in space (first index) and time (second index) and the field
    ! The third array index takes values 1 or 2, indicating down or up spin
	write(6,*) 'space coupling', mc%par%V*mc%dtau, 'time coupling', mc%K_cl
    do j=1,5
        do k=1,3
            mc%weights(j,k,1)=exp(-dble(4*(k-2))*mc%K_cl)*exp(-dble(4*(j-3))*mc%par%V*mc%dtau)*exp(2.0d0*mc%par%hz*mc%dtau)
            mc%weights(j,k,2)=exp(-dble(4*(k-2))*mc%K_cl)*exp(-dble(4*(j-3))*mc%par%V*mc%dtau)*exp(-2.0d0*mc%par%hz*mc%dtau)
    !        write(6,*) mc%weights(j,k,1),mc%weights(j,k,2)

        end do
    end do


!~     if ((mc%bcx.ne.'periodic').and.(mc%bcx.ne.'antiperiodic')) stop '****  bcx must be periodic or antiperiodic  *******'
!~     if ((mc%bcy.ne.'periodic').and.(mc%bcy.ne.'antiperiodic')) stop '****  bcy must be periodic or antiperiodic  *******'
!~

	mc%measflags%is_measure_G = is_measure_G
	mc%measflags%is_measure_Gk = is_measure_Gk
	mc%measflags%is_measure_ETCF = is_measure_ETCF
	mc%measflags%is_measure_ETPC = is_measure_ETPC
	mc%measflags%is_measure_tdgf = is_measure_tdgf
	mc%measflags%is_measure_tdgf_K = is_measure_tdgf_K
	mc%measflags%is_measure_tdcf = is_measure_tdcf
	mc%measflags%is_measure_cccf = is_measure_cccf
    mc%measflags%is_measure_chi = is_measure_chi
	mc%measflags%is_measure_tdpc = is_measure_tdpc
	mc%measflags%is_measure_A = is_measure_A
	mc%measflags%is_overwrite= is_overwrite
	mc%measflags%is_measure_E= is_measure_E
    mc%measflags%is_measure_vertex= is_measure_vertex
	mc%measflags%n_write_meas = n_write_meas

    !!! open files
    if (purpose.eq.'dqmc') then
        if (is_overwrite.eq.1) then
            open(13,file=trim(adjustl(mc%run_name))//'_cluster.out', status='replace',buffered='yes')
            open(14,file=trim(adjustl(mc%run_name))//'.out',    status='replace',buffered='yes')
            open(16,file=trim(adjustl(mc%run_name))//'_eta.out', status='replace', access='stream', form='unformatted')
        else
            open(13,file=trim(adjustl(mc%run_name))//'_cluster.out', status='unknown', access='append',buffered='yes')
            open(14,file=trim(adjustl(mc%run_name))//'.out',status='unknown', access='append',buffered='yes')
            open(16,file=trim(adjustl(mc%run_name))//'_eta.out',status='unknown', buffered='yes', access='stream', form='unformatted', position='append')
        endif
    elseif (purpose.eq.'measure') then
        open(14,file=trim(adjustl(mc%run_name))//'.out',status='unknown', access='append',buffered='yes')
        open(16,file=trim(adjustl(mc%run_name))//'_eta.out',status='old', buffered='yes', access='stream', form='unformatted', position='rewind')
    endif

    !!! write header for run_name.out file:
    	write(14,*) '# iteration, bosonic action, avg fm , ferro susc.,'
    	write(14,*) '#xy bosonic action, tau bosonic action, signed magnetization'

    ! Initialize kinetic energy matrices
    call allocate_kin_mat(1, Lx*Ly/2, mc)
    call allocate_kin_mat(2, Lx*Ly/2, mc)
    call allocate_kin_mat(3, Lx*Ly/2, mc)
    call allocate_kin_mat(4, Lx*Ly/2, mc)

    ! Initialize nematic OP NEWSL
    call allocate_o_p(1, Lx*Ly/2, Ntau, mc)
    call allocate_o_p(2, Lx*Ly/2, Ntau, mc)
    call allocate_o_p(3, Lx*Ly/2, Ntau, mc)
    call allocate_o_p(4, Lx*Ly/2, Ntau, mc)

    call assign_kin_mat(1, HORIZONTAL, 0, thop, mc)
    call assign_kin_mat(2, HORIZONTAL, 1, thop, mc)
    call assign_kin_mat(3, VERTICAL  , 0, thop, mc)
    call assign_kin_mat(4, VERTICAL  , 1, thop, mc)


    call assign_o_p(1, mc)
    call assign_o_p(2, mc)
    call assign_o_p(3, mc)
    call assign_o_p(4, mc)


    !!! adding an action matrix for the bosons (for a local interaction)
    if (is_local_int.eq.1) then
        call assign_bosonic_action(mc)
    endif
    !!! set the initial dphi
    mc%delta = 0.1
    mc%opt%num_check = 10
    mc%opt%num_adjust_mu = 10

    !!!!! Seed Random Number Generator
    call srand(mc%seed)


    if (is_local_int.eq.1) then
        write(6,*) 'running local interaction'
    else
        write(6,*) 'using the usual, transverse field ising interaction'
    endif

    if (is_complex.eq.1) then
        write(6,*) 'using complex numbers'
    else
        write(6,*) 'using real numbers'
    endif


    !!!! Initialize eta NEWSL
    if ((is_cont_run.eq.1).and.(.not.(config_file_name.eq.'none'))) then
        ! Read initial configuration from file
        open(15,file=config_file_name,status='old')
        write(6,*) 'Reading initial configuration from: ', config_file_name
 		do jmat=1,mc%num_kin_mat
			do jt=1,mc%Ntau
			    do j=1,mc%o_p(jmat)%num_pairs
	            	read(15,*) mc%o_p(jmat)%eta(j,jt)
	            end do
            end do
        end do
        close(15)
    else
	! Random initialization of nematic OP NEWSL
    	do jmat=1,4
		    do jt = 1,mc%Ntau
		        do j = 1, Lx*Ly/2
			if(is_rand_initial.eq.0) then
                mc%o_p(jmat)%eta(j,jt)=1
            else
                mc%o_p(jmat)%eta(j,jt)=2*int(2.0d0*rand())-1
			endif
                end do
	        end do
    	end do
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Test Green's function    !!!!!!
!    call center_G(mc)
!    call savemat(mc%lat_dim, mc%G%m, 'G.dat')
!    stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! tests !!!!!!!!!!!!!!!!!!!!!!
!	call test_rdet(mc)
!    call test_rotateG(mc)
!    call test_global_update(mc)
!   call test_update_G(mc)
!	call test_loop_rotation(mc)
!	call test_operateT(mc)
!    call test_reconstructG(mc)
!    call test_S_phi(mc)

	!mc%is_verbose = 1
end subroutine init_par


!!! copies the parameters needed for allocating a copy of a mc config:
subroutine copy_param(mc_new,mc)
	use global
	implicit none
	type(mcconfig) mc, mc_new
	mc_new%Ntau = mc%Ntau
    mc_new%Ntau_mul = mc%Ntau_mul
    mc_new%Lx = mc%Lx
    mc_new%Ly = mc%Ly
    mc_new%num_kin_mat = mc%num_kin_mat
end subroutine copy_param
