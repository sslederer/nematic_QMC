program dqmc
    use global !module has global variables
    implicit none

    type(McConfig) mc ! defines a structure in global of type McConfig; almost everything is in there; defined in module_global
    real clock_start, clock_finish
    integer(8) sclock_start,sclock_finish,rate
    call cpu_time(clock_start)
    call system_clock(sclock_start,rate)
    write(6,*) 'Begin DQMC' ! output to terminal
    ! write version number
	write(6,*) 'Running version 18'
    call random_seed !seed is input to the program; start various instances by choosing appropriate seed
    call init_par(mc) !this is within mc_config.f90

    if (mc%purpose.eq.'dqmc') then
        call qmc_loop(mc)
    elseif(mc%purpose.eq.'measure') then
        call compute_means(mc)
    else
        stop '**** Variable purpose must take one of the values: dqmc, measure ****'
    endif

    call deallocate_mc_config(mc) !important to do this or else comp will eventually crash

    !!! close files
    close(14); close(16); close(17);close(13)



    write(6,*) 'End DQMC'
    call cpu_time(clock_finish)
    call system_clock(sclock_finish)
!    write(6,*) 'OperateT_left  took:', mc%prof%t(1), 'sec'
!    write(6,*) 'OperateT_right took:', mc%prof%t(2), 'sec'
!    write(6,*) 'reconstruct_G(left)  took:', mc%prof%t(3), 'sec'
!    write(6,*) 'reconstruct_G(right) took:', mc%prof%t(4), 'sec'
!    write(6,*) 'update_G(left)  took:', mc%prof%t(5), 'sec'
!    write(6,*) 'update_G(right) took:', mc%prof%t(6), 'sec'
!    write(6,*) 'AdvanceSVD(left)  took:', mc%prof%t(7), 'sec'
!    write(6,*) 'AdvanceSVD(right) took:', mc%prof%t(8), 'sec'

!    write(6,*) 'Total CPU time:', (clock_finish-clock_start)/60.d0, 'min'

    write(6,*) 'OperateT_left  took:', mc%prof_s%t(1), 'sec'
    write(6,*) 'OperateT_right took:', mc%prof_s%t(2), 'sec'
    write(6,*) 'reconstruct_G(left)  took:', mc%prof_s%t(3), 'sec'
    write(6,*) 'reconstruct_G(right) took:', mc%prof_s%t(4), 'sec'
    write(6,*) 'update_G(left)  took:', mc%prof_s%t(5), 'sec'
    write(6,*) 'update_G(right) took:', mc%prof_s%t(6), 'sec'
    write(6,*) 'Each local spacetime sweep took:', mc%prof_s%t(9)*dble(2*mc%Lx*mc%Ly*mc%Ntau)/dble(mc%opt%num_trial), 'sec'
    write(6,*) 'Each cluster update took', mc%prof_s%t(10)/mc%opt%num_trial_global, 'sec'

    write(6,*) 'Total time:', dble(sclock_finish-sclock_start)/(dble(rate)*60.d0), 'min'

    write(6,*) 'Tried',mc%opt%num_trial,', accepted',dble(mc%opt%num_accept)/dble(mc%opt%num_trial)
    write(6,*) 'Tried',mc%opt%num_trial_global,'global, accepted',dble(mc%opt%num_accept_global)/dble(mc%opt%num_trial_global)
    write(6,*) 'Spins flipped per global update', dble(mc%opt%size_global_tot)/dble(mc%opt%num_trial_global)
    write(6,*) 'Final ht=',mc%par%ht,', final Vt=',mc%par%Vt
end program dqmc

subroutine profiling(message)
    character(*) :: message
    real(4) T1, TA(2)

    T1=DTIME(TA)
    write(6,*) message, ' took: ', T1, 'sec'
    write(6,*) ''
end subroutine


!!!! Main QMC loop
subroutine qmc_loop(mc)
    use global
    implicit none
    type(McConfig) mc
    integer(4) j, dir, nt, jmat,k, is_accept
    real(8) fmr,fmsr
    character(*), parameter :: message = 'QMC iteration'

	fmr=0.d0 !running averages of the fm order parameter, and the square thereof.
	fmsr= 0.0d0

	mc%opt%size_global=0

    dir = LEFT
	if (mc%is_run_fermions.eq.1) then
		call initG(mc)
	endif


	do j=1,mc%max_mc_it
		write(6,*) 'starting iteration', j
		if (mc%is_run_fermions.ne.1) then
			if(mod(j,2).eq.0) then
				do nt=1,mc%Ntau
					mc%n=mod(nt-2+mc%Ntau,mc%Ntau)+1
					call sweep_timeslice(mc,RIGHT)
				end do
			else
				do nt=1,mc%Ntau
					mc%n=mc%Ntau+1-nt
					call sweep_timeslice(mc,LEFT)
				end do
			end if
		else
			do nt=1,(mc%Ntau)
				if (mod(mc%n+dir,mc%Ntau_mul).ne.0) then
					call sweep_timeslice(mc,dir)
                    if (j.le.mc%n_eq) call adjust_mu_eta(mc)
				else !! This can be optimized by incorporating the rotateG into stab_reconstructG.
					call sweep_timeslice(mc,dir)
                    if (j.le.mc%n_eq) call adjust_mu_eta(mc)
					call rotateG(dir,mc)
					call stab_reconstructG(-dir,mc)
				endif
			end do
		endif
		if(mc%is_run_fermions.ne.1) then
			if (mc%num_global_update.eq.0) then
				call cluster_update(mc)
			elseif (mod(j,mc%num_global_update).eq.0) then
				call cluster_update(mc)
			endif
		else !!! When the fermions are around, global updates should only happen for j=odd.
			if ((mod(j,2).eq.0).and.(mc%num_global_update.ne.0)) then
				do k=1,mc%num_global_update
					call cluster_update(mc)
				end do
			endif
		endif
    	call write_out(j,mc,fmr,fmsr)
    	call profiling(message)

	    dir = -dir

   	end do
!	write(6,*) 'fm',fmr/dble(2*mc%Lx*mc%Ly*mc%Ntau),'chi',mc%dtau*(fmsr-fmr**2)/dble(2*mc%Lx*mc%Ly*mc%Ntau)
end subroutine qmc_loop


subroutine sweep_timeslice(mc, dir)
    use global
    use random_numbers
    implicit none
    type(McConfig) mc
    integer(4) nt, jmat, j, is_accept, dir
    integer(4) jmat_i,jmat_f, jmat_dir,is_even, is_inv
    integer(8) sclock_start,sclock_finish,rate
    real(8) dSe, action_ratio,W ,rdet
    $op_type$ eta_new, eta_old

	call system_clock(sclock_start,rate)
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
		if(mc%is_verbose.eq.1) then
			write(6,*) 'sweeping timeslice, nt, jmat'
			write(6,*) nt,jmat
		endif
		do j=1,mc%Lx*mc%Ly/2
            eta_old = mc%o_p(jmat)%eta(j,nt)
            if(is_local_int.eq.1) then
                eta_new = eta_old + (rand()-0.5d0) * mc%delta
                call delta_S_eta_new(dSe, j, nt, mc, jmat, eta_new)
                action_ratio = exp(-dSe)
                !
                ! call calc_S_eta_new(temp1, mc)
                ! mc%o_p(jmat)%eta(j,nt) = eta_new
                ! call calc_S_eta_new(temp2, mc)
                ! mc%o_p(jmat)%eta(j,nt) = eta_old
                ! write(6,*) 'action difference', dSe / (temp2 - temp1)

            else 		!proposed move is flipping the spin on the given link
                eta_new = -eta_old
                call weight_flip(action_ratio,j,nt,mc,jmat)
            endif

			if (mc%is_run_fermions.eq.1) then
				call det_ratio(rdet, mc, j, jmat, eta_old, eta_new)
				action_ratio=action_ratio*rdet
			endif

			!!decide whether to accept
			is_accept=1;
			if(action_ratio.lt.1.d0) then
				if (rand().gt.action_ratio) is_accept=0
			endif
			if (is_accept.eq.1) then
				if(mc%is_run_fermions.eq.1) then
					call update_G(dir, mc, j, jmat, eta_old, eta_new)
					mc%det_change = mc%det_change*rdet
				end if
				mc%o_p(jmat)%eta(j,nt)=eta_new
				mc%opt%num_accept=mc%opt%num_accept+1
			end if
			mc%opt%num_trial=mc%opt%num_trial+1
		end do
		! rotate the current kinetic energy matrix to the other side
		if(mc%is_run_fermions.eq.1) then
			if (mc%is_verbose.eq.1) then
				write(6,*) 'finished timeslice, rotating the current kinetic energy matrix to other side'
			endif
			call OperateT(nt, LEFT, jmat, is_inv, mc%dtau, mc%G, mc)
			call OperateT(nt, RIGHT, jmat, -is_inv, mc%dtau, mc%G, mc)
		end if
	end do
	!now change mc%n:
	if(dir.eq.RIGHT) mc%n = mod(mc%n,mc%Ntau)+1
	if(dir.eq.LEFT)  mc%n = mod(mc%n-2+mc%Ntau, mc%Ntau)+1

	call system_clock(sclock_finish)
	mc%prof_s%t(9)=mc%prof_s%t(9)+dble(sclock_finish-sclock_start)/dble(rate)
end subroutine sweep_timeslice


! a global update which flips a cluster of like spins using the Wolff algorithm
!A cluster is formed probabilistically according to the spin part of the action, and the cluster is flipped according
!To the change in the fermion component

subroutine cluster_update(mc)
    use global
    use random_numbers
    implicit none
    type(McConfig) mc,mc_new
    integer(4) k, nt,j, n,initial_spin, cluster_size,is_accepted,m, attempts
    real(8) action_ratio, H0s,Hfs,H0t,Hft,det_ratio, bare_boson_ratio
    integer(8) sclock_start, sclock_finish,rate
    character(*), parameter :: message = 'Populating cluster'

    if(is_local_int.eq.1) then
        write(6,*) 'global updates have not been implemented for the local interaction model.'
        return
    endif
    call system_clock(sclock_start,rate)

!! Create another copy of mc, on which a cluster will be built and proposed:
    call copy_param(mc_new,mc)
    call allocate_mc_config(mc_new)
    mc_new = mc


!! Calculate the initial bosonic action:
	if((mc%is_run_fermions.eq.1).or.(mc%par%hz.ne.0)) then
		call calc_S_eta(H0s,mc_new,1)
		call calc_S_eta(H0t,mc_new,2)
	end if


!! Randomly pick a seed for the cluster:
	k=int(4.d0*rand())+1
	nt= int(dble(mc%Ntau)*rand())+1
	j= int(dble(mc%o_p(k)%num_pairs)*rand())+1
	cluster_size=1
	initial_spin=mc_new%o_p(k)%eta(j,nt)

!! Populate the cluster
	call expand_cluster(mc_new,k,nt,j,initial_spin,cluster_size)

!! If the cluster isn't big enough keep trying, up to 10000 times
	attempts=1
	do while ((cluster_size.lt.mc%opt%min_cluster_size).and.(attempts.le.10000))
		k=int(4.d0*rand())+1
		nt= int(dble(mc%Ntau)*rand())+1
		j= int(dble(mc%o_p(k)%num_pairs)*rand())+1
		cluster_size=1
		mc_new=mc
		initial_spin=mc_new%o_p(k)%eta(j,nt)
		call expand_cluster(mc_new,k,nt,j,initial_spin,cluster_size)
		attempts=attempts+1
	end do

	!! default values for the case is_run_fermions=0 and hz=0
	action_ratio=1.d0
	bare_boson_ratio=1.d0
	det_ratio=1.d0
	call profiling(message)

	is_accepted=0
	if(cluster_size.ge.mc%opt%min_cluster_size) then
		mc_new%opt%num_trial_global=mc_new%opt%num_trial_global+1
		mc%opt%num_trial_global=mc%opt%num_trial_global+1
	!! Calculate bosonic part of the action after the flip:
		if((mc%is_run_fermions.eq.1).or.(mc%par%hz.ne.0)) then
			call calc_S_eta(Hfs,mc_new,1)
			call calc_S_eta(Hft,mc_new,2)
			action_ratio=exp(-2*mc%par%hz*dble(initial_spin*cluster_size)*mc%dtau-(Hfs-H0s)*(mc%par%V-mc%par%Vt)/mc%par%V-(Hft-H0t)*(mc%K_cl-mc%K_clt)/mc%K_cl)
			bare_boson_ratio=exp(-2*mc%par%hz*dble(initial_spin*cluster_size)*mc%dtau-(Hfs-H0s)-(Hft-H0t))

		!! Calculate fermionic Green's function, determinant ratio:

			if (mc%is_run_fermions.eq.1) then
				det_ratio=1.d0/mc%det_change
				call initG(mc_new)
				do j = 1,mc%lat_dim
					det_ratio = det_ratio*(mc_new%S(j)/mc%S(j))**2
				end do
			endif
			action_ratio=action_ratio*det_ratio
			if (action_ratio.ge.rand()) is_accepted=1

		else
			is_accepted=1
		end if

		write(13,*) mc%opt%num_trial_global
		write(13,*) cluster_size
		write(13,*) action_ratio/det_ratio
		write(13,*) bare_boson_ratio
		write(13,*) det_ratio
		write(13,*) action_ratio
		write(13,*) is_accepted

	else
		write(6,*) 'cluster too small'
	end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!Profiling the cluster updates:
	m=mc%opt%num_trial_global/mc%opt%num_check_global
	if (mod(m,2).eq.0) then
		mc_new%opt%size_global(2)=mc_new%opt%size_global(2)+cluster_size
		mc_new%opt%size_global_tot=mc_new%opt%size_global_tot+cluster_size
	else
		mc_new%opt%size_global(4)=mc_new%opt%size_global(4)+cluster_size
		mc_new%opt%size_global_tot=mc_new%opt%size_global_tot+cluster_size
	end if


	if (is_accepted.eq.1) then
		mc = mc_new
		mc%opt%num_accept_global=mc%opt%num_accept_global+1
	end if




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Attempt to rescale Ht, Vt
	!!! I think this should be placed here... but it doesn't really matter. Yoni
	call change_ht_Vt(mc)

    call deallocate_mc_config(mc_new)
    call system_clock(sclock_finish)
    mc%prof_s%t(10)=mc%prof_s%t(10)+dble(sclock_finish-sclock_start)/dble(rate)
end subroutine cluster_update

subroutine change_ht_Vt(mc)
	use global
	implicit none
	type(mcconfig) mc
	integer(4) m

    m=mc%opt%num_trial_global/mc%opt%num_check_global
    !adjust fictitious ising couplings every num_check_global global updates
    if(mod(mc%opt%num_trial_global,mc%opt%num_check_global).eq.0) then
		if (mod(m,2).eq.0) then !adjust Vt every other check
			write(6,*) 'last time flipped', dble(mc%opt%size_global(2))/dble(mc%opt%num_check_global),'per trial'
			if(mc%opt%size_global(2).lt.mc%opt%size_global(1)) mc%par%dVt=-mc%par%dVt
			write(6,*) 'changing Vt from', mc%par%Vt,'to',mc%par%Vt+mc%par%dVt
			mc%par%Vt=mc%par%Vt+mc%par%dVt
			mc%opt%size_global(1)=mc%opt%size_global(2) !
			mc%opt%size_global(2)=0
		else		!adjust ht every other check
			write(6,*) 'last time flipped', dble(mc%opt%size_global(4))/dble(mc%opt%num_check_global),'per trial'
			if(mc%opt%size_global(4).lt.mc%opt%size_global(3)) mc%par%dht=-mc%par%dht
			write(6,*) 'changing ht from', mc%par%ht,'to',mc%par%ht+mc%par%dht
			mc%par%ht=mc%par%ht+mc%par%dht
			mc%opt%size_global(3)=mc%opt%size_global(4)
			mc%opt%size_global(4)=0
		end if
		!recalculate the fictitious probabilities with the new values of ht and Vt
		mc%K_clt=-0.5d0*log(tanh(mc%par%ht*mc%dtau))
		write(6,*) 'Old ps',mc%ps,', oldpt',mc%pt
		mc%ps=1.0d0-exp(-2.0d0*mc%par%Vt*mc%dtau)
		mc%pt=1.0d0-exp(-2.0d0*mc%K_clt)
		write(6,*) 'new ps',mc%ps,',new pt',mc%pt
    end if

end subroutine change_ht_Vt


recursive subroutine expand_cluster(mc,k,nt,j,initial_spin,cluster_size) !NEWSL; this assumes the system prefers nematic spin ferromagnetism;
    use global
    use random_numbers
    implicit none
    type(McConfig) mc
    integer(4) k, nt,j, n,initial_spin, q1, q2, ntp,ntm,cluster_size

	ntp=mod(nt,mc%Ntau)+1
	ntm=mod(nt-2+mc%Ntau,mc%Ntau)+1

	mc%o_p(k)%eta(j,nt)=-initial_spin
	do n=1,4  !Try to expand the cluster along spatial neighbors
		q1=mc%o_p(k)%Q(j,n,1)
		q2=mc%o_p(k)%Q(j,n,2)
		if (mc%o_p(q1)%eta(q2,nt).eq.initial_spin) then
			if (rand().lt.mc%ps) then
				cluster_size=cluster_size+1
				call expand_cluster(mc,q1,nt,q2,initial_spin,cluster_size)
			end if
		end if
	end do

!Try to expand along time neighbors

	if (mc%o_p(k)%eta(j,ntp).eq.initial_spin) then
		if (rand().lt.mc%pt) then
			cluster_size=cluster_size+1
			call expand_cluster(mc,k,ntp,j,initial_spin,cluster_size)
		end if
	end if

	if (mc%o_p(k)%eta(j,ntm).eq.initial_spin) then
		if (rand().lt.mc%pt) then
			cluster_size=cluster_size+1
			call expand_cluster(mc,k,ntm,j,initial_spin,cluster_size)
		end if
	end if

end subroutine expand_cluster


!!!! write output configuration !!!!!
subroutine write_out(it, mc,fmr, fmsr)
    use global
    implicit none
    type(McConfig) mc
    integer(4) it, nt, jx,jy,jt,k,j, jmat,s1,s2,c
    real(8) fmr,fmsr,E,Ss,St, signedM, moment_size,density

    !!!! For each sweep, write running average of both fm and afm order parameter to <run_name>.out file
	call calc_S_eta(Ss,mc,1)
	call calc_S_eta(St,mc,2)

    signedM = 0.d0
    do jmat=1,4
        do jt=1,mc%Ntau
            do j=1,mc%o_p(jmat)%num_pairs
                signedM = signedM + mc%o_p(jmat)%eta(j,jt)
            enddo
        enddo
    enddo

	fmr= (fmr*dble(it-1)+dble(abs(signedM)))/dble(it)
	fmsr= (fmsr*dble(it-1)+dble(signedM)**2)/dble(it)

    moment_size = 0.d0
    do j=1,4
        moment_size = moment_size + sum(mc%o_p(j)%eta**2)
    enddo
    moment_size = sqrt(moment_size /dble(2 * mc%Lx * mc%Ly  * mc%ntau))

    density = 0.d0
    do j=1,mc%lat_dim
        density = density + 1.d0 - dble(mc%G%m(j,j))
    enddo
    density = density / mc%lat_dim

	write(14,*)  it 													!iteration
	write(14,*) E/dble(2*mc%Lx*mc%Ly*mc%Ntau)							!action
	write(14,*) fmr/dble(2*mc%Lx*mc%Ly*mc%Ntau) 						!'avg fm #'
	write(14,*) mc%dtau*(fmsr-fmr**2)/dble(2*mc%Lx*mc%Ly*mc%Ntau)		!ferro susceptibility
	write(14,*) Ss
	write(14,*) St
	write(14,*) dble(signedM)/dble(2*mc%Lx*mc%Ly*mc%Ntau)				!signed magnetization
    write(14,*) moment_size
    write(14,*) density



!!!Write eta configuration to <run_name>_config.out (overwrite previous configs) !!!
    if (mod(it,mc%n_write_config).eq.0) then
        open(12,file=trim(adjustl(mc%run_name))//'_config.out',status='replace',buffered='yes')
        call system_clock(s1,c)
		do jmat=1,mc%num_kin_mat
			do jt=1,mc%Ntau
			    do j=1,mc%o_p(jmat)%num_pairs
	            	write(12,*) mc%o_p(jmat)%eta(j,jt)
	            end do
            end do
        end do
        close(12)
    endif


    !!!! write eta configuration to <run_name>_eta.out file (saves previous configs) !!!!!
    if (mod(it,mc%n_write_out).eq.0) then
        call system_clock(s1,c)
		do jmat=1,mc%num_kin_mat
			do jt=1,mc%Ntau
			    do j=1,mc%o_p(jmat)%num_pairs
	            	write(16) mc%o_p(jmat)%eta(j,jt)
	            end do
            end do
        end do
    endif

end subroutine write_out


!!!!!Calculate eta field action from eta configuration NEWSL
subroutine calc_S_eta(Se,mc,a)
    use global
    implicit none
    type(McConfig) mc
    integer(4) :: j, jt, jmat,jtp,k,Q1,Q2,a
    real(8) :: Se
	Se=0.d0
	do jt=1,mc%Ntau
		jtp = mod(jt,mc%Ntau)+1
		do jmat=1,4
			do j=1,mc%o_p(jmat)%num_pairs

				if((a.eq.2).or.(a.eq.0)) then
	!time derivative term (transverse field) and symm. breaking field hz
					Se=Se-mc%par%hz*dble(mc%o_p(jmat)%eta(j,jt))*mc%dtau
					Se= Se - mc%K_cl*dble( mc%o_p(jmat)%eta(j,jt)*(mc%o_p(jmat)%eta(j,jtp)))
	!ising term. I have opted to make a positive V mean ferromagnetic ordering of the spins. Will need to incorporate a relative
	!minus sign between horizontal and vertical bonds in the coupling to fermions
				end if
				if((a.eq.1).or.(a.eq.0)) then
					do k=1,2
						Q1=mc%o_p(jmat)%Q(j,k,1)
						Q2=mc%o_p(jmat)%Q(j,k,2)

						Se= Se-mc%par%V*dble(mc%o_p(jmat)%eta(j,jt)*mc%o_p(Q1)%eta(Q2,jt))*mc%dtau
					end do
				end if
			end do
		end do
	end do

end subroutine calc_S_eta

!calculates the acceptance probability directly by lookup, avoiding recalculating exponentials
subroutine weight_flip(W, j, jt, mc,jmat)
    use global
    implicit none
    type(McConfig) mc
    integer(4) :: j, jt, jmat,jtp, jtm,k,Q1,Q2,nt,ns
    real(8) :: W

    jtp = mod(jt          ,mc%Ntau)+1
    jtm = mod(jt-2+mc%Ntau,mc%Ntau)+1
    nt=mc%o_p(jmat)%eta(j,jt)*(mc%o_p(jmat)%eta(j,jtm) + mc%o_p(jmat)%eta(j,jtp))/2+2
    ns=0
!   ising term. I have opted to make a positive V mean ferromagnetic ordering of the spins. Will need to incorporate a relative
!minus sign between horizontal and vertical bonds in the coupling to fermions
    do k=1,4
        Q1=mc%o_p(jmat)%Q(j,k,1)
        Q2=mc%o_p(jmat)%Q(j,k,2)

        ns=ns+mc%o_p(jmat)%eta(j,jt)*mc%o_p(Q1)%eta(Q2,jt)
    end do
ns=ns/2+3

W=mc%weights(ns,nt,(mc%o_p(jmat)%eta(j,jt)+1)/2+1)

end subroutine weight_flip

!!!!! Calculate change in S_eta due to change in eta NEWSL
subroutine delta_S_eta(dSe, j, jt, mc,jmat,a)
	use global
	implicit none
	type(McConfig) mc
	integer(4) :: j, jt, jmat,jtp, jtm,k,Q1,Q2,a
	real(8) :: dSe,Si,Sf

	dSe=0.d0

	jtp = mod(jt          ,mc%Ntau)+1
	jtm = mod(jt-2+mc%Ntau,mc%Ntau)+1

	if((a.eq.2) .or. (a.eq.0)) then
	!	time derivative term (transverse field) and symm. breaking field hz
		dSe=dSe + 2.d0*mc%par%hz*dble(mc%o_p(jmat)%eta(j,jt))*mc%dtau
		dSe= dSe + 2.d0*mc%K_cl*dble(mc%o_p(jmat)%eta(j,jt)* ( mc%o_p(jmat)%eta(j,jtm) + mc%o_p(jmat)%eta(j,jtp) ) )
	end if

	if((a.eq.1) .or. (a.eq.0)) then
!	ising term
	do k=1,4
		Q1=mc%o_p(jmat)%Q(j,k,1)
		Q2=mc%o_p(jmat)%Q(j,k,2)

		dSe= dSe + 2.d0*mc%par%V*dble(mc%o_p(jmat)%eta(j,jt)*mc%o_p(Q1)%eta(Q2,jt))*mc%dtau
	end do
	end if

end subroutine delta_S_eta

subroutine delta_S_eta_new(dSe, jp, jt, mc, jmat, eta_new)
    use global
    implicit none
    type(McConfig) mc
    integer(4) :: jt, jp, jp2, jmat, jmat2
    real(8) :: dSe
    $op_type$:: eta_new, delta

    dSe=0.d0
    delta = eta_new - mc%o_p(jmat)%eta(jp, jt)
    do jmat2 =1,4
        do jp2 = 1, mc%kin_mat(jmat2)%num_pairs
            dSe = dSe + 2.0 * mc%dtau * delta * mc%S0(jp, jmat, jp2, jmat2) * mc%o_p(jmat2)%eta(jp2,jt)
        enddo
    enddo
    dSe = dSe + mc%dtau * mc%S0(jp,jmat,jp,jmat) * delta **2

end subroutine delta_S_eta_new

subroutine calc_S_eta_new(S, mc)
    use global
    implicit none
    type(McConfig) mc
    integer(4) :: jt, jmat1, jmat2, jp1, jp2
    real(8) :: S

    S=0.d0
    do jt = 1,mc%Ntau
        do jmat1=1,4; do jp1 = 1,mc%o_p(jmat1)%num_pairs
            do jmat2 = 1,4; do jp2 = 1,mc%o_p(jmat2)%num_pairs
                S = S +  mc%dtau * mc%o_p(jmat1)%eta(jp1,jt) * mc%S0(jp1, jmat1, jp2, jmat2) * mc%o_p(jmat2)%eta(jp2,jt)
            enddo; enddo
        enddo; enddo
    enddo
end subroutine calc_S_eta_new

subroutine adjust_mu_eta(mc)
    use global
    implicit none
    type(mcconfig) mc
    real(8) ratio, rho, dmu
    real(8) delta_max, delta_min
    integer(4) j

    delta_min = 0.000001
    delta_max = 5.0
    ! delta_max = 10.0
    if(is_local_int.eq.1) then
        if (mc%opt%num_trial.gt.mc%opt%num_check) then
            ratio = dble(mc%opt%num_accept)/dble(mc%opt%num_trial)
            if (ratio.gt.0.6d0) then
                mc%delta = min(mc%delta*1.05d0, delta_max)
    			write(6,'("Curent acceptance ratio ",(F6.3),". Increasing dphi to ", (F))'), ratio, mc%delta
            elseif (ratio.lt.0.4d0) then
                mc%delta = max(mc%delta/1.05d0, delta_min)
    			write(6,'("Curent acceptance ratio ",(F6.3),". Decreasing dphi to ", (F))'), ratio, mc%delta
            endif
            mc%opt%num_trial = 0; mc%opt%num_accept = 0
        endif
    endif

!!! todo: add adjustment of mu to fix the density
    ! dmu = 0.00
    ! rho = 0.d0
    ! mc%par%target_rho = 0.6
    ! do j=1,mc%lat_dim
    !     rho = rho + 1.0 - dble(mc%G%m(j,j))
    ! enddo
    ! rho = rho / mc%lat_dim
    !
    ! mc%opt%rho = mc%opt%rho + rho
    ! mc%opt%num_add_rho = mc%opt%num_add_rho +1
    !
    ! if(mc%opt%num_add_rho.ge.mc%opt%num_adjust_mu) then
    !     rho = mc%opt%rho / mc%opt%num_add_rho
    !     if(rho.gt.mc%par%target_rho) then
    !         write(6,*) 'rho=', rho, 'target_rho=', mc%par%target_rho, 'adjusting mu=', mc%par%mu
    !         mc%par%mu = mc%par%mu -dmu
    !     elseif(rho.lt.mc%par%target_rho) then
    !         mc%par%mu = mc%par%mu + dmu
    !     endif
    !     mc%opt%rho = 0.d0
    !     mc%opt%num_add_rho = 0
    ! endif
    !
end subroutine adjust_mu_eta
