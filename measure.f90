subroutine compute_means(mc)
    use global
    implicit none
    type(McConfig) mc
    type(McMeasure) meas
	integer (4)  nt, dir, readstat

	call initmeas(meas,mc)

	if (meas%flags%is_measure_A.eq.1) call boson_meas(mc)

	if (meas%flags%is_measure_fermions.eq.1) then
		!!! Read phi file until end of file
		do
			call readeta(mc,readstat)
			if (readstat /= 0) exit

			meas%n_meas = meas%n_meas + 1
			write(6,*) 'Measurement number', meas%n_meas

		    !!! Construct Green's function for all time slices.
			call get_equal_time_G(meas%G,meas%Gh,mc)

			!!time displaced G measurement
			if (meas%flags%is_measure_tdgf.eq.1) then
				call scalettar_tdgf_both_sides(meas%tdgf,mc)
				call add_mean_tdgf(meas%tdgf,meas%mean_tdgf,mc)
                ! call get_tdgf_k1k2(meas,mc)
                ! call calc_TDPC2(meas, mc)
			end if




			!!time displaced correlation function (must be after tdgf)
		    if (meas%flags%is_measure_tdcf.eq.1) then
		    	call add_tdcf(meas,mc)
		    end if

		    !! time displaced current-current correlation function (must be after tdgf)
		    if (meas%flags%is_measure_cccf.eq.1) then
		    	call calc_cccf(meas,mc)
		    end if

		    !! time displaced pairing correlations (must be after tdgf)
		    if (meas%flags%is_measure_tdpc.eq.1) then
		    	call add_tdpc(meas,mc)
		    end if

            !! fermionic nematic correlations (must be after tdgf)
            if (meas%flags%is_measure_chi.eq.1) then
                call calc_nematic_corr(meas,mc)
            end if

			!!equal time measurements

			if(meas%flags%is_measure_E.eq.1) then
				call add_E(meas,mc)
			end if

			do nt=1,(mc%Ntau)
				!!! Measure, Update mean measured values
				mc%G%m=meas%G(nt)%m
				if ((meas%flags%is_measure_G.eq.1).or.(meas%flags%is_measure_Gk.eq.1)) then
					 meas%meanG%m =  meas%meanG%m + mc%G%m
				end if
				if (meas%flags%is_measure_ETCF.eq.1) then
					call add_ETCF(meas,mc)
				end if
				if (meas%flags%is_measure_ETPC.eq.1) then
					call add_ETPC(meas,mc)
				end if

			end do

			call profiling('Equal time measurements')

			if(mod(meas%n_meas, meas%n_write_meas).eq.0) then
				call write_measure(meas, mc)
			end if
		end do
	end if
	call end_measure(meas,mc)

end subroutine compute_means

!!!!! Write mean Green's function (averaged over all time slices)
subroutine write_mean_G(meanG,mc)
    use global
    implicit none
    $G_type$ meanG
    type(mcconfig) mc
    integer(4) j1, j2

    do j2 = 1, meanG%d
        do j1 = 1, meanG%d
            write(17) dble(meanG%m(j1, j2))
        end do
    end do
end subroutine write_mean_G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Transform Green's function to k-space !!!!!!!!!
subroutine transform_Gk(G, mc, Gk)
    use global
    implicit none
    $G_type$ G
    type(McConfig) mc
    complex(8), allocatable :: Gxy(:,:), Vk(:), Kx(:), Ky(:)
    complex(8) Gk(mc%Lx, mc%Ly), coef, theta1, theta2
    integer(4) jx, jy, ix1, ix2, iy1, iy2, ix, iy, Nkx, Nky, j1, j2

    Nkx = mc%Lx; Nky = mc%Ly
    allocate(Kx(Nkx), Ky(Nky), Gxy(mc%Lx*mc%Ly,mc%Lx*mc%Ly), Vk(mc%Lx*mc%Ly))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! define Kx, Ky, taking care of different boundary conditions !!!!!!
    do jx = 1, Nkx
        Kx(jx) = cmplx(dble(jx-1 + mc%Fx)*(2.d0*PI)/dble(mc%Lx))
    end do
    do jy = 1, Nky
        Ky(jy) = cmplx(dble(jy-1 + mc%Fy)*(2.d0*PI)/dble(mc%Ly))
    end do

    Gxy = (0.d0, 0.d0)
    do ix1 = 1, mc%Lx; do iy1 = 1, mc%Ly
        do ix2 = 1, mc%Lx; do iy2 = 1, mc%Ly
            Gxy(ix1+(iy1-1)*mc%Lx, ix2+(iy2-1)*mc%Lx) = G%m(mc%i(ix1,iy1), mc%i(ix2,iy2))
        end do; end do
	end do; end do

	do jy = 1, Nkx; do jx = 1, Nky
		do ix = 1, mc%Lx; do iy = 1, mc%Ly
			Vk(ix+(iy-1)*mc%Lx) = exp((0.d0,1.d0)*(cmplx(ix)*Kx(jx)+cmplx(iy)*Ky(jy)) )
		end do; end do
        Gk(jx, jy) = dot_product(Vk, matmul(Gxy, Vk))
    end do; end do

    Gk = Gk/cmplx(mc%Lx*mc%Ly)  !!! Normalization
    deallocate(Kx, Ky, Gxy, Vk)
end subroutine transform_Gk


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Transform Green's function to k-space, and write to file fileid
subroutine write_Gk(G, mc,fileid)
    use global
    implicit none
    $G_type$ G
    type(McConfig) mc
    complex(8), allocatable :: Gk(:,:)
    integer(4) jx, jy
	integer(4) fileid

    allocate(Gk(mc%Lx,mc%Ly))
    call transform_Gk(G, mc, Gk)

    do jy = 1, mc%Ly; do jx = 1, mc%Lx
        write(fileid) dble(Gk(jx, jy))
        write(fileid) aimag(Gk(jx, jy))
    end do; end do
    deallocate(Gk)
end subroutine write_Gk


subroutine readeta(mc,readstat)
	use global
	implicit none
	type(McConfig) mc
	integer (4) readstat
    integer(4) jt,k,j
    $op_type$ temp_eta

	readstat=0

	do k=1,mc%num_kin_mat
		do jt=1,mc%Ntau
			do j=1,mc%o_p(k)%num_pairs
		    	read(16, iostat=readstat) temp_eta
		    	mc%o_p(k)%eta(j,jt) = temp_eta
		    end do
    	end do
	end do

end subroutine

!! Computes the equal time, local correlation function ETCF(s4,s3,s2,s1,d34,d12,deltar)=\sum_r <c^dag(r,s1)c(r+d12,s2)c^dag(r+dr,s3)c(r+dr+d34,s4)> for the current timeslice.
subroutine add_ETCF(meas,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(mcMeasure) meas
	$thop_type$ :: G23,G41,G21,G43
	! $thop_type$ :: dETCF
	integer(4) :: d(5,2)
	integer(4) :: x1,y1,d12,d34,j1,j2,j3,j4
	integer(4) :: deltax,deltay
	integer(4) :: s1,s2,s3,s4

	d=meas%d

	do deltax=0,(mc%Lx-1); do deltay=0,(mc%Ly-1)
		do x1=1,mc%Lx; do y1=1,mc%Ly
			do d12=1,5; do d34=1,5
				call get_r_index(j1,x1,y1,mc)
				call get_r_index(j2,x1+d(d12,1),y1+d(d12,2),mc)
				call get_r_index(j3,x1+deltax,y1+deltay,mc)
				call get_r_index(j4,x1+deltax+d(d34,1),y1+deltay+d(d34,2),mc)

				G43=-mc%G%m(j4,j3)
				if (j4.eq.j3) G43 = 1.d0 + G43
				G21=-mc%G%m(j2,j1)
				if (j2.eq.j1) G21 = 1.d0 + G21
				G23=mc%G%m(j2,j3)
				G41=-mc%G%m(j4,j1)
				if (j4.eq.j1) G41 = 1.d0 + G41


				meas%ETCF(1,1,1,1,d34,d12,deltay+1,deltax+1)=meas%ETCF(1,1,1,1,d34,d12,deltay+1,deltax+1) + G43 * G21 + G23 * G41
				meas%ETCF(2,2,2,2,d34,d12,deltay+1,deltax+1)=meas%ETCF(2,2,2,2,d34,d12,deltay+1,deltax+1) + myconj(G43) * myconj(G21) + myconj(G23) * myconj(G41)

				meas%ETCF(1,1,2,2,d34,d12,deltay+1,deltax+1)=meas%ETCF(1,1,2,2,d34,d12,deltay+1,deltax+1) + G43 * myconj(G21)
				meas%ETCF(2,2,1,1,d34,d12,deltay+1,deltax+1)=meas%ETCF(2,2,1,1,d34,d12,deltay+1,deltax+1) + myconj(G43) * G21

				meas%ETCF(1,2,2,1,d34,d12,deltay+1,deltax+1)=meas%ETCF(1,2,2,1,d34,d12,deltay+1,deltax+1) + myconj(G23) * G41
				meas%ETCF(1,2,2,1,d34,d12,deltay+1,deltax+1)=meas%ETCF(1,2,2,1,d34,d12,deltay+1,deltax+1) + G23 * myconj(G41)


!~ 				do s1 =1,2; do s2=1,2; do s3=1,2; do s4=1,2
!~ 					dETCF=0.d0
!~ 					if((s1.eq.s2).and.(s3.eq.s4)) dETCF = dETCF + G43 * G21
!~ 					if((s2.eq.s3).and.(s4.eq.s1)) dETCF = dETCF + G23 * G41
!~ 					meas%ETCF(s4,s3,s2,s1,d34,d12,deltay+1,deltax+1)=meas%ETCF(s4,s3,s2,s1,d34,d12,deltay+1,deltax+1) +dETCF
!~ 				end do; end do ; end do ; end do

			end do; end do
		end do; end do;
	end do; end do

end subroutine add_ETCF

subroutine write_CF(CF,unitnum,mc)
	use global
	implicit none
	type(mcconfig) mc
	$thop_type$ :: CF(2,2,2,2,5,5,mc%Lx,mc%Ly)
    $thop_type$ c
	integer(4) :: deltax,deltay, d12,d34,s1,s2,s3,s4, unitnum

		do deltax=0,(mc%Lx-1); do deltay=0,(mc%Ly-1)
			do d12=1,5; do d34=1,5
				do s1 =1,2; do s2=1,2; do s3=1,2; do s4=1,2
                    c = CF(s4,s3,s2,s1,d34,d12,deltay+1,deltax+1)
					write(unitnum) dble(c)
                    if (is_complex.eq.1) then
                        write(unitnum) aimag(c + cmplx(0.d0))
                    endif
				end do; end do; end do; end do;
			end do; end do
		end do; end do
        ! do deltax=1,mc%Lx; do deltay=1,mc%Ly
        !     write(6,*) CF(1,1,1,1,1,1,deltax,deltay)
        ! enddo; enddo;
end subroutine write_CF


subroutine write_CCCF(meas, mc)
	use global
	implicit none
	type(mcconfig) mc
	type(mcMeasure) meas
	integer(4) jdir, jx, jy, jt, jqx, jqy

	do jdir=1,2
		do jx=1,mc%Lx; do jy=1,mc%Ly
	        write(28) meas%cccf(jdir, jx, jy)
		end do; end do
	end do

	do jdir=1,2
        do jqx=1,meas%cccf2_N_q
            do jqy=1,meas%cccf2_N_q
        		do jt=1,mc%Ntau
        	        write(31) meas%cccf2(jdir, jqx, jqy, jt)
                enddo
            enddo
		end do
	end do

end subroutine

!! Computes the equal time, local pairing correlation ETCF(s4,s3,s2,s1,d34,d12,deltar)\sum_r <c^dag(r,s1)c^dag(r+d12,s2)c(r+dr,s3)c(r+dr+d34,s4)> for the current timeslice.
subroutine add_ETPC(meas,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(mcMeasure) meas
	$thop_type$ :: G32,G41,G42,G31
	! $G_type$ :: dETPC
	integer(4) :: d(5,2)
	integer(4) :: x1,y1,d12,d34,j1,j2,j3,j4
	integer(4) :: deltax,deltay
	integer(4) :: s1,s2,s3,s4

	d=meas%d

	do deltax=0,(mc%Lx-1); do deltay=0,(mc%Ly-1)
		do x1=1,mc%Lx; do y1=1,mc%Ly
			do d12=1,5; do d34=1,5
				call get_r_index(j1,x1,y1,mc)
				call get_r_index(j2,x1+d(d12,1),y1+d(d12,2),mc)
				call get_r_index(j3,x1+deltax,y1+deltay,mc)
				call get_r_index(j4,x1+deltax+d(d34,1),y1+deltay+d(d34,2),mc)

				G32=-mc%G%m(j3,j2)
				if (j3.eq.j2) G32 = 1.d0 + G32
				G41=-mc%G%m(j4,j1)
				if (j4.eq.j1) G41 = 1.d0 + G41
				G42=-mc%G%m(j4,j2)
				if (j4.eq.j2) G42 = 1.d0 + G42
				G31=-mc%G%m(j3,j1)
				if (j3.eq.j1) G31 = 1.d0 + G31

				meas%ETPC(1,1,1,1,d34,d12,deltay+1,deltax+1)=meas%ETPC(1,1,1,1,d34,d12,deltay+1,deltax+1) + G32 * G41 - G42 * G31
				meas%ETPC(2,2,2,2,d34,d12,deltay+1,deltax+1)=meas%ETPC(2,2,2,2,d34,d12,deltay+1,deltax+1) + myconj(G32) * myconj(G41) - myconj(G42) * myconj(G31)

				meas%ETPC(2,1,1,2,d34,d12,deltay+1,deltax+1)=meas%ETPC(1,1,2,2,d34,d12,deltay+1,deltax+1) + G32 * myconj(G41)
				meas%ETPC(1,2,2,1,d34,d12,deltay+1,deltax+1)=meas%ETPC(2,2,1,1,d34,d12,deltay+1,deltax+1) + myconj(G32) * G41

				meas%ETPC(1,2,1,2,d34,d12,deltay+1,deltax+1)=meas%ETPC(1,2,1,2,d34,d12,deltay+1,deltax+1) - G42 * myconj(G31)
				meas%ETPC(2,1,2,1,d34,d12,deltay+1,deltax+1)=meas%ETPC(2,1,2,1,d34,d12,deltay+1,deltax+1) - myconj(G42) * G31


!~ 				do s1 =1,2; do s2=1,2; do s3=1,2; do s4=1,2
!~ 					dETPC=0.d0
!~ 					if((s3.eq.s2).and.(s4.eq.s1)) dETPC = dETPC + G32 * G41
!~ 					if((s4.eq.s2).and.(s3.eq.s1)) dETPC = dETPC - G42 * G31
!~ 					meas%ETPC(s4,s3,s2,s1,d34,d12,deltay+1,deltax+1)=meas%ETPC(s4,s3,s2,s1,d34,d12,deltay+1,deltax+1) +dETPC
!~ 				end do; end do ; end do ; end do

			end do; end do
		end do; end do;
	end do; end do

end subroutine add_ETPC

subroutine get_r_index(j,x,y,mc)
	use global
	implicit none
	integer(4):: j,x,y
	type(mcconfig) mc
	j=mc%i(modulo(x-1,mc%Lx)+1,modulo(y-1,mc%Ly)+1)
end subroutine



subroutine allocate_tdgf(tdgf, mc)
	use global
	implicit none
	type(mcconfig) mc
	type(meas_tdgf) tdgf
	integer(4) j,Nbackwards

	allocate(tdgf%Gt0(mc%Ntau))
	allocate(tdgf%G0t(mc%Ntau))

    if (is_complex.eq.1) then
    	do j=1,mc%Ntau
    		call allocate_cmatrix(tdgf%Gt0(j),mc%lat_dim)
    		tdgf%Gt0(j)%m=0.d0
    		call allocate_cmatrix(tdgf%G0t(j),mc%lat_dim)
    		tdgf%G0t(j)%m=0.d0
    	end do
    else
        do j=1,mc%Ntau
            call allocate_rmatrix(tdgf%Gt0(j),mc%lat_dim)
            tdgf%Gt0(j)%m=0.d0
            call allocate_rmatrix(tdgf%G0t(j),mc%lat_dim)
            tdgf%G0t(j)%m=0.d0
        end do
    endif

end subroutine

subroutine deallocate_tdgf(tdgf,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(meas_tdgf) tdgf
	integer(4) j

    if (is_complex.eq.1) then
    	do j=1,mc%Ntau
    		call deallocate_cmatrix(tdgf%Gt0(j))
    		call deallocate_cmatrix(tdgf%G0t(j))
    	end do
    else
        do j=1,mc%Ntau
            call deallocate_rmatrix(tdgf%Gt0(j))
            call deallocate_rmatrix(tdgf%G0t(j))
        end do
    endif

	deallocate(tdgf%Gt0)
	deallocate(tdgf%G0t)

end subroutine

subroutine add_mean_tdgf(tdgf,mean_tdgf,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(meas_tdgf) tdgf,mean_tdgf
	integer(4) j

	do j=1,mc%Ntau
		mean_tdgf%Gt0(j)%m=	mean_tdgf%Gt0(j)%m + tdgf%Gt0(j)%m
		mean_tdgf%G0t(j)%m=	mean_tdgf%G0t(j)%m + tdgf%G0t(j)%m
	end do
end subroutine


!!!!! Write mean time displaced Green's function in K space:
subroutine write_tdgf_K(tdgf,mc)
    use global
    implicit none
    type(mcconfig) mc
    type(meas_tdgf) tdgf
    integer(4) jt

	do jt=1,mc%Ntau
		call write_Gk(tdgf%Gt0(jt),mc,25)
	end do

	do jt=1,mc%Ntau
		call write_GK(tdgf%G0t(jt),mc,25)
	end do

end subroutine write_tdgf_K

subroutine write_tdgf(tdgf,mc)
    use global
    implicit none
    type(mcconfig) mc
    type(meas_tdgf) tdgf
    $thop_type$ g
    integer(4) j1, j2, jt

		do jt=1, mc%Ntau
			do j2 = 1, tdgf%Gt0(jt)%d
			    do j1 = 1, tdgf%Gt0(jt)%d
                    g = tdgf%Gt0(jt)%m(j1, j2)
			        write(25) dble(g)
                    if (is_complex.eq.1) then
			            write(25) aimag(g + cmplx(0.d0))
                    endif
			    end do
			end do
		end do

		do jt=1, mc%Ntau
			do j2 = 1, tdgf%G0t(jt)%d
			    do j1 = 1, tdgf%G0t(jt)%d
                    g = tdgf%G0t(jt)%m(j1, j2)
                    write(25) dble(g)
                    if (is_complex.eq.1) then
			             write(25) aimag(g + cmplx(0.d0))
                    endif
			    end do
			end do
		end do

end subroutine write_tdgf

subroutine initMeas(meas,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(mcmeasure) meas
	integer(4) d(5,2)

    meas%cccf2_N_q = 4 !number of momenta to store in cccf2
	meas%flags=mc%measflags
	call init_meas_logic(meas%flags)
	call open_meas_files(meas,mc)
	call allocate_meas(meas,mc)

!! matrix of displacements: the first index is just a name (1:5) , the second indicates the direction (x or y). the value is the displacement: 0,1, or -1.
	meas%d(:,1)=(/0,1,-1,0,0/);
	meas%d(:,2)=(/0,0,0,1,-1/);
	meas%minusd = (/ 1,3,2,5,4 /) !For a given d, this is the index of -d.

	meas%n_meas = 0
	meas%n_write_meas = mc%measflags%n_write_meas


end subroutine initMeas


subroutine open_meas_files(meas,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(mcmeasure) meas

	if (meas%flags%is_overwrite.eq.1) then
		if (meas%flags%is_measure_G.eq.1)    open(17,file=trim(adjustl(mc%run_name))//'_G.out', status='replace', form='unformatted', access='stream')
		if (meas%flags%is_measure_Gk.eq.1)   open(19,file=trim(adjustl(mc%run_name))//'_Gk.out', status='replace', form='unformatted', access='stream')
		if (meas%flags%is_measure_ETCF.eq.1)   open(23,file=trim(adjustl(mc%run_name))//'_ETCF.out', status='replace', form='unformatted', access='stream')
		if (meas%flags%is_measure_ETPC.eq.1)   open(24,file=trim(adjustl(mc%run_name))//'_ETPC.out', status='replace', form='unformatted', access='stream')
		if (meas%flags%is_measure_TDGF.eq.1)   open(25,file=trim(adjustl(mc%run_name))//'_TDGF.out', status='replace', form='unformatted', access='stream')
		if (meas%flags%is_measure_TDCF.eq.1)   open(26,file=trim(adjustl(mc%run_name))//'_TDCF.out', status='replace', form='unformatted', access='stream')
		if (meas%flags%is_measure_TDPC.eq.1)   open(27,file=trim(adjustl(mc%run_name))//'_TDPC.out', status='replace', form='unformatted', access='stream')
		if (meas%flags%is_measure_CCCF.eq.1)   open(28,file=trim(adjustl(mc%run_name))//'_CCCF.out', status='replace', form='unformatted', access='stream')
		if (meas%flags%is_measure_E.eq.1)   open(29,file=trim(adjustl(mc%run_name))//'_E.out', status='replace', form='unformatted', access='stream')
		if (meas%flags%is_measure_A.eq.1)   open(30,file=trim(adjustl(mc%run_name))//'_A.out', status='replace', form='unformatted', access='stream')
		if (meas%flags%is_measure_CCCF.eq.1)   open(31,file=trim(adjustl(mc%run_name))//'_CCCF2.out', status='replace', form='unformatted', access='stream')
        if (meas%flags%is_measure_chi.eq.1)   open(32,file=trim(adjustl(mc%run_name))//'_chi.out', status='replace', form='unformatted', access='stream')
        if (meas%flags%is_measure_chi.eq.1)   open(34,file=trim(adjustl(mc%run_name))//'_chi2.out', status='replace', form='unformatted', access='stream')
        if (meas%flags%is_measure_vertex.eq.1)   open(33,file=trim(adjustl(mc%run_name))//'_vertex.out', status='replace', form='unformatted', access='stream')
	else
		if (meas%flags%is_measure_G.eq.1)    open(17,file=trim(adjustl(mc%run_name))//'_G.out',status='unknown', form='unformatted', access='stream', position='append')
		if (meas%flags%is_measure_Gk.eq.1)   open(19,file=trim(adjustl(mc%run_name))//'_Gk.out',status='unknown', form='unformatted', access='stream', position='append')
		if (meas%flags%is_measure_ETCF.eq.1)   open(23,file=trim(adjustl(mc%run_name))//'_ETCF.out',status='unknown', form='unformatted', access='stream', position='append')
		if (meas%flags%is_measure_ETPC.eq.1)   open(24,file=trim(adjustl(mc%run_name))//'_ETPC.out',status='unknown', form='unformatted', access='stream', position='append')
		if (meas%flags%is_measure_TDGF.eq.1)   open(25,file=trim(adjustl(mc%run_name))//'_TDGF.out',status='unknown', form='unformatted', access='stream', position='append')
		if (meas%flags%is_measure_TDCF.eq.1)   open(26,file=trim(adjustl(mc%run_name))//'_TDCF.out',status='unknown', form='unformatted', access='stream', position='append')
		if (meas%flags%is_measure_TDPC.eq.1)   open(27,file=trim(adjustl(mc%run_name))//'_TDPC.out',status='unknown', form='unformatted', access='stream', position='append')
		if (meas%flags%is_measure_CCCF.eq.1)   open(28,file=trim(adjustl(mc%run_name))//'_CCCF.out',status='unknown', form='unformatted', access='stream', position='append')
		if (meas%flags%is_measure_E.eq.1)   open(29,file=trim(adjustl(mc%run_name))//'_E.out',status='unknown', form='unformatted', access='stream', position='append')
		if (meas%flags%is_measure_A.eq.1)   open(30,file=trim(adjustl(mc%run_name))//'_A.out',status='unknown', form='unformatted', access='stream', position='append')
		if (meas%flags%is_measure_CCCF.eq.1)   open(31,file=trim(adjustl(mc%run_name))//'_CCCF2.out',status='unknown', form='unformatted', access='stream', position='append')
        if (meas%flags%is_measure_chi.eq.1)   open(32,file=trim(adjustl(mc%run_name))//'_chi.out',status='unknown', form='unformatted', access='stream', position='append')
        if (meas%flags%is_measure_chi.eq.1)   open(34,file=trim(adjustl(mc%run_name))//'_chi2.out',status='unknown', form='unformatted', access='stream', position='append')
        if (meas%flags%is_measure_vertex.eq.1)   open(33,file=trim(adjustl(mc%run_name))//'_vertex.out',status='unknown', form='unformatted', access='stream', position='append')
	endif
end subroutine open_meas_files


subroutine init_meas_logic(flags)
	use global
	implicit none
	type(MCMeasureFlags) flags

		if ((flags%is_measure_Gk.eq.1).or.(flags%is_measure_G.eq.1).or.(flags%is_measure_tdcf.eq.1)&
			&.or.(flags%is_measure_etcf.eq.1).or.(flags%is_measure_etpc.eq.1).or.(flags%is_measure_etcf.eq.1).or.(flags%is_measure_etpc.eq.1)) then !we need the etgf for the tdcf, etcf
			flags%is_measure_etgf=1
		else
			flags%is_measure_etgf=0
		end if

		if ((flags%is_measure_tdgf_K.eq.1)) then
			flags%is_measure_tdgf=1
		end if

		if (((flags%is_measure_tdcf.eq.1).or.(flags%is_measure_tdpc.eq.1).or.&
        &(flags%is_measure_cccf.eq.1).or.(flags%is_measure_chi.eq.1).or.(flags%is_measure_vertex.eq.1))&
        &.and.(flags%is_measure_tdgf.ne.1)) then !we need the tdgf
			flags%is_measure_tdgf=1
			flags%is_measure_tdgf_K =1   !just so it won't write down too much data...
		end if

	flags%is_measure_fermions =0
	if((flags%is_measure_etgf.eq.1).or.(flags%is_measure_tdgf.eq.1).or.(flags%is_measure_E.eq.1)) flags%is_measure_fermions = 1

end subroutine init_meas_logic

subroutine allocate_meas(meas,mc)
	use global
	implicit none
	type(mcMeasure) meas
	type(mcconfig) mc
	integer(4) j
	if(meas%flags%is_measure_fermions.eq.1) then
		allocate(meas%G(mc%Ntau))
		allocate(meas%Gh(mc%Ntau))
        if (is_complex.eq.1) then
    		do j=1,mc%Ntau
    			call allocate_cmatrix(meas%G(j),mc%lat_dim)
    			call allocate_cmatrix(meas%Gh(j),mc%lat_dim)
    		enddo
        else
            do j=1,mc%Ntau
    			call allocate_rmatrix(meas%G(j),mc%lat_dim)
    			call allocate_rmatrix(meas%Gh(j),mc%lat_dim)
    		enddo
        endif

	endif
		!!! Initialize mean G matrix
	if ((meas%flags%is_measure_G.eq.1).or.(meas%flags%is_measure_Gk.eq.1)) then
        if (is_complex.eq.1) then
            call allocate_cmatrix(meas%meanG,mc%lat_dim)
        else
            call allocate_rmatrix(meas%meanG,mc%lat_dim)
        endif
		meas%meanG%m=0.d0
	    if (meas%flags%is_measure_G.eq.1) write(6,*) 'Writing Greens function measurements to *_G.out'
	    if (meas%flags%is_measure_Gk.eq.1) then
	        write(6,*) 'Writing Fourier transformed Greens function measurements to *_Gk.out'
	    endif
	end if


	if (meas%flags%is_measure_ETCF.eq.1) then
		allocate(meas%ETCF(2,2,2,2,5,5,mc%Ly,mc%Lx))
		meas%ETCF=0.d0
		write(6,*) 'Writing equal time correlation function to *_ETCF.out'
	end if

	if (meas%flags%is_measure_ETPC.eq.1) then
		allocate(meas%ETPC(2,2,2,2,5,5,mc%Ly,mc%Lx))
		meas%ETPC=0.d0
		write(6,*) 'Writing equal time pairing correlations to *_ETPC.out'
	end if

	!!! Initialize time-displaced Greens function.
	if (meas%flags%is_measure_tdgf.eq.1) then
		call allocate_tdgf(meas%tdgf, mc) !! start calculating backwards from tau=dtau*Nsvd/2.
		call allocate_tdgf(meas%mean_tdgf, mc)

		if(meas%flags%is_measure_tdgf_K.eq.1) then
			write(6,*) 'Writing time displaced Gk to *_tdgf.out'
		else
			write(6,*) 'Writing time displaced G to *_tdgf.out'
		end if
	end if

    if (meas%flags%is_measure_vertex.eq.1) then
        allocate(meas%Gk1k2t0(mc%Ntau))
        allocate(meas%Gk1k20t(mc%Ntau))
        do j=1,mc%Ntau
            call allocate_cmatrix(meas%Gk1k2t0(j), mc%lat_dim)
            call allocate_cmatrix(meas%Gk1k20t(j), mc%lat_dim)
        enddo
        allocate(meas%P_vertex(mc%Lx, mc%Ly, mc%Lx, mc%Ly))
    end if

	if (meas%flags%is_measure_TDPC.eq.1) then
		allocate(meas%TDPC(2,2,2,2,5,5,mc%Ly,mc%Lx))
		meas%TDPC=0.d0
		write(6,*) 'Writing time displaced pairing correlations to *_TDPC.out'
	end if

	if (meas%flags%is_measure_TDCF.eq.1) then
		allocate(meas%TDCF(2,2,2,2,5,5,mc%Ly,mc%Lx))
		meas%TDCF=0.d0
		write(6,*) 'Writing time displaced correlation function to *_TDCF.out'
	end if

	if (meas%flags%is_measure_CCCF.eq.1) then
		allocate(meas%CCCF(2,mc%Lx,mc%Ly))
		meas%CCCF=0.d0

		allocate(meas%CCCF2(2,meas%cccf2_N_q,meas%cccf2_N_q,mc%Ntau))
		meas%CCCF2=0.d0
		write(6,*) 'Writing time displaced current-current correlation function to *_CCCF.out and *_CCCF2.out'
	end if
    if (meas%flags%is_measure_chi.eq.1) then
        ! allocate(meas%chi(2,2,mc%Lx, mc%Ly, mc%ntau))
        allocate(meas%chi(2, mc%Lx, mc%Ly, mc%ntau))
        meas%chi = 0.d0
        write(6,*) ' Writing fermionic nematic correlation function to *_chi.out'
    endif

	end subroutine allocate_meas

subroutine end_measure(meas,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(mcmeasure) meas
	integer(4) j

	if(meas%flags%is_measure_fermions.eq.1) then
        if (is_complex.eq.1) then
            do j=1,mc%Ntau
    			call deallocate_cmatrix(meas%G(j))
    			call deallocate_cmatrix(meas%Gh(j))
    		enddo
        else
            do j=1,mc%Ntau
    			call deallocate_rmatrix(meas%G(j))
    			call deallocate_rmatrix(meas%Gh(j))
    		enddo
        endif

		deallocate(meas%G)
		deallocate(meas%Gh)
	endif

	if ((meas%flags%is_measure_G.eq.1).or.(meas%flags%is_measure_Gk.eq.1)) then
        if(is_complex.eq.1) then
            call deallocate_cmatrix(meas%meanG)
        else
            call deallocate_rmatrix(meas%meanG)
        endif
    endif
    if (meas%flags%is_measure_tdgf.eq.1) call deallocate_tdgf(meas%tdgf,mc)

    if (meas%flags%is_measure_vertex.eq.1) then
        do j=1,mc%Ntau
            call deallocate_cmatrix(meas%Gk1k2t0(j))
            call deallocate_cmatrix(meas%Gk1k20t(j))
        enddo
        deallocate(meas%Gk1k2t0)
        deallocate(meas%Gk1k20t)
        deallocate(meas%P_vertex)
    end if

	if (meas%flags%is_measure_ETCF.eq.1) deallocate(meas%ETCF)
	if (meas%flags%is_measure_ETPC.eq.1) deallocate(meas%ETPC)
	if (meas%flags%is_measure_TDCF.eq.1) deallocate(meas%TDCF)
	if (meas%flags%is_measure_TDPC.eq.1) deallocate(meas%TDPC)
	if (meas%flags%is_measure_CCCF.eq.1) deallocate(meas%CCCF,meas%CCCF2)
    if (meas%flags%is_measure_chi.eq.1) deallocate(meas%chi)

!	close files:
	if (meas%flags%is_measure_G.eq.1)    close(17)
	if (meas%flags%is_measure_Gk.eq.1)   close(19)
	if (meas%flags%is_measure_ETCF.eq.1)   close(23)
	if (meas%flags%is_measure_ETPC.eq.1)   close(24)
	if (meas%flags%is_measure_TDGF.eq.1)   close(25)
	if (meas%flags%is_measure_TDCF.eq.1)   close(26)
	if (meas%flags%is_measure_TDPC.eq.1)   close(27)
	if (meas%flags%is_measure_CCCF.eq.1)   then
		close(28)
		close(31)
	endif
	if (meas%flags%is_measure_E.eq.1)   close(29)
    if (meas%flags%is_measure_chi.eq.1)   close(32)
    if (meas%flags%is_measure_vertex.eq.1)   close(33)

end subroutine end_measure


subroutine write_measure(meas,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(mcmeasure) meas
	integer(4) j

	write(6,*) 'Writing out the measurements'
	if(meas%n_meas.eq.0) then
		write(6,*) 'No measurements taken.'
	else
		!!!! Write mean equal time G. averaged over time and measurements
		if ((meas%flags%is_measure_G.eq.1).or.(meas%flags%is_measure_Gk.eq.1)) then
			meas%meanG%m = meas%meanG%m/dble(mc%Ntau*meas%n_write_meas)
			if (meas%flags%is_measure_G.eq.1)  call write_mean_G(meas%meanG,mc)
			if (meas%flags%is_measure_Gk.eq.1) call write_Gk(meas%meanG, mc,19)
			meas%meanG%m = 0.d0
		endif

		!!!! Write mean equal time ETCF. averaged over time and measurements
		if (meas%flags%is_measure_ETCF.eq.1) then
			meas%ETCF=meas%ETCF/dble(mc%Lx*mc%Ly*mc%Ntau*meas%n_write_meas) !Divide by spacetime volume beacause it was averaged over.
			call write_CF(meas%ETCF,23,mc)
			meas%ETCF = 0.d0
		end if

		!!!! Write mean equal time ETPC. averaged over time and measurements
		if (meas%flags%is_measure_ETPC.eq.1) then
			meas%ETPC=meas%ETPC/dble(mc%Lx*mc%Ly*mc%Ntau*meas%n_write_meas)
			call write_CF(meas%ETPC,24,mc)
			meas%ETPC = 0.d0
		end if

		!!! Write time displaced Green's function. Will never write both tdgf_K and tdgf at the same run.
		if (meas%flags%is_measure_tdgf_K.eq.1) then
			do j=1,mc%Ntau
				meas%mean_tdgf%Gt0(j)%m=meas%mean_tdgf%Gt0(j)%m/dble(meas%n_write_meas)
				meas%mean_tdgf%G0t(j)%m=meas%mean_tdgf%G0t(j)%m/dble(meas%n_write_meas)
			end do
			call write_tdgf_K(meas%mean_tdgf,mc)
			do j=1,mc%Ntau
				meas%mean_tdgf%Gt0(j)%m=0.d0
				meas%mean_tdgf%G0t(j)%m=0.d0
			end do
		else if (meas%flags%is_measure_tdgf.eq.1) then
			do j=1,mc%Ntau
				meas%mean_tdgf%Gt0(j)%m=meas%mean_tdgf%Gt0(j)%m/dble(meas%n_write_meas)
				meas%mean_tdgf%G0t(j)%m=meas%mean_tdgf%G0t(j)%m/dble(meas%n_write_meas)
			end do
			call write_tdgf(meas%mean_tdgf,mc)
			do j=1,mc%Ntau
				meas%mean_tdgf%Gt0(j)%m=0.d0
				meas%mean_tdgf%G0t(j)%m=0.d0
			end do
		end if

		if (meas%flags%is_measure_TDCF.eq.1) then
			meas%TDCF=meas%TDCF/dble(meas%n_write_meas*mc%Lx*mc%Ly)*dble(mc%dtau)   !divide by LxLy because there's an average over the initial position, multiply by dtau for the omega=0 transform
			call write_CF(meas%TDCF,26,mc)
			meas%TDCF = 0.d0
		end if

		!!!! Write mean TDPC. averaged over initial time and measurements
		if (meas%flags%is_measure_TDPC.eq.1) then
			meas%TDPC=meas%TDPC/dble(meas%n_write_meas*mc%Lx*mc%Ly)*dble(mc%dtau)   !divide by LxLy because there's an average over the initial position, multiply by dtau for the omega=0 transform
			call write_CF(meas%TDPC,27,mc)
			meas%TDPC = 0.d0
		end if

		!!!! Write mean CCCF. Averaged over measurements
		if (meas%flags%is_measure_CCCF.eq.1) then
			meas%CCCF=meas%CCCF/dble(meas%n_write_meas*mc%Lx*mc%Ly)*dble(mc%dtau)   !divide by LxLy because there's an average over the initial position, multiply by dtau for the omega=0 transform
			meas%CCCF2=meas%CCCF2/dble(meas%n_write_meas*mc%Lx*mc%Ly)  !divide by LxLy because there's an average over the initial position.
			call write_CCCF(meas,mc)
			meas%CCCF = 0.d0
			meas%CCCF2 = 0.d0
		end if

		!!! Write mean energy
		if (meas%flags%is_measure_E.eq.1) then
			meas%E=meas%E/dble(meas%n_write_meas)
			do j=1,7; write(29) meas%E(j); enddo
			meas%E  = 0.d0
		end if

        !!! Write mean chi
		if (meas%flags%is_measure_chi.eq.1) then
			meas%chi=meas%chi/dble(meas%n_write_meas)
            call write_nematic_corr(mc, meas)
			meas%chi = 0.d0
		end if

        !!! Write mean chi
		if (meas%flags%is_measure_vertex.eq.1) then
			meas%P_vertex=meas%P_vertex/dble(meas%n_write_meas)
            call write_vertex(mc, meas)
			meas%P_vertex = 0.d0
		end if

	end if

end subroutine write_measure


!!! Calculate the complete zero frequency expectation value \sum_\tau \sum_r <cdagger1(r,tau)c2(r+d12,tau)cdagger3(r+deltar,0)c4(r+deltar+d34,0)>
!! The index convention tdcf(j4,j3,j2,j1,d12,d34,deltar_y+1,deltar_x+1). The j's are ordered as (up,down), d's are (0,+x,-x,+y,-y), deltar starts at 0.
!!! I think this works only for periodic boundary conditions
subroutine add_tdcf(meas,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(mcMeasure) meas
	integer(4) :: x1,y1,d12,d34,jl1,jl2,jl3,jl4,jt
	integer(4) :: deltax,deltay
	integer(4) :: j1,j2,j3,j4
	$thop_type$ :: G23, G41, G21,G43, Dcorr, Duncorr
	integer(4) :: d(5,2)

	d=meas%d

	do jt=1,mc%Ntau
		do deltax=0,(mc%Lx-1); do deltay=0,(mc%Ly-1)
			do x1=1,mc%Lx; do y1=1,mc%Ly
				do d12=1,5; do d34=1,5
					call get_r_index(j1,x1,y1,mc)
					call get_r_index(j2,x1+d(d12,1),y1+d(d12,2),mc)
					call get_r_index(j3,x1+deltax,y1+deltay,mc)
					call get_r_index(j4,x1+deltax+d(d34,1),y1+deltay+d(d34,2),mc)

					G23=meas%tdgf%gt0(jt)%m(j2,j3)
					G41=meas%tdgf%g0t(jt)%m(j4,j1)

					G21=meas%Gh(jt)%m(j2,j1)
					G43=meas%Gh(1)%m(j4,j3)

					Dcorr=-G23*G41
					Duncorr=G21*G43


					meas%tdcf(1,1,1,1,d34,d12,deltay+1,deltax+1)=meas%tdcf(1,1,1,1,d34,d12,deltay+1,deltax+1) - G23 * G41 + G21 * G43
					meas%tdcf(2,2,2,2,d34,d12,deltay+1,deltax+1)=meas%tdcf(2,2,2,2,d34,d12,deltay+1,deltax+1) - myconj(G23) * myconj(G41) + myconj(G21) * myconj(G43)

					meas%tdcf(1,2,2,1,d34,d12,deltay+1,deltax+1)=meas%tdcf(1,2,2,1,d34,d12,deltay+1,deltax+1) - myconj(G23) * G41
					meas%tdcf(2,1,1,2,d34,d12,deltay+1,deltax+1)=meas%tdcf(2,1,1,2,d34,d12,deltay+1,deltax+1) - G23 * myconj(G41)

					meas%tdcf(1,1,2,2,d34,d12,deltay+1,deltax+1)=meas%tdcf(1,1,2,2,d34,d12,deltay+1,deltax+1) + myconj(G21) * G43
					meas%tdcf(2,2,1,1,d34,d12,deltay+1,deltax+1)=meas%tdcf(2,2,1,1,d34,d12,deltay+1,deltax+1) + G21 * myconj(G43)

!~ 					do jl1 =1,2; do jl2=1,2;
!~ 						meas%tdcf(jl1,jl2,jl2,jl1,d34,d12,deltay+1,deltax+1)=meas%tdcf(jl1,jl2,jl2,jl1,d34,d12,deltay+1,deltax+1) +Dcorr
!~ 					end do; end do
!~
!~ 					do jl1 =1,2; do jl3=1,2;
!~ 						meas%tdcf(jl3,jl3,jl1,jl1,d34,d12,deltay+1,deltax+1)=meas%tdcf(jl3,jl3,jl1,jl1,d34,d12,deltay+1,deltax+1) +Duncorr
!~ 					end do; end do

				end do; end do
			end do; end do;
		end do; end do
	end do
	call profiling('TDCF calculation')
end subroutine add_tdcf

!! creates the effective thop(dir, site_index, time_index) where dir=x,y, for each site.
!! The site is taken to be the 'first' site in the bond. Hopefully that means the site with smaller x/y (depending on the direction of the bond).
!! Shouldn't work when thop is complex!!
subroutine get_thop(thop, mc)
	use global
	implicit none
	type(mcconfig) mc
	real(8) thop(2,mc%lat_dim,mc%Ntau)
	integer(4) ind, jtau, jmat,j

	do jmat=1,2
		do j = 1, mc%kin_mat(jmat)%num_pairs
			ind = mc%kin_mat(jmat)%ind(j,1)
			do jtau=1,mc%Ntau
	    		thop(1,ind,jtau) = mc%par%thop*(1.d0-mc%par%alpha*mc%o_p(jmat)%eta(j,jtau))
	    	end do
	    end do
    end do
	do jmat=3,4
		do j = 1, mc%kin_mat(jmat)%num_pairs
			ind = mc%kin_mat(jmat)%ind(j,1)
			do jtau=1,mc%Ntau
	    		thop(2,ind,jtau) = mc%par%thop*(1.d0+mc%par%alpha*mc%o_p(jmat)%eta(j,jtau))
	    	end do
	    end do
    end do

end subroutine get_thop

!!! calculating <J(q,tau) J(-q, 0)> = 4 tr(G(0,0+) I(-q, 0)) tr(G(tau,tau+) I(q,tau)) -2 Tr(G(tau,0) I( -q, 0) G(0, tau) I( q, tau))
!!!!!!!! Calculate the current-current correlation function (\sum_r_0 \int \sum_\tau jx(r_0,tau)jx(r+r_0,0) and the same for jyjy).
!!!!!!!! To get the right correlation function I later divide by dtau.
!!!!!!!! For the Fz case, we define j = (jup, jdn)
!!!!!!!! Here I use:
!!!!!!!! <jup jdn >_corr = 0
!!!!!!!! <jdn jdn >_corr = <jup jup>_corr ^*.
!!!!!!!! <j1 j2 >_uncorr = <j1><j2> = (<jup1> + <jdn1>)(<jup2> + <jdn2>) = 4 Im(<jup1>) Im(<jup2>)
subroutine calc_cccf(meas,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(mcMeasure) meas
	integer :: jt, jdir, jmat1, jmat3, jp1, jp3, deltax, deltay, ind1, ind2, ind3, ind4, nqx, nqy
    integer :: N_q_max
	$thop_type$ :: thop1, thop3, Dcorr, Duncorr, D
	real(8) :: alpha1, alpha3, qx, qy

    N_q_max = meas%cccf2_N_q

	do jdir = 0,1 !Jx, Jy
		do jmat1=(1+2*jdir),(2+2*jdir); do jp1 = 1, mc%kin_mat(jmat1)%num_pairs
			do jmat3 = (1+2*jdir),(2+2*jdir); do jp3 = 1, mc%kin_mat(jmat3)%num_pairs
				alpha1 =mc%par%alpha
				if((jmat1.eq.1).or.(jmat1.eq.2)) alpha1=-alpha1

				alpha3 =mc%par%alpha
				if((jmat3.eq.1).or.(jmat3.eq.2)) alpha3=-alpha3

				Dcorr = cmplx(0.d0)
				Duncorr = cmplx(0.d0)
				do jt = 1,mc%Ntau
					thop1 = mc%kin_mat(jmat1)%thop(jp1) * (1.d0+alpha1*mc%o_p(jmat1)%eta(jp1,jt))
					thop3 = mc%kin_mat(jmat3)%thop(jp3) * (1.d0+alpha3*mc%o_p(jmat3)%eta(jp3,1))
					ind1 = mc%kin_mat(jmat1)%ind(jp1,1)
					ind2 = mc%kin_mat(jmat1)%ind(jp1,2)
					ind3 = mc%kin_mat(jmat3)%ind(jp3,1)
					ind4 = mc%kin_mat(jmat3)%ind(jp3,2)

					deltax = mc%ix(ind1)-mc%ix(ind3)
					if(deltax.lt.0) deltax = deltax + mc%Lx
					deltay = mc%iy(ind1)-mc%iy(ind3)
					if(deltay.lt.0) deltay = deltay + mc%Ly

					!correlated part
					D = 	thop1		 * thop3 		* meas%tdgf%gt0(jt)%m(ind2,ind3) * meas%tdgf%g0t(jt)%m(ind4,ind1)

					D = D - thop1 		 * myconj(thop3) * meas%tdgf%gt0(jt)%m(ind2,ind4) * meas%tdgf%g0t(jt)%m(ind3,ind1)
					D = D - myconj(thop1) * thop3 		* meas%tdgf%gt0(jt)%m(ind1,ind3) * meas%tdgf%g0t(jt)%m(ind4,ind2)
					D = D + myconj(thop1) * myconj(thop3) * meas%tdgf%gt0(jt)%m(ind1,ind4) * meas%tdgf%g0t(jt)%m(ind3,ind2)
					D = D + myconj(D)
					!uncorrelated part
					D = D - 4.d0 * dble(thop1 * meas%Gh(jt)%m(ind1,ind2) - myconj(thop1) * meas%Gh(jt)%m(ind2,ind1)) * dble(thop3 * meas%Gh(1)%m(ind3,ind4) - myconj(thop3) * meas%Gh(1)%m(ind4,ind3))
					meas%cccf(jdir+1,deltax+1,deltay+1) = meas%cccf(jdir+1,deltax+1,deltay+1) + dble(D)

                    do nqx=1,N_q_max
                        qx = (nqx -1) * 2 * PI / dble(mc%Lx)
                        do nqy=1,N_q_max
                            qy = (nqy -1) * 2 * PI / dble(mc%Ly)
                            meas%cccf2(jdir+1,nqx,nqy,jt) = meas%cccf2(jdir+1,nqx,nqy,jt) + dble(D) * cos(qx * deltax) * cos(qy * deltay) !imposing reflection symmetries.
                        enddo
                    enddo

				enddo
			enddo;enddo
		enddo;enddo
	enddo
	call profiling('CCCF calculation')
end subroutine calc_cccf

! Calculates the expectation value of the Hamiltonian and adds to meas%E.
subroutine add_E(meas, mc)
	use global
	implicit none
	type(mcconfig) mc
	type(mcmeasure) meas
	integer jdir, jmat, jp, ind1, ind2, jt,j
	$thop_type$ thop
	real(8) H_spatial, H_temporal, H_x, H_y, H_mu, H_const, alpha

	! bosonic spatial part
	H_spatial = 0.d0
	call calc_S_eta(H_spatial,mc,1) ! a=1: only the spatial part will be calculated.
	H_spatial = H_spatial / mc%dtau / dble(mc%Ntau)


	! bosonic temporal part
	H_temporal = 0.d0
	call calc_S_eta(H_temporal,mc,2) ! a=2: only the temporal (for h_z =0) part will be calculated. The difference between the action and the hamiltonian amounts to a multiplicative constant.
	H_temporal = -H_temporal / mc%K_cl * mc%par%h/sinh(2.d0 * mc%dtau * mc%par%h) /dble(mc%Ntau)


	! Constant term from the temporal part
	H_const = -2.0 * mc%Lx * mc%Ly * mc%par%h * 0.5 * (tanh(mc%dtau * mc%par%h) + 1.d0 / tanh(mc%dtau * mc%par%h))

	! kinetic energy:
	H_x = 0.d0
	H_y = 0.d0
	do jmat=1,4;
		alpha = mc%par%alpha
		if((jmat.eq.1).or.(jmat.eq.2)) alpha = -alpha
		do jp = 1, mc%kin_mat(jmat)%num_pairs
			ind1 = mc%kin_mat(jmat)%ind(jp,1)
			ind2 = mc%kin_mat(jmat)%ind(jp,2)
			do jt = 1, mc%Ntau
				thop = mc%kin_mat(jmat)%thop(jp) * (1.d0+alpha*mc%o_p(jmat)%eta(jp,jt))
				if((jmat.eq.1).or.(jmat.eq.2)) then
					H_x = H_x - 2.d0 * dble(thop * meas%Gh(jt)%m(ind2,ind1) + myconj(thop) * meas%Gh(jt)%m(ind1,ind2)) / dble(mc%Ntau) !looks complicated due to the spin structure.
				elseif ((jmat.eq.3).or.(jmat.eq.4)) then
					H_y = H_y - 2.d0 * dble(thop * meas%Gh(jt)%m(ind2,ind1) + myconj(thop) * meas%Gh(jt)%m(ind1,ind2)) / dble(mc%Ntau) !looks complicated due to the spin structure.
				endif
			enddo
		enddo
	enddo

	!chemical potential:
	H_mu = 0.d0
	do jt=1, mc%Ntau
		do j=1, mc%lat_dim
			H_mu = H_mu -2 * mc%par%mu * dble(meas%Gh(jt)%m(j,j)) / dble(mc%Ntau) ! factor of 2 for spin is included.
		enddo
	enddo

	meas%E(1) = meas%E(1) + H_spatial + H_temporal + H_const + H_x + H_y + H_mu
	meas%E(2) = meas%E(2) + H_spatial
	meas%E(3) = meas%E(3) + H_temporal
	meas%E(4) = meas%E(4) + H_const
	meas%E(5) = meas%E(5) + H_x
	meas%E(6) = meas%E(6) + H_y
	meas%E(7) = meas%E(7) + H_mu

end subroutine


!! computes tdpc(j4,j3,j2,j1,d34,d12,deltar) = sum_r sum_tau <c1(r,tau)c2(r+d12,tau)c3^dag(r+deltar,0) c4^dag(r+deltar+d34, 0)>
subroutine add_tdpc(meas,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(meas_tdgf) tdgf
	type(mcMeasure) meas
	integer(4) :: x1,y1,d12,d34,jl1,jl2,jl3,jl4,jt
	integer(4) :: deltax,deltay
	integer(4) :: j1,j2,j3,j4
	$thop_type$ :: G23, G14, G24, G13
	integer(4) :: d(5,2)

	d=meas%d

	do jt=1,mc%Ntau
		do deltax=0,(mc%Lx-1); do deltay=0,(mc%Ly-1)
			do x1=1,mc%Lx; do y1=1,mc%Ly
				do d12=1,5; do d34=1,5
					call get_r_index(j1,x1,y1,mc)
					call get_r_index(j2,x1+d(d12,1),y1+d(d12,2),mc)

					call get_r_index(j3,x1+deltax,y1+deltay,mc)
					call get_r_index(j4,x1+deltax+d(d34,1),y1+deltay+d(d34,2),mc)

					G23=meas%tdgf%gt0(jt)%m(j2,j3)
					G14=meas%tdgf%gt0(jt)%m(j1,j4)

					G24=meas%tdgf%gt0(jt)%m(j2,j4)
					G13=meas%tdgf%gt0(jt)%m(j1,j3)

					meas%tdpc(1,1,1,1,d34,d12,deltay+1,deltax+1)=meas%tdpc(1,1,1,1,d34,d12,deltay+1,deltax+1) + G23 * G14 - G24 * G13
					meas%tdpc(2,2,2,2,d34,d12,deltay+1,deltax+1)=meas%tdpc(2,2,2,2,d34,d12,deltay+1,deltax+1) + myconj(G23) * myconj(G14) - myconj(G24) * myconj(G13)

					meas%tdpc(1,2,2,1,d34,d12,deltay+1,deltax+1)=meas%tdpc(1,2,2,1,d34,d12,deltay+1,deltax+1) + myconj(G23) * G14
					meas%tdpc(2,1,1,2,d34,d12,deltay+1,deltax+1)=meas%tdpc(2,1,1,2,d34,d12,deltay+1,deltax+1) + G23 * myconj(G14)

					meas%tdpc(2,1,2,1,d34,d12,deltay+1,deltax+1)=meas%tdpc(2,1,2,1,d34,d12,deltay+1,deltax+1) - myconj(G24) * G13
					meas%tdpc(1,2,1,2,d34,d12,deltay+1,deltax+1)=meas%tdpc(1,2,1,2,d34,d12,deltay+1,deltax+1) - G24 * myconj(G13)


!~ 					do jl1 =1,2; do jl2=1,2;
!~ 						meas%tdpc(jl1,jl2,jl2,jl1,d34,d12,deltay+1,deltax+1)=meas%tdpc(jl1,jl2,jl2,jl1,d34,d12,deltay+1,deltax+1) +G23*G14
!~ 						meas%tdpc(jl2,jl1,jl2,jl1,d34,d12,deltay+1,deltax+1)=meas%tdpc(jl2,jl1,jl2,jl1,d34,d12,deltay+1,deltax+1) -G24*G13
!~ 					end do; end do;
				end do; end do
			end do; end do;
		end do; end do
	end do

	call profiling('TDPC calculation')
end subroutine add_tdpc

subroutine scalettar_tdgf_both_sides(tdgf,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(meas_tdgf) tdgf
	integer(4) taun, j, n, d
	type(svd), allocatable :: BT0(:), BBetaT(:), BT0Inv(:), BBetaTInv(:)
	type(svd) USV

	d= mc%lat_dim
	allocate(BT0(mc%Nsvd),BBetaT(mc%Nsvd),BT0Inv(mc%Nsvd),BBetaTInv(mc%Nsvd))
	do j=1,mc%Nsvd
		call allocate_usv(BT0(j),d)
		call allocate_usv(BBetaT(j),d)
		call allocate_usv(BT0Inv(j),d)
		call allocate_usv(BBetaTInv(j),d)
	enddo

	call allocate_usv(USV,d)

	call calc_B_USV(mc,BT0,LEFT,NO_INV) ! this is not necessary if initG was run before.
	call calc_B_USV(mc,BBetaT,RIGHT,NO_INV)

	call calc_B_USV(mc,BBetaTInv,LEFT,INV)
	call calc_B_USV(mc,BT0Inv,RIGHT,INV)


	do n=1,mc%Nsvd
		taun=min((n-1)*mc%Ntau_mul+1,mc%Ntau)
		if (n.eq.1) then
			call stable_invert_B(USV, BBetaT(n), NO_INV)
			call USV_to_mat(tdgf%Gt0(taun),USV, NO_INV)
			call stable_invert_B(USV, BBetaTInv(n), NO_INV)
			call USV_to_mat(tdgf%G0t(taun),USV, NO_INV)
		else
			call stable_inverse_sum(USV, BT0Inv(n-1), BBetaT(n), NO_INV, NO_INV)
			call USV_to_mat(tdgf%Gt0(taun),USV, NO_INV)
			call stable_inverse_sum(USV, BT0(n-1), BBetaTInv(n), NO_INV, NO_INV)
			call USV_to_mat(tdgf%G0t(taun),USV, NO_INV)
		endif
		tdgf%G0t(taun)%m=-tdgf%G0t(taun)%m

	end do

	call fill_tdgf(tdgf,mc)

	do j=1,mc%Nsvd
		call deallocate_usv(BT0(j))
		call deallocate_usv(BBetaT(j))
		call deallocate_usv(BT0Inv(j))
		call deallocate_usv(BBetaTInv(j))
	enddo

	deallocate(BT0,BBetaT,BT0Inv,BBetaTInv)
	call profiling('TDGF calculation')
end subroutine scalettar_tdgf_both_sides

!Given tdgf at every n*Ntau_mul, fill in between the n*Ntau_mul's:
subroutine fill_tdgf(tdgf,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(meas_tdgf) tdgf
	integer(4) taun, j, n

	do n=1,mc%Nsvd
		taun=min((n-1)*mc%Ntau_mul+1,mc%Ntau)
		do j=1, mc%Ntau_mul-1
			if (taun+j>mc%Ntau) exit
			tdgf%Gt0(taun+j)%m = tdgf%Gt0(taun+j-1)%m
			tdgf%G0t(taun+j)%m = tdgf%G0t(taun+j-1)%m
			call OperateB(LEFT, NO_INV, tdgf%Gt0(taun+j) , taun+j-1, taun+j-1,mc)
			call OperateB(RIGHT, INV, tdgf%G0t(taun+j) , taun+j-1, taun+j-1,mc)
		end do
	end do

end subroutine fill_tdgf

subroutine get_equal_time_G(G,Gh,mc)
	use global
	implicit none
	type(mcconfig) mc
	$G_type$ :: G(mc%Ntau), Gh(mc%Ntau)
	integer(4) nt, dir, j

	call initG(mc)
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 	open(94,file='epsilon.out', status='replace', form='unformatted', access='stream')
!~ 	write(94)  log(mc%USV(1)%S)/(mc%dtau * mc%Ntau_mul)
!~ 	close(94)
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	dir = LEFT
	do nt=1,(mc%Ntau)
		if(nt.gt.1) then
	    	!!! Advance Green's function
	    	if (mod(mc%n+dir,mc%Ntau_mul).ne.0) then
	        	call rotateG(-dir, mc)
	    	else
	        	call stab_reconstructG(-dir,mc)
	    	endif
	    end if
	    G(mc%n)%m = mc%G%m
	    Gh(mc%n)%m = -mc%G%m
	    do j=1,mc%lat_dim
	    	Gh(mc%n)%m(j,j)=Gh(mc%n)%m(j,j)+1.d0
    	enddo
   enddo

end subroutine

! !calculates the full, electronic bond-density bond-density correlation function.
! subroutine calc_nematic_corr_old(meas,mc)
! 	use global
! 	implicit none
! 	type(mcconfig) mc
! 	type(mcMeasure) meas
! 	integer jt, jdir1, jdir3, jmat1, jmat3, jp1, jp3, deltax, deltay, ind1, ind2, ind3, ind4
! 	complex(8) thop1, thop3, thop1_full, thop3_full, Dcorr, Duncorr, D
! 	real(8) alpha1, alpha3
!
! 	do jmat1=1,4; do jp1 = 1, mc%kin_mat(jmat1)%num_pairs
! 		do jmat3 = 1,4; do jp3 = 1, mc%kin_mat(jmat3)%num_pairs
! 			alpha1 = mc%par%alpha
!             jdir1 = 1
! 			if((jmat1.eq.1).or.(jmat1.eq.2)) then
!                 alpha1=-alpha1
!                 jdir1 = 2
!             endif
!
! 			alpha3 = mc%par%alpha
!             jdir3 = 1
! 			if((jmat3.eq.1).or.(jmat3.eq.2)) then
!                 alpha3=-alpha3
!                 jdir3 = 2
!             endif
!
! 			Dcorr = cmplx(0.d0)
! 			Duncorr = cmplx(0.d0)
! 			do jt = 1,mc%Ntau
! 				thop1 = mc%kin_mat(jmat1)%thop(jp1) !* (1.d0+alpha1*mc%o_p(jmat1)%eta(jp1,jt))
! 				thop3 = mc%kin_mat(jmat3)%thop(jp3) !* (1.d0+alpha3*mc%o_p(jmat3)%eta(jp3,1))
! 				ind1 = mc%kin_mat(jmat1)%ind(jp1,1)
! 				ind2 = mc%kin_mat(jmat1)%ind(jp1,2)
! 				ind3 = mc%kin_mat(jmat3)%ind(jp3,1)
! 				ind4 = mc%kin_mat(jmat3)%ind(jp3,2)
!
! 				deltax = mc%ix(ind1)-mc%ix(ind3)
! 				if(deltax.lt.0) deltax = deltax + mc%Lx
! 				deltay = mc%iy(ind1)-mc%iy(ind3)
! 				if(deltay.lt.0) deltay = deltay + mc%Ly
!
! 				!correlated part
! 				D = 	thop1		 * thop3 		* meas%tdgf%gt0(jt)%m(ind2,ind3) * meas%tdgf%g0t(jt)%m(ind4,ind1)
! 				D = D + thop1 		 * conjg(thop3) * meas%tdgf%gt0(jt)%m(ind2,ind4) * meas%tdgf%g0t(jt)%m(ind3,ind1)
! 				D = D + conjg(thop1) * thop3 		* meas%tdgf%gt0(jt)%m(ind1,ind3) * meas%tdgf%g0t(jt)%m(ind4,ind2)
! 				D = D + conjg(thop1) * conjg(thop3) * meas%tdgf%gt0(jt)%m(ind1,ind4) * meas%tdgf%g0t(jt)%m(ind3,ind2)
! 				D = D + conjg(D)
! 				!uncorrelated part
! 				! D = D - 4.d0 * dble(thop1 * meas%Gh(jt)%m(ind1,ind2) + conjg(thop1) * meas%Gh(jt)%m(ind2,ind1)) * dble(thop3 * meas%Gh(1)%m(ind3,ind4) + conjg(thop3) * meas%Gh(1)%m(ind4,ind3))
!
!                 meas%chi(jdir1, jdir3, jt, deltax+1,deltay+1) = meas%chi(jdir1, jdir3, jt, deltax+1,deltay+1) + dble(D)
!
! 			enddo
! 		enddo;enddo
! 	enddo;enddo
! end subroutine calc_nematic_corr_old

!calculates the electronic nematic correlation function and the electronic-bosonic cross-correlator.
! The order parameters are defined on site.
subroutine calc_nematic_corr(meas,mc)
	use global
	implicit none
	type(mcconfig) mc
	type(mcMeasure) meas
	integer :: j1, j3, jt, jdir1, jdir3, jmat1, jmat3, jp1, jp3, deltax, deltay, ind1, ind2, ind3, ind4
	$thop_type$ :: thop1, thop3, thop1_full, thop3_full, D
    integer :: sgn1, sgn3

	do jmat1=1,4; do j1 = 1, mc%lat_dim
		do jmat3 = 1,4; do j3 = 1, mc%lat_dim
            jdir1 = 2
            sgn1 = 1
			if((jmat1.eq.1).or.(jmat1.eq.2)) then
                jdir1 = 1
                sgn1 = -1
            endif
            jp1 = mc%kin_mat(jmat1)%pair_ind(j1,2)
            ind1 = mc%kin_mat(jmat1)%ind(jp1,1)
            ind2 = mc%kin_mat(jmat1)%ind(jp1,2)

            jdir3 = 2
            sgn3 = 1
			if((jmat3.eq.1).or.(jmat3.eq.2)) then
                jdir3 = 1
                sgn3 = -1
            endif
            jp3 = mc%kin_mat(jmat3)%pair_ind(j3,2)
            ind3 = mc%kin_mat(jmat3)%ind(jp3,1)
            ind4 = mc%kin_mat(jmat3)%ind(jp3,2)

            deltax = mc%ix(j1)-mc%ix(j3)
            if(deltax.lt.0) deltax = deltax + mc%Lx
            deltay = mc%iy(j1)-mc%iy(j3)
            if(deltay.lt.0) deltay = deltay + mc%Ly

            ! write(6,*) jmat1, sgn1, jp1, j1, ind1, ind2
            ! write(6,*) jmat3, sgn3, jp3, j3, ind3, ind4
            ! write(6,*) '!!!'

			do jt = 1,mc%Ntau
				thop1 = sgn1 * mc%kin_mat(jmat1)%thop(jp1) !* (1.d0 + mc%par%alpha * sgn1 * mc%o_p(jmat1)%eta(jp1,jt))
				thop3 = sgn3 * mc%kin_mat(jmat3)%thop(jp3) !* (1.d0 + mc%par%alpha * sgn3 * mc%o_p(jmat3)%eta(jp3,1))

				!correlated part
				D = 	-thop1		 * thop3 		* meas%tdgf%gt0(jt)%m(ind2,ind3) * meas%tdgf%g0t(jt)%m(ind4,ind1)
				D = D - thop1 		 * myconj(thop3) * meas%tdgf%gt0(jt)%m(ind2,ind4) * meas%tdgf%g0t(jt)%m(ind3,ind1)
				D = D - myconj(thop1) * thop3 		* meas%tdgf%gt0(jt)%m(ind1,ind3) * meas%tdgf%g0t(jt)%m(ind4,ind2)
				D = D - myconj(thop1) * myconj(thop3) * meas%tdgf%gt0(jt)%m(ind1,ind4) * meas%tdgf%g0t(jt)%m(ind3,ind2)
				D = D + myconj(D)
				! uncorrelated part
				D = D + 4.d0 * dble(thop1 * meas%Gh(jt)%m(ind1,ind2) + myconj(thop1) * meas%Gh(jt)%m(ind2,ind1)) * dble(thop3 * meas%Gh(1)%m(ind3,ind4) + myconj(thop3) * meas%Gh(1)%m(ind4,ind3))

                meas%chi(1,deltax+1,deltay+1,jt) = meas%chi(1,deltax+1,deltay+1,jt) + dble(D)

                D = (thop1 * meas%Gh(1)%m(ind1,ind2) +myconj(thop1) * meas%Gh(jt)%m(ind2,ind1) ) * mc%o_p(jmat3)%eta(jp3,jt)
                meas%chi(2,deltax+1,deltay+1,jt) = meas%chi(2,deltax+1,deltay+1,jt) + dble(D)
			enddo
		enddo;enddo
	enddo;enddo
end subroutine calc_nematic_corr


subroutine write_nematic_corr(mc, meas)
    use global
    implicit none
    type(mcconfig) mc
    type(mcmeasure) meas
    real(8) chi, q_max, omega_max, omega0, q0, qx, qy, omega
    integer(4) jo, jqx, jqy, n_omega_max, n_q_max
    integer(4) jx, jy, jt, jcross
    integer(4) unitnum1, unitnum2
    unitnum1 = 32
    unitnum2 = 34
    q_max = PI
    omega_max = 10.0

    q0 = 2 * PI / mc%Lx
    omega0 = 2 * PI / mc%ntau
    n_omega_max = int(omega_max / omega0 * mc%dtau)
    n_q_max = int(q_max /q0)


    do jcross=1,2
        do jqx=1,n_q_max+1; do jqy=1,n_q_max+1; do jo=1,n_omega_max+1
            chi = 0.d0
            do jx = 1,mc%Lx; do jy = 1,mc%Ly; do jt=1,mc%ntau
                qx = q0 * (jqx - 1)
                qy = q0 * (jqy - 1)
                omega = omega0 * (jo-1)
                chi = chi + cos(qx * (jx-1) + qy * (jy-1) + omega * (jt-1)) * meas%chi(jcross,jx, jy,jt) * mc%dtau / dble(8 * mc%Lx * mc%Ly)
                ! divided by the volume because of the average over the COM. multiplied by dtau for the freq. transform. The 8 is just a convention
            enddo; enddo; enddo
            write(unitnum1) chi
        enddo; enddo; enddo
    enddo

    do jcross=1,2
        do jx=1,mc%Lx; do jy=1,mc%Ly
            write(unitnum2) meas%chi(jcross,jx,jy,1)
            write(unitnum2) meas%chi(jcross,jx,jy,mc%ntau/2)
        enddo; enddo
    enddo

end subroutine

subroutine transform_Gk1k2(G, Gk, mc)
    use global
    implicit none
    type(mcconfig) mc
    $G_type$ G
    type(cmatrix) Gk, U, G_temp
    integer(4) d, jk, jr
    real(8) dkx, dky
    complex(8) i

    d = mc%lat_dim

    call allocate_cmatrix(U,d)
    call allocate_cmatrix(G_temp,d)
    dkx = 2 * PI /dble(mc%Lx)
    dkx = 2 * PI /dble(mc%Ly)

    i = (0.d0, 1.d0)
    do jr=1,d
        do jk=1,d
            U%m(jk,jr) = exp( i * dkx * (mc%ix(jk) +mc%Fx) * mc%ix(jr) + i * dky * (mc%iy(jk) + mc%Fy) * mc%iy(jr) )
        enddo
    enddo
    U%m = U%m / dble(mc%Lx * mc%Ly)

    call zsq_matmul(G_temp,U,G,'N', 'N')
    call zsq_matmul(Gk,G_temp,U,'N', 'C')


    call deallocate_cmatrix(U)
    call deallocate_cmatrix(G_temp)

end subroutine transform_Gk1k2


subroutine calc_P_vertex(meas, mc)
    use global
    implicit none
    type(mcconfig) mc
    type(mcMeasure) meas
    type(cmatrix) Gk, G_temp
    integer(4) d, jt, j1, j2, j3, j4, unitnum
    integer(4) jkx1, jky1, jkx2, jky2, jkx3, jky3, jkx4, jky4
    real(8) P

    ! d = mc%lat_dim
    ! call allocate_cmatrix(Gk, d)
    ! call allocate_cmatrix(G_temp, d)
    ! allocate(tdpc(mc%Lx, mc%Ly, mc%Lx, mc%Ly))
    ! tdpc = 0.d0

    do jt=1,mc%ntau
        do j1 = 1, mc%lat_dim
            jkx1= mc%ix(j1)
            jky1 = mc%iy(j1)

            jkx2 = mc%Lx - jkx1
            if(jkx2.eq.0) jkx2 = mc%Lx
            jky2 = mc%Ly - jky1
            if(jky2.eq.0) jky2 = mc%Ly
            j2 = mc%i(jkx2,jky2)

            do j3 = 1, mc%lat_dim
                jkx3= mc%ix(j3)
                jky3 = mc%iy(j3)

                jkx4 = mc%Lx - jkx3
                if(jkx4.eq.0) jkx4 = mc%Lx
                jky4 = mc%Ly - jky3
                if(jky4.eq.0) jky4 = mc%Ly
                j4 = mc%i(jkx4,jky4)

                P =     abs(meas%Gk1k2t0(jt)%m(j1, j3))**2 + abs(meas%Gk1k2t0(jt)%m(j2, j3))**2 &
                + abs(meas%Gk1k2t0(jt)%m(j1, j4))**2 + abs(meas%Gk1k2t0(jt)%m(j2, j4))**2

                meas%P_vertex(jkx1, jky1, jkx3, jky3) = meas%P_vertex(jkx1, jky1, jkx3, jky3) + P
            enddo
        enddo
    enddo

end subroutine calc_P_vertex


subroutine get_tdgf_k1k2(meas,mc)
    use global
    implicit none
    type(mcconfig) mc
    type(mcMeasure) meas
    type(cmatrix) G_temp
    integer(4) d, jt

    d = mc%lat_dim
    call allocate_cmatrix(G_temp, d)

    do jt=1,mc%ntau
        G_temp%m = meas%tdgf%gt0(jt)%m
        call transform_Gk1k2(meas%Gk1k2t0(jt), G_temp ,mc)
        G_temp%m = meas%tdgf%g0t(jt)%m
        call transform_Gk1k2(meas%Gk1k20t(jt), G_temp ,mc)
    enddo

    call deallocate_cmatrix(G_temp)
    call profiling('Fourier transforming G ')
end subroutine


! subroutine write_nematic_corr_old(mc, meas)
!     use global
!     implicit none
!     type(mcconfig) mc
!     type(mcmeasure) meas
!     real(8) chi, q_max, omega_max, omega0, q0, arg, qx, qy, f
!     integer(4) jo, jqx, jqy, n_omega_max, n_q_max, jdir1, jdir2, ex1, ey1, ex2, ey2
!     integer(4) jx, jy, jt
!     integer(4) unitnum
!     unitnum = 32
!     q_max = PI / 2.0
!     omega_max = 10.0
!
!     q0 = 2 * PI / mc%Lx
!     omega0 = 2 * PI / mc%ntau
!     n_omega_max = int(omega_max / omega0 * mc%dtau)
!     n_q_max = int(q_max /q0)
!
!
!     do jqx=1,n_q_max+1; do jqy=1,n_q_max+1; do jo=1,n_omega_max+1
!         chi = 0.d0
!         do jdir1=1,2;
!             if (jdir1.eq.1) then
!                 ex1 = 1
!                 ey1 = 0
!             else
!                 ex1 = 0
!                 ey1 = 1
!             endif
!             do jdir2=1,2
!                 if (jdir2.eq.1) then
!                     ex2 = 1
!                     ey2 = 0
!                 else
!                     ex2 = 0
!                     ey2 = 1
!                 endif
!                 do jx = 1,mc%Lx; do jy = 1,mc%Ly; do jt=1,mc%ntau
!                     qx = q0 * (jqx - 1)
!                     qy = q0 * (jqy - 1)
!                     arg =        qx * (dble(jx - 1) + dble(ex1 - ex2)/2.d0)
!                     arg = arg + qy * (dble(jy - 1) + dble(ey1 - ey2)/2.d0)
!                     f = cos(arg) * cos((qx * ex1 + qy * ey1)/2.d0) * cos((qx * ex2 + qy * ey2)/2.d0) * cos(omega0 * (jo -1) * (jt-1))
!                     chi = chi + f * meas%chi(jdir1, jdir2, jt, jx, jy) * mc%dtau / dble(mc%Lx * mc%Ly)
!                     ! divided by the volume because of the average over the COM. multiplied by dtau for the freq. transform
!                 enddo; enddo; enddo
!             enddo
!         enddo
!         write(unitnum) chi
!         write(6,*) q0, omega0, jqx-1, jqy-1, jo-1, chi
!     enddo; enddo; enddo
!
! end subroutine

!~ subroutine add_ebonds(meas, mc)
!~ 	use global
!~ 	implicit none
!~ 	type(mcconfig) mc
!~ 	type(mcmeasure) meas
!~ 	integer jdir, jmat, jp, ind1, ind2, jt
!~ 	complex(8) thop
!~ 	real alpha
!~
!~ 	do jdir = 0,1 !Jx, Jy
!~ 		do jmat=(1+2*jdir),(2+2*jdir)
!~ 			alpha = mc%par%alpha
!~ 			if((jmat.eq.1).or.(jmat.eq.2)) alpha=-alpha
!~
!~ 			do jp = 1, mc%kin_mat(jmat)%num_pairs
!~ 				thop = mc%kin_mat(jmat)%thop(jp) * (1.d0 + alpha * mc%o_p(jmat)%eta(j,jtau))
!~
!~ 				ind1 = mc%kin_mat(jmat)%ind(jp,1)
!~ 				ind2 = mc%kin_mat(jmat)%ind(jp,2)
!~ 				do jt = 1, mc%Ntau
!~ 					meas%Ebonds(jdir+1) = meas%Ebonds(jdir+1)- 2.d0 * dble(thop * meas%Gh(jt)%m(ind2,ind1) + conjg(thop) * meas%Gh(jt)%m(ind1,ind2))
!~ 				enddo
!~ 			enddo
!~ 		enddo
!~ 	enddo
!~
!~ end subroutine

subroutine write_vertex(mc, meas)
    use global
    implicit none
    type(mcconfig) mc
    type(McMeasure) meas
    integer(4) jkx1, jkx2, jky1, jky2

    do jkx1 = 1,mc%Lx; do jky1 = 1,mc%Ly
        do jkx2 = 1,mc%Lx; do jky2 = 1,mc%Ly
            write(33) meas%P_vertex(jkx1, jky1, jkx2, jky2)
        enddo; enddo
    enddo; enddo

end subroutine write_vertex
