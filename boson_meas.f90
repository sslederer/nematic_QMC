!!measurements done on the eta part only, not involving the fermions.
subroutine boson_meas(mc)
	use global
	implicit none
	type(mcconfig) mc
	type(Rmatrix) A(2,2), meanA(2,2) !A is the boson correlator, which for now has a 2x2 matrix structure in "orbital" space.
	integer (4)  nt, dir, n_meas, stat,k1,k2
		
		!!! Initialize mean G matrix and A matrix
	
	do k1=1,2
		do k2=1,2
			call allocate_Rmatrix(meanA(k1,k2),mc%lat_dim*mc%Ntau)
			meanA(k1,k2)%m=0.d0
			call allocate_Rmatrix(A(k1,k2),mc%lat_dim*mc%Ntau)
		end do
	end do
	
	!!!Skip the first n_eq iterations worth of data
k1=0
	do
		call readeta(mc,stat)
		k1=k1+1
		if((stat /= 0).or.(k1.ge.mc%n_eq)) exit
	end do
	!!! Read eta file until end of file
	

	n_meas = 0
	do
		call readeta(mc,stat)
		if (stat /= 0) exit        

		k2=0
		do k1=1,2
		k2=k2+sum(mc%o_p(k1)%eta)
		end do
		write(6,*) dble(k2)**2
	
	    n_meas = n_meas + 1
		write(6,*) 'Measurement number', n_meas
	
	!!!!! construct boson green's function
		do k1=1,2
			do k2=1,2
				call boson_correlations(A(k1,k2)%m,mc,k1,k2)
			end do
		end do
		write(6,*) sum(A(1,1)%m)

		do k1=1,2
			do k2=1,2
				meanA(k1,k2)%m=meanA(k1,k2)%m+A(k1,k2)%m
			end do
		end do

	end do

	if(n_meas.eq.0) then
		write(6,*) 'No boson measurements taken.'
	else
		!!! Write mean fourier transformed A avergaed over measurements
		do k1=1,2
			do k2=1,2
				meanA(k1,k2)%m=meanA(k1,k2)%m/dble(n_meas)
				call write_Ak(meanA(k1,k2), mc)
			end do
		end do
	end if	

	!!! Deallocate meanA and A
	do k1=1,2
		do k2=1,2
		    deallocate(meanA(k1,k2)%m)
		    deallocate(A(k1,k2)%m)
		    A(k1,k2)%d=0
		    meanA(k1,k2)%d=0
		end do
	end do	
	
		!! rewind the *_eta.out file so that fermion measurements can take place.
		rewind(16)
		
end subroutine boson_meas

subroutine write_Ak(A, mc)
    use global
    implicit none 
    type(CMatrix) A
    type(McConfig) mc
    complex(8), allocatable :: Ak(:,:,:)
    integer(4) jx, jy,jz

    allocate(Ak(mc%Lx,mc%Ly,mc%Ntau))
    call transform_Ak(A, mc, Ak)

    do jy = 1, mc%Ly; do jx = 1, mc%Lx; do jz= 1,mc%Ntau
        write(30) dble(Ak(jx, jy,jz))
        write(30) aimag(Ak(jx, jy,jz))
    end do; end do; end do
    deallocate(Ak)
end subroutine write_Ak
	
subroutine boson_correlations(A,mc,k1,k2)  !Calculates boson correlations for a given monte carlo configuration
					! k1 or k2 =1 is horizontal, 2 vertical.
	use global
	implicit none
	type(McConfig) mc
	integer(4) k1,k2,jmat1,nt1,j1,nt1m,jmat2,nt2,j2,nt2m
	real(8) A(mc%Lx*mc%Ly*mc%Ntau,mc%Lx*mc%Ly*mc%Ntau), delta, kron
	A=0.0
	
	if((k1.eq.1).and.(k2.eq.1)) then
		do jmat1=1,2
			do jmat2=1,2
				do j1=1,mc%Lx*mc%Ly/2
					do j2=1,mc%Lx*mc%Ly/2
						do nt1=1,mc%Ntau
							nt1m=mod(nt1-2+mc%Ntau,mc%Ntau)+1
							do nt2=1,mc%Ntau
								nt2m=mod(nt2-2+mc%Ntau,mc%Ntau)+1
								A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)=dble(mc%o_p(jmat1)%eta(j1,nt1)*mc%o_p(jmat2)%eta(j2,nt2))
								delta=kron(mc%o_p(jmat2)%eta(j2,nt2),mc%o_p(jmat2)%eta(j2,nt2m))*kron(mc%o_p(jmat1)%eta(j1,nt1),mc%o_p(jmat1)%eta(j1,nt1m))
								A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)=A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)*delta/cosh(mc%par%h*mc%dtau)**2
								if((jmat1.eq.jmat2).and.(j1.eq.j2).and.(nt1.eq.nt2)) then
									A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)=1.0
								endif
							end do
						end do
					end do
				end do
			end do
		end do
		
	else if((k1.eq.1).and.(k2.eq.2)) then
		do jmat1=1,2
			do jmat2=3,4
				do j1=1,mc%Lx*mc%Ly/2
					do j2=1,mc%Lx*mc%Ly/2
						do nt1=1,mc%Ntau
							nt1m=mod(nt1-2+mc%Ntau,mc%Ntau)+1
							do nt2=1,mc%Ntau
								nt2m=mod(nt2-2+mc%Ntau,mc%Ntau)+1
								A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)=dble(mc%o_p(jmat1)%eta(j1,nt1)*mc%o_p(jmat2)%eta(j2,nt2))
								delta=kron(mc%o_p(jmat2)%eta(j2,nt2),mc%o_p(jmat2)%eta(j2,nt2m))*kron(mc%o_p(jmat1)%eta(j1,nt1),mc%o_p(jmat1)%eta(j1,nt1m))
								A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)=A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)*delta/cosh(mc%par%h*mc%dtau)**2
								if((jmat1.eq.jmat2).and.(j1.eq.j2).and.(nt1.eq.nt2)) then
									A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)=1.0
								endif
							end do
						end do
					end do
				end do
			end do
		end do
		
	else if((k1.eq.2).and.(k2.eq.1)) then
		do jmat1=3,4
			do jmat2=1,2
				do j1=1,mc%Lx*mc%Ly/2
					do j2=1,mc%Lx*mc%Ly/2
						do nt1=1,mc%Ntau
							nt1m=mod(nt1-2+mc%Ntau,mc%Ntau)+1
							do nt2=1,mc%Ntau
								nt2m=mod(nt2-2+mc%Ntau,mc%Ntau)+1
								A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)=dble(mc%o_p(jmat1)%eta(j1,nt1)*mc%o_p(jmat2)%eta(j2,nt2))
								delta=kron(mc%o_p(jmat2)%eta(j2,nt2),mc%o_p(jmat2)%eta(j2,nt2m))*kron(mc%o_p(jmat1)%eta(j1,nt1),mc%o_p(jmat1)%eta(j1,nt1m))
								A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)=A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)*delta/cosh(mc%par%h*mc%dtau)**2
								if((jmat1.eq.jmat2).and.(j1.eq.j2).and.(nt1.eq.nt2)) then
									A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)=1.0
								endif
							end do
						end do
					end do
				end do
			end do
		end do
	else if((k1.eq.2).and.(k2.eq.2)) then
		do jmat1=3,4
			do jmat2=3,4
				do j1=1,mc%Lx*mc%Ly/2
					do j2=1,mc%Lx*mc%Ly/2
						do nt1=1,mc%Ntau
							nt1m=mod(nt1-2+mc%Ntau,mc%Ntau)+1
							do nt2=1,mc%Ntau
								nt2m=mod(nt2-2+mc%Ntau,mc%Ntau)+1
								A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)=dble(mc%o_p(jmat1)%eta(j1,nt1)*mc%o_p(jmat2)%eta(j2,nt2))
								delta=kron(mc%o_p(jmat2)%eta(j2,nt2),mc%o_p(jmat2)%eta(j2,nt2m))*kron(mc%o_p(jmat1)%eta(j1,nt1),mc%o_p(jmat1)%eta(j1,nt1m))
								A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)=A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)*delta/cosh(mc%par%h*mc%dtau)**2
								if((jmat1.eq.jmat2).and.(j1.eq.j2).and.(nt1.eq.nt2)) then
									A((mc%o_p(jmat1)%ind(j1,1)-1)*mc%Ntau+nt1,(mc%o_p(jmat2)%ind(j2,1)-1)*mc%Ntau+nt2)=1.0
								endif
							end do
						end do
					end do
				end do
			end do
		end do
	endif

end subroutine boson_correlations

real function kron(a,b)
    integer(4) a,b
    if(a.eq.b) then
    	kron=1
    else
    	kron=0
    endif
    return
    end


subroutine transform_Ak(A, mc, Ak)
    use global
    implicit none 
    type(RMatrix) A
    type(McConfig) mc
    complex(8), allocatable :: Axyz(:,:), Vk(:), Kx(:), Ky(:),Kz(:)
    complex(8) Ak(mc%Lx, mc%Ly,mc%Ntau), coef, theta1, theta2
    integer(4) jx, jy,jz, ix1, ix2, iy1, iy2, iz1,iz2,ix, iy,iz, Nkx, Nky,Nkz

    Nkx = mc%Lx; Nky = mc%Ly; Nkz=mc%Ntau
    allocate(Kx(Nkx), Ky(Nky), Kz(Nkz),Axyz(mc%Lx*mc%Ly*mc%Ntau,mc%Lx*mc%Ly*mc%Ntau), Vk(mc%Lx*mc%Ly*mc%Ntau))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! define Kx, Ky, Kz (ie omega)
    do jx = 1, Nkx
        Kx(jx) = cmplx(dble(jx-1)*(2.d0*PI)/dble(mc%Lx))
    end do
    do jy = 1, Nky
        Ky(jy) = cmplx(dble(jy-1)*(2.d0*PI)/dble(mc%Ly))
    end do
    do jz = 1, Nkz
        Kz(jz) = cmplx(dble(jz-1)*(2.d0*PI)/dble(mc%Ntau))
    end do

    Axyz = (0.d0, 0.d0)
    do ix1 = 1, mc%Lx; do iy1 = 1, mc%Ly; do iz1 = 1, mc%Ntau
        do ix2 = 1, mc%Lx; do iy2 = 1, mc%Ly; do iz2=1,mc%Ntau
            Axyz(ix1+(iy1-1)*mc%Lx+(iz1-1)*mc%Lx*mc%Ly, ix2+(iy2-1)*mc%Lx+(iz2-1)*mc%Lx*mc%Ly) = cmplx(A%m((mc%i(ix1,iy1)-1)*mc%Ntau+iz1, (mc%i(ix2,iy2)-1)*mc%Ntau+iz2))
        end do; end do; end do
	end do; end do; end do

	do jy = 1, Nkx; do jx = 1, Nky; do jz=1,Nkz
		do ix = 1, mc%Lx; do iy = 1, mc%Ly; do iz=1,mc%Ntau
			Vk(ix+(iy-1)*mc%Lx+(iz-1)*mc%Lx*mc%Ly) = exp((0.d0,1.d0)*(cmplx(ix)*Kx(jx)+cmplx(iy)*Ky(jy)+cmplx(iz)*Kz(jz)))
		end do; end do; end do
        Ak(jx, jy,jz) = dot_product(Vk, matmul(Axyz, Vk)) 
    end do; end do; end do

    Ak = Ak/dble(mc%Lx*mc%Ly*mc%Ntau)**2  !!! Don't normalize for now


end subroutine transform_Ak


