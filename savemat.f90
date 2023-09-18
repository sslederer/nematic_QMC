subroutine savemat(dim,mat,unitnum)
	implicit none
	integer(4) :: dim
	integer ios,unitnum
	real(8) :: mat(dim,dim)

	integer(4) k1, k2

!	OPEN(5, IOSTAT=ios, FILE=fname, STATUS='unknown')

	do k1=1, dim
	   do k2=1,dim
		  write(unitnum,*) mat(k2,k1)
	   end do
	end do

	!do k1=1, dim
	!   do k2=1,dim
	!      write(5,*) aimag(mat(k2,k1))
	!   end do
	!end do

!	CLOSE(5, IOSTAT=ios, STATUS='keep')

end subroutine savemat


subroutine save_real_vec(dim,vec,fname)
implicit none
integer(4) :: dim
integer ios
real(8) :: vec(dim)
character :: fname*20 

integer(4) k1

OPEN(5, IOSTAT=ios, FILE=fname, STATUS='unknown')

do k1=1, dim
   write(5,*) vec(k1)
end do

CLOSE(5, IOSTAT=ios, STATUS='keep')

end subroutine save_real_vec

