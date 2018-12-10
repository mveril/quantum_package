subroutine build_residual(r,t2,vvoo_db_spin_int,u,v,oovv_Delta)

  BEGIN_DOC
  ! Build the CCD residual
  END_DOC
  
  implicit none
  
! Output variables 
  double precision,intent(in)  :: abij_antispinint(n_spin_virt,n_spin_virt,n_spin_occ,n_spin_occ)
  double precision,intent(in)  :: u(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in)  :: v(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in)  :: t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in)  :: oovv_Delta(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

! Local variables
  integer                      :: i,j,a,b

! Output variables 
  double precision,intent(out) :: r(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

! Initialize residual

  r(:,:,:,:) = 0d0

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          r(i,j,a,b) = vvoo_db_spin_int(a,b,i,j) + oovv_Delta(i,j,a,b)*t2(i,j,a,b) + u(i,j,a,b) + v(i,j,a,b)
        end do
      end do
    end do
  end do

end subroutine build_residual
