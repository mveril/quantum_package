subroutine Build_Residual(residual,t2,abij_antispinint,u,v,Delta)
  BEGIN_DOC
  ! residual
  END_DOC
  
  implicit none
  
  integer :: i,j,a,b
  double precision,intent(out),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: residual
  double precision,intent(in) :: abij_antispinint(n_spin_virt,n_spin_virt,n_spin_occ,n_spin_occ)
  double precision,intent(in),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: u,v,t2,Delta
  residual(:,:,:,:) = 0d0
  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          residual(i,j,a,b) = abij_antispinint(a,b,i,j) + Delta(i,j,a,b)*t2(i,j,a,b)+u(i,j,a,b)+v(i,j,a,b)
        end do
      end do
    end do
  end do
end subroutine
