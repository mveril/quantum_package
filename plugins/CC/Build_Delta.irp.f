subroutine Build_Delta(Delta)  
  implicit none
  BEGIN_DOC
  ! Delta
  END_DOC
  double precision,intent(out) :: Delta(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  integer :: i,j,a,b
  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          Delta(i,j,a,b) = eigenvalues_fock_matrix_ao((a+1+n_spin_occ)/2)+eigenvalues_fock_matrix_ao((b+1+n_spin_occ)/2)-eigenvalues_fock_matrix_ao((i+1)/2)-eigenvalues_fock_matrix_ao((j+1)/2)
        enddo
      enddo
    enddo
  enddo
end subroutine
