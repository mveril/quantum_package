subroutine Build_Delta(Delta)
  implicit none

  BEGIN_DOC
  ! Delta
  END_DOC

  double precision,intent(out) :: Delta(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  integer :: i, j, a, b, Conv_Spin_Index
  Delta(:,:,:,:) = 0d0
  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          Delta(i, j, a, b) = eigenvalues_fock_matrix_ao(Conv_Spin_Index(a+n_spin_occ))+eigenvalues_fock_matrix_ao(Conv_Spin_Index(b+n_spin_occ))-eigenvalues_fock_matrix_ao(Conv_Spin_Index(i))-eigenvalues_fock_matrix_ao(Conv_Spin_Index(j))
        enddo
      enddo
    enddo
  enddo

end subroutine Build_Delta
