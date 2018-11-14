BEGIN_PROVIDER [ double precision, ijab_D, (n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)]
  implicit none
  BEGIN_DOC
  ! ijab_D ijab
  END_DOC
  
  integer :: i,j,a,b
  ! At start Δ is the MP2 try vector
  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt

        ijab_D(i,j,a,b)=eigenvalues_fock_matrix_ao((a+1+n_spin_occ)/2)+eigenvalues_fock_matrix_ao((b+1+n_spin_occ)/2)-eigenvalues_fock_matrix_ao((i+1)/2)-eigenvalues_fock_matrix_ao((j+1)/2)

        enddo
      enddo
    enddo
  enddo
END_PROVIDER
