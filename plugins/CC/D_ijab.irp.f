
double precision function D_ijab(i,j,a,b)
  
  BEGIN_DOC
  ! Delta ijab
  END_DOC
  
  implicit none
  integer,intent(in) :: i,j,a,b
  D_ijab=eigenvalues_fock_matrix_ao((a+1+n_spin_occ)/2)+eigenvalues_fock_matrix_ao((b+1+n_spin_occ)/2)-eigenvalues_fock_matrix_ao((i+1)/2)-eigenvalues_fock_matrix_ao((j+1)/2)
end function
