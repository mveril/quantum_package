
double precision function D_ijab(i,j,a,b)
  
  BEGIN_DOC
  ! Delta ijab
  END_DOC
  
  implicit none
  integer,intent(in) :: i,j,a,b
  D_ijab=eigenvalues_fock_matrix_ao((a-n_spin_occ)/2+1)+eigenvalues_fock_matrix_ao((b-n_spin_occ)/2+1)-eigenvalues_fock_matrix_ao(i)-eigenvalues_fock_matrix_ao(j)
end function
