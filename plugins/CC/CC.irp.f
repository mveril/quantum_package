program CC
  
  implicit none
  
  integer :: i,Conv_Spin_Index,a
  BEGIN_DOC
! Coupled Cluster module
  END_DOC
  call Write_Double(6,HF_energy,"Hartree Fock energy")
  ! At this time CC can only perform CCD
  if (Debug) then
    print *, "Occupeid"
    do i=1,n_spin_occ
      print *, fock_matrix_mo(Conv_Spin_Index(i),Conv_Spin_Index(i))
    end do
    print *, "Virtuals"
    do a=1,n_spin_virt
      print *,fock_matrix_mo(Conv_Spin_Index(a+n_spin_occ),Conv_Spin_Index(a+n_spin_occ))
      !call write_double(6,eigenvalues_fock_matrix_ao(Conv_Spin_Index(a+n_spin_occ)))
    end do
  end if
  call CCD()
end
