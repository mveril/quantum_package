program CC
  
  implicit none
  
  BEGIN_DOC
!   Coupled Cluster module
  END_DOC
  
  integer :: i,Conv_Spin_Index,a

  call Write_Double(6,HF_energy,"Hartree Fock energy")
  ! At this time CC can only perform CCD

  if (Debug) then
    print *, "Occupied"
    do i=1,n_spin_occ
      print *, fock_matrix_mo(Conv_Spin_Index(i),Conv_Spin_Index(i))
    end do
    print *, "Virtuals"
    do a=1,n_spin_virt
      print *,fock_matrix_mo(Conv_Spin_Index(a+n_spin_occ),Conv_Spin_Index(a+n_spin_occ))
    end do
  end if

  call CCD()
end
