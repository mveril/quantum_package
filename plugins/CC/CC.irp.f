program CC
  
  implicit none
  
  BEGIN_DOC
! Coupled Cluster module
  END_DOC
  call write_double(6,HF_energy,"Hartree Fock energy")
  ! At this time CC can only perform CCD
  call CCD()
end
