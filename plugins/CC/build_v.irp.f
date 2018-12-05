subroutine Build_v(v,ijab_antispinint,t2)
  
  BEGIN_DOC
  ! V matrix 
  END_DOC
  

  implicit none
double precision,intent(out) :: V(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: ijab_antispinint,t2

! Initialise

  V(:,:,:,:) = 0d0

! Accumulate

  call Accumulate_X1(V,ijab_antispinint,t2) 

  call Accumulate_X2(V,ijab_antispinint,t2) 

  call Accumulate_X3(V,ijab_antispinint,t2) 

  call Accumulate_X4(V,ijab_antispinint,t2) 

end subroutine