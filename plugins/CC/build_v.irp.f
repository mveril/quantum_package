subroutine Build_V(v,oovv_db_spin_int,t2)
  
  BEGIN_DOC
  ! V matrix (defined on equation 21 of pople paper 
  END_DOC

  implicit none

  double precision,intent(out) :: v(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

  double precision,intent(in)  :: oovv_db_spin_int(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in)  :: t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

! Initialise

  v(:,:,:,:) = 0d0

! Accumulate
  call accumulate_X1(v,oovv_db_spin_int,t2) 

  call accumulate_X2(v,oovv_db_spin_int,t2) 

  call accumulate_X3(v,oovv_db_spin_int,t2) 

  call accumulate_X4(v,oovv_db_spin_int,t2) 

end subroutine
