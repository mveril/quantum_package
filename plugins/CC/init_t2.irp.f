  subroutine init_t2(t2,oovv_db_spin_int,oovv_Delta)

  BEGIN_DOC
  ! Initialize CCD amplitudes with MP2 amplitudes
  END_DOC

  implicit none

  double precision,intent(in)  :: ijab_antispinint(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in)  :: oovv_Delta(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

  integer::i,j,a,b
  double precision:: Ec_MP2

  double precision,intent(out) :: t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

  ! Initialize variables 

  Ec_MP2 = 0d0
  t2(:,:,:,:) = 0d0

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          t2(i,j,a,b) = -oovv_db_spin_int(i,j,a,b)/oovv_Delta(i,j,a,b)
          Ec_MP2 += 0.25d0*oovv_db_spin_int(i,j,a,b)*t2(i,j,a,b)
        enddo
      enddo
    enddo
  enddo

  call write_double(6,Ec_MP2,"MP2 correction: ")

end subroutine init_t2
