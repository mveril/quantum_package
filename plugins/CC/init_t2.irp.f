subroutine Init_t2(t2,ijab_antispinint,inv_Delta)
implicit none

BEGIN_DOC
! initialise t for CCD 
END_DOC

double precision,intent(out) :: t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
double precision,intent(in),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: ijab_antispinint,inv_Delta
integer::i,j,a,b
double precision::MP2_corr

  t2(:,:,:,:) = 0d0
  MP2_corr = 0d0
  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          t2(i,j,a,b) = ijab_antispinint(i,j,a,b)*(-inv_Delta(i, j, a, b))
          MP2_corr += 0.25d0*ijab_antispinint(i, j, a , b)*t2(i, j, a, b)
        enddo
      enddo
    enddo
  enddo

  call Write_Double(6,MP2_corr,"MP2 correction")
end subroutine
