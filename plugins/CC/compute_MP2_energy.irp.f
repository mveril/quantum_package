subroutine compute_MP2_energy(E_HF,oovv_Delta,oovv_db_spin_int,Ec_MP2,E_MP2)

BEGIN_DOC
! Compute MP2 energy
END_DOC

  implicit none

! Input variable

  double precision,intent(in)  :: E_HF
  double precision,intent(in)  :: oovv_db_spin_int(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in)  :: oovv_Delta(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

! Local variable

  integer                      :: i,j,a,b

! Ouput variable

  double precision,intent(out) :: Ec_MP2,E_MP2

  Ec_MP2 = 0d0
  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          Ec_MP2 = Ec_MP2 + oovv_db_spin_int(i,j,a,b)*oovv_db_spin_int(i,j,a,b)/oovv_Delta(i,j,a,b)
        end do
      end do
    end do
  end do
  Ec_MP2 = 0.25d0*Ec_MP2

  E_MP2 = E_HF + Ec_MP2

  print*, 'MP2 correlation energy = ',Ec_MP2
  print*, 'MP2 total       energy = ',E_MP2

end subroutine compute_MP2_energy
