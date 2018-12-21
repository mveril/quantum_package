subroutine compute_HF_energy(oooo_db_spin_int,E_HF)

BEGIN_DOC
! Compute Hartree-Fock energy
END_DOC

  implicit none

! Input variable

  double precision,intent(in)  :: oooo_db_spin_int(n_spin_occ,n_spin_occ,n_spin_occ,n_spin_occ)

! Local variable

  integer                      :: i,j
  double precision             :: E_1e,E_2e

! Output variable

  double precision,intent(out) :: E_HF

! Compute one-electron part of the HF energy

  E_1e = 0d0
  do i=1,elec_alpha_num
    E_1e = E_1e + mo_kinetic_integral(i,i) + mo_nucl_elec_integral(i,i)
  end do
  do i=1,elec_beta_num
    E_1e = E_1e + mo_kinetic_integral(i,i) + mo_nucl_elec_integral(i,i)
  end do

! Compute two-electron part of the HF energy

  E_2e = 0d0
  do i=1,n_spin_occ
    do j=1,n_spin_occ
      E_2e = E_2e + oooo_db_spin_int(i,j,i,j)
    end do
  end do
  E_2e = 0.5d0*E_2e

  E_HF = E_1e + E_2e
  print*, 'HF 1e     = ',E_1e
  print*, 'HF 2e     = ',E_2e
  print*, 'HF energy = ',E_HF

end subroutine compute_HF_energy
