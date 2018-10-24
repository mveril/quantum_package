subroutine getTe(Te)
  BEGIN_DOC
  ! Kinetic energy
  END_DOC
  implicit none
  PROVIDE mo_bielec_integrals_in_map 
  integer::i
  double precision,intent(out)::Te
  Te=0d0
  do i=1,mo_tot_num
    Te += mo_occ(i)*mo_kinetic_integral(i,i)
  enddo
  print *, '******'
  print *, 'Kinetic energy', Te
  print *, '******'
end subroutine
