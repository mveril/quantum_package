subroutine getVe(Ve)
  BEGIN_DOC
  ! Electrons-nucleous repulsion energy
  END_DOC
  implicit none
  PROVIDE mo_bielec_integrals_in_map 
  integer::i
  double precision,intent(out)::Ve
  Ve=0d0
  do i=1,mo_tot_num
    Ve += mo_occ(i)*mo_nucl_elec_integral(i,i)
  enddo
  print *, '******'
  print *, 'Nuclear potential', Ve
  print *, '******'
end subroutine
