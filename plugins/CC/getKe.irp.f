subroutine getKe(Ke)
  BEGIN_DOC
  ! Kinetic energy
  END_DOC
  implicit none
  PROVIDE mo_bielec_integrals_in_map 
  integer::i,j
  double precision,allocatable::orb(:),Kint(:,:)
  double precision::get_mo_bielec_integral
  double precision,intent(out)::Ke
  allocate(Kint(mo_tot_num,mo_tot_num),orb(mo_tot_num))
  orb(:)=0.5d0*mo_occ(:)
  do j=1, mo_tot_num
    do i=1, mo_tot_num
      Kint(i,j)=get_mo_bielec_integral(i,j,j,i,mo_integrals_map)
    enddo
  enddo
  Ke=dot_product(orb,matmul(Kint,orb))
  print *, '******'
  print *, 'exchange energy', Ke
  print *, '******'
end subroutine
