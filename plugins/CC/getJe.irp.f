subroutine getJe(Je)
  BEGIN_DOC
  ! Kinetic energy
  END_DOC
  implicit none
  PROVIDE mo_bielec_integrals_in_map 
  integer::i,j
  double precision,allocatable::orb(:),jint(:,:)
  double precision::get_mo_bielec_integral
  double precision,intent(out)::Je
  allocate(Jint(mo_tot_num,mo_tot_num),orb(mo_tot_num))
  orb(:)=0.5d0*mo_occ(:)
  do j=1, mo_tot_num
    do i=1, mo_tot_num
      Jint(i,j)=get_mo_bielec_integral(i,j,i,j,mo_integrals_map)
    enddo
  enddo
  Je=dot_product(orb,matmul(Jint,orb))
  print *, '******'
  print *, 'Coulonb energy', Je
  print *, '******'
end subroutine
