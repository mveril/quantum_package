BEGIN_PROVIDER [ double precision, t_coeff, (elec_num,elec_num,mo_tot_num*2-elec_num,mo_tot_num*2-elec_num)]
  implicit none

  BEGIN_DOC
 ! provide t coefficents 
  END_DOC

  integer::i,j,a,b
  integer::n_spin_virt
  double precision :: ijab_antispinint,D_ijab

  n_spin_virt=mo_tot_num*2-elec_num
  do i=1,elec_num
    do j=1,elec_num
      do a=1,n_spin_virt
        do b=1,n_spin_virt
        t_coeff(i,j,a,b)=ijab_antispinint(i,j,a,b)/(-D_ijab(i,j,a,b))
        enddo
      enddo
    enddo
  enddo
END_PROVIDER
