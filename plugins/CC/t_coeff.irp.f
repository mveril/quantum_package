BEGIN_PROVIDER [ double precision, t_coeff, (n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)]
implicit none

BEGIN_DOC
! provide t coefficents 
END_DOC

integer::i,j,a,b
double precision :: ijab_antispinint
  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          t_coeff(i,j,a,b)=ijab_antispinint(i,j,a,b)/(-ijab_D(i,j,a,b))
        enddo
      enddo
    enddo
  enddo
END_PROVIDER
