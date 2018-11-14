BEGIN_PROVIDER [ double precision, t_ampli, (n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)]
implicit none

BEGIN_DOC
! provide t coefficents 
END_DOC

integer::i,j,a,b
double precision::MP2_corr
double precision :: ijab_antispinint
  MP2_corr=0d0
  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          t_ampli(i,j,a,b)=ijab_antispinint(i,j,a,b)/(-ijab_D(i,j,a,b))
          MP2_corr+=0.25d0*ijab_antispinint(i,j,a,b)*t_ampli(i,j,a,b)
        enddo
      enddo
    enddo
  enddo
  call write_double(6,MP2_corr,"MP2 correction")
END_PROVIDER
