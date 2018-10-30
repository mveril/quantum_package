BEGIN_PROVIDER [ integer, n_spin_virt]
implicit none

BEGIN_DOC
! number of virtual spinorbitals
END_DOC

n_spin_virt=n_spin_tot-n_spin_occ

END_PROVIDER

BEGIN_PROVIDER [ integer, n_spin_occ]
implicit none

BEGIN_DOC
! number of occuped spinorbitals
END_DOC

n_spin_occ=elec_num

END_PROVIDER

BEGIN_PROVIDER [ integer, n_spin_tot]
implicit none

BEGIN_DOC
! number of spinorbitals
END_DOC

n_spin_tot=mo_tot_num*2

END_PROVIDER
