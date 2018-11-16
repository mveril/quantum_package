double precision function spinint(p,q,r,s)
  
  BEGIN_DOC
  ! spinorbital integrals 
  END_DOC
  
  implicit none
  
  PROVIDE mo_bielec_integrals_in_map 
  
  integer,intent(in) :: p,q,r,s
  double precision :: get_mo_bielec_integral
  integer :: c,d
  spinint=0d0
  !Non zero only if same spin
  if (mod(p,2) == mod(r,2) .and. mod(q,2)==mod(s,2))  then
    spinint = get_mo_bielec_integral((p+1)/2,(q+1)/2,(r+1)/2,(s+1)/2,mo_integrals_map)
  end if
end function

double precision function ijab_spinint(i,j,a,b)
  
  BEGIN_DOC
  ! ijab spinorbital integral
  END_DOC
  
  implicit none
  
  double precision :: spinint
  
  integer,intent(in) :: i,j,a,b
  
  integer :: aoff,boff 
  
  aoff=a+n_spin_occ
  boff=b+n_spin_occ
  ijab_spinint=spinint(i,j,aoff,boff)
end function

double precision function abcd_spinint(a,b,c,d)
  
  BEGIN_DOC
  ! abcd spinorbital integral
  END_DOC
  
  implicit none
  double precision :: spinint
  
  integer,intent(in) :: a,b,c,d
  
  integer :: aoff,boff,coff,doff 
  aoff=a+n_spin_occ
  boff=b+n_spin_occ
  coff=c+n_spin_occ
  doff=d+n_spin_occ
  abcd_spinint=spinint(aoff,boff,coff,doff)
end function
