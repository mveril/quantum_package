double precision function spinint(p,q,r,s)
  
  BEGIN_DOC
  ! spinorbital integrals 
  END_DOC
  
  implicit none
  
  PROVIDE mo_bielec_integrals_in_map 
  
  integer,intent(in) :: p,q,r,s
  double precision :: get_mo_bielec_integral
  integer :: c,d,Conv_spin_Index
  spinint = 0d0
  !Non zero only if same spin
  if (mod(p,2) == mod(r,2) .and. mod(q,2)==mod(s,2))  then
    spinint = get_mo_bielec_integral(Conv_Spin_Index(p),Conv_Spin_Index(q),Conv_Spin_Index(r),Conv_Spin_Index(s),mo_integrals_map)
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
  
  aoff = a + n_spin_occ
  boff = b + n_spin_occ
  ijab_spinint = spinint(i,j,aoff,boff)
end function

double precision function abcd_spinint(a,b,c,d)
  
  BEGIN_DOC
  ! abcd spinorbital integral
  END_DOC
  
  implicit none
  double precision :: spinint
  
  integer,intent(in) :: a,b,c,d
  
  integer :: aoff,boff,coff,doff 
  aoff= a + n_spin_occ
  boff = b + n_spin_occ
  coff = c + n_spin_occ
  doff = d + n_spin_occ
  abcd_spinint = spinint(aoff,boff,coff,doff)
end function

double precision function iajb_spinint(i,a,j,b)
  
  BEGIN_DOC
  ! iajb spinorbital integral
  END_DOC
  
  implicit none
  double precision :: spinint
  
  integer,intent(in) :: i,a,j,b
  
  integer :: aoff,boff
  aoff = a + n_spin_occ
  boff = b + n_spin_occ
  iajb_spinint = spinint(i,aoff,j,boff)
end function

double precision function iabj_spinint(i,a,b,j)

  
  BEGIN_DOC
  ! iajb spinorbital integral
  END_DOC
  
  implicit none
  double precision :: spinint
  
  integer,intent(in) :: i,a,b,j
  
  integer :: aoff,boff
  aoff = a + n_spin_occ
  boff = b + n_spin_occ
  iabj_spinint=spinint(i,aoff,boff,j)
end function

double precision function abij_spinint(a,b,i,j)
  
  BEGIN_DOC
  ! abij spinorbital integral
  END_DOC
  
  implicit none
  double precision :: spinint
  
  integer,intent(in) :: a,b,i,j
  
  integer :: aoff,boff
  aoff = a + n_spin_occ
  boff = b + n_spin_occ
  abij_spinint = spinint(aoff,boff,i,j)
end function
