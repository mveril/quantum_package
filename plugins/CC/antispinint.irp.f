double precision function antispinint(p,q,r,s)
  
  BEGIN_DOC
  ! antisymetric spinorbital integral
  END_DOC
  
  implicit none
  
  integer,intent(in) :: p,q,r,s
  
  double precision :: spinint 
 
  antispinint=spinint(p,q,r,s)-spinint(p,q,s,r)
end function

double precision function ijab_antispinint(i,j,a,b)
  
  BEGIN_DOC
  ! ijab antisymetric spinorbital integral
  END_DOC
  
  implicit none
  
  integer,intent(in) :: i,j,a,b
  integer :: aoff,boff
  double precision :: antispinint 
  
  aoff=a+n_spin_occ
  boff=b+n_spin_occ 
  ijab_antispinint=antispinint(i,j,aoff,boff)
end function


double precision function abcd_antispinint(a,b,c,d)
  
  BEGIN_DOC
  ! abcd antisymetric spinorbital integral
  END_DOC
  
  implicit none
  
  integer,intent(in) :: a,b,c,d

  integer :: aoff,boff,coff,doff 
  double precision :: antispinint 
  
  aoff=a+n_spin_occ
  boff=b+n_spin_occ 
  coff=c+n_spin_occ
  doff=d+n_spin_occ
  abcd_antispinint=antispinint(aoff,boff,coff,doff)
end function

double precision function iajb_antispinint(i,a,j,b)
  
  BEGIN_DOC
  ! abcd antisymetric spinorbital integral
  END_DOC
  
  implicit none
  
  integer,intent(in) :: i,a,j,b

  integer :: aoff,boff
  double precision :: antispinint 
  
  aoff=a+n_spin_occ
  boff=b+n_spin_occ 
  iajb_antispinint=antispinint(i,aoff,j,boff)
end function

double precision function abij_antispinint(a,b,i,j)
  
  BEGIN_DOC
  ! abcd antisymetric spinorbital integral
  END_DOC
  
  implicit none
  
  integer,intent(in) :: a,b,i,j

  integer :: aoff,boff
  double precision :: antispinint 
  
  aoff=a+n_spin_occ
  boff=b+n_spin_occ
  abij_antispinint=antispinint(aoff,boff,i,j)
end function
