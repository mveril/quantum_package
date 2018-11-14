double precision function antispinint(p,q,r,s)
  
  BEGIN_DOC
  ! antisymetric spinorbital integral
  END_DOC
  
  implicit none
  
  integer,intent(in) :: p,q,r,s
  
  double precision :: spinint 
 
  antispinint=spinint(p,q,r,s)-spinint(p,q,s,r)
end function

double precision function offset_antispinint(p,q,r,s,isvirtualp,isvirtualq,isvirtualr,isvirtuals)
  
  BEGIN_DOC
  ! antisymetric spinorbital integral with offset
  END_DOC
  
  implicit none
  
  integer,intent(in) :: p,q,r,s
  logical,intent(in) :: isvirtualp,isvirtualq,isvirtualr,isvirtuals
  
  double precision :: offset_spinint
 
  offset_antispinint=offset_spinint(p,q,r,s,isvirtualp,isvirtualq,isvirtualr,isvirtuals)-offset_spinint(p,q,s,r,isvirtualp,isvirtualq,isvirtuals,isvirtualr)
end function


double precision function ijab_antispinint(i,j,a,b)
  
  BEGIN_DOC
  ! ijab antisymetric spinorbital integral
  END_DOC
  
  implicit none
  
  integer,intent(in) :: i,j,a,b
  
  double precision :: offset_antispinint 
  ijab_antispinint=offset_antispinint(i,j,a,b,.False.,.False.,.True.,.True.)
end function


double precision function abcd_antispinint(a,b,c,d)
  
  BEGIN_DOC
  ! abcd antisymetric spinorbital integral
  END_DOC
  
  implicit none
  
  integer,intent(in) :: a,b,c,d
  
  double precision :: offset_antispinint 
  
  abcd_antispinint=offset_antispinint(a,b,c,d,.True.,.True.,.True.,.True.)
end function

double precision function iajb_antispinint(i,a,j,b)
  
  BEGIN_DOC
  ! abcd antisymetric spinorbital integral
  END_DOC
  
  implicit none
  
  integer,intent(in) :: i,a,j,b
  
  double precision :: offset_antispinint 
  
  iajb_antispinint=offset_antispinint(i,a,j,b,.False.,.True.,.False.,.True.)
end function

double precision function abij_antispinint(a,b,i,j)
  
  BEGIN_DOC
  ! abcd antisymetric spinorbital integral
  END_DOC
  
  implicit none
  
  integer,intent(in) :: a,b,i,j
  
  double precision :: offset_antispinint 
  
  abij_antispinint=offset_antispinint(a,b,i,j,.True.,.True.,.False.,.False.)
end function
