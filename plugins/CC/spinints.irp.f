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
  if (mod(p,2) == mod(r,2)) then
    spinint = get_mo_bielec_integral(p/2+1,q/2+1,r/2+1,s/2+1,mo_integrals_map)
  end if
end function

double precision function offset_spinint(p,q,r,s,isvirtp,isvirtq,isvirtr,isvirts)

  BEGIN_DOC
  ! ijab spinorbital integrals with virtual offset 
  END_DOC
  implicit none
  
  double precision :: spinint
  integer,intent(in) :: p,q,r,s
  logical,intent(in) :: isvirtp,isvirtq,isvirtr,isvirts
  integer::poff,qoff,roff,soff
  poff=p
  qoff=s
  roff=r
  soff=s

  if (isvirtp) then
    poff+=elec_num
  end if
  if (isvirtq) then
    qoff+=elec_num
  end if
  if (isvirtr) then
    roff+=elec_num
  end if
  if (isvirts) then
    soff+=elec_num
  end if
  offset_spinint=spinint(poff,qoff,roff,soff)
end function
double precision function ijab_spinint(i,j,a,b)
  
  BEGIN_DOC
  ! ijab spinorbital integral
  END_DOC
  
  implicit none
  
  double precision :: offset_spinint
  integer,intent(in) :: i,j,a,b
  ijab_spinint=offset_spinint(i,j,a,b,.False.,.False.,.True.,.True.) 
end function

double precision function abcd_spinint(a,b,c,d)
  
  BEGIN_DOC
  ! abcd spinorbital integral
  END_DOC
  
  implicit none
  double precision :: offset_spinint
  integer,intent(in) :: a,b,c,d
  abcd_spinint=offset_spinint(a,b,c,d,.True.,.True.,.True.,.True.) 
end function
