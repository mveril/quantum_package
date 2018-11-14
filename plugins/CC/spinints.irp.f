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
  qoff=q
  roff=r
  soff=s
  !Apply ofset if is virtual
  if (isvirtp) then
    poff+=n_spin_occ
  end if
  if (isvirtq) then
    qoff+=n_spin_occ
  end if
  if (isvirtr) then
    roff+=n_spin_occ
  end if
  if (isvirts) then
    soff+=n_spin_occ
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
