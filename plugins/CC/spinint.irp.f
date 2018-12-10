double precision function spinint(p,q,r,s)
  
  BEGIN_DOC
  ! spinorbital integrals 
  END_DOC
  
  implicit none
  
  PROVIDE mo_bielec_integrals_in_map 
  
  integer,intent(in) :: p,q,r,s
  double precision :: get_mo_bielec_integral
  integer :: c,d,Conv_Spin_Index

  ! Initialize variables 

  spinint = 0d0

  ! Non zero only if same spin

  if (mod(p,2) == mod(r,2) .and. mod(q,2) == mod(s,2)) &
    spinint = get_mo_bielec_integral(Conv_Spin_Index(p),Conv_Spin_Index(q),Conv_Spin_Index(r),Conv_Spin_Index(s),mo_integrals_map)

end function
