double precision function get_u_el(i,j,a,b)
  
  BEGIN_DOC
  ! u element 
  END_DOC
  
  implicit none
  
  
  integer,intent(in) :: i,j,a,b
  double precision :: abcd_antispinint,antispinint,iajb_antispinint
  double precision :: pcd, pkl, pkc 
  integer :: k,l,c,d,n_spinvirt

  pcd=0d0
  pkl=0d0
  pkc=0d0
  n_spinvirt=mo_tot_num*2-elec_num
  do c = 1,n_spinvirt
    do d = 1,n_spinvirt
  !cd loop
      pcd += abcd_antispinint(a,b,c,d)*t_coeff(i,j,c,d)
    end do
  end do
  do k = 1,elec_num
    do l= 1 ,elec_num
    !Kl loop
      pkl += antispinint(k,l,i,j)*t_coeff(k,l,a,b)
    end do
    do c = 1,n_spinvirt
      !KC only part
      pkc -=iajb_antispinint(k,b,j,c)*t_coeff(i,k,a,c)
      pkc +=iajb_antispinint(k,a,j,c)*t_coeff(i,k,b,c)
      pkc -=iajb_antispinint(k,a,i,c)*t_coeff(j,k,b,c)
      pkc +=iajb_antispinint(k,b,i,c)*t_coeff(j,k,a,c)
    end do
  end do
  get_u_el=0.5d0*pcd+0.5d0*pkl+pkc
end function
