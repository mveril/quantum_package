double precision function get_u_el(i,j,a,b)
  
  BEGIN_DOC
  ! u element 
  END_DOC
  
  implicit none
  
  
  integer,intent(in) :: i,j,a,b
  double precision :: abcd_antispinint,antispinint,iajb_antispinint
  double precision :: pcd, pkl, pkc 
  integer :: k,l,c,d

  pcd=0d0
  pkl=0d0
  pkc=0d0
  do c = 1,n_spin_virt
    do d = 1,n_spin_virt
  !cd loop
      pcd += abcd_antispinint(a,b,c,d)*t_ampli(i,j,c,d)
    end do
  end do
  do k = 1,n_spin_occ
    do l= 1 ,n_spin_occ
    !Kl loop
      pkl += antispinint(k,l,i,j)*t_ampli(k,l,a,b)
    end do
    do c = 1,n_spin_virt
      !KC only part
      pkc -=iajb_antispinint(k,b,j,c)*t_ampli(i,k,a,c)
      pkc +=iajb_antispinint(k,a,j,c)*t_ampli(i,k,b,c)
      pkc -=iajb_antispinint(k,a,i,c)*t_ampli(j,k,b,c)
      pkc +=iajb_antispinint(k,b,i,c)*t_ampli(j,k,a,c)
    end do
  end do
  get_u_el=0.5d0*pcd+0.5d0*pkl+pkc
end function
