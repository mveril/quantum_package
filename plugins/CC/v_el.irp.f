double precision function get_v_el(i,j,a,b)
  
  BEGIN_DOC
  ! V element 
  END_DOC
  
  implicit none
  
  
  integer,intent(in) :: i,j,a,b
  double precision :: get_X1_el,get_X2_el,get_X3_el,get_X4_el
  double precision :: pkl,pkc,pk,pc 
  integer :: k,l,c,n_spinvirt
  pkl=0d0
  pkc=0d0
  pk=0d0
  pc=0d0
  n_spinvirt=mo_tot_num*2-elec_num
  do k = 1,elec_num
    do l = 1,elec_num
  !KL loop
      pkl += get_X1_el(k,l,i,j)*t_coeff(k,l,a,b)
    end do
    do c = 1,n_spinvirt
  !KC loop
      pkc += get_X4_el(i,k,a,c)*t_coeff(j,k,b,c)+get_X4_el(i,k,b,c)*t_coeff(k,j,a,c)
    end do
    !K only part
    pk +=get_X3_el(k,j)*t_coeff(i,k,a,b)+get_X3_el(k,i)*t_coeff(k,j,a,b)
  end do
  do c = 1,n_spinvirt
  !c only part
    pc+=get_X2_el(b,c)*t_coeff(i,j,a,c)+get_X2_el(a,c)*t_coeff(i,j,c,b)
  end do
  get_v_el=0.25d0*pkl-0.5d0*pc-0.5d0*pk+pkc
end function

