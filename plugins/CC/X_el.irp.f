double precision function get_X1_el(k,l,i,j)
  
  BEGIN_DOC
  ! X1 intermediate Element 
  END_DOC
  
  implicit none
  
  
  integer,intent(in) :: k,l,i,j
  double precision :: ijab_antispinint 
  integer :: c,d,n_spinvirt
  get_X1_el=0d0
  n_spinvirt=mo_tot_num*2-elec_num
  do c = 1,n_spinvirt
    do d = 1,n_spinvirt
      get_X1_el += ijab_antispinint(k,l,c,d)*t_coeff(k,l,c,d)
    end do
  end do
end function

double precision function get_X2_el(b,c)
  
  BEGIN_DOC
  ! X2 intermediate Element 
  END_DOC

  implicit none
  
  integer,intent(in) :: b,c
  double precision :: ijab_antispinint
  integer::k,l,d,n_spinvirt
  get_X2_el = 0d0
  n_spinvirt=mo_tot_num*2-elec_num
  do k = 1,elec_num
    do l = 1,elec_num
      do d = 1,n_spinvirt
        get_X2_el += ijab_antispinint(k,l,c,d)*t_coeff(k,l,b,d)
      enddo
    enddo
  enddo
end function
double precision function get_x3_el(k,j)

  BEGIN_DOC
  ! x3 intermediate element 
  END_DOC

  implicit none

  integer,intent(in) :: k,j
  double precision :: ijab_antispinint
  integer :: l,c,d,n_spinvirt
  get_X3_el=0d0
  n_spinvirt=mo_tot_num*2-elec_num
  do l = 1, elec_num
    do c = 1, n_spinvirt
      do d =1, n_spinvirt 
        get_X3_el += ijab_antispinint(k,l,c,d)*t_coeff(j,l,c,d)
      enddo
    enddo
  enddo
end function

double precision function get_x4_el(i,l,a,d)
  
  BEGIN_DOC
  ! x4 intermediate element 
  END_DOC

  implicit none
  
  integer,intent(in)::i,l,a,d
  double precision:: ijab_antispinint
  integer::k,c,n_spinvirt
  get_X4_el=0d0
  n_spinvirt=mo_tot_num*2-elec_num
  do k = 1,elec_num
    do c = 1,n_spinvirt
      get_X4_el += ijab_antispinint(k,l,c,d)*t_coeff(i,k,a,c)
    end do
  end do
end function
