double precision function get_X1_el(k,l,i,j)
  
  BEGIN_DOC
  ! X1 intermediate Element 
  END_DOC
  
  implicit none
  
  
  integer,intent(in) :: k,l,i,j
  double precision :: ijab_antispinint 
  integer :: c,d
  get_X1_el=0d0
  do c = 1,n_spin_virt
    do d = 1,n_spin_virt
      get_X1_el += ijab_antispinint(k,l,c,d)*t_ampli(i,j,c,d)
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
  integer::k,l,d
  get_X2_el = 0d0
  do k = 1,n_spin_occ
    do l = 1,n_spin_occ
      do d = 1,n_spin_virt
        get_X2_el += ijab_antispinint(k,l,c,d)*t_ampli(k,l,b,d)
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
  integer :: l,c,d
  get_X3_el=0d0
  do l = 1, n_spin_occ
    do c = 1, n_spin_virt
      do d =1, n_spin_virt 
        get_X3_el += ijab_antispinint(k,l,c,d)*t_ampli(j,l,c,d)
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
  integer::k,c
  get_X4_el=0d0
  do k = 1,n_spin_occ
    do c = 1,n_spin_virt
      get_X4_el += ijab_antispinint(k,l,c,d)*t_ampli(i,k,a,c)
    end do
  end do
end function
