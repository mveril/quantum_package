double precision function get_conv_max(e_corr)
  BEGIN_DOC
  ! u element 
  END_DOC
  
  implicit none
  
  
  integer :: i,j,a,b
  double precision :: abij_antispinint,get_u_el,old_conv_max
  double precision,intent(in) :: e_corr
  old_conv_max=0d0
  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          get_conv_max=abij_antispinint(a,b,i,j)+(ijab_D(i,j,a,b)+HF_energy-e_corr)*t_coeff(i,j,a,b)+get_u_el(i,j,a,b)
          get_conv_max=max(get_conv_max,old_conv_max)
          old_conv_max=get_conv_max
        end do
      end do
    end do
  end do
end function
