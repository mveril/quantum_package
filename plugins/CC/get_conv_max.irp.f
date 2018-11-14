double precision function get_conv_max(e_corr)
  BEGIN_DOC
  ! u element 
  END_DOC
  
  implicit none
  
  
  integer :: i,j,a,b
  double precision :: abij_antispinint,get_u_el,get_v_el,old_conv_max
  double precision,intent(in) :: e_corr
  old_conv_max=0d0
  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          ! Calculate one element
          get_conv_max=abij_antispinint(a,b,i,j)+ijab_D(i,j,a,b)*t_ampli(i,j,a,b)+get_u_el(i,j,a,b)+get_v_el(i,j,a,b)
          ! Keep the biggest element 
          get_conv_max=max(dabs(get_conv_max),dabs(old_conv_max))
          ! Continue to loop
          old_conv_max=get_conv_max
  
        end do
      end do
    end do
  end do
end function
