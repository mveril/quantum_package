subroutine CCD_Pople()

BEGIN_DOC
! Pople CCD algorithm 
END_DOC
  implicit none

  double precision::abij_antispinint,ijab_antispinint,get_u_el,get_v_el,D_ijab,CCD_energy
  integer :: iteration_CCD,i,j,a,b
  iteration_CCD=0
  call write_time(6)
  do while(iteration_CCD < 1000)
      iteration_CCD += 1

      call write_int(6,iteration_CCD,'Current CCD iteration')
      CCD_energy=0d0

      do i=1,size(t_coeff,1)
        do j=1,size(t_coeff,2)
          do a=1,size(t_coeff,3)
            do b=1,size(t_coeff,4)
              t_coeff(i,j,a,b)=-(abij_antispinint(a,b,i,j)+get_u_el(i,j,a,b)+get_v_el(i,j,a,b))/D_ijab(i,j,a,b)
               CCD_energy+=ijab_antispinint(i,j,a,b)*t_coeff(i,j,a,b)
            end do
          end do
        end do
      end do

      CCD_energy = 0.25d0*CCD_energy
      call write_double(6,HF_energy+CCD_energy,"CCD corrected energy")
      TOUCH t_coeff 
  end do
  call write_time(6)
  call write_double(6,hf_energy+CCD_energy,"Final CCD corrected energy")
end subroutine
