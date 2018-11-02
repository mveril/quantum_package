subroutine CCD_Pople()

BEGIN_DOC
! Pople CCD algorithm 
END_DOC
  implicit none

  double precision::abij_antispinint,ijab_antispinint,get_u_el,get_v_el,CCD_corr
  integer :: iteration_CCD,i,j,a,b
  double precision :: conv,get_conv_max
  iteration_CCD=0
  conv=thresh_CC+1d0
  call write_time(6)
  do while(iteration_CCD < n_it_CC_max  .and. conv>thresh_CC)
      iteration_CCD += 1
      call write_int(6,iteration_CCD,'Current CCD iteration')
      CCD_corr=0d0
      do i=1,size(t_ampli,1)
        do j=1,size(t_ampli,2)
          do a=1,size(t_ampli,3)
            do b=1,size(t_ampli,4)
              t_ampli(i,j,a,b)=-(abij_antispinint(a,b,i,j)+get_u_el(i,j,a,b)+get_v_el(i,j,a,b))/ijab_D(i,j,a,b)
               CCD_corr+=ijab_antispinint(i,j,a,b)*t_ampli(i,j,a,b)
            end do
          end do
        end do
      end do

      CCD_corr*=0.25d0
      call write_double(6,CCD_corr,"Current CCD correction energy")
      call write_double(6,HF_energy+CCD_corr,"Current CCD corrected energy")
      conv=get_conv_max(CCD_corr)
      call write_double(6,conv,"Current convergence")
      TOUCH t_ampli 
  end do
  call write_time(6)
  call write_int(6,iteration_CCD,'Number of CCD iteration')
  call write_double(6,conv,"Final convergence value")
  call write_double(6,CCD_corr,"Final CCD correction")
  call write_double(6,hf_energy+CCD_corr,"Final CCD corrected energy")
end subroutine
