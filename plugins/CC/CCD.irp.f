subroutine CCD

BEGIN_DOC
! CCD algorithm 
END_DOC

  implicit none

  double precision:: CCD_corr
  integer :: iteration_CCD,i,j,a,b
  double precision :: conv
  double precision,allocatable,dimension(:,:,:,:) :: ijkl_antispinint,ijab_antispinint,abij_antispinint,iajb_antispinint,abcd_antispinint,u,v,residual,t2
  iteration_CCD=0
  conv=thresh_CC+1d0
  call write_time(6)

! Allocate All constant arrays

  allocate(ijkl_antispinint(n_spin_occ,n_spin_occ,n_spin_occ,n_spin_occ))
  allocate(ijab_antispinint(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  allocate(abij_antispinint(n_spin_virt,n_spin_virt,n_spin_occ,n_spin_occ))
  allocate(iajb_antispinint(n_spin_occ,n_spin_virt,n_spin_occ,n_spin_virt))
  allocate(abcd_antispinint(n_spin_virt,n_spin_virt,n_spin_virt,n_spin_virt)) 
  allocate(t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))

! Build all constant array  

  call Build_ijkl_antispinint(ijkl_antispinint)
  call Build_ijab_antispinint(ijab_antispinint)
  call Build_abij_antispinint(abij_antispinint)
  call Build_iajb_antispinint(iajb_antispinint)
  call Build_abcd_antispinint(abcd_antispinint)
  call init_t2(t2,ijab_antispinint)
! Allocate all loop dependent arrays

  allocate(u(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  allocate(v(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  allocate(residual(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))

! Main loop

  do while(iteration_CCD < n_it_CC_max  .and. conv>thresh_CC)
    iteration_CCD += 1
    call write_int(6,iteration_CCD,'Current CCD iteration')

!   Build u and v

    call build_u(u,ijkl_antispinint,iajb_antispinint,abcd_antispinint,t2)
    call build_v(v,ijab_antispinint,t2)

!   Update residual

    call build_residual(residual,t2,abij_antispinint,u,v)

!   Update amplitudes and calculate current CCD correlation energy

    CCD_corr=0d0
    do i=1,n_spin_occ
      do j=1,n_spin_occ
        do a=1,n_spin_virt
          do b=1,n_spin_virt
            t2(i,j,a,b) -= residual(i,j,a,b)/ijab_D(i,j,a,b)
            CCD_corr += ijab_antispinint(i,j,a,b)*t2(i,j,a,b)
          end do
        end do
      end do
    end do    
    CCD_corr *= 0.25d0

!   Convergence criteria

    conv = maxval(dabs(residual))

!   Dump current results

    call write_double(6,CCD_corr,"Current CCD correction energy")
    call write_double(6,HF_energy+CCD_corr,"Current CCD corrected energy")
    call write_double(6,conv,"Current convergence")

  end do

! Deallocate arrays

  deallocate(u,v,ijkl_antispinint,ijab_antispinint,abij_antispinint,iajb_antispinint,abcd_antispinint)

! Dump final results
  
  call write_time(6)
  call write_int(6,iteration_CCD,'Number of CCD iteration')
  call write_double(6,conv,"Final convergence value")
  call write_double(6,CCD_corr,"Final CCD correction")
  call write_double(6,hf_energy+CCD_corr,"Final CCD corrected energy")

end subroutine CCD
