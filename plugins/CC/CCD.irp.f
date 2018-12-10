subroutine CCD

BEGIN_DOC
! CCD algorithm 
END_DOC

  implicit none

  double precision:: CCD_corr
  integer :: iteration_CCD,i,j,a,b
  double precision :: conv
  double precision,allocatable,dimension(:,:,:,:) :: ijkl_antispinint,ijab_antispinint,abij_antispinint,iajb_antispinint,abcd_antispinint,u,v,residual,t2,Delta,inv_Delta
  iteration_CCD=0
  conv = thresh_CC+1d0
  call Write_Time(6)

! Allocate All constant arrays

  allocate(ijkl_antispinint(n_spin_occ,n_spin_occ,n_spin_occ,n_spin_occ))
  allocate(ijab_antispinint(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  allocate(abij_antispinint(n_spin_virt,n_spin_virt,n_spin_occ,n_spin_occ))
  allocate(iajb_antispinint(n_spin_occ,n_spin_virt,n_spin_occ,n_spin_virt))
  allocate(abcd_antispinint(n_spin_virt,n_spin_virt,n_spin_virt,n_spin_virt)) 
  allocate(t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  allocate(Delta(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  allocate(inv_Delta(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))

! Build all constant arrays

  call Build_ijkl_antispinint(ijkl_antispinint)
  if (Debug) then
    print *, "OOOO"
    call matout(n_spin_occ*2,n_spin_occ*2,ijkl_antispinint)
  end if
  call Build_ijab_antispinint(ijab_antispinint)
  if (Debug) then
    print *, "OOVV"
    call matout(n_spin_occ*2,n_spin_virt*2,ijab_antispinint)
  end if
  call Build_abij_antispinint(abij_antispinint)
  if (Debug) then
    print *, "VVOO"
    call matout(n_spin_virt*2,n_spin_occ*2,abij_antispinint)
  end if
  call Build_iajb_antispinint(iajb_antispinint)
  if (Debug) then
    print *, "OVOV"
    call matout(n_spin_occ*2,n_spin_virt*2,iajb_antispinint)
  end if
  call Build_abcd_antispinint(abcd_antispinint)
  if (Debug) then
    print *, "VVVV"
    call matout(n_spin_virt*2,n_spin_virt*2,abcd_antispinint)
  end if
  call Build_Delta(Delta)
  if (Debug) then
    print *, "Delta"
    call matout(n_spin_occ*2,n_spin_virt*2,Delta)
  end if
  inv_Delta(:,:,:,:) = 1d0/Delta(:,:,:,:)
  if (Debug) then
    print *, "inv_delta"
    call matout(n_spin_occ*2,n_spin_occ*2,inv_delta)
  end if
  call Init_t2(t2, ijab_antispinint, inv_Delta)
  if (Debug) then
    print *, "t2"
    call matout(n_spin_occ*2,n_spin_virt*2,t2)
  end if


! Allocate all loop dependent arrays

  allocate(u(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  allocate(v(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  allocate(residual(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))

! Main loop

  do while(iteration_CCD < n_it_CC_max  .and. conv>thresh_CC)
    iteration_CCD += 1
    call Write_Int(6,iteration_CCD,'Current CCD iteration')

!   Build u and v

    call Build_U(u,ijkl_antispinint,iajb_antispinint,abcd_antispinint,t2)
    if (Debug) then
      print *, "U"
      call matout(n_spin_occ*2,n_spin_virt*2,U)
    end if
    call Build_V(v,ijab_antispinint,t2)
    if (Debug) then
      print *, "V"
      call matout(n_spin_occ*2,n_spin_virt*2,V)
    end if

!   Update residual

    call Build_Residual(residual, t2, abij_antispinint, u, v, Delta)
    if (Debug) then
      print *, "residual"
      call matout(n_spin_occ*2,n_spin_virt*2,residual)
    end if
!   Update amplitudes and calculate current CCD correlation energy

    CCD_corr = 0d0
    do i=1,n_spin_occ
      do j=1,n_spin_occ
        do a=1,n_spin_virt
          do b=1,n_spin_virt
            t2(i,j,a,b) -= residual(i, j, a, b)*inv_Delta(i, j, a, b)
            CCD_corr += ijab_antispinint(i, j, a, b)*t2(i, j, a, b)
          end do
        end do
      end do
    end do    
    if (Debug) then
      print *, "t2"
      call matout(n_spin_occ*2,n_spin_virt*2,t2)
    end if
    CCD_corr *= 0.25d0

!   Convergence criteria

    conv = maxval(dabs(residual))

!   Dump current results

    call Write_Double(6, CCD_corr, "Current CCD correction energy")
    call Write_Double(6, HF_energy+CCD_corr, "Current CCD corrected energy")
    call Write_Double(6, conv, "Current convergence")

  end do

! Deallocate arrays

  deallocate(u,v,ijkl_antispinint,ijab_antispinint,abij_antispinint,iajb_antispinint,abcd_antispinint)

! Dump final results
  
  call write_time(6)
  call write_int(6, iteration_CCD, 'Number of CCD iteration')
  call write_double(6, conv, "Final convergence value")
  call write_double(6, CCD_corr, "Final CCD correction")
  call write_double(6, hf_energy+CCD_corr, "Final CCD corrected energy")

end subroutine CCD
