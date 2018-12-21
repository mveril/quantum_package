subroutine CCD

BEGIN_DOC
! Coupled cluster doubles (CCD) algorithm.
! See Pople et al. IJQC 14 545 (1978) for more details
END_DOC

  implicit none

  integer                      :: it
  double precision             :: E_HF,E_MP2,Ec_MP2,E_CCD,Ec_CCD
  double precision             :: conv,get_mo_bielec_integral,Conv_Spin_Index

  double precision,allocatable :: oooo_db_spin_int(:,:,:,:)
  double precision,allocatable :: oovv_db_spin_int(:,:,:,:)
  double precision,allocatable :: vvoo_db_spin_int(:,:,:,:)
  double precision,allocatable :: ovov_db_spin_int(:,:,:,:)
  double precision,allocatable :: vvvv_db_spin_int(:,:,:,:)
  double precision,allocatable :: u(:,:,:,:)
  double precision,allocatable :: v(:,:,:,:)
  double precision,allocatable :: r(:,:,:,:)
  double precision,allocatable :: t2(:,:,:,:)
  double precision,allocatable :: oovv_Delta(:,:,:,:)

  integer                      :: i,j,a,b

! Initialize variables

  it   = 0
  conv = 1d0

! Timing

  call write_time(6)

! Allocate arrays

  allocate(oooo_db_spin_int(n_spin_occ,n_spin_occ,n_spin_occ,n_spin_occ))
  allocate(oovv_db_spin_int(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  allocate(vvoo_db_spin_int(n_spin_virt,n_spin_virt,n_spin_occ,n_spin_occ))
  allocate(ovov_db_spin_int(n_spin_occ,n_spin_virt,n_spin_occ,n_spin_virt))
  allocate(vvvv_db_spin_int(n_spin_virt,n_spin_virt,n_spin_virt,n_spin_virt)) 
  allocate(t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  allocate(oovv_Delta(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))

! Build integral arrays

  call build_oooo_db_spin_int(oooo_db_spin_int)

! Compute Hartree-Fock energy

  call compute_HF_energy(oooo_db_spin_int,E_HF)

  call build_oovv_db_spin_int(oovv_db_spin_int)
  if (DEBUG) then
  ! do i=1,n_spin_occ
  !   do j=1,n_spin_occ
  !     do a=1,n_spin_virt
  !       do b=1,n_spin_virt
  !         print *, i,j,a,b,  oovv_db_spin_int(i,j,a,b)
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  !  stop
  end if

  call build_vvoo_db_spin_int(vvoo_db_spin_int)
  call build_ovov_db_spin_int(ovov_db_spin_int)
  call build_vvvv_db_spin_int(vvvv_db_spin_int)

! Compute denominators

  call build_oovv_Delta(oovv_Delta)

  if (Debug)
!   Compute Hartree-Fock energy
    call compute_MP2_energy(E_HF,oovv_Delta,oovv_db_spin_int,Ec_MP2,E_MP2)
  End if
! Initialize amplitudes

  call init_t2(t2,oovv_db_spin_int,oovv_Delta)

! Allocate all loop dependent arrays

  allocate(u(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  allocate(v(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  allocate(r(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))

! Main loop

  do while(it < max_it_CC  .and. conv > thresh_CC) 

!   Increment 

    it += 1

!   Timing

    call write_int(6,it,'CCD iteration n. ')

!   Build linear array u

    call build_u(u,oooo_db_spin_int,ovov_db_spin_int,vvvv_db_spin_int,t2)

!   Build quadratic array v
    call build_v(v,oovv_db_spin_int,t2)

!   Update residual

    call build_residual(r,t2,vvoo_db_spin_int,u,v,oovv_Delta)

!   Update amplitudes and calculate current CCD correlation energy

    Ec_CCD = 0d0
    do i=1,n_spin_occ
      do j=1,n_spin_occ
        do a=1,n_spin_virt
          do b=1,n_spin_virt
            t2(i,j,a,b) -= r(i,j,a,b)/oovv_Delta(i,j,a,b)
            Ec_CCD      += oovv_db_spin_int(i,j,a,b)*t2(i,j,a,b)
          end do
        end do
      end do
    end do    
    Ec_CCD *= 0.25d0

!   Convergence criteria

    conv = maxval(abs(r))

!   Dump current results

    call write_double(6,Ec_CCD,            "Current CCD correction energy")
    call write_double(6,HF_energy + Ec_CCD,"Current CCD corrected energy" )
    call write_double(6,conv,              "Current convergence"          )

  end do

! Deallocate arrays

  deallocate(u,v,oooo_db_spin_int,oovv_db_spin_int,vvoo_db_spin_int,ovov_db_spin_int,vvvv_db_spin_int)

! Dump final results
  
  call write_time  (6)
  call write_int   (6,it,              "Number of CCD iteration"   )
  call write_double(6,conv,            "Final convergence value"   )
  call write_double(6,Ec_CCD,          "Final CCD correction"      )
  call write_double(6,HF_energy+Ec_CCD,"Final CCD corrected energy")

end subroutine CCD
